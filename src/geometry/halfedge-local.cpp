
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
void printHalfedgeData(Halfedge_Mesh::HalfedgeRef h);
void printVertexPositions(Halfedge_Mesh::FaceRef f) {
  Halfedge_Mesh::HalfedgeRef h = f->halfedge; // get the first halfedge of the face
  std::cout<<std::endl<<std::endl;
  do {
    printHalfedgeData(h);
    h = h->next;               // move to the next halfedge around the face
  } while (h != f->halfedge);  // keep going until we're back at the beginning
}
void printHalfedgeData(Halfedge_Mesh::HalfedgeRef h) {
    Halfedge_Mesh::VertexRef v = h->vertex;  
    Halfedge_Mesh::HalfedgeRef nextHalfedge = h->next;
    Halfedge_Mesh::FaceRef face = h->face;

    std::cout<<"------------------ Current Halfedge Data -----------------"<<std::endl;
    std::cout << "Vertex ID: " << v->id << " Position: " << v->position
              << " | Halfedge ID: " << h->id
              << " | Next Halfedge ID: " << nextHalfedge->id
              << " | Face ID: " << face->id
              << " | Edge ID: " << h->edge->id <<std::endl;
    
    std::cout<< "Twin Halfedge ID: " << h->twin->id 
             << " | Twin Vertex ID: " << h->twin->vertex->id 
             << " | Twin Face ID: " << h->twin->face->id <<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
}
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge
  HalfedgeRef initialHalf = e->halfedge;
  HalfedgeRef initialTwinHalf = e->halfedge->twin;

  HalfedgeRef in = initialHalf->next;
  VertexRef iv = initialHalf->vertex;

  HalfedgeRef itn = initialTwinHalf->next;
  VertexRef itv = initialTwinHalf->vertex;
  HalfedgeRef itnn = itn->next;

  HalfedgeRef itp = getPrev(e->halfedge->twin);

  Vec3 vPos(e->center());
  VertexRef v = emplace_vertex();
  v->position = vPos;
  //old face left bt
  EdgeRef ed1 = emplace_edge();
  HalfedgeRef h1 = emplace_halfedge();
  bool onBoundary = e->halfedge->face->id == e->halfedge->twin->next->twin->face->id ||
                    e->halfedge->twin->face->id == e->halfedge->next->twin->face->id;
  if(onBoundary)
  {
    //debug all initial value, for initialHalf face, initialTwinHalf face, next halfedge of both
    printVertexPositions(initialHalf->face);
    printVertexPositions(initialTwinHalf->face);
  }
 

  v->halfedge = h1;
  h1->vertex = v;
  h1->next = in->next;
  h1->edge = ed1;
  h1->face = initialHalf->face;
  h1->face->halfedge = h1;
  ed1->halfedge = h1;

  initialHalf->next = h1;


  //new face left tp
  FaceRef f1 = emplace_face();
  EdgeRef ed2 = emplace_edge();
  HalfedgeRef h2 = emplace_halfedge();
  HalfedgeRef h3 = emplace_halfedge();

  f1->halfedge = h2;
  ed2->halfedge = h2;


  h2->vertex = v;
  h2->edge = ed2;
  h2->next = in;
  h2->face = f1;

  h3->vertex = in->next->vertex;
  h3->next = h2;
  h3->edge = ed1;
  h3->face = f1;
  h3->twin = h1;
  h1->twin = h3;

  in->next = h3;
  in->face = f1;

  if(onBoundary)
  {
    
    initialTwinHalf->vertex = v;
    HalfedgeRef h5 = emplace_halfedge();
    h5->vertex = itv;
    h5->edge = ed2;
    h5->twin = h2;
    
    //id 33 eh o twin do h2
    h5->next = initialTwinHalf;
    getPrev(initialTwinHalf)->next = h5;
    h5->face = initialTwinHalf->face;
    h2->twin = h5;

    
    interpolate_data({iv, itv}, v);
    return v;
  }
  //old face rigth bt
  EdgeRef ed3 = emplace_edge();
  HalfedgeRef h4 = emplace_halfedge();

  h4->vertex = itp->vertex;
  h4->next = initialTwinHalf;
  
  h4->face = initialTwinHalf->face;
  h4->face->halfedge = h4;
  h4->edge = ed3;
  initialTwinHalf->vertex = v;

  itn->next = h4;
  ed3->halfedge = h4;

  //new face rigth tp
  FaceRef f2 = emplace_face();
  HalfedgeRef h5 = emplace_halfedge();
  HalfedgeRef h6 = emplace_halfedge();

  f2->halfedge = h5;


  h5->vertex = itv;
  h5->edge = ed2;
  h5->next = h6;
  h5->face = f2;
  h5->twin = h2;
  h2->twin = h5;

  h6->vertex = v;
  h6->next = itnn;
  h6->edge = ed3;
  h6->face = f2;
  h6->twin = h4;
  h4->twin = h6;

  itnn->face = f2;
  itnn->next = h5;
  interpolate_data({iv, itv}, v);

  return v;

}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	(void)f;
    return std::nullopt;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 */


std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
  
	//A2L1: Flip Edge
  if(e->halfedge->face->id == e->halfedge->twin->next->twin->face->id || 
     e->halfedge->twin->face->id == e->halfedge->next->twin->face->id)
  {
    return std::nullopt;
  }
  HalfedgeRef initialHalf = e->halfedge;
  HalfedgeRef initialTwinHalf = e->halfedge->twin;

  HalfedgeRef in = initialHalf->next;
  VertexRef iv = initialHalf->vertex;

  HalfedgeRef itn = initialTwinHalf->next;
  VertexRef itv = initialTwinHalf->vertex;

  HalfedgeRef ip = getPrev(e->halfedge);
  HalfedgeRef itp = getPrev(e->halfedge->twin);

  initialHalf->vertex = initialTwinHalf->next->next->vertex;
  initialTwinHalf->vertex = initialHalf->next->next->vertex;

  initialHalf->next = in->next;
  initialTwinHalf->next = itn->next;

  in->next = initialTwinHalf;
  itn->next = initialHalf;
  
  in->face = initialTwinHalf->face;
  itn->face = initialHalf->face;
  
  itp->next = in;
  ip->next = itn;

  
  iv->halfedge = ip->twin;
  itv->halfedge = itp->twin;
  
  e->halfedge->face->halfedge = ip;
  e->halfedge->twin->face->halfedge = itp;
  return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */

int countFaceSides(Halfedge_Mesh::FaceRef f)
{
  Halfedge_Mesh::HalfedgeRef h = f->halfedge->next;
  int counter = 1;
  while(f->halfedge != h)
  {
    counter ++;
    h = h->next;
  }
  return counter;
}
bool possibleEdgeCollapse(Halfedge_Mesh::EdgeRef e)
{
  if(countFaceSides(e->halfedge->face) > 3 || countFaceSides(e->halfedge->twin->face) > 3)
    return true;
  
  return false;
}
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edgea

  //if(!possibleEdgeCollapse(e))
    //return e->halfedge->vertex;
  std::set<EdgeRef> updateEdges;

  VertexRef newVertex = emplace_vertex();
  newVertex->position = e->center();

  VertexRef oldVertex1 = e->halfedge->vertex;
  VertexRef oldVertex2 = e->halfedge->twin->vertex;

  HalfedgeRef h1 = e->halfedge;
  HalfedgeRef h2 = h1->twin;

  FaceRef h1F = h1->face;
  FaceRef h2F = h2->face;

  HalfedgeRef h1Prev = getPrev(h1);
  HalfedgeRef h2Prev = getPrev(h2);

  HalfedgeRef hIdx = h1;
  bool initial = true;

  while(hIdx != h1 || initial)
  {
    updateCurrentNextAndPrevVertex(hIdx, oldVertex1, oldVertex2, newVertex, updateEdges);
    //std::cout<<"updating halfedge "<<hIdx->id<<std::endl;
    hIdx = getPrev(hIdx)->twin;
    initial = false;
  }
  initial = true;
  hIdx = h2;
  while(hIdx != h2 || initial)
  {
    updateCurrentNextAndPrevVertex(hIdx, oldVertex1, oldVertex2, newVertex, updateEdges);
    //std::cout<<"updating halfedge "<<hIdx->id<<std::endl;
    hIdx = getPrev(hIdx)->twin;
    initial = false;
  }
  

  h1Prev->next = h1->next;
  h2Prev->next = h2->next;

  std::set<EdgeRef> alreadyErased;
  for(auto edge: updateEdges)
  {
    for(auto edge1: updateEdges)
    {
      if(edge == edge1 || alreadyErased.count(edge) || alreadyErased.count(edge1)) 
        continue;

      bool overlapping = edge->id != edge1->id && 
      ((edge->halfedge->vertex == edge1->halfedge->vertex && edge->halfedge->twin->vertex == edge1->halfedge->twin->vertex) ||
      (edge->halfedge->vertex == edge1->halfedge->twin->vertex && edge->halfedge->twin->vertex == edge1->halfedge->vertex));

      if(overlapping)
      {
        HalfedgeRef h;
        HalfedgeRef h0;
      
        if(edge->halfedge->face->id != h1F->id && edge->halfedge->face->id != h2F->id )
        {
          h = edge->halfedge;
        }
        else
        {
          h = edge->halfedge->twin;
        }
        edge->halfedge = h;

        if(edge1->halfedge->face->id != h1F->id && edge1->halfedge->face->id != h2F->id )
        {
          h0 = edge1->halfedge;
        }
        else
        {
          h0 = edge1->halfedge->twin;
        }
        // Atualmente,quando o objeto colapsa para 2D, temos problemas ao subdividir ainda mais. O problema eh
        // que as faces colapsam exatamente no mesmo plano, e eu não sei como atribuir os halfedges de um plano para um objeto 3D
        // Posso talvez deletar uma face e para o twin desses halfedges, atribuir os halfedges da face mantida como seus twins
        if(h0->twin->face != h->twin->face)
          continue;
        alreadyErased.insert(edge1);
        alreadyErased.insert(edge);
          
        h0->twin->vertex->halfedge = h0->twin->next;
        getPrev(h0->twin)->next = h0->twin->next;
        printHalfedgeData(h0->twin);
        if(h->twin->vertex == h->next->vertex)
        {
          h->next->vertex->halfedge = h->next;
        }
        else if(h->twin->vertex == getPrev(h)->vertex)
        {
          getPrev(h)->vertex->halfedge = getPrev(h);
        }
        
        getPrev(h->twin)->next = h->twin->next;
        erase_face(h->twin->face);
        erase_halfedge(h->twin);
        h0->twin->face->halfedge = getPrev(h0->twin);
        
        erase_halfedge(h0->twin);
        h->twin = h0;
        h0->twin = h;
        h0->edge = edge;
        h->edge = edge;
        h->vertex->halfedge = h;
        h0->vertex->halfedge = h0;
        erase_edge(edge1);
      }
    }
  }
  
  if(h1->face->id < 1073741888)
    h1->face->halfedge = h1Prev;
  if(h2->face->id < 1073741888)
    h2->face->halfedge = h2Prev;
  erase_edge(e);
  erase_vertex(oldVertex1);
  erase_vertex(oldVertex2);
  erase_halfedge(h1);
  erase_halfedge(h2);

   

  return newVertex;
	
  
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move
	
}


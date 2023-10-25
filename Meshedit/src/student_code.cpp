#include "student_code.h"
#include "halfEdgeMesh.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  {
      std::vector<Vector2D> inter = std::vector<Vector2D>();
      for (int i = 0; i < points.size() - 1; i ++) {
          Vector2D p = lerp2D(points[i], points[i+1], t);
          inter.push_back(p);
      }
    // TODO Part 1.
      return inter;
  }

Vector2D BezierCurve::lerp2D(Vector2D p1, Vector2D p2, float t) {
    return ((1 - t) * p1) + (t * p2);
}

Vector3D BezierPatch::lerp3D(Vector3D p1, Vector3D p2, float t) const {
    return ((1 - t) * p1) + (t * p2);
}


  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
      std::vector<Vector3D> inter = std::vector<Vector3D>();
      for (int i = 0; i < points.size() - 1; i ++) {
          Vector3D p = lerp3D(points[i], points[i+1], t);
          inter.push_back(p);
      }
    // TODO Part 2.
      return inter;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
      std::vector<Vector3D> v = points;
      for (int i = 0; i < points.size() - 1; i++) {
          v = evaluateStep(v, t);
      }
    // TODO Part 2.
    return v[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const
  {
      std::vector<Vector3D> vec;
      for (int i = 0; i < controlPoints.size(); i ++) {
          vec.push_back(evaluate1D(controlPoints[i], u));
      }
      Vector3D value = evaluate1D(vec, v);
    // TODO Part 2.
      return value;
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      HalfedgeCIter half = halfedge();
      Vector3D sum;
      do {
          Vector3D u = (half->twin()->vertex()->position) - position;
          Vector3D v = (half->twin()->next()->twin()->vertex()->position) - position;
          sum += cross(u, v);
          half = half->twin()->next();
      }
      while (half != halfedge());

      sum.normalize();
      return -sum;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
      if (!e0->isBoundary()) {
          HalfedgeIter h0 = e0->halfedge();
          HalfedgeIter h1 = h0->next();
          HalfedgeIter h2 = h1->next();
          HalfedgeIter h3 = h0->twin();
          HalfedgeIter h4 = h3->next();
          HalfedgeIter h5 = h4->next();
          HalfedgeIter h6 = h1->twin();
          HalfedgeIter h7 = h2->twin();
          HalfedgeIter h8 = h4->twin();
          HalfedgeIter h9 = h5->twin();

          VertexIter v0 = h0->vertex();
          VertexIter v1 = h3->vertex();
          VertexIter v2 = h5->vertex();
          VertexIter v3 = h2->vertex();

          EdgeIter e1 = h5->edge();
          EdgeIter e2 = h4->edge();
          EdgeIter e3 = h2->edge();
          EdgeIter e4 = h1->edge();

          FaceIter f0 = h0->face();
          FaceIter f1 = h3->face();

          // For each vertex, edge, and face, set its halfedge pointer.
          v0->halfedge() = h4;
          v1->halfedge() = h1;
          v2->halfedge() = h0;
          v3->halfedge() = h3;

          e0->halfedge() = h0;
          e1->halfedge() = h5;
          e2->halfedge() = h4;
          e3->halfedge() = h2;
          e4->halfedge() = h1;

          f0->halfedge() = h2;
          f1->halfedge() = h1;

          // Setting all halfedges
          h0->setNeighbors(h2, h3, v2, e0, f0);
          h1->setNeighbors(h3, h6, v1, e4, f1);
          h2->setNeighbors(h4, h7, v3, e3, f0);
          h3->setNeighbors(h5, h0, v3, e0, f1);
          h4->setNeighbors(h0, h8, v0, e2, f0);
          h5->setNeighbors(h1, h9, v2, e1, f1);

          return e0;
      } else {
          return e0;
      }
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();

      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h5->vertex();
      VertexIter v3 = h2->vertex();

      EdgeIter e1 = h5->edge();
      EdgeIter e2 = h4->edge();
      EdgeIter e3 = h2->edge();
      EdgeIter e4 = h1->edge();

      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();

      // new elements
      VertexIter vm = newVertex();
      vm->isNew = true;

      HalfedgeIter h10 = newHalfedge();
      HalfedgeIter h11 = newHalfedge();
      HalfedgeIter h12 = newHalfedge();
      HalfedgeIter h13 = newHalfedge();
      HalfedgeIter h14 = newHalfedge();
      HalfedgeIter h15 = newHalfedge();

      EdgeIter e5 = newEdge();
      EdgeIter e6 = newEdge();
      EdgeIter e7 = newEdge();
      e5->isNew = true;
      e6->isNew = false; // since the original spline split into edge 6,
                         // technically it is not "new"
      e7->isNew = true;

      FaceIter f2 = newFace();
      FaceIter f3 = newFace();

      // e0, e1, e2, e3 don't change
      e5->halfedge() = h10;
      e6->halfedge() = h12;
      e7->halfedge() = h14;

      f0->halfedge() = h1;
      f1->halfedge() = h3;
      f2->halfedge() = h2;
      f3->halfedge() = h4;

      v0->halfedge() = h7;
      v1->halfedge() = h9;
      v2->halfedge() = h8;
      v3->halfedge() = h6;
      vm->halfedge() = h0;

      vm->position = (v0->position + v1->position)/2;

      h10->twin() = h11;
      h11->twin() = h10;
      h12->twin() = h13;
      h13->twin() = h12;
      h14->twin() = h15;
      h15->twin() = h14;

      h1->next() = h10;
      h10->next() = h0;
      h3->next() = h15;
      h15->next() = h5;
      h2->next() = h12;
      h12->next() = h11;
      h11->next() = h2;
      h4->next() = h14;
      h14->next() = h13;
      h13->next() = h4;

      h0->vertex() = vm;
      h11->vertex() = vm;
      h13->vertex() = vm;
      h15->vertex() = vm;
      h10->vertex() = v3;
      h12->vertex() = v0;
      h14->vertex() = v2;

      h10->edge() = e5;
      h11->edge() = e5;
      h12->edge() = e6;
      h13->edge() = e6;
      h14->edge() = e7;
      h15->edge() = e7;

      h10->face() = f0;
      h15->face() = f1;
      h2->face() = f2;
      h11->face() = f2;
      h12->face() = f2;
      h4->face() = f3;
      h13->face() = f3;
      h14->face() = f3;

      return vm;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)

    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.


      // Step A: Compute the positions of both new and old vertices using the original mesh. We want to perform these computations before subdivision because traversing a coarse mesh is much easier than traversing a subdivided mesh with more elements

      //Loop through all the vertices and find n and new position
      for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          HalfedgeIter half = v->halfedge();
          int n = 0;
          Vector3D sum = 0;
          float u = 0;

          do {
              Vector3D u = (half->twin()->vertex()->position);
              half = half->twin()->next();
              sum += u;
              n += 1;
          }
          while (half != v->halfedge());

          if (n == 3) {
              u = 3.0/16.0;
          } else {
              u = 3.0/(8.0*n);
          }

          Vector3D new_position = (1 - n * u) * v->position + u * sum;
          v->newPosition = new_position;
          v->isNew = false;
      }

      for(EdgeIter edges = mesh.edgesBegin(); edges != mesh.edgesEnd(); edges++) {
          HalfedgeIter half = edges->halfedge();
          VertexIter A = half->vertex();
          VertexIter B = half->twin()->vertex();
          VertexIter C = half->next()->next()->vertex();
          VertexIter D = half->twin()->next()->next()->vertex();

          Vector3D vertex_position = 3.0/8.0 * (A->position + B->position) + 1.0/8.0 * (C->position + D->position);

          edges->newPosition = vertex_position;
          edges->isNew = false;
      }

      // Step B: Subdivide the original mesh via edge splits
      for (EdgeIter edges = mesh.edgesBegin(); edges != mesh.edgesEnd(); edges++) {
          // If both vertices on the edge are old, then the edge is old and we can split it
          // checking edges->isNew may be weird due to the split edge function. The ta told me that
          // one of the edges is not technically new, and so we should not be splitting it. So instead
          // just check if the two vertices are old to make sure edge is old.
          if (!edges->halfedge()->vertex()->isNew && !edges->halfedge()->twin()->vertex()->isNew) {
            VertexIter v = mesh.splitEdge(edges);
            v->newPosition = edges->newPosition;
          }
      }

      // Now flip
      for (EdgeIter edges = mesh.edgesBegin(); edges != mesh.edgesEnd(); edges++) {
          if (edges->isNew) {
              if (edges->halfedge()->vertex()->isNew != edges->halfedge()->twin()->vertex()->isNew) {
                mesh.flipEdge(edges);
              }
          }
      }

      // Step C: Copy the new vertex positions into final Vertex::position.
      for(VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++ ) {
        v->position = v->newPosition;
      }
  }
}

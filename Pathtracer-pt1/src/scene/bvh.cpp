#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.
  // https://cs184.eecs.berkeley.edu/sp22/lecture/10-62/ray-tracing-acceleration
    
    BBox bbox;
    int prims;
    BBox centroids;
    for (auto p = start; p != end; p++) {
        BBox bb = (*p)->get_bbox();
        bbox.expand(bb);
        centroids.expand(bb.centroid());
        prims += 1;
    }
    
    if (prims <= max_leaf_size) {
        BVHNode *node = new BVHNode(bbox);
        node->start = start;
        node->end = end;
        return node;
    }
    
    // construct a bounding box of all the centroids
    // bbox.extend gives a vector beginning at min corner and ending at max corner
    // take largest of these coordinates, this is the longest axis we can split on
    int axis;
    Vector3D longest = centroids.extent;
    // setting the longest axis
    if (longest.x > longest.y && longest.x > longest.z) {
        axis = 0;
    } else if (longest.y > longest.x && longest.y > longest.z) {
        axis = 1;
    } else {
        axis = 2;
    }
    
    double splitPoint = (centroids.max[axis] + centroids.min[axis]) / 2;
    
    vector<Primitive *> right;
    vector<Primitive *> left;
    
    vector<Primitive *>::iterator ordered;
    ordered = partition(start, end, [&axis, &splitPoint](Primitive *p) {
        return p->get_bbox().centroid()[axis] < splitPoint;
    });
    
    if (ordered == start || ordered == end) {
        BVHNode *node = new BVHNode(bbox);
        node->start = start;
        node->end = end;
        node->l = NULL;
        node->r = NULL;
        return node;
    }
    
    BVHNode *node = new BVHNode(bbox);
    node->l = construct_bvh(start, ordered, max_leaf_size);
    node->r = construct_bvh(ordered, end, max_leaf_size);
    
  return node;


}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
    double t0;
    double t1;
    BBox bb = node->bb;
    if (!(bb.intersect(ray, t0, t1))) {
        return false;
    }
    if (node->isLeaf()) {
        for (auto p = node->start; p != node->end; p++) {
          total_isects++;
          if ((*p)->has_intersection(ray))
            return true;
        }
    }
    bool i1 = has_intersection(ray, node->l);
    bool i2 = has_intersection(ray, node->r);
    
  return i1 || i2;
}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
    bool hit = false;
    BBox bb = node->bb;
    double t0;
    double t1;
    if (!(bb.intersect(ray, t0, t1))) {
        return hit;
    }
    if (node->isLeaf()) {
        for (auto p = node->start; p != node->end; p++) {
          total_isects++;
          hit = (*p)->intersect(ray, i) || hit;
        }
        return hit;
    }
    bool hit1;
    bool hit2;
    hit1 = intersect(ray, i, node->l);
    hit2 = intersect(ray, i, node->r);

    return hit1 || hit2;
    //return smallest hit
}

} // namespace SceneObjects
} // namespace CGL

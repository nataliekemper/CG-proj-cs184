#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
    
    double tmin_x = (min.x - r.o.x) / r.d.x;
    double tmin_y = (min.y - r.o.y) / r.d.y;
    double tmin_z = (min.z - r.o.z) / r.d.z;

    double tmax_x = (max.x - r.o.x) / r.d.x;
    double tmax_y = (max.y - r.o.y) / r.d.y;
    double tmax_z = (max.z - r.o.z) / r.d.z;
    
    if (tmin_x > tmax_x) {
        std::swap(tmin_x, tmax_x);
    }
    if (tmin_y > tmax_y) {
        std::swap(tmin_y, tmax_y);
    }
    if (tmin_z > tmax_z) {
        std::swap(tmin_z, tmax_z);
    }
    // max of the mins
    double temp1 = std::max(tmin_x, tmin_y);
    double temp2 = std::max(temp1, tmin_z);
    double tmin = std::max(temp2, r.min_t);
    //min of the maxs
    double temp3 = std::min(tmax_x, tmax_y);
    double temp4 = std::fmin(temp3, tmax_z);
    double tmax = std::min(temp4, r.max_t);

    //double tmin = fmax(t0, min1);
    //double tmax = fmin(t1, max1);

    if (tmax < tmin || r.max_t < tmin || r.min_t > tmax) {
        return false;
    } else {
        t0 = tmin;
        t1 = tmax;
      return true;
    }
}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL

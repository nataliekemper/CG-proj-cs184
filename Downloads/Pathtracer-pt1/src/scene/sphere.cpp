#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.




  return true;

}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
    Vector3D center = this->o;
    int R = this->r;
    Vector3D o = r.o;
    Vector3D d = r.d;
    float a = dot(d, d);
    float b = dot(2 * (o - center), d);
    float c = dot((o - center), (o - center)) - R * R;
    float t1 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    float t2 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);

    return (r.min_t <= t1 && t1 <= r.max_t)
      || (r.min_t <= t2 && t2 <= r.max_t);
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
    Vector3D center = this->o;
    float R = this->r;
    Vector3D o = r.o;
    Vector3D d = r.d;
    float a = dot(d, d);
    float b = dot(2 * (o - center), d);
    float c = dot((o - center), (o - center)) - R * R;
    float t1 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
    float t2 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);

    if (r.min_t <= t1 && t1 <= r.max_t) {
      r.max_t = t1;
      i->t = t1;
      i->n = (o + t1 * d - center).unit();
      i->primitive = this;
      i->bsdf = get_bsdf();
      return true;
    } else if (r.min_t <= t2 && t2 <= r.max_t) {
      r.max_t = t2;
      i->t = t2;
      i->n = (o + t2 * d - center).unit();
      i->primitive = this;
      i->bsdf = get_bsdf();
      return true;
    }
    return false;
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL

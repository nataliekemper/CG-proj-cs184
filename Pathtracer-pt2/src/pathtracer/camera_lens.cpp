#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Project 3-2: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.

// Transforms normalized image coordinates to campera space coordinates
  Vector2D translated = Vector2D(x - 0.5, y - 0.5);
  Vector3D camera_coord = Vector3D(translated.x * tan(0.5 * radians(hFov)) * 2, translated.y * tan(0.5 * radians(vFov)) * 2, -1);

  // Convert a ray from camera space to world space

  Vector3D pLens = Vector3D(lensRadius * sqrt(rndR) * cos(rndTheta), lensRadius * sqrt(rndR) * sin(rndTheta), 0);

  // follow direction of camera_coord -> focus center to get red segment's direction
  // use similar triangles. x / 1 = -pfocus_x / focalDistance and same for y
  Vector3D pFocus = Vector3D(camera_coord.x * focalDistance, camera_coord.y * focalDistance, -focalDistance);

  // make a ray that originates at pLens and goes to pFocus
  Vector3D direction = (pFocus - pLens).unit();

  Ray worldRay = Ray(c2w * pLens + pos, c2w * direction);
  worldRay.min_t = nClip;
  worldRay.max_t = fClip;
  return worldRay;
}


} // namespace CGL

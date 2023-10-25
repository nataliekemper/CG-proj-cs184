#include "pathtracer.h"

#include "pathtracer/intersection.h"
#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"
#include "util/random_util.h"
#include "vector3D.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading

  for (int i = 0; i < num_samples; i++) {
    Vector3D wi = hemisphereSampler->get_sample();
    Ray ir = Ray(hit_p, o2w * wi);
    ir.min_t = EPS_F;
    Intersection iisect;
    if (bvh->intersect(ir, &iisect)) {
      L_out += isect.bsdf->f(w_out, wi) * iisect.bsdf->get_emission() * 2 * PI * wi.z;
    }
  }

  L_out /= num_samples;

  return L_out;

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  int num_samples = ns_area_light;

  for (auto light : scene->lights) {
    for (int i = 0; i < num_samples; i++) {
      Vector3D wi;
      double distToLight, pdf;
      Intersection iisect;

      Vector3D emission = light->sample_L(hit_p, &wi, &distToLight, &pdf);
      Ray between = Ray(hit_p, wi);
      between.min_t = EPS_F;
      between.max_t = distToLight - EPS_F;

      Vector3D wi_local = w2o * wi;
      if (wi_local.z > 0 && !bvh->intersect(between, &iisect)) {
        L_out += isect.bsdf->f(w_out, wi) * emission * wi_local.z / pdf;
      }

      if (light->is_delta_light()) break;
    }
    if (!light->is_delta_light()) L_out /= num_samples;
  }

  return L_out;

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light

  return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  if (direct_hemisphere_sample) {
    return estimate_direct_lighting_hemisphere(r, isect);
  } else {
    return estimate_direct_lighting_importance(r, isect);
  }
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.
  Vector3D L_out = one_bounce_radiance(r, isect);

  double cpdf = 0.3; //continuation probability
  if (r.depth == max_ray_depth || (r.depth + 1 > 2 && !coin_flip(cpdf))) {
    return L_out; // Random termination with probability 1 - cpdf
  }

  Vector3D wi;
  double pdf;
  Vector3D sample = isect.bsdf->sample_f(w_out, &wi, &pdf);

  Ray next = Ray(hit_p, o2w * wi);
  next.min_t = EPS_F;
  next.depth = r.depth + 1;
  Intersection iisect;

  if (bvh->intersect(next, &iisect)) {
    Vector3D rad = at_least_one_bounce_radiance(next, iisect) * sample * wi.z / pdf;
    if (r.depth + 1 > 2) rad /= cpdf;
    L_out += rad;
  }

  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.

  // if (!bvh->intersect(r, &isect))
    // return envLight ? envLight->sample_dir(r) : L_out;

  if (bvh->intersect(r, &isect)) {
    L_out = zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect);
  } else {
    L_out = Vector3D(0, 0, 0);
  }
  // L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"


    int num_samples = ns_aa;          // total samples to evaluate
    int adaptive_samples = 0;
    double illuminance;
    double s1;
    double s2;

    Vector3D rad_sum;

    while (adaptive_samples < num_samples) {
        Vector2D sample = gridSampler->get_sample();
        Ray r = camera->generate_ray((sample.x + x) / sampleBuffer.w, (sample.y + y) / sampleBuffer.h);
        Vector3D rad = est_radiance_global_illumination(r);
        rad_sum += rad;

        // implement adaptive sampling
        illuminance = rad.illum();
        s1 += illuminance;
        s2 += illuminance * illuminance;
        adaptive_samples++;
        if (adaptive_samples % samplesPerBatch == 0) {
            double mean = s1/adaptive_samples;
            double variance = (1.0/(adaptive_samples - 1)) * (s2 - ((s1 * s1) / adaptive_samples));
            double sd = std::sqrt(variance);
            double convergence = 1.96 * (sd / std::sqrt(adaptive_samples));
            if (convergence <= maxTolerance * mean) {
                // printf("Break at %d samples\n", adaptive_samples);
                break;
            }
        }

    }
    // if (adaptive_samples == num_samples) printf("No break\n");
    Vector3D avg = rad_sum / adaptive_samples;
    sampleBuffer.update_pixel(avg, x, y);
    // updating sampleCountBuffer
    sampleCountBuffer[x + y * sampleBuffer.w] = adaptive_samples;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL

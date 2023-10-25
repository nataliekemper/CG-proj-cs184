#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
    // Building Mass Grid
    for (int i = 0; i < num_height_points; i++) {
        for (int j = 0; j < num_width_points; j++) {
            bool pin = false;
            double w = width / num_width_points;
            double h = height / num_height_points;
            for (int p = 0; p < pinned.size(); p++) {
                if (pinned[p][0] == j && pinned[p][1] == i) {
                    pin = true;
                    break;
                }
            }
            
            if (orientation == HORIZONTAL) {
                Vector3D pos = Vector3D(j * w, 1.0, i * h);
                PointMass mass = PointMass(pos, pin);
                point_masses.emplace_back(mass);
            }
            if (orientation == VERTICAL) {
                double z = ((double)rand() / (double)RAND_MAX) * (1.0/1000 - (-1.0/1000)) + (-1.0/1000);
                Vector3D pos = Vector3D(j * w, i * h, z);
                PointMass mass = PointMass(pos, pin);
                point_masses.emplace_back(mass);
            }
        }
    }
    
    // Building Spring Grid
    for (int h = 0; h < num_height_points; h++) {
        for (int w = 0; w < num_width_points; w++) {
            PointMass* mass = &point_masses[h * num_width_points + w];
            // Structural
            if (h > 0) {
                Spring structural = Spring(mass, mass - num_width_points, STRUCTURAL); // pointmass above
                springs.emplace_back(structural);
            }
            if (w > 0) {
                Spring structural = Spring(mass, mass - 1, STRUCTURAL); // pointmass to the left
                springs.emplace_back(structural);
            }
            // Shearing
            if (h > 0 && w < num_width_points - 1) {
                Spring shear = Spring(mass, mass - num_width_points + 1, SHEARING); // diagonal upper right
                springs.emplace_back(shear);
            }
            if (h > 0 && w > 0) {
                Spring shear = Spring(mass, mass - num_width_points - 1, SHEARING); // diagonal upper left
                springs.emplace_back(shear);
            }
            // Bending
            if (w > 1) {
                Spring bend = Spring(mass, mass - 2, BENDING); // two to the left
                springs.emplace_back(bend);
            }
            if (h > 1) {
                Spring bend = Spring(mass, mass - 2 * num_width_points, BENDING); // directly above
                springs.emplace_back(bend);
            }
        }
    }
    
  // TODO (Part 1): Build a grid of masses and springs.

}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
  for (auto& pm : point_masses) {
    pm.forces = 0;
    for (auto acceleration : external_accelerations) {
      pm.forces += acceleration * mass;
    }
  }

  for (auto& spring : springs) {
    if (!cp->enable_structural_constraints && spring.spring_type == STRUCTURAL) continue;
    if (!cp->enable_shearing_constraints && spring.spring_type == SHEARING) continue;
    if (!cp->enable_bending_constraints && spring.spring_type == BENDING) continue;

    Vector3D displacement = spring.pm_a->position - spring.pm_b->position;
    double Fs = cp->ks * (displacement.norm() - spring.rest_length);

    spring.pm_a->forces += Fs * -displacement;
    spring.pm_b->forces += Fs * displacement;
  }

  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  for (auto& pm : point_masses) {
    if (pm.pinned) continue;

    Vector3D new_position = pm.position + (1 - cp->damping/100) * (pm.position - pm.last_position) + (pm.forces / mass) * delta_t * delta_t;
    pm.last_position = pm.position;
    pm.position = new_position;
  }

  // TODO (Part 4): Handle self-collisions.
  build_spatial_map();
  for (auto& pm : point_masses) {
    self_collide(pm, simulation_steps);
  }


  // TODO (Part 3): Handle collisions with other primitives.
    for (auto& pm : point_masses) {
        for (auto *co : *collision_objects) {
            co->collide(pm);
        }
    }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for (auto& spring : springs) {
    Vector3D displacement = spring.pm_a->position - spring.pm_b->position;
    if (displacement.norm() > spring.rest_length * 1.1) {
      if (spring.pm_a->pinned && spring.pm_b->pinned) // Do nothing if both pinned
         continue;
      else if (spring.pm_a->pinned) // Correct spring b
        spring.pm_b->position = spring.pm_a->position - spring.rest_length * 1.1 * displacement.unit();
      else if (spring.pm_b->pinned) // Correct spring a
        spring.pm_a->position = spring.pm_b->position + spring.rest_length * 1.1 * displacement.unit();
      else { // Correct both
        Vector3D center = (spring.pm_a->position + spring.pm_b->position) / 2;
        spring.pm_a->position = center + spring.rest_length * 0.55 * displacement.unit();
        spring.pm_b->position = center - spring.rest_length * 0.55 * displacement.unit();
      }
    }
  }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
  for (auto& pm : point_masses) {
    float hash = hash_position(pm.position);
    auto vec = map.find(hash); // Find entry with the same hash as the point
    if (vec == map.end()) // Create entry if it doesn't exist
      vec = map.emplace(hash, new vector<PointMass *>()).first;
    vec->second->push_back(&pm); // Append to entry
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
  Vector3D correction = Vector3D(0.0);
  int num_corrections = 0;
  
  for (auto& candidate : *map[hash_position(pm.position)]) {
    if (&pm == candidate) continue; // Don't collide the mass with itself
    Vector3D displacement = candidate->position - pm.position;
    if (displacement.norm() < 2 * thickness) {
      correction += displacement - displacement.unit() * 2 * thickness;
      num_corrections++;
    }
  }

  if (num_corrections == 0) return;

  Vector3D total_correction = correction / num_corrections / simulation_steps;
  pm.position += total_correction;
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
  float w = 3 * width / num_width_points;
  float h = 3 * height / num_height_points;
  float t = max(w, h);

  float x = floor(pos.x / w);
  float y = floor(pos.y / h);
  float z = floor(pos.z / t);

  // https://stackoverflow.com/a/2634715
  return ((53 + x) * 53 + y) * 53 + z; 
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}

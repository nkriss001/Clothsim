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
  // TODO (Part 1.1): Build a grid of masses.
  double w_disp = width/num_width_points;
  double h_disp = height/num_height_points;
  for (int h = 0; h < num_height_points; h++) {
    for (int w = 0; w < num_width_points; w++) {
      Vector3D pos;
      if (orientation == HORIZONTAL) {
        pos = Vector3D(w_disp * w, 1, h_disp * h);
      } else {
        double offset = (double)rand() / (500.0*RAND_MAX) - (1.0/1000.0);
        pos = Vector3D(w_disp * w, h_disp * h, offset);
      }
      point_masses.emplace_back(PointMass(pos, false));
    }
  }
  for (int i = 0; i < pinned.size(); i++) {
    vector<int> xy = pinned[i];
    int index = xy[0] + xy[1] * num_width_points;
    point_masses[index].pinned = true;
  }

  // TODO (Part 1.2): Add springs 
  for (int j = 0; j < point_masses.size(); j++) {
    int x = j % num_width_points;
    int y = (j - x) / num_width_points;
    if (x > 0) {
      springs.emplace_back(Spring(&point_masses[j], &point_masses[x - 1 + y * num_width_points], STRUCTURAL));
    }
    if (y < num_height_points - 1) {
      springs.emplace_back(Spring(&point_masses[j], &point_masses[x + (y + 1) * num_width_points], STRUCTURAL));
    }
    if (x > 0 && y < num_height_points - 1) {
      springs.emplace_back(Spring(&point_masses[j], &point_masses[x - 1 + (y + 1) * num_width_points], SHEARING));
    }
    if (x < num_width_points - 1 && y < num_height_points - 1) {
      springs.emplace_back(Spring(&point_masses[j], &point_masses[x + 1 + (y + 1) * num_width_points], SHEARING));
    }
    if (x > 1) {
      springs.emplace_back(Spring(&point_masses[j], &point_masses[x - 2 + y * num_width_points], BENDING));
    }
    if (y < num_height_points - 2) {
      springs.emplace_back(Spring(&point_masses[j], &point_masses[x + (y + 2) * num_width_points], BENDING));
    }
  }

}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2.1): Compute total force acting on each point mass.
  Vector3D externalForce;
  for (Vector3D &accel : external_accelerations) {
    externalForce += accel;
  }
  externalForce *= mass;
  for (PointMass &pointMass : point_masses) {
    pointMass.forces = externalForce;
  }
  for (Spring &spring : springs) {
    if ((spring.spring_type == STRUCTURAL && !cp->enable_structural_constraints)
      || (spring.spring_type == SHEARING && !cp->enable_shearing_constraints)
      || (spring.spring_type == BENDING && !cp->enable_bending_constraints)) {
      continue;
    } else {
      Vector3D connect = spring.pm_b->position - spring.pm_a->position;
      double springMag = cp->ks * (connect.norm() - spring.rest_length);
      if (spring.spring_type == BENDING) {
        springMag *= 0.2;
      }
      Vector3D springForce = connect * springMag / connect.norm();
      spring.pm_a->forces += springForce;
      spring.pm_b->forces -= springForce;
    }
  }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass positions
  for (PointMass &pointMass : point_masses) {
    if (!pointMass.pinned) {
      Vector3D v = (1 - cp->damping / 100) * (pointMass.position - pointMass.last_position);
      pointMass.last_position = pointMass.position;
      pointMass.position = pointMass.position + v + pointMass.forces * delta_t * delta_t / mass;
    }
  }

  // This won't do anything until you complete Part 4.
  build_spatial_map();
  for (PointMass &pm : point_masses) {
    self_collide(pm, simulation_steps);
  }

  // This won't do anything until you complete Part 3.
  for (PointMass &pm : point_masses) {
    for (CollisionObject *co : *collision_objects) {
      co->collide(pm);
    }
  }


  // TODO (Part 2.3): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for (Spring &spring : springs) {
    Vector3D connect = spring.pm_b->position - spring.pm_a->position;
    double length = connect.norm();
    if (length > 1.1 * spring.rest_length) {
      if (!spring.pm_a->pinned && !spring.pm_b->pinned) {
        double fraction = 0.5 * (1 - 1.1 * spring.rest_length / length);
        spring.pm_a->position += fraction * connect;
        spring.pm_b->position -= fraction * connect;
      } else if (spring.pm_a->pinned && !spring.pm_b->pinned) {
        double fraction = 1 - 1.1 * spring.rest_length / length;
        spring.pm_b->position -= fraction * connect;
      } else if (!spring.pm_a->pinned && spring.pm_b->pinned) {
        double fraction = 1 - 1.1 * spring.rest_length / length;
        spring.pm_a->position += fraction * connect;
      }
    }
  }

}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4.2): Build a spatial map out of all of the point masses.
  for (PointMass &pointMass: point_masses) {
    double hash = hash_position(pointMass.position);
    if (map[hash] == NULL) {
      map[hash] = new vector<PointMass *>;
    }
    map[hash]->emplace_back(&pointMass);
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4.3): Handle self-collision for a given point mass.
  Vector3D correction;
  for (PointMass *pointMass : *map[hash_position(pm.position)]) {
    double dist = (pm.position - pointMass->position).norm();
    if (dist != 0 && dist < 2 * thickness) {
      correction += ((pm.position - pointMass->position) * ((2 * thickness / dist) - 1));
    }
  }
  pm.position += (correction / simulation_steps);
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.
  double w = 3 * width / num_width_points;
  double h = 3 * height / num_height_points;
  double t = max(w, h);
  Vector3D box((pos.x - fmod(pos.x, w)) / w, (pos.y - fmod(pos.y, h)) / h, (pos.z - fmod(pos.z, t)) / t);
  return box.x + box.y * num_width_points + box.z * num_height_points * num_width_points;
  return 0.f;
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
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm, pm + num_width_points, pm + 1));
      triangles.push_back(new Triangle(pm + 1, pm + num_width_points,
                                       pm + num_width_points + 1));
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

#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3.1): Handle collisions with spheres.
	if ((pm.position - origin).norm() <= radius) {
		Vector3D ray = pm.position - origin;
		ray.normalize();
		Vector3D tangent = ray * radius + origin;
		Vector3D correction = tangent - pm.last_position;
		pm.position = pm.last_position + (1 - friction) * correction;
	}
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  Misc::draw_sphere(shader, origin, radius * 0.92);
}

#ifndef CONSTRUCTORS_H
#define CONSTRUCTORS_H
#include "phys_types.h"
#include "data.h"
#include <string.h>

Mass mk_Point(Data *data, double *pos, double *vel, double m, double g_mult, bool mobile, int parent_i, double cr, double cg, double cb) {
  Mass point;
  point.pos = malloc(DIMS * sizeof(double));
  point.vel = malloc(DIMS * sizeof(double));
  memcpy(point.pos, pos, DIMS * sizeof(double));
  memcpy(point.vel, vel, DIMS * sizeof(double));
  point.g_mult = g_mult;
  point.m = m;
  point.a = calloc(DIMS, sizeof(double));
  //memset(point.a, 0, DIMS * sizeof(double));
  point.parent_i = parent_i;
  point.mobile = mobile;
  point.r = 0;
  point.mu = 0;
  point.i = data->mass_i++;
  point.type = 0;
  point.cr = cr; point.cg = cg; point.cb = cb;
  return point;
}

Mass mk_s_Point(Data *data, double *pos, double m, bool mobile, int parent_i, double cr, double cg, double cb) {
  return mk_Point(data, pos, (double[]){0,0,0}, m, (double) mobile, mobile, parent_i, cr, cg, cb);
}

Spring mk_Spring(Data *data, Mass *p1, Mass *p2, double l, double k, double f) {
  Spring spring;
  spring.p1 = p1;
  spring.p2 = p2;
  spring.k = k;
  spring.l = l;
  spring.f = f;
  spring.inert = !p1->mobile && !p2->mobile;
  spring.i = data->spring_i++;
  return spring;
}

Spring mk_eq_Spring(Data *data, int pi1, int pi2, double k, double f) {
  Mass *p1 = data->masses[pi1], *p2 = data->masses[pi2];
  double diff[DIMS];
  vec_op(diff, p1->pos, p2->pos, -1, DIMS);
  double dist = vec_norm(diff, DIMS);
  return mk_Spring(data, p1, p2, dist, k, f);
}
Spring mk_eq_Spring_fr_mass(Data *data, int pi1, int pi2, double k, double f, Mass *masses) {
  Mass *p1 = masses + pi1, *p2 = masses + pi2;
  double diff[DIMS];
  vec_op(diff, p1->pos, p2->pos, -1, DIMS);
  double dist = vec_norm(diff, DIMS);
  return mk_Spring(data, p1, p2, dist, k, f);
}

Mass mk_Sphere(Data *data, double *pos, double *vel, double r, double m, double mu) {
  Mass sphere;
  sphere.pos = malloc(DIMS * sizeof(double));
  sphere.vel = malloc(DIMS * sizeof(double));
  memcpy(sphere.pos, pos, DIMS * sizeof(*pos));
  memcpy(sphere.vel, vel, DIMS * sizeof(*vel));
  sphere.r = r;
  sphere.m = m;
  sphere.a = calloc(DIMS, sizeof(double));
  sphere.g_mult = 1;
  sphere.mobile = true;
  sphere.mu = mu;
  sphere.i = data->mass_i++;
  sphere.type = 3;
  return sphere;
}

Object mk_Box(Data *data,
    double x1, double y1, double z1, double x2, double y2, double z2,
    int x_s, int y_s, int z_s,
    double total_m, double k, double f, double r,
    Faces fixed, Faces collide,
    double cr, double cg, double cb) {

  Object box;
  box.i = data->obj_i++;
  box.init_mass_i = data->mass_i;
  box.springs_free = false;
  box.masses_free = false;
  box.n_masses = x_s * y_s * z_s;
  double m = total_m/box.n_masses;
  box.n_springs = (x_s * y_s * 7 - (x_s + y_s - 1) * 5 - 1) * (z_s);
  box.masses = malloc(box.n_masses * sizeof(Mass));
  box.collide = malloc(box.n_masses * sizeof(int));
  box.n_collide = 0;
  box.springs = malloc(box.n_springs * sizeof(Spring));
  box.x_s = x_s; box.y_s = y_s; box.z_s = z_s;
#define BOX_IND(x, y, z) (x) * y_s * z_s + (y) * z_s + (z)
  for (int x = 0; x < x_s; x++) {
    for (int y = 0; y < y_s; y++) {
      for (int z = 0; z < z_s; z++) {
        bool mobile = !(
          (x == 0 && fixed.nx) || (x == x_s - 1 && fixed.px) ||
          (y == 0 && fixed.ny) || (y == y_s - 1 && fixed.py) ||
          (z == 0 && fixed.nz) || (z == z_s - 1 && fixed.pz));
        bool has_collision = (
          (x == 0 && collide.nx) || (x == x_s - 1 && collide.px) ||
          (y == 0 && collide.ny) || (y == y_s - 1 && collide.py) ||
          (z == 0 && collide.nz) || (z == z_s - 1 && collide.pz));
        box.masses[BOX_IND(x, y, z)] = mk_s_Point(data, (double[]){
          x1 + (x2 - x1) * (double) x / (x_s - 1),
          y1 + (y2 - y1) * (double) y / (y_s - 1),
          z1 + (z2 - z1) * (double) z / (z_s - 1)}, m, mobile, box.i,
          cr * (1.0 + mobile)/2, cg * (1.0 + mobile)/2, cb * (1.0 + mobile)/2);
        box.masses[BOX_IND(x, y, z)].r = r;
        if (has_collision) {
          box.collide[box.n_collide++] = BOX_IND(x, y, z);
        }
      }
    }
  }
  box.collide = realloc(box.collide, box.n_collide * sizeof(int));
  int i = 0;
  for (int z = 0; z < z_s; z++) {
    for (int x = 0; x < x_s; x++) {
      for (int y = 0; y < y_s; y++) {
        if (z) {
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x, y, z    ),
            BOX_IND(x, y, z - 1),
            k, f, box.masses);
        }
        if (x) {
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x    , y, z),
            BOX_IND(x - 1, y, z),
            k, f, box.masses);
        }
        if (y) {
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x, y    , z),
            BOX_IND(x, y - 1, z),
            k, f, box.masses);
        }
        if (x && y && z) {
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x    , y    , z    ),
            BOX_IND(x - 1, y - 1, z - 1),
            k, f, box.masses);
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x    , y    , z - 1),
            BOX_IND(x - 1, y - 1, z    ),
            k, f, box.masses);
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x    , y - 1, z    ),
            BOX_IND(x - 1, y    , z - 1),
            k, f, box.masses);
          box.springs[i++] = mk_eq_Spring_fr_mass(data,
            BOX_IND(x - 1, y    , z    ),
            BOX_IND(x    , y - 1, z - 1),
            k, f, box.masses);
        }
      }
    }
  }
  box.n_springs = i;
  return box;
}

/*void fix_mass(Mass *mass) {
    mass->mobile = false;
    mass->g_mult = 0;
    mass->type = 0;
}

int fix_face(Object *box, int face) {
  if (box->shape != Box) {
    return 1;
  }
  int size;
  switch (face) {
    case nx:
    case px:
      size = box->y_s * box->z_s;
      break;
    case ny:
    case py:
      size = box->x_s * box->z_s;
      break;
    case nz:
    case pz:
      size = box->x_s * box->z_s;
      break;
    default:
      return 1;
  }
  for (int i = 0; i < size; i++) {
    fix_mass(box->masses + box->faces[face][i]);
  }
  return 0;
}*/

#endif
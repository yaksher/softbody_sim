#ifndef PHYSICS_TYPES_H
#define PHYSICS_TYPES_H
#include <stdlib.h>

typedef struct Mass {
  double *pos, *vel, g_mult, m, *a, r, mu, cr, cg, cb;
  bool mobile;
  size_t parent_i;
  size_t i;
  // Codes: 0 = fixed point, 1 = fixed sphere, 2 = mobile point, 3 = mobile sphere
  int type;
} Mass;

typedef struct Spring {
  Mass *p1, *p2;
  double k, l, f;
  bool inert;
  size_t i;
} Spring;

void free_mass(Mass mass) {
  free(mass.pos);
  free(mass.vel);
  free(mass.a);
}

#define LIST(type) typedef struct type ## List {\
  size_t i, size;\
  type *a;\
} type ## List;\
\
type ## List mk_ ## type ## List() {\
  type ## List list;\
  list.a = malloc(sizeof(type));\
  list.size = 1;\
  list.i = 0;\
  return list;\
}\
\
void app_ ## type ## List (type ## List list, type elem) {\
  if (list.i >= list.size) {\
    list.a = realloc(list.a, list.size * 2 * sizeof(type));\
  }\
  list.a[list.i++] = elem;\
}\
\
type ## List require_space_ ## type ## List (type ## List list, int space) {\
  if (list.size - list.i <= space) {\
    list.a = realloc(list.a, (list.i + space) * sizeof(type));\
  }\
  return list;\
}

//LIST(Mass)
//LIST(Spring)

typedef enum Shape{Ball, Box} Shape;
typedef struct Faces {
  bool nx, px, ny, py, nz, pz;
} Faces;

typedef struct Object {
  size_t n_masses;
  Mass *masses;
  size_t n_springs;
  Spring *springs;
  void *kd;
  Shape shape;
  size_t init_mass_i;
  size_t i;
  int x_s, y_s, z_s;
  size_t *collide;
  size_t n_collide;
  bool springs_free;
  bool masses_free;
  //int **faces;//*nx, *px, *ny, *py, *nz, *pz;
} Object;

void free_obj(Object obj) {
  if (!obj.masses_free) {
    for (size_t i = 0; i < obj.n_masses; i++) {
      free_mass(obj.masses[i]);
    }
    free(obj.masses);
  }
  if (!obj.springs_free) {
    free(obj.springs);
  }
  free(obj.collide);
  /*if (obj.shape == Box) {
    for (int s = 0; s < 6; s++) {
      free(obj.faces[s]);
    }
    free(obj.faces);
  }*/
}

#endif
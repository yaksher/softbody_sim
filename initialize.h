#ifndef PHYS_INITIALIZE_H
#define PHYS_INITIALIZE_H
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "data.h"
#include "constructors.h"

void flatten_objects(Data *data) {
  data->n_masses = 0;
  data->n_springs = 0;
  data->n_collides = 0;
  for (size_t i = 0; i < data->n_objects; i++) {
    data->n_masses += data->objects[i].n_masses;
    data->n_springs += data->objects[i].n_springs;
    data->n_collides += data->objects[i].n_collide;
  }
  data->masses = malloc(data->n_masses * sizeof(Mass*));
  data->springs = malloc(data->n_springs * sizeof(Spring));
  data->collides = malloc(data->n_collides * sizeof(int));
  size_t cur_mass = 0;
  size_t cur_spring = 0;
  size_t cur_collide = 0;
  for (size_t i = 0; i < data->n_objects; i++) {
    //memcpy(masses + cur_mass, objects[i].masses, objects[i].n_masses * sizeof(Mass));
    for (size_t j = 0; j < data->objects[i].n_masses; j++) {
      data->masses[cur_mass + j] = &data->objects[i].masses[j];
    }
    memcpy(data->springs + cur_spring, data->objects[i].springs, data->objects[i].n_springs * sizeof(Spring));
    for (size_t j = cur_collide; j < cur_collide + data->objects[i].n_collide; j++) {
      data->collides[j] = data->objects[i].collide[j - cur_collide] + cur_mass;
    }
    //free(objects[i].masses);
    free(data->objects[i].springs);
    //objects[i].masses = masses + cur_mass;
    data->objects[i].springs = data->springs + cur_spring;
    //objects[i].masses_free = true;
    data->objects[i].springs_free = true;
    cur_mass += data->objects[i].n_masses;
    cur_spring += data->objects[i].n_springs;
    cur_collide += data->objects[i].n_collide;
  }
}
void flatten_springs(Data *data) {
  data->n_springs = 0;
  for (size_t i = 0; i < data->n_objects; i++) {
    data->n_springs += data->objects[i].n_springs;
  }
  data->springs = malloc(data->n_springs * sizeof(Spring));
  size_t cur_spring = 0;
  for (size_t i = 0; i < data->n_objects; i++) {
    memcpy(data->springs + cur_spring, data->objects[i].springs, data->objects[i].n_springs * sizeof(Spring));
    free(data->objects[i].springs);
    data->objects[i].springs = data->springs + cur_spring;
    data->objects[i].springs_free = true;
    cur_spring += data->objects[i].n_springs;
  }
}

//#define SPRING_BOARD
//#define BLOCK_DROP
//#define POOL
#define BOX_ON_BOX
#if defined(BLOCK_DROP)
#define NO_OBJECTS
int populate_objects(Data *data) {
  // Jelly/trampoline ball drop
  int x_s = 9, y_s = 9, z_s = 5;
#define POINT_IND(x, y, z) (x) * y_s * z_s + (y) * z_s + (z)
  n_points = x_s * y_s * z_s;
  //n_springs = 1924;
  n_springs = (x_s * y_s * 7 - (x_s + y_s - 1) * 5 - 1) * (z_s - 1);
  n_spheres = 1;
  n_masses = n_points + n_spheres;
  points = malloc(n_points * sizeof(Mass));
  springs = malloc(n_springs * sizeof(Spring));
  spheres = malloc(n_spheres * sizeof(Mass));
  double m = -1.00;
  // Jelly ball drop
  for (int x = 0; x < x_s; x++) {
    for (int y = 0; y < y_s; y++) {
      points[POINT_IND(x, y, 0)] = mk_s_Point((double[]){x - x_s/2, y - y_s/2, 0}, m, false);
    }
  }
  for (int z = 1; z < z_s; z++) {
    for (int x = 0; x < x_s; x++) {
      for (int y = 0; y < y_s; y++) {
        points[POINT_IND(x, y, z)] = mk_s_Point((double[]){x - x_s/2, y - y_s/2, z}, m, true);
      }
    }
  }
  // Trampoline ball drop
  /*for (int z = 0; z < z_s; z++) {
    for (int x = 0; x < x_s; x++) {
      for (int y = 0; y < y_s; y++) {
        if ((x != 0 && x != x_s - 1) && (y != 0 && y != y_s - 1)) {
          points[POINT_IND(x, y, z)] = mk_s_Point((double[]){x - x_s/2, y - y_s/2, z}, m, true);
          types[POINT_IND(x, y, z)] = 2;
        } else {
          points[POINT_IND(x, y, z)] = mk_s_Point((double[]){x - x_s/2, y - y_s/2, z}, m, false);
          types[POINT_IND(x, y, z)] = 0;
        }
      }
    }
  }*/
  double k = 1;
  double f = 0.05;
  int i = 0;
  for (int z = 1; z < z_s; z++) {
    for (int x = 0; x < x_s; x++) {
      for (int y = 0; y < y_s; y++) {
        springs[i++] = mk_eq_Spring(
          POINT_IND(x, y, z    ),
          POINT_IND(x, y, z - 1),
          k, f);
        if (x) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y, z),
            POINT_IND(x - 1, y, z),
            k, f);
        }
        if (y) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x, y    , z),
            POINT_IND(x, y - 1, z),
            k, f);
        }
        if (x && y) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y    , z    ),
            POINT_IND(x - 1, y - 1, z - 1),
            k, f);
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y    , z - 1),
            POINT_IND(x - 1, y - 1, z    ),
            k, f);
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y - 1, z    ),
            POINT_IND(x - 1, y    , z - 1),
            k, f);
          springs[i++] = mk_eq_Spring(
            POINT_IND(x - 1, y    , z    ),
            POINT_IND(x    , y - 1, z - 1),
            k, f);
        }
      }
    }
  }
  spheres[0] = mk_Sphere((double[]){0,0,15},(double[]){0,0,-8}, 2, 8, 0);
  printf("populated objects\n");
  return n_masses * DIMS * 2;
}
#elif defined(SPRING_BOARD)
#define NO_OBJECTS
int populate_objects() {
  int x_s = 21, y_s = 9, z_s = 3;
#define POINT_IND(x, y, z) (x) * y_s * z_s + (y) * z_s + (z)
  n_points = x_s * y_s * z_s;
  n_springs = (x_s * y_s * 7 - (x_s + y_s - 1) * 5 - 1) * (z_s);// - x_s * y_s;
  n_spheres = 1;
  n_masses = n_points + n_spheres;
  masses = malloc(n_masses * sizeof(Mass));
  points = masses;
  spheres = masses + n_points;
  springs = malloc(n_springs * sizeof(Spring));
  double m = 0.001;
  double k = 3000;
  double f = 0.05;
  for (int z = 0; z < z_s; z++) {
    for (int x = 0; x < x_s; x++) {
      for (int y = 0; y < y_s; y++) {
        if (x != x_s - 1) {
          points[POINT_IND(x, y, z)] = mk_s_Point((double[]){x - 4, y - y_s/2, z}, m, true);
        } else {
          points[POINT_IND(x, y, z)] = mk_s_Point((double[]){x - 4, y - y_s/2, z}, m, false);
        }
      }
    }
  }
  int i = 0;
  for (int z = 0; z < z_s; z++) {
    for (int x = 0; x < x_s; x++) {
      for (int y = 0; y < y_s; y++) {
        if (z) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x, y, z    ),
            POINT_IND(x, y, z - 1),
            k, f);
        }
        if (x) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y, z),
            POINT_IND(x - 1, y, z),
            k, f);
        }
        if (y) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x, y    , z),
            POINT_IND(x, y - 1, z),
            k, f);
        }
        if (x && y && z) {
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y    , z    ),
            POINT_IND(x - 1, y - 1, z - 1),
            k, f);
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y    , z - 1),
            POINT_IND(x - 1, y - 1, z    ),
            k, f);
          springs[i++] = mk_eq_Spring(
            POINT_IND(x    , y - 1, z    ),
            POINT_IND(x - 1, y    , z - 1),
            k, f);
          springs[i++] = mk_eq_Spring(
            POINT_IND(x - 1, y    , z    ),
            POINT_IND(x    , y - 1, z - 1),
            k, f);
        }
      }
    }
  }
  n_springs = i;
  double sphere_pos[] = {0,0,15}, sphere_vel[] = {0,0,-8},
    sphere_r = 2, sphere_m = 8, sphere_mu = 0;
  spheres[0] = mk_Sphere(sphere_pos,sphere_vel, sphere_r, sphere_m, sphere_mu);
  printf("populated objects\n");
  return n_masses * DIMS * 2;
}
#elif defined(BOX_ON_BOX)
size_t populate_objects(Data *data) {
  printf("Creating objects\n");
  double k = 3000;
  double f = 0.05;
  double r = 0.7;
  data->n_objects = 2;
  data->objects = malloc(data->n_objects * sizeof(Object));
  data->objects[0] = mk_Box(data,
    -4, -4, 0, 16, 4, 4,
    21, 9, 5,
    5, k, f, r,
    (Faces) {false, true, false, false, false, false},
    (Faces) {false, false, false, false, false, true},
    0, 0, 1
  );
  data->objects[1] = mk_Box(data,
    -2, -2, 15, 2, 2, 19,
    5, 5, 5,
    10, k*2, f, r,
    (Faces) {false, false, false, false, false, false},
    (Faces) {false, false, false, false, true, false},
    0.5, 0.5, 0
  );
  data->n_masses = data->objects[0].n_masses + data->objects[1].n_masses;
  data->n_springs = data->objects[0].n_springs + data->objects[1].n_springs;
  data->n_points = data->n_masses;
  printf("populated objects\n");
  //flatten_springs();
  flatten_objects(data);
  return data->n_masses * DIMS * 2;
}
#elif defined(POOL)
#define NO_OBJECTS
int populate_objects() {
  n_points = 11;
  n_spheres = 0;
  n_springs = 0;
  n_masses = 11;
  double m = 0.001;
  masses = malloc(n_masses * sizeof(Mass));
  points = masses;
  spheres = masses + n_points;
  springs = malloc(n_springs * sizeof(Spring));
  double coords[10][3] = {
    {  0,   0, 0},
    {-.2, -.1, 0}, {-.2,  .1, 0},
    {-.4, -.2, 0}, {-.4,   0, 0}, {-.4, .2, 0},
    {-.6, -.3, 0}, {-.6, -.1, 0}, {-.6, .1, 0}, {-.6, .3, 0}
  };
  for (int i = 0; i < 10; i++) {
    points[i] = mk_s_Point(coords[i], m, true);
  }
  points[10] = mk_Point((double[]){1, 0, 0}, (double[]){-1, 0, 0}, m, 1, true, 2);
  return n_masses * DIMS * 2;
}
#endif
#endif
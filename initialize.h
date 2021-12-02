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

// #define BOX_ON_BOX
#define ONE_BOX
#if defined(BOX_ON_BOX)
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
#elif defined(ONE_BOX)
size_t populate_objects(Data *data) {
  printf("Creating objects\n");
  double k = 20;
  double f = 0.05;
  double r = 0.7;
  data->n_objects = 1;
  data->objects = malloc(data->n_objects * sizeof(Object));
  data->objects[0] = mk_Box(data,
    -32, -32, 0, 32, 32, 2,
    65, 65, 3,
    500, k, f, r,
    (Faces) {true, true, true, true, false, false},
    (Faces) {false, false, false, false, false, false},
    0, 0, 1
  );
  data->n_masses = data->objects[0].n_masses;
  data->n_springs = data->objects[0].n_springs;
  data->n_points = data->n_masses;
  printf("populated objects\n");
  //flatten_springs();
  flatten_objects(data);
  return data->n_masses * DIMS * 2;
}
#endif

#endif // PHYS_INITIALIZE_H
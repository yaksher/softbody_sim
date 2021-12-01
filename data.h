#ifndef PHYS_GLOBALS_H
#define PHYS_GLOBALS_H
#include <pthread.h>

#include "phys_types.h"
const size_t DIMS = 3; // Number of dimension the simulation exists in.
const size_t NUM_THREADS = 16;
// double g[DIMS] = {0, 0, 0};
// double g[DIMS] = {0, 0, -10};

// int mass_i = 0;
// int obj_i = 0;
// int spring_i = 0;

// int n_points = 0;
// int n_springs = 0;
// int n_spheres = 0;
// int n_masses = 0;
// int n_objects = 0;
// int n_collides = 0;
// Mass *points = NULL;
// Spring *springs = NULL;
// Mass *spheres = NULL;
// Mass **masses = NULL;
// Object *objects = NULL;
// int *collides = NULL;

typedef struct Data {
    double g[DIMS];

    double *accelerations;

    size_t mass_i;
    size_t obj_i;
    size_t spring_i;

    size_t n_points;
    size_t n_springs;
    size_t n_spheres;
    size_t n_masses;
    size_t n_objects;
    size_t n_collides;
    Mass *points;
    Spring *springs;
    Mass *spheres;
    Mass **masses;
    Object *objects;
    int *collides;
} Data;
#endif
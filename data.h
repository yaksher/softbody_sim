#ifndef PHYS_GLOBALS_H
#define PHYS_GLOBALS_H

#include "phys_types.h"
#include "threading.h"
const size_t DIMS = 3; // Number of dimension the simulation exists in.

typedef struct {
    double g[DIMS];

    double *accelerations;

    thread_handler_t thread_handler;

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
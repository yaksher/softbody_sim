#ifndef VEC_FUNCS_H
#define VEC_FUNCS_H
#include <stdlib.h>

void vec_op(double *out_v, const double *v1, const double *v2, double c, size_t dims);

void vec_mult(double *out_v, const double *v, double c, size_t dims);

double vec_norm(const double *v, size_t dims);

double dot(const double *v1, const double *v2, size_t size);

void mat_vec_prod(const double **m, const double *v, size_t h, size_t w, double *out_v);
#endif
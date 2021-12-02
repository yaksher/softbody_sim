#include "vec_funcs.h"
#include <math.h>

void vec_op(double* out_v, const double* v1, const double* v2, double c, size_t dims) {
    for (size_t i = 0; i < dims; i++) {
        out_v[i] = v1[i] + v2[i] * c;
    }
}

void vec_mult(double* out_v, const double* v, double c, size_t dims) {
    for (size_t i = 0; i < dims; i++) {
        out_v[i] = v[i] * c;
    }
}

double vec_norm(const double* v, size_t dims) {
    double norm_sq = 0;
    for (size_t i = 0; i < dims; i++) {
        norm_sq += v[i] * v[i];
    }
    return sqrt(norm_sq);
}

double dot(const double* v1, const double* v2, size_t size) {
    double dot_prod = 0;
    for (size_t i = 0; i < size; i++) {
        dot_prod += v1[i] * v2[i];
    }
    return dot_prod;
}
void mat_vec_prod(const double** m, const double* v, size_t h, size_t w, double* out_v) {
    for (size_t i = 0; i < h; i++) {
        out_v[i] = 0;
        for (size_t j = 0; j < w; j++) {
            out_v[i] += m[j][i] * v[j];
        }
    }
}
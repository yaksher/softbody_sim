#include <stdbool.h>
#include <stdlib.h>

#include "vec_funcs.h"

int sign(double x);
double max(double a, double b);
double min(double a, double b);
void print_arr(double* arr, int size);
void print_arr2d(double** arr, int s1, int s2);

typedef int err_code_t;

typedef struct {
    double *y0;
    size_t y_size;
    double t0;
    double tf;
} ivp_params_t;

typedef struct {
    size_t save_limit;
    size_t save_count;
    double **saved;
    bool progress;
} solver_params_t;

typedef err_code_t (*deriv_t)(double t, double *y, size_t y_size, double *out, void *data);

err_code_t solve_ivp(deriv_t deriv, ivp_params_t *ivp, solver_params_t *solver, void *data);

#include <stdbool.h>
#include <stdlib.h>

#include "vec_funcs.h"

int sign(double x);
double max(double a, double b);
double min(double a, double b);
void print_arr(double* arr, int size);
void print_arr2d(double** arr, int s1, int s2);

typedef int err_code_t;

typedef err_code_t (*deriv_t)(double t, double *y, size_t y_size, double *out, void *data);

typedef struct {
    deriv_t deriv;
    void *data;
    double *y0;
    size_t y_size;
    double t0;
    double tf;
} ivp_t;

typedef struct {
    size_t save_limit;
    size_t save_count;
    double **saved;
    bool progress;
} solver_params_t;

err_code_t solve_ivp(ivp_t *ivp, solver_params_t *solver);

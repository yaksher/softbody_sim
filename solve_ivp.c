#include "solve_ivp.h"
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "vec_funcs.h"

const double MIN_FACTOR = 0.2;
const double MAX_FACTOR = 10;
const double MAX_STEP = 0.1;
const double MIN_STEP = 0.0009765625;
const double MIN_STEP_FPE_MARGIN = 1 + 1e-3;
const double SAFETY = 0.9;
const double ORDER = 4;
#define ERROR_ORDER 5
const double ERROR_EXP = -1.0 / (ERROR_ORDER + 1);


err_code_t solve_ivp(ivp_t *ivp, solver_params_t *solver) {
// int solve_ivp(deriv_t deriv, double* y0, int y_size, double t0, double tf, int save_limit, double** saved, bool progress, void *data) {
    /* Roughly equivalent to scipy.integrate.solve_ivp using rk45 solver, but
     * with significantly less customization, debug, and robustness.
     * Returns an array of y value arrays concatenated with time_stamp for every
     * step taken. Truncates at save_limit.
     */
	// unpack ivp
    deriv_t deriv = ivp->deriv;
    void *data = ivp->data;
	double t0 = ivp->t0, tf = ivp->tf, *y0 = ivp->y0;
	size_t y_size = ivp->y_size;

	// unpack solver params
	size_t save_limit = solver->save_limit;
	double **saved = solver->saved;
	bool progress = solver->progress;

    double save_interval = fabs(tf - t0) / save_limit;
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    double rtol = 1E-3;
    double atol = 1E-6;
    int direction = sign(tf - t0);

    double f0[y_size];
    err_code_t err_code;
    if ((err_code = deriv(t0, y0, y_size, f0, data))) {
        return err_code;
    }
    double scale[y_size];
    double d0 = 0;
    double d1 = 0;
    for (size_t i = 0; i < y_size; i++) {
        scale[i] = atol + fabs(y0[i]) * rtol;
        double y0_scale = y0[i] / scale[i];
        d0 += y0_scale * y0_scale;
        double f0_scale = f0[i] / scale[i];
        d1 += f0_scale * f0_scale;
    }
    d0 = sqrt(d0 / y_size);
    d1 = sqrt(d1 / y_size);

    double h0;
    if (d0 < 1E-5 || d1 < 1E-5) {
        h0 = 1E-6;
    } else {
        h0 = 0.01 * d0 / d1;
    }

    double y1[y_size];
    for (size_t i = 0; i < y_size; i++) {
        y1[i] = y0[i] + h0 * direction * f0[i];
    }
    double f1[y_size];
    if ((err_code = deriv(t0 + h0 * direction, y1, y_size, f1, data))) {
        return err_code;
    }
    double d2 = 0;
    for (size_t i = 0; i < y_size; i++) {
        double fdif_scale = (f1[i] - f0[i]) / scale[i];
        d2 += fdif_scale * fdif_scale;
    }
    d2 = sqrt(d2 / y_size) / h0;
    double h1;
    if (d1 <= pow(10, -15) && d2 <= pow(10, -15)) {
        h1 = max(1E-6, h0 * 1E-3);
    } else {
        h1 = pow(0.01 / max(d1, d2), 1 / (ERROR_ORDER + 1));
    }
    double h_abs = min(100 * h0, h1);

    double t = t0;
    double* y = y0;
    double* f = f0;
    double t_new;
    double y_new[y_size];
    double f_new[y_size];
    size_t step_count = 1;
    size_t real_steps = 1;
    for (size_t i = 0; i < y_size; i++) {
        saved[0][i] = y0[i];
    }
    saved[0][y_size] = t0;

    const size_t n_stages = 6;
    double* K[n_stages + 1];
    K[0] = f;
    for (size_t i = 1; i < n_stages; i++) {
        K[i] = malloc(y_size * sizeof(double));
    }
    K[n_stages] = f_new;
    const double A[6][5] = {
                {0, 0, 0, 0, 0},
                {1.0/5, 0, 0, 0, 0},
                {3.0/40, 9.0/40, 0, 0, 0},
                {44.0/45, -56.0/15, 32.0/9, 0, 0},
                {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729, 0},
                {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656}};
    const double B[] = {35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84};
    const double C[] = {0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1};
    const double E[] = {-71.0/57600, 0, 71.0/16695, -71.0/1920, 17253.0/339200, -22.0/525, 1.0/40};
    const double min_step = 1E-25;

    while (direction * (t - tf) < 0) {
        h_abs = max(min_step, h_abs);
        if (progress) {
            clock_gettime(CLOCK_MONOTONIC, &finish);
            double elapsed = (finish.tv_sec - start.tv_sec);
            elapsed += (finish.tv_nsec - start.tv_nsec) / 1E9;
            printf("\r%.2f / %.2f | run : %.2fs | left : %.2fs | step : %.2e | %lu", t, tf, elapsed, (tf/t - 1) * elapsed, h_abs, real_steps);
            fflush(stdout);
        }
        bool step_accepted = false;
        bool step_rejected = false;
        while (!step_accepted) {
            if (h_abs < min_step) {
                fprintf(stderr, "Step %lu too small\n", step_count);
                for (size_t i = 1; i < n_stages; i++) {
                    free(K[i]);
                }
                return -1;
            }
            double h = max(MIN_STEP, min(h_abs, MAX_STEP)) * direction;
            t_new = t + h;

            if (direction * (t_new - tf) > 0) {
                t_new = tf;
            }
            h = t_new - t;
            h_abs = fabs(h);
            for (size_t s = 1; s < n_stages; s++) {
                const double *a = A[s];
                double c = C[s];
                double dy[y_size];
                mat_vec_prod((const double **) K, a, y_size, s, dy);
                double ydy[y_size];
                for (size_t i = 0; i < y_size; i++) {
                    ydy[i] = y[i] + dy[i] * h;
                }
                if ((err_code = deriv(t + c * h, ydy, y_size, K[s], data))) {
                    for (size_t i = 1; i < n_stages; i++) {
                        free(K[i]);
                    }
                    return err_code;
                }
            }
            double y_delta[y_size];
            mat_vec_prod((const double **) K, B, y_size, n_stages, y_delta);
            for (size_t i = 0; i < y_size; i++) {
                y_new[i] = y[i] + h * y_delta[i];
            }
            if ((err_code = deriv(t + h, y_new, y_size, f_new, data))) {
                for (size_t i = 1; i < n_stages; i++) {
                    free(K[i]);
                }
                return err_code;
            }
            double K_err[y_size];
            mat_vec_prod((const double **) K, E, y_size, n_stages + 1, K_err);
            double error_norm = 0;
            double scale[y_size];
            for (size_t i = 0; i < y_size; i++) {
                K_err[i] *= h;
                scale[i] = atol + max(fabs(y[i]), fabs(y_new[i])) * rtol;
                double K_err_scale = K_err[i] / scale[i];
                error_norm += K_err_scale * K_err_scale;
            }
            error_norm = sqrt(error_norm / y_size);
            if (isnan(error_norm)) {
                fprintf(stderr, "\nERROR. NAN ENCOUNTERED\n");
                for (size_t i = 1; i < n_stages; i++) {
                    free(K[i]);
                }
                return 3;
            }
            if (error_norm < 1 || h_abs <= MIN_STEP * MIN_STEP_FPE_MARGIN) {
                double factor;
                if (error_norm == 0) {
                    factor = MAX_FACTOR;
                } else {
                    factor = min(MAX_FACTOR, SAFETY * pow(error_norm, ERROR_EXP));
                }

                if (step_rejected) {
                    factor = min(1, factor);
                }
                h_abs *= factor;
                step_accepted = true;
            } else {
                h_abs *= max(MIN_FACTOR, SAFETY * pow(error_norm, ERROR_EXP));
                step_rejected = true;
            }
        }
        if (step_count < save_limit && direction * (t_new - saved[step_count-1][y_size]) > save_interval) {
            for (size_t i = 0; i < y_size; i++) {
                saved[step_count][i] = y_new[i];
            }
            saved[step_count][y_size] = t_new;
            step_count++;
        }
        real_steps++;


        t = t_new;
        memcpy(y, y_new, y_size * sizeof(double));
        memcpy(f, f_new, y_size * sizeof(double));
        //y = y_new;
        //f = f_new;
    }
    for (size_t i = 1; i < n_stages; i++) {
        free(K[i]);
    }
    // *saved[save_limit] = (double) step_count;
	solver->save_count = step_count;
    if (!progress) {
        clock_gettime(CLOCK_MONOTONIC, &finish);
        double elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1E9;
        fprintf(stderr, "Solve_IVP completed in %f seconds", elapsed);
    }
    return 0;
}

int sign(double x) {
    return x > 0.0 - x < 0.0;
}
double max(double a, double b) {
    return a > b ? a : b;
}
double min(double a, double b) {
    return a < b ? a : b;
}
void print_arr(double* arr, int size) {
    printf("size : %d | ", size);
    for (int i = 0; i < size; i++) {
        printf("%.2f,", arr[i]);
    }
    printf("\n");
}
void print_arr2d(double** arr, int s1, int s2) {
    for (int i = 0; i < s1; i++) {
        for (int j = 0; j < s2; j++) {
            printf("%.2f, ", arr[i][j]);
        }
        printf("\n");
    }
}

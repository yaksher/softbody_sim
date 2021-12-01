#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#include "vec_funcs.h"
#include "phys_types.h"
#include "constructors.h"
#include "solve_ivp.h"
#include "sim.h"
#include "kdtree.h"
#include "data.h"
#include "initialize.h"

const double RENDER_FRAMERATE = 50.0; // At most 10 * timeinterval * framerate frames will be saved by the code for interpolation.
const double COLLISION_THRESHOLD = 1.0/5.0;

int apply_force(Spring spring, Data *data) {
    if (spring.inert) {
        return 0;
    }
    Mass *p1 = spring.p1, *p2 = spring.p2;
    if (p1->parent_i != p2->parent_i) {
        fprintf(stderr, "\nspring %lu has points %lu and %lu in different objects %lu and %lu\n", spring.i, p1->i, p2->i, p1->parent_i, p2->parent_i);
        return 1;
    }
    double diff[DIMS];
    //#define PROPER_FRICTION
    #ifdef PROPER_FRICTION
    double r_dot_v = 0, r_dot_self = 0;
    for (size_t d = 0; d < DIMS; d++) {
        double v_diff_d = p2->vel[d] - p1->vel[d];
        diff[d] = p2->pos[d] - p1->pos[d];
        r_dot_v += v_diff_d * diff[d];
        r_dot_self += diff[d] * diff[d];
    }
    double dist = sqrt(r_dot_self);
    double deflection = dist - spring.l;

    for (size_t d = 0; d < DIMS; d++) {
        double unit_d = diff[d] / dist;
        double v_r_proj_d = diff[d] * r_dot_v / r_dot_self;
        double v_diff_d = p2->vel[d] - p1->vel[d];
        double F_d = - unit_d * deflection * spring.k - v_r_proj_d * spring.f;
        p1->a[d] -= F_d * p1->mobile / p1->m;
        p2->a[d] += F_d * p2->mobile / p2->m;
    }
    #else
    double r_dot_self = 0;
    for (size_t d = 0; d < DIMS; d++) {
        diff[d] = p2->pos[d] - p1->pos[d];
        r_dot_self += diff[d] * diff[d];
    }
    double dist = sqrt(r_dot_self);
    double deflection = dist - spring.l;

    for (size_t d = 0; d < DIMS; d++) {
        double unit_d = diff[d] / dist;
        double v_diff_d = p2->vel[d] - p1->vel[d];
        double F_d = - unit_d * deflection * spring.k - v_diff_d * spring.f;
        data->accelerations[p1->i * DIMS + d] -= F_d * p1->mobile / p1->m;
        data->accelerations[p2->i * DIMS + d] += F_d * p2->mobile / p2->m;
    }
    #endif
    if (isnan(p1->a[0])) {
        fprintf(stderr, "\nspring: %lu %d %.2f %.2f %.2f\n", p1->i, p1->type, p1->a[0], p1->a[1], p1->a[2]);
        return 3;
    }
    return 0;
}

void collide(Mass *self, Mass *others, size_t n_others) {
    #define COLLISION_FRICTION
    for (size_t i = 0; i < n_others; i++) {
        if (self == others + i) {
            continue;
        }
        double r_vec[DIMS];
        double v_diff[DIMS];
        double r_dot_self = 0;
        double r_dot_v = 0;
        for (size_t d = 0; d < DIMS; d++) {
            r_vec[d] = others[i].pos[d] - self->pos[d];
            r_dot_self += r_vec[d] * r_vec[d];
            #ifdef COLLISION_FRICTION
            v_diff[d] = others[i].vel[d] - self->vel[d];
            r_dot_v += r_vec[d] * v_diff[d];
            #endif
        }
        double dist = sqrt(r_dot_self);

        #ifdef COLLISION_FRICTION
        double v_orth[DIMS];
        double v_orth_dot_self = 0;
        for (size_t d = 0; d < DIMS; d++) {
            double v_r_proj_d = r_vec[d] * r_dot_v / r_dot_self;
            v_orth[d] = v_diff[d] - v_r_proj_d;
            v_orth_dot_self += v_orth[d] * v_orth[d];
        }
        double v_orth_size = sqrt(v_orth_dot_self);
        if (v_orth_size == 0) {
            memset(v_orth, 0, DIMS * sizeof(double));
            v_orth_size = 1;
        }
        #endif

        double r_diff = max(dist - self->r - others[i].r, 0.001);
        if (r_diff < COLLISION_THRESHOLD) {
            double F = 1 / (5 * r_diff) - 1;
            double mu = self->mu + others[i].mu;
            for (size_t d = 0; d < DIMS; d++) {
                #ifdef COLLISION_FRICTION
                self->a[d] -= F * (r_vec[d] / dist - mu * v_orth[d] / v_orth_size) / self->m;
                others[i].a[d] += others[i].mobile * F * (r_vec[d] / dist - mu * v_orth[d] / v_orth_size) / others[i].m;
                #else
                self->a[d] -= F * (r_vec[d] / dist) / self->m;
                others[i].a[d] += others[i].mobile * F * (r_vec[d] / dist) / others[i].m;
                #endif
            }
        }
    }
    if (isnan(self->a[0])) {
        fprintf(stderr, "\n%lu %d %.2f %.2f %.2f", self->i, self->type, self->a[0], self->a[1], self->a[2]);
    }
}

void collide_set(Mass *self, void *others, double collision_threshold, Data *data) {
    double other_pos[DIMS];
    while (!kd_res_end(others)) {
        Mass *other = (Mass *)kd_res_item(others, other_pos);
        kd_res_next(others);
        if (self == other) {
            continue;
        }
        double r_vec[DIMS];
        double v_diff[DIMS];
        double r_dot_self = 0;
        double r_dot_v = 0;
        for (size_t d = 0; d < DIMS; d++) {
            r_vec[d] = other_pos[d] - self->pos[d];
            r_dot_self += r_vec[d] * r_vec[d];
            v_diff[d] = other->vel[d] - self->vel[d];
            r_dot_v += r_vec[d] * v_diff[d];
        }
        double dist = sqrt(r_dot_self);

        double v_orth[DIMS];
        double v_orth_dot_self = 0;
        for (size_t d = 0; d < DIMS; d++) {
            double v_r_proj_d = r_vec[d] * r_dot_v / r_dot_self;
            if (isnan(v_r_proj_d)) {
                v_r_proj_d = 0;
            }
            v_orth[d] = v_diff[d] - v_r_proj_d;
            v_orth_dot_self += v_orth[d] * v_orth[d];
        }
        double v_orth_size = sqrt(v_orth_dot_self);
        if (isnan(v_orth_size)) {
            v_orth_size = 0;
        }
        if (v_orth_size == 0) {
            for (size_t d = 0; d < DIMS; d++) {
                v_orth[d] = 0;
            }
            v_orth_size = 1;
        }
        double r_diff = max(dist - self->r - other->r, 0.001);
        if (r_diff < collision_threshold) {
            double F = collision_threshold / r_diff - 1;
            double mu = self->mu + other->mu;
            for (size_t d = 0; d < DIMS; d++) {
                data->accelerations[self->i * DIMS + d] -= self->mobile * F * (r_vec[d] / dist - mu * v_orth[d] / v_orth_size) / self->m;
                data->accelerations[other->i * DIMS + d] += other->mobile * F * (r_vec[d] / dist - mu * v_orth[d] / v_orth_size) / other->m;
            }
        }
    }
    kd_res_free(others);
}

void initial_state(double* out, int y_size, Data *data) {
    #if false//!defined(NO_OBJECTS)
    size_t writer = 0;
    for (size_t o = 0; o < n_objects; o++) {
        for (size_t i = 0; i < objects[o].n_masses; i++) {
            for (size_t d = 0; d < DIMS; d++) {
                out[writer] = objects[o].masses[i].pos[d];
                out[n_masses * DIMS + writer++] = objects[o].masses[i].vel[d];
            }
        }
    }
    #else
    (void) y_size;
    for (size_t i = 0; i < data->n_masses; i++) {
        for (size_t d = 0; d < DIMS; d++) {
            out[i * DIMS + d] = data->masses[i]->pos[d];
            out[(data->n_masses + i) * DIMS + d] = data->masses[i]->vel[d];
        }
    }
    #endif
    fprintf(stderr, "compiled initial state\n");
}

err_code_t deriv(double t, double *y, size_t y_size, double *out, void *_data) {
    Data *data = _data;
    (void) t; // system is time invariant

    int vel_offset = y_size / 2;
    double *vels = y + vel_offset;
    data->accelerations = out + vel_offset;
    for (size_t i = 0; i < data->n_masses; i++) {
        for (size_t d = 0; d < DIMS; d++) {
            data->masses[i]->pos[d] = y[i * DIMS + d];
            data->masses[i]->vel[d] = vels[i * DIMS + d];
            //masses[i].a[d] = g[d] * masses[i].g_mult;
            data->accelerations[i * DIMS + d] = data->g[d] * data->masses[i]->g_mult;
        }
    }
    for (size_t o = 0; o < data->n_objects; o++) {
        Object *obj = &data->objects[o];
        if (data->n_objects > 1) {
            obj->kd = kd_create(DIMS);
            for (size_t i = 0; i < obj->n_collide; i++) {
                kd_insert(obj->kd, obj->masses[obj->collide[i]].pos, obj->masses + obj->collide[i]);
            }
        }
    }
    int spring_err;
    for (size_t i = 0; i < data->n_springs; i++) {
        if ((spring_err = apply_force(data->springs[i], data))) {
            return spring_err;
        }
    }
    for (size_t o1 = 0; o1 < data->n_objects; o1++) {
        Object *obj1 = &data->objects[o1];
        /*for (size_t i = 0; i < obj1->n_springs; i++) {
            if ((spring_err = apply_force(obj1->springs[i]))) {
                return spring_err;
            }
        }*/
        for (size_t o2 = o1 + 1; o2 < data->n_objects; o2++) {
            for (size_t i = 0; i < obj1->n_collide; i++) {
                Object *obj2 = &data->objects[o2];
                collide_set(obj1->masses + obj1->collide[i],
                    kd_nearest_range(obj2->kd,
                        obj1->masses[obj1->collide[i]].pos,
                        COLLISION_THRESHOLD + obj1->masses[obj1->collide[i]].r),
                    COLLISION_THRESHOLD, data);
            }
        }
    }
    if (data->n_objects > 1) {
        for (size_t o = 0; o < data->n_objects; o++) {
            kd_free(data->objects[o].kd);
        }
    }
    memcpy(out, vels, vel_offset * sizeof(double));
    return 0;
}

void save_data(double** solved, size_t frames, size_t y_size, double delta_t, Data *data) {
    FILE *out = fopen("sim_data.txt", "w+");

    (void) delta_t;
    fprintf(out, "%lu\n%lu\n", DIMS, data->n_masses);
    for (size_t o = 0; o < data->n_objects; o++) {
        Object *obj = &data->objects[o];
        for (size_t i = 0; i < obj->n_masses; i++) {
            if (o != data->n_objects - 1 || i != obj->n_masses - 1) {
                fprintf(out, "%f,%f,%f,", obj->masses[i].cr, obj->masses[i].cg, obj->masses[i].cb);
            }
        }
    }
    Mass *last = &data->objects[data->n_objects - 1].masses[data->objects[data->n_objects - 1].n_masses - 1];
    fprintf(out, "%f,%f,%f\n", last->cr, last->cg, last->cb);
    fprintf(out, "\n");
    for (size_t i = 0; i < frames; i++) {
        for (size_t j = 0; j < y_size / 2; j ++) {
            fprintf(out, "%f,", solved[i][j]);
        }
        fprintf(out, "%f\n", solved[i][y_size]);
    }

    fclose(out);
}

int main() {
    Data data;

    memset(&data, 0, sizeof(data));
    data.g[2] = -10.0;

    const double t0 = 0, tf = 5;
    const double delta_t = tf - t0;
    size_t y_size = populate_objects(&data);
    const size_t save_limit = RENDER_FRAMERATE * delta_t * 10;

    double* state_0 = malloc(y_size * sizeof(double));
    initial_state(state_0, y_size, &data);

    // fprintf(stderr, "n_objects: %lu | n_masses: %lu | n_points: %lu | n_spheres: %lu | n_springs: %lu\n",
    //                                    data.n_objects,data.n_masses,data.n_points,data.n_spheres,data.n_springs);

    double** solved = malloc(save_limit * sizeof(double*));
    for (size_t i = 0; i < save_limit; i++) {
        solved[i] = malloc((y_size + 1) * sizeof(double));
    }

    fprintf(stderr, "beginning integration\n");

    ivp_params_t ivp = (ivp_params_t) {
        .y0 = state_0,
        .y_size = y_size,
        .t0 = t0,
        .tf = tf
    };
    solver_params_t solver = (solver_params_t) {
        .save_limit = save_limit,
        .save_count = 0,
        .saved = solved,
        .progress = false
    };

    // run the solver
    if (solve_ivp(deriv, &ivp, &solver, &data)) {
        fprintf(stderr, "solve_ivp failed\n");
        return 1;
    }
    fprintf(stderr, "\ncomputation complete\n");

    // save the data
    size_t frames = solver.save_count;
    fprintf(stderr, "frames: %lu\n", frames);

    save_data(solved, frames, y_size, delta_t, &data);
    fprintf(stderr, "saved\n");

    // free resources
    for (size_t i = 0; i < save_limit; i++) {
        free(solved[i]);
    }
    free(solved);

    if (data.objects) {
        for (size_t o = 0; o < data.n_objects; o++) {
            free_obj(data.objects[o]);
        }
        free(data.objects);
    }
    if (data.masses) {
        free(data.masses);
    }
    if (data.springs) {
        free(data.springs);
    }
    return 0;
}
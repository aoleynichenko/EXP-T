/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

/*
 * Direct Inversion in the Iterative Subspace for Coupled Cluster.
 * See also:
 * - G. E. Scuseria, T. J. Lee, H. F. Schaefer III.
 *   Chem. Phys. Lett. 130. 236-239 (1986)
 * But the actually implemented algorithm is taken from:
 * https://github.com/psi4/psi4numpy/tree/master/Coupled-Cluster/Spin_Orbitals/CCSD
 */

#include "diis.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "error.h"
#include "engine.h"
#include "linalg.h"
#include "memory.h"
#include "options.h"
#include "timer.h"


/**
 * allocates memory, initializes the empty diis_queue_t structure and
 * returns pointer to it.
 *
 * @param do_t1 extrapolate singles amplitudes (T1)
 * @param do_t2 extrapolate doubles amplitudes (T2)
 * @param do_t3 extrapolate doubles amplitudes (T3)
 */
diis_queue_t *new_diis_queue(int do_t1, int do_t2, int do_t3)
{
    diis_queue_t *q = (diis_queue_t *) cc_malloc(sizeof(diis_queue_t) * 1);
    q->err_matrix = (double *) cc_malloc(sizeof(double) * DIIS_MAX * DIIS_MAX);
    q->n = 0;
    q->do_t1 = do_t1;
    q->do_t2 = do_t2;
    q->do_t3 = do_t3;
    return q;
}


/**
 * destructor for the diis_queue_t object
 */
void delete_diis_queue(diis_queue_t *q)
{
    timer_new_entry("diis", "DIIS extrapolation");
    timer_start("diis");

    for (int i = 0; i < q->n; i++) {
        diagram_stack_erase(q->e1[i]);
        diagram_stack_erase(q->e2[i]);
        diagram_stack_erase(q->e3[i]);
        diagram_stack_erase(q->t1[i]);
        diagram_stack_erase(q->t2[i]);
        diagram_stack_erase(q->t3[i]);
    }
    cc_free(q->err_matrix);
    cc_free(q);

    timer_stop("diis");
}


/**
 * adds amplitudes to the DIIS queue. calculates corresponding error vectors
 *
 * @param diag_t1new T1 amplitudes obtained at the current CC iteration
 * @param diag_t1old old T1 amplitudes, in fact, previous DIIS extrapolant
 * @param diag_t2new see above, for T2
 * @param diag_t2old see above, for T2
 * @param diag_t3new see above, for T3
 * @param diag_t3old see above, for T3
 * @param iter number of current iteration, for unique names
 *
 * @todo 'iter' is unnecessary, we can count iteration inside the queue object
 */
void diis_put(diis_queue_t *q,
              char *diag_t1new, char *diag_t1old,
              char *diag_t2new, char *diag_t2old,
              char *diag_t3new, char *diag_t3old,
              int iter)
{
    char buf[CC_DIAGRAM_MAX_NAME];
    char buf_err[CC_DIAGRAM_MAX_NAME];

    timer_new_entry("diis", "DIIS extrapolation");
    timer_start("diis");

    if (q->n >= DIIS_MAX) {
        errquit("Too many vectors (%d) in DIIS! Please, increase DIIS_MAX (src/methods/diis.h)", q->n);
    }

    if (q->do_t1) {
        // store new T1 amplitudes
        sprintf(buf, "t1_%s_%d", diag_t1new, iter);
        copy(diag_t1new, buf);
        strcpy(q->t1[q->n], buf);

        // store T1 error vector
        sprintf(buf_err, "e1_%s_%d", diag_t1new, iter);
        copy(diag_t1new, buf_err);
        update(buf_err, -1.0, diag_t1old);
        strcpy(q->e1[q->n], buf_err);
    }

    if (q->do_t2) {
        // store new T2 amplitudes
        sprintf(buf, "t2_%s_%d", diag_t2new, iter);
        copy(diag_t2new, buf);
        strcpy(q->t2[q->n], buf);

        // store T2 error vector
        sprintf(buf_err, "e2_%s_%d", diag_t2new, iter);
        copy(diag_t2new, buf_err);
        update(buf_err, -1.0, diag_t2old);
        strcpy(q->e2[q->n], buf_err);
    }

    if (q->do_t3) {
        // store new T2 amplitudes
        sprintf(buf, "t3_%s_%d", diag_t3new, iter);
        copy(diag_t3new, buf);
        strcpy(q->t3[q->n], buf);

        // store T2 error vector
        sprintf(buf_err, "e3_%s_%d", diag_t3new, iter);
        copy(diag_t3new, buf_err);
        update(buf_err, -1.0, diag_t3old);
        strcpy(q->e3[q->n], buf_err);
    }

    // calculate scalar products with itself and all previous error vectors
    int n = q->n;
    for (int i = 0; i <= n; i++) {

        double s1 = 0.0, s2 = 0.0, s3 = 0.0;

        if (q->do_t1) {
            double complex z1 = scalar_product("C", "N", q->e1[i], q->e1[n]);
            s1 = creal(z1);
        }
        if (q->do_t2) {
            double complex z2 = scalar_product("C", "N", q->e2[i], q->e2[n]);
            s2 = creal(z2);
        }
        if (q->do_t3) {
            double complex z3 = scalar_product("C", "N", q->e3[i], q->e3[n]);
            s3 = creal(z3);
        }

        q->err_matrix[i * DIIS_MAX + n] = s1 + s2 + s3;
        q->err_matrix[n * DIIS_MAX + i] = s1 + s2 + s3;
    }

    q->n++;

    timer_stop("diis");
}


/**
 * truncates the DIIS queue object (new length = len).
 * amplutides and error vectors are removed from the bottom of the stack (queue).
 */
void diis_truncate(diis_queue_t *q, int len)
{
    timer_new_entry("diis", "DIIS extrapolation");
    timer_start("diis");

    if (q->n <= len) { return; }

    int offs = q->n - len;

    for (int i = 0; i < offs; i++) {
        if (q->do_t1) {
            diagram_stack_erase(q->e1[i]);
            diagram_stack_erase(q->t1[i]);
        }
        if (q->do_t2) {
            diagram_stack_erase(q->e2[i]);
            diagram_stack_erase(q->t2[i]);
        }
        if (q->do_t3) {
            diagram_stack_erase(q->e3[i]);
            diagram_stack_erase(q->t3[i]);
        }
    }

    for (int i = 0; i < len; i++) {
        if (q->do_t1) {
            strcpy(q->t1[i], q->t1[i + offs]);
            strcpy(q->e1[i], q->e1[i + offs]);
        }
        if (q->do_t2) {
            strcpy(q->t2[i], q->t2[i + offs]);
            strcpy(q->e2[i], q->e2[i + offs]);
        }
        if (q->do_t3) {
            strcpy(q->t3[i], q->t3[i + offs]);
            strcpy(q->e3[i], q->e3[i + offs]);
        }
    }

    // shift the submatrix of the whole error matrix to the upper-left corner
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            q->err_matrix[i * DIIS_MAX + j] = q->err_matrix[(i + offs) * DIIS_MAX + (j + offs)];
        }
    }

    q->n = len;

    timer_stop("diis");
}


/**
 * Performs DIIS extrapolation.
 *
 * @note for complex-valued amplitudes only real parts of <e_i|e_j> dot products
 * are used, and it works, I don't know, why
 */
void diis_extrapolate(diis_queue_t *q, char *extrap_t1, char *extrap_t2, char *extrap_t3)
{
    timer_new_entry("diis", "DIIS extrapolation");
    timer_start("diis");

    int dim = q->n;
    int bdim = dim + 1;

    // Pulay matrix
    double *B = (double *) cc_malloc(sizeof(double) * bdim * bdim);
    double *right = (double *) cc_malloc(sizeof(double) * bdim);
    double *diis_coeffs = (double *) cc_malloc(bdim * sizeof(double));

    // Pulay matrix construction
    // common for both Singles and Doubles amplitudes and error vectors
    for (int i = 0; i < dim; i++) {
        int j;
        for (j = i; j < dim; j++) {
            /*double s1 = 0.0, s2 = 0.0, s3 = 0.0;
            if (q->do_t1) {
                double complex z1 = scalar_product("C", "N", q->e1[i], q->e1[j]);
                s1 = creal(z1);
            }
            if (q->do_t2) {
                double complex z2 = scalar_product("C", "N", q->e2[i], q->e2[j]);
                s2 = creal(z2);
            }
            if (q->do_t3) {
                double complex z3 = scalar_product("C", "N", q->e3[i], q->e3[j]);
                s3 = creal(z3);
            }
            B[i * bdim + j] = s1 + s2 + s3;
            B[j * bdim + i] = s1 + s2 + s3;*/

            B[i * bdim + j] = q->err_matrix[i * DIIS_MAX + j];
            B[j * bdim + i] = B[i * bdim + j];
        }
        B[i * bdim + j] = -1.0;
    }

    // -1 -1 -1 ... 0
    for (int i = 0; i < bdim - 1; i++) {
        B[(bdim - 1) * bdim + i] = -1.0;
    }
    B[bdim * bdim - 1] = 0.0;

    // normalize Pulay matrix by abs max
    double absmax = 0.0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (fabs(B[i * bdim + j]) > absmax) {
                absmax = fabs(B[i * bdim + j]);
            }
        }
    }
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            B[i * bdim + j] /= absmax;
        }
    }

    // r.h.s. of linear equations
    for (int i = 0; i < bdim; i++) {
        right[i] = 0.0;
    }
    right[bdim - 1] = -1;

    // solve system of linear equations for weights
    int info = linsys(bdim, B, right, diis_coeffs);
    // Check for the exact singularity
    if (info != 0) {
        printf(" DIIS:\n");
        printf(" The diagonal element of the triangular factor of A, "
               "U(%i,%i) is zero, so that A is singular;\n", info, info);
        printf(" the solution could not be computed."
               " DIIS will be turned off\n");
        cc_opts->diis_enabled = 0;
        goto finalize;
    }

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        printf(" DIIS coeff-s:\n");
        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
            printf("%12.3e\n", diis_coeffs[i]);
            sum += diis_coeffs[i];
        }
        printf(" sum of DIIS coeff-s = %12.6e\n", sum);
    }

    // construct extrapolations for T1 and T2
    if (q->do_t1) {
        clear(extrap_t1);
        for (int i = 0; i < dim; i++) {
            update(extrap_t1, diis_coeffs[i], q->t1[i]);
        }
    }
    if (q->do_t2) {
        clear(extrap_t2);
        for (int i = 0; i < dim; i++) {
            update(extrap_t2, diis_coeffs[i], q->t2[i]);
        }
    }
    if (q->do_t3) {
        clear(extrap_t3);
        for (int i = 0; i < dim; i++) {
            update(extrap_t3, diis_coeffs[i], q->t3[i]);
        }
    }

    finalize:
    cc_free(B);
    cc_free(right);
    cc_free(diis_coeffs);

    timer_stop("diis");
}

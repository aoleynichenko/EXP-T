/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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

#include "linalg.h"

#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "memory.h"
#include "cblas_lapacke.h"


void ascrea(double complex *arr, int32_t *order, size_t n);


/**
 * Solves eigenproblem for complex non-hermitian matrix. See documentation for
 * the 'zgeev' routine from LAPACK for details.
 *
 * @param n     matrix size
 * @param A     n x n matrix to be diagonalized (on exit, A has been overwritten)
 * @return ev   eigenvalues
 * @return vl   left eigenvectors (stored by rows, in the same order as their eigenvalues)
 * @return vr   right eigenvectors (stored by rows, in the same order as their eigenvalues)
 * @return 0 if success, 1 if fail
 */
int eigensolver(int32_t n, double complex *A, double complex *ev, double complex *vl, double complex *vr)
{
    /*
     * solve eigenproblem using LAPACK
     */
    int32_t lda = n;
    int32_t ldvl = n;
    int32_t ldvr = n;
    int32_t info = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'V', 'V', n, A, lda, ev, vl, ldvl, vr, ldvr);
    if (info > 0) {
        printf("The algorithm failed to compute eigenvalues (zgeev).\n");
        return 1;
    }

    /*
     * now eigenvectors are biorthogonal, but eigenvalues are not ordered
     * in ascending order.
     * we reorder eigenvalues and corresponding eigenvectors:
     */
    int32_t *order = (int32_t *) cc_malloc(n * sizeof(int32_t));
    double complex *vectmp = (double complex *) cc_malloc(n * sizeof(double complex));
    ascrea(ev, order, n);

    for (size_t i = 0; i < n; i++) {
        // reorder left eigenvectors
        for (size_t j = 0; j < n; j++) {
            vectmp[j] = vl[n * order[j] + i];
        }
        for (size_t j = 0; j < n; j++) {
            vl[n * j + i] = vectmp[j];
        }
        // reorder right eigenvectors
        for (size_t j = 0; j < n; j++) {
            vectmp[j] = vr[n * order[j] + i];
        }
        for (size_t j = 0; j < n; j++) {
            vr[n * j + i] = vectmp[j];
        }
    }

    /*
     * normalize vectors to get a biorthonormal set
     */
    for (size_t i = 0; i < n; i++) {
        double complex c = 0.0 + 0.0 * I;
        for (size_t j = 0; j < n; j++) {
            c += conj(vl[i * n + j]) * vr[i * n + j];
        }
        for (size_t j = 0; j < n; j++) {
            vl[i * n + j] /= conj(c);
        }
    }

    cc_free(order);
    cc_free(vectmp);

    return 0;
}


/**
 * computes eigenvalues eva /eigenvectors a  of complex hermitian matrix A.
 * matrix A is destructed.
 * eigenvalues in descending order.
 * eigenvectors are stored by rows.
 *
 * @note uses eigensolver() = zgeev, in fact
 */
void eigher(double complex *A, double *eva, int32_t n)
{

    double complex *eva_buf;
    double complex *A_buf;
    double complex *vl_buf;

    A_buf = z_zeros(n, n);
    vl_buf = z_zeros(n, n);
    eva_buf = z_zeros(n, 1);
    memcpy(A_buf, A, sizeof(double complex) * n * n);

    eigensolver(n, A_buf, eva_buf, vl_buf, A);

    for (int i = 0; i < n; i++) {
        eva[i] = creal(eva_buf[i]);
    }

    cc_free(A_buf);
    cc_free(vl_buf);
    cc_free(eva_buf);
}


/**
 * Sorts an array of complex numbers in ascending order, only real parts are
 * compared (ASCending REAl parts).
 *
 * @param arr  array of complex numbers to be sorted
 * @param ord  resulting permutation of the "arr" elements
 * @param n    vector length
 */
void ascrea(double complex *arr, int32_t *order, size_t n)
{
    double complex ztmp;
    size_t i, j;
    int32_t itmp;

    for (i = 0; i < n; i++) {
        order[i] = i;
    }

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (creal(arr[i]) <= creal(arr[j])) {
                continue;
            }
            // swap!
            ztmp = arr[i];
            arr[i] = arr[j];
            arr[j] = ztmp;
            itmp = order[i];
            order[i] = order[j];
            order[j] = itmp;
        }
    }
}

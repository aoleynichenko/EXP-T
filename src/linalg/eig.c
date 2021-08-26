/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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

#include "memory.h"

// ZGEEV prototype -- is required by the GNU compiler.
// Note: zgeev_ declaration seems to be OK only in case LAPACK
// is built with 4-byte integers (?)
// maybe it is better simply to ignore compiler's warnings
#if defined(__GNUC__) && !defined(BLAS_MKL)
extern void zgeev_(char *jobvl, char *jobvr, int *n, double complex *a,
        int *lda, double complex *w, double complex *vl, int *ldvl,
        double complex *vr, int *ldvr, double complex *work, int *lwork,
        double *rwork, int *info);
extern void LAPACKE_zgeev(int matrix_layout, char jobvl, char jobvr, int n,
        complex double *a, int lda, complex double *w, complex double *vl,
        int ldvl, complex double *vr, int ldvr);
#endif /* COMPILER_GNU */

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
int eig(int32_t n, double complex *A, double complex *ev, double complex *vl, double complex *vr)
{
    int32_t lda, ldvl, ldvr, info, lwork;
    double complex wkopt;
    double complex *work;
    double *rwork;
    int32_t *order;
    double complex *vectmp;
    size_t i, j;

    lda = n;
    ldvl = n;
    ldvr = n;
    rwork = (double *) cc_malloc(2 * n * sizeof(double));

    // Query and allocate the optimal workspace
    lwork = -1;
#ifdef BLAS_MKL // in order to avoid warnings during compilation
    zgeev_("Vectors", "Vectors", &n, (MKL_Complex16 *) A, &lda,
           (MKL_Complex16 *) ev, (MKL_Complex16 *) vl, &ldvl, (MKL_Complex16 *) vr, &ldvr,
           (MKL_Complex16 *) &wkopt, &lwork, rwork, &info);
#else
    zgeev_("Vectors", "Vectors", &n, A, &lda, ev, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, rwork, &info);
#endif
    lwork = (int) creal(wkopt);

    work = (double complex *) cc_malloc(lwork * sizeof(double complex));

    // Solve eigenproblem
#ifdef BLAS_MKL
    // what is the size of MKL_INT? 4 or 8? Seems to be 8
    LAPACKE_zgeev(CblasColMajor, 'V', 'V',
                  n, (MKL_Complex16 *) A, lda, (MKL_Complex16 *) ev,
                  (MKL_Complex16 *) vl, ldvl, (MKL_Complex16 *) vr, ldvr);
#elif defined BLAS_OPENBLAS
    LAPACKE_zgeev(CblasColMajor, 'V', 'V', n, A, lda, ev, vl, ldvl, vr, ldvr);
#else
    zgeev_("Vectors", "Vectors", &n, A, &lda, ev, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
#endif

    // Check for convergence
    if (info > 0) {
        printf("The algorithm failed to compute eigenvalues (zgeev).\n");
        return 1;
    }

    cc_free(work);
    cc_free(rwork);

    // now eigenvectors are biorthogonal, but eigenvalues are not ordered
    // in ascending order.
    // we reorder eigenvalues and corresponding eigenvectors:
    order = (int32_t *) cc_malloc(n * sizeof(int32_t));
    vectmp = (double complex *) cc_malloc(n * sizeof(double complex));
    ascrea(ev, order, n);

    for (i = 0; i < n; i++) {
        // reorder left eigenvectors
        for (j = 0; j < n; j++) {
            vectmp[j] = vl[n * order[j] + i];
        }
        for (j = 0; j < n; j++) {
            vl[n * j + i] = vectmp[j];
        }
        // reorder right eigenvectors
        for (j = 0; j < n; j++) {
            vectmp[j] = vr[n * order[j] + i];
        }
        for (j = 0; j < n; j++) {
            vr[n * j + i] = vectmp[j];
        }
    }

    // normalize vectors to get biorthonormal set
    for (i = 0; i < n; i++) {
        double complex c = 0.0 + 0.0 * I;
        for (j = 0; j < n; j++) {
            c += conj(vl[i * n + j]) * vr[i * n + j];
        }
        for (j = 0; j < n; j++) {
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
 * @note uses eig() = zgeev, in fact
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

    eig(n, A_buf, eva_buf, vl_buf, A);

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

/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
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

#include <math.h>
#include <stdint.h>

#include "cblas_lapacke.h"
#include "memory.h"

void sort_vectors(int n, double complex *ev, double complex *vl, double complex *vr,
                  int (*cmp)(double complex, double complex));

void sort_record_order(double complex *arr, int32_t *order, size_t n,
                       int (*cmp)(double complex, double complex));

int complex_cmp_descending(double complex a, double complex b);


/**
 * Singular value decomposition.
 *    A = U * lambda * V^H
 * U, V -- unitary matrices
 * lambda -- diagonal matrix; non-negative REAL numbers
 * H -- hermitian conjugated
 *
 * Algorithm notes:
 * U -- eigenvectors of the A*A^H matrix, lambda^2 = its eigenvalues
 * V -- the same for the A^H*A matrix
 * SVD is implemented via two diagonalizations, without SVD subroutine from BLAS
 *
 * Input:
 *   square complex-valued matrices
 *   pre-allocated arrays A, U, V (of size n x n) and 'lambda' (vector dim = n)
 *   NOTE: this code can be used only for square matrices
 * Returns:
 *   lambda -- singular values, in descending order
 *   U      -- left-singular vectors (stored by rows)
 *   V      -- right-singular vectors (stored by rows)
 *
 */
void svd(int n, double complex *A, double *lambda, double complex *U, double complex *V)
{
    size_t nbytes = n * n * sizeof(double complex);
    double complex *AAH = cc_calloc(sizeof(double complex), n * n);
    double complex *AHA = cc_calloc(sizeof(double complex), n * n);
    double complex *lambda2 = cc_calloc(sizeof(double complex), 1 * n);
    double complex alpha = 1.0 + 0.0 * I;
    double complex beta = 0.0 + 0.0 * I;

    // A*A^H
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, n, n, n, &alpha, A, n, A, n, &beta, AAH, n);
    eigensolver(n, AAH, lambda2, AHA /* used as buffer */, U);
    sort_vectors(n, lambda2, AHA /* buffer */, U, complex_cmp_descending);

    // A^H*A
    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, n, n, n, &alpha, A, n, A, n, &beta, AHA, n);
    eigensolver(n, AHA, lambda2, AAH /* buffer */, V);
    sort_vectors(n, lambda2, AAH /* buffer */, V, complex_cmp_descending);

    // calculate sqrt of eigenvalues in order to obtain singular values
    for (int i = 0; i < n; i++) {
        if (creal(lambda2[i]) < 1e-15) {
            lambda[i] = 0.0;
        }
        else {
            lambda[i] = sqrt(creal(lambda2[i]));
        }
    }

    cc_free(AAH);
    cc_free(AHA);
    cc_free(lambda2);
}


/*******************************************************************************
 * sort_vectors
 *
 * Sort array of eigenvalues and corresponding eigenvectors(both left and right)
 * using some function-comparator.
 ******************************************************************************/
void sort_vectors(int n, double complex *ev, double complex *vl, double complex *vr,
                  int (*cmp)(double complex, double complex))
{
    int32_t *order;
    double complex *vectmp;
    size_t i, j;

    order = (int32_t *) cc_malloc(n * sizeof(int32_t));
    vectmp = (double complex *) cc_malloc(n * sizeof(double complex));
    sort_record_order(ev, order, n, cmp);

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

    cc_free(order);
    cc_free(vectmp);
}


/*******************************************************************************
 * subroutine sort_record_order
 *
 * Sorts an array of complex numbers using comparator function.
 * New order will be recorded.
 *
 * Arguments:
 *   arr  --  array of complex numbers to be sorted
 *   ord  --  resulting permutation of the "arr" elements
 *   n    --  vector size
 *   cmp  --  comparator function
 ******************************************************************************/
void sort_record_order(double complex *arr, int32_t *order, size_t n,
                       int (*cmp)(double complex, double complex))
{
    double complex ztmp;
    size_t i, j;
    int32_t itmp;

    for (i = 0; i < n; i++) {
        order[i] = i;
    }

    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (cmp(arr[i], arr[j]) < 0) {
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


// helper function
int complex_cmp_descending(double complex a, double complex b)
{
    double re_a = creal(a);
    double re_b = creal(b);
    return re_b - re_a;
}

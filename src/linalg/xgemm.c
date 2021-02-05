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

#include <assert.h>
#include <string.h>


/**
 * Matrix-matrix multiplication. The choice of the actual BLAS
 * subroutine depends on the data type used (the first parameter, int data_type)
 *
 * @param data_type CC_DOUBLE or CC_DOUBLE_COMPLEX
 * @param trans_a "N" -- do nothing, "T" -- transpose, "C" -- transpose and complex conjugate
 * @param trans_b see trans_a
 * @param other arguments - see BLAS or MKL documentation
 * @return C = alpha * A^trans_a * B^trans_b + beta * C
 * @note matrices must be stored in the row-major (C) format!
 */
void xgemm(
        data_type_t data_type,
        char *trans_a, char *trans_b,    // "N", "T" or "C"
        int m, int n, int k,
        void *alpha, void *A, int lda, void *B, int ldb,
        void *beta, void *C, int ldc
)
{
    CBLAS_TRANSPOSE trans_a_op, trans_b_op;

    assert(data_type == CC_DOUBLE || data_type == CC_DOUBLE_COMPLEX);
    assert(strcmp(trans_a, "N") == 0 || strcmp(trans_a, "T") == 0 || strcmp(trans_a, "C") == 0);
    assert(strcmp(trans_b, "N") == 0 || strcmp(trans_b, "T") == 0 || strcmp(trans_b, "C") == 0);

    if (strcmp(trans_a, "N") == 0) {
        trans_a_op = CblasNoTrans;
    }
    else if (strcmp(trans_a, "T") == 0) {
        trans_a_op = CblasTrans;
    }
    else {
        trans_a_op = CblasConjTrans;
    }

    if (strcmp(trans_b, "N") == 0) {
        trans_b_op = CblasNoTrans;
    }
    else if (strcmp(trans_b, "T") == 0) {
        trans_b_op = CblasTrans;
    }
    else {
        trans_b_op = CblasConjTrans;
    }

    if (data_type == CC_DOUBLE) {
        double dalpha = *((double *) alpha);
        double dbeta = *((double *) beta);
        cblas_dgemm(CblasRowMajor, trans_a_op, trans_b_op,
                    m, n, k, dalpha, A, lda, B, ldb, dbeta, C, ldc);
    }
    else { // CC_DOUBLE_COMPLEX
        double complex zalpha = *((double complex *) alpha);
        double complex zbeta = *((double complex *) beta);
        cblas_zgemm(CblasRowMajor, trans_a_op, trans_b_op,
                    m, n, k, &zalpha, A, lda, B, ldb, &zbeta, C, ldc);
    }
}

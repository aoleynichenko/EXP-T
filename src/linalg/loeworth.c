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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "memory.h"

void eigher(double complex *A, double *eva, int32_t n);


/**
 * Symmetric (Loewdin) orthogonalization
 *
 * @param n      number of vectors to be orthogonalized
 * @param C      input vectors (stored row-wise). Will be overwritten
 * @param C_orth output vectors (orthogonalized, row-wise)
 * @param print  print level (0 - nothing, 1 - resulting overlap, > 1 - all)
 * @return C_orth
 *
 * @author A. Zaitsevskii (as a part of the HEFFMAN code)
 * @author A. Oleynichenko
 */
void loewdin_orth(size_t n, double complex *C, double complex *C_orth, int print)
{
    double complex *S, *U, *S12;
    double *s_eigval;
    double complex alpha = 1.0 + 0.0 * I;
    double complex beta = 0.0 + 0.0 * I;

    size_t nbytes = sizeof(double complex) * n * n;
    S = cc_malloc(nbytes);
    U = cc_malloc(nbytes);
    S12 = cc_malloc(nbytes);
    s_eigval = cc_malloc(sizeof(double) * n);

    if (print >= 2) {
        printf(" *** LOEWDIN ORTHOGONALIZATION ***\n");
        xprimat(CC_COMPLEX, C, n, n, "Input non-orthogonal vectors");
    }

    overlap(n, n, C, C, S);
    if (print >= 2) {
        xprimat(CC_COMPLEX, S, n, n, "OVERLAP");
    }

    // diagonalize S (S is hermitian matrix)
    // S_diag = U^H * S * U
    // U -- eigenvectors, in mathematical notation they are stored columnwise
    eigher(S, s_eigval, n);
    if (print >= 2) {
        xprimat(CC_COMPLEX, S, n, n, "EIGHER eigenvectors");
    }
    memcpy(U, S, n * n * sizeof(double complex));

    // S^{-1/2}_diag = 1/sqrt(S_diag)
    // inverse transformation: S = U * S12_diag * U^H
    // but we store C in "transposed" manner => S = U^H * S_12 * U (why ?!)
    // "undiagonalization"
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double complex zzz = 0.0 + 0.0 * I;
            for (int k = 0; k < n; k++) {
                zzz += U[k * n + i] * conj(U[k * n + j]) / sqrt(s_eigval[k]);
            }
            S12[j * n + i] = zzz;
        }
    }
    if (print >= 2) {
        xprimat(CC_COMPLEX, S12, n, n, "S^{-1/2}");
    }

    // C' = S^{-1/2} * C
    // But in fact: C' = C^T S^{-1/2}
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double complex zzz = 0.0 + 0.0 * I;
            for (int k = 0; k < n; k++) {
                zzz = zzz + C[k * n + i] * S12[k * n + j];
            }
            U[j * n + i] = zzz;
        }
    }
    if (print >= 2) {
        xprimat(CC_COMPLEX, U, n, n, "COEF ORTHOGONALIZED");
    }
    memcpy(C_orth, U, nbytes);

    // check orthogonalization
    overlap(n, n, C_orth, C_orth, S);
    if (print >= 1) {
        xprimat(CC_COMPLEX, S, n, n, "NEW OVERLAP");
    }

    cc_free(S);
    cc_free(U);
    cc_free(s_eigval);
    cc_free(S12);
}

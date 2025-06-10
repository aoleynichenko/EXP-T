/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
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
#include "memory.h"

/**
 * Calculates overlap matrix S = <C1|C2>, where C1, C2 -- n x n square matrices
 * (matrix version of the dot product operation)
 *
 * @param n  dimension of all square matrices
 * @param C1 coeff matrix 1 (vectors are stored row-wise)
 * @param C2 coeff matrix 2
 * @param S  (output) overlap matrix
 *
 * @note S matrix must be pre-allocated
 * @note vectors must be stored row-wise
 *
 * @author 2017-2019 Andrei Zaitsevskii
 * @author 2019 Alexander Oleynichenko
 */
void overlap(size_t dim, size_t n_states, double complex *C1, double complex *C2, double complex *S)
{
    // Naive implementation:

    size_t i, j, k;

    for (i = 0; i < n_states; i++) {
        for (j = 0; j < n_states; j++) {
            double complex s = 0.0 + 0.0 * I;
            for (k = 0; k < dim; k++) {
                s = s + conj(C1[i * dim + k]) * C2[j * dim + k];
            }
            S[i * n_states + j] = s;
        }
    }


    // LAPACK-based implementation (for square matrices):

    /*
    double complex alpha = 1.0 + 0.0*I;
    double complex beta  = 0.0 + 0.0*I;

    // if vectors in C1 and C2 are stored row-wise,
    // S* = C1 x C2^T
    // S = conj(C1 x C2^T)
    xgemm(CC_COMPLEX, "N", "C", n, n, n, &alpha, C1, n, C2, n, &beta, S, n);
    for (size_t i = 0; i < n * n; i++) {
        S[i] = conj(S[i]);
    }
    */
}

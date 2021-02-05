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
void overlap(size_t n, double complex *C1, double complex *C2, double complex *S)
{
    size_t i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double complex s = 0.0 + 0.0 * I;
            for (k = 0; k < n; k++) {
                s = s + conj(C1[i * n + k]) * C2[j * n + k];
            }
            S[i * n + j] = s;
        }
    }
}

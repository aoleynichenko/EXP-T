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


/*
 * Construct projector onto non-orthogonal basis in the basis of Slater determinants.
 * Basis vectors are stored row-wise in the 'C' matrix.
 *
 * dim      - dimension of the basis of Slater determinants
 * n_states - number of states in the basis
 * C        - matrix of basis vectors (model vectors), stored row-wise
 *            dimensions: 'n_states' rows, 'dim' columns
 * P        - (target) square matrix of the projector of size dim x dim
 */
void construct_projector(int dim, int n_states, double complex *C, double complex *P)
{
    if (n_states == 0) {
        return;
    }
    double complex *S = z_zeros(n_states, n_states);
    double complex *S_inv = z_zeros(n_states, n_states);

    overlap(dim, n_states, C, C, S);
    inverse_matrix(n_states, S, S_inv);

    for (int a = 0; a < dim; a++) {
        for (int b = 0; b < dim; b++) {
            double complex P_ab = 0.0 + 0.0*I;

            for (int mu = 0; mu < n_states; mu++) {
                for (int nu = 0; nu < n_states; nu++) {
                    P_ab += C[mu*dim+a] * S_inv[mu*n_states+nu] * conj(C[nu*dim+b]);
                }
            }

            P[a*dim + b] = P_ab;
        }
    }

    cc_free(S);
    cc_free(S_inv);
}



/*
 * Construct projector onto the basis of bi-orthogonal model vectors.
 * in the basis of Slater determinants.
 * Basis vectors are stored row-wise in the 'C_left' and 'C_right' matrix.
 *
 * dim      - dimension of the basis of Slater determinants
 * n_states - number of states in the basis
 * C_left   - matrix of basis vectors (left model vectors), stored row-wise
 *            dimensions: 'n_states' rows, 'dim' columns
 * C_right  - the same, but for the right model vectors
 * P        - (target) square matrix of the projector of size dim x dim
 */
void construct_biorth_projector(int dim, int n_states, double complex *C_left, double complex *C_right, double complex *P)
{
    if (n_states == 0) {
        return;
    }
    double complex *S = z_zeros(n_states, n_states);
    double complex *S_inv = z_zeros(n_states, n_states);

    overlap(dim, n_states, C_left, C_right, S);

    for (int a = 0; a < dim; a++) {
        for (int b = 0; b < dim; b++) {
            double complex P_ab = 0.0 + 0.0*I;

            for (int mu = 0; mu < n_states; mu++) {
                P_ab += C_right[mu*dim+a] * conj(C_left[mu*dim+b]);
            }

            P[a*dim + b] = P_ab;
        }
    }

    cc_free(S);
    cc_free(S_inv);
}

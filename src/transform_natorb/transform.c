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

#include "transform.h"

#include <complex.h>
#include <stdlib.h>
#include <string.h>

void zaxpy(int n, double complex a, const double complex *x, double complex *y);
void zaxpy_conj(int n, double complex a, const double complex *x, double complex *y);


/**
 * transforms natural spinors
 * from the basis of molecular spinors to the atomic orbital (AO) basis.
 */
void transform(
        int n_nat_spinors, int n_act_spinors, int basis_dim,
        int *map_act_to_all, int *kramers_conj,
        // natural orbitals expanded in the basis of (active) molecular spinors
        // NOs are stored by rows
        // n_nat_spinors x n_act_spinors
        double complex *nat_coeffs_spinor_basis,
        // molecular spinors expanded in the AO basis
        // n_spinors x basis_dim
        double complex *cmo_alpha, double complex *cmo_beta,
        // (quasi) natural spinors expanded in the AO basis
        // n_nat_spinors x basis_dim
        double complex *cno_alpha, double complex *cno_beta
)
{
    for (int i_nat = 0; i_nat < n_nat_spinors; i_nat++) {

        double complex *nat_coef_alpha = cno_alpha + i_nat * basis_dim;
        double complex *nat_coef_beta = cno_beta + i_nat * basis_dim;

        memset(nat_coef_alpha, 0, sizeof(double complex) * basis_dim);
        memset(nat_coef_beta, 0, sizeof(double complex) * basis_dim);

        for (int i_act = 0; i_act < n_act_spinors; i_act++) {

            int spinor_idx = map_act_to_all[i_act];

            double complex coef = nat_coeffs_spinor_basis[i_nat * n_nat_spinors + i_act];
            double complex *spinor_alpha = cmo_alpha + spinor_idx * basis_dim;
            double complex *spinor_beta = cmo_beta + spinor_idx * basis_dim;

            if (kramers_conj[i_act] == 0) {
                zaxpy(basis_dim, coef, spinor_alpha, nat_coef_alpha);
                zaxpy(basis_dim, coef, spinor_beta, nat_coef_beta);
            }
            else { // use Kramers partner
                zaxpy_conj(basis_dim, - coef, spinor_beta, nat_coef_alpha);
                zaxpy_conj(basis_dim, coef, spinor_alpha, nat_coef_beta);
            }
        }

    }
}


/**
 * constant times a vector plus a vector:
 * y = a * x + y
 */
void zaxpy(int n, double complex a, const double complex *x, double complex *y)
{
    for (int i = 0; i < n; i++) {
        y[i] = a * x[i] + y[i];
    }
}


/**
 * constant times a (complex conjugated) vector plus a vector:
 * y = a * conj(x) + y
 */
void zaxpy_conj(int n, double complex a, const double complex *x, double complex *y)
{
    for (int i = 0; i < n; i++) {
        y[i] = a * conj(x[i]) + y[i];
    }
}




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


#include "print_vectors.h"

#include <stdio.h>
#include <math.h>
#include <stdio.h>

void print_vectors(FILE *out_file,
                   int n_vectors, int basis_dim,
                   const double *eigenvalues,
                   const double complex *coeffs_alpha,
                   const double complex *coeffs_beta,
                   char **basis_fun_labels)
{
    const double PRINT_THRESH = 1e-10;

    printf(" > vectors (DIRAC output format):\n");
    fprintf(out_file, "\n");
    fprintf(out_file, "  Number of large orbitals in each symmetry:     %d\n", basis_dim);
    fprintf(out_file, "No. of positive energy orbitals (NESH):    %d\n", n_vectors);
    fprintf(out_file, "\n");

    for (int ivec = 0; ivec < n_vectors; ivec++) {

        fprintf(out_file, "* Electronic eigenvalue no.%4d: %.12f\n", ivec + 1, eigenvalues[ivec]);
        fprintf(out_file, "====================================================\n");

        for (int ibas = 0; ibas < basis_dim; ibas++) {
            double complex coef_alpha = coeffs_alpha[ivec * basis_dim + ibas];
            double complex coef_beta = coeffs_beta[ivec * basis_dim + ibas];
            double alpha_re = creal(coef_alpha);
            double alpha_im = cimag(coef_alpha);
            double beta_re = creal(coef_beta);
            double beta_im = cimag(coef_beta);

            double abs_sum = fabs(alpha_re) + fabs(alpha_im) + fabs(beta_re) + fabs(beta_im);

            if (abs_sum < PRINT_THRESH) {
                continue;
            }

            char *label = basis_fun_labels[ibas];
            printf("%8d  %-15s%14.10f%14.10f%14.10f%14.10f\n",
                   ibas + 1, (label != NULL) ? label : "", alpha_re, alpha_im, -beta_re, beta_im);
        }

        printf("\n");
    }
}

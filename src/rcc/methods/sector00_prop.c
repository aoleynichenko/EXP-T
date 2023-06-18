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


#include "cc_properies.h"
#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ccutils.h"
#include "engine.h"
#include "datamodel.h"
#include "diis.h"
#include "heff.h"
#include "linalg.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void save_density_matrix(char *path, int nspinors, double complex *dm);

int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

void guess_operator_symmetry(int nspinors, double complex *prop_matrix);

double complex calculate_vacuum_expectation_value(int nspinors, double complex *prop_matrix);

double complex sector_0h0p_calculate_overlap(int nspinors, double complex *dm);

void sector_0h0p_analytic_density_matrix_expectation(int pt_order);


/**
 * Contracts density matrix in the 0h0p sector with property matrix.
 */
void sector_0h0p_calculate_properties()
{
    if (cc_opts->n_analyt_prop == 0) {
        return;
    }

    /*
     * read density matrix from file
     */
    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);
    read_prop_single_file(nspinors, "density_0h0p.txt", dm);

    /*
     * loop over property matrices
     */
    for (int iprop = 0; iprop < cc_opts->n_analyt_prop; iprop++) {

        char *prop_file_name = cc_opts->analyt_prop_files[iprop];
        printf("\n Property: %s\n", prop_file_name);
        printf(" -----------------------------------\n\n");

        /*
         * read matrix of the operator in the spinor basis
         */
        int nspinors = get_num_spinors();
        double complex *prop_matrix = z_zeros(nspinors, nspinors);
        int status = read_prop_single_file(nspinors, prop_file_name, prop_matrix);
        if (status == EXIT_FAILURE) {
            cc_free(prop_matrix);
            continue;
        }

        guess_operator_symmetry(nspinors, prop_matrix);

        /*
         * calculate correlation contribution to the expectation value
         */
        double complex expect_value = 0.0;
        for (int p = 0; p < nspinors; p++) {
            for (int q = 0; q < nspinors; q++) {
                expect_value += dm[p * nspinors + q] * prop_matrix[p * nspinors + q];
            }
        }

        /*
         * print expectation value
         */
        double complex scf_contrib = calculate_vacuum_expectation_value(nspinors, prop_matrix);
        double complex ccsd_contrib = expect_value - scf_contrib;

        printf(" SCF contribution          %20.12f%20.12f\n", creal(scf_contrib), cimag(scf_contrib));
        printf(" Correlation contribution  %20.12f%20.12f\n", creal(ccsd_contrib), cimag(ccsd_contrib));
        printf(" Final expectation value   %20.12f%20.12f\n", creal(expect_value), cimag(expect_value));

        cc_free(prop_matrix);
    }

    /*
     * finalize
     */
    cc_free(dm);
}


/**
 * Overlap integral <psi|psi> for the CC wave function |psi> = e^T |0>.
 *
 * Only "correlation" part of the density matrix is used.
 * Basic idea: contraction of the identity matrix with the density matrix.
 * Diagonal matrix elements of the identity matrix should be multiplied by -1 for holes.
 */
double complex sector_0h0p_calculate_overlap(int nspinors, double complex *dm)
{
    double complex *delta_matrix = z_zeros(nspinors, nspinors);

    for (int i = 0; i < nspinors; i++) {
        delta_matrix[i * nspinors + i] = is_hole(i) ? -1 : +1;
    }

    double complex overlap = 1.0 + 0.0 * I;
    for (int p = 0; p < nspinors; p++) {
        for (int q = 0; q < nspinors; q++) {
            double complex d_pq = dm[p * nspinors + q];

            // get "correlation part" of the DM
            if (p == q && is_hole(p)) {
                d_pq -= 1.0;
            }

            overlap += d_pq * delta_matrix[p * nspinors + q];
        }
    }

    cc_free(delta_matrix);

    return overlap;
}


/**
 * Natural orbitals (spinors) and their occupation numbers in the 0h0p sector.
 */
void sector_0h0p_calculate_natural_spinors()
{
    /*
     * read density matrix from file
     */
    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);
    read_prop_single_file(nspinors, "density_0h0p.txt", dm);

    /*
     * DM diagonalization => occupation numbers & natural spinors.
     * natural spinors are stored by rows.
     */
    double complex *natorb_left = z_zeros(nspinors, nspinors);
    double complex *natorb_right = z_zeros(nspinors, nspinors);
    double complex *occ_numbers = z_zeros(nspinors, 1);

    eigensolver(nspinors, dm, occ_numbers, natorb_left, natorb_right);

    /*
     * print beautiful table
     */
    printf("\n");
    printf(" Occupation numbers of natural spinors:\n");
    printf(" --------------------------------------\n");
    printf("\n");

    for (int i = nspinors - 1; i >= 0; i--) {
        double occ_num = creal(occ_numbers[i]);

        // find spinor with max weight => determine NO symmetry
        double max_coef = 0.0;
        int max_index = 0;
        for (int j = 0; j < nspinors; j++) {
            double complex coef = natorb_right[nspinors * i + j];
            if (cabs(coef) > max_coef) {
                max_coef = cabs(coef);
                max_index = j;
            }
        }

        char *natorb_sym = get_irrep_name(get_spinor_irrep(max_index));
        int seq_no = nspinors - i;

        printf("%4d %-6s%10.7f    ", seq_no, natorb_sym, occ_num);
        if (seq_no % 5 == 0) {
            printf("\n");
        }
    }
    printf("\n\n");

    /*
     * clean up
     */
    cc_free(natorb_left);
    cc_free(natorb_right);
    cc_free(occ_numbers);
    cc_free(dm);
}

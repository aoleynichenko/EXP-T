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

#include "natural_spinors.h"

#include <math.h>
#include <stdio.h>

#include "engine.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "finite_order_overlap.h"


int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);


/**
 * Calculates natural spinors (NS) and natural occupation numbers for the pure electronic state.
 *
 * Natural spinors will be written to the formatted txt file 'natural_spinors_*.txt'
 */
void calculate_natural_spinors_and_occ_numbers(char *path_to_dm_file, char *path_to_nat_spinors_file)
{
    /*
     * read density matrix from file
     */
    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);
    read_prop_single_file(nspinors, path_to_dm_file, dm);

    /*
     * DM diagonalization => occupation numbers & natural spinors.
     * natural spinors are stored by rows.
     */
    double complex *natorb_left = z_zeros(nspinors, nspinors);
    double complex *natorb_right = z_zeros(nspinors, nspinors);
    double complex *occ_numbers = z_zeros(nspinors, 1);

    eigensolver(nspinors, dm, occ_numbers, natorb_left, natorb_right);

    /*
     * print beautiful table with occupation numbers
     */
    printf("\n");
    printf(" occupation numbers of natural spinors:\n");
    printf(" --------------------------------------\n");
    printf("\n");

    double occ_num_sum = 0.0;
    for (int i = nspinors - 1; i >= 0; i--) {
        double occ_num = creal(occ_numbers[i]);
        occ_num_sum += occ_num;

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
    if (nspinors % 5 != 0) {
        printf("\n");
    }
    printf("\n");
    printf(" sum of occupation numbers: %.7f\n", occ_num_sum);
    printf("\n");

    /*
     * print brief information about expansions of natural spinors to stdout
     */
    double const OCC_THRESH = 1e-2;    // threshold for printing natural orbitals
    double const COEF_THRESH = 1e-1;   // threshold for printing model vec-s coeff-s

    printf("\n");
    printf(" natural spinors - expansion coefficients\n");
    printf(" ----------------------------------------\n");
    printf("\n");
    printf(" threshold for occ numbers : %.2e\n", OCC_THRESH);
    printf(" threshold for coefficients: %.2e\n", COEF_THRESH);
    printf("\n");

    for (int i = nspinors - 1; i >= 0; i--) {
        double complex occ_num = occ_numbers[i];
        double occ_num_re = creal(occ_num);
        double occ_num_im = cimag(occ_num);

        if (fabs(creal(occ_num)) < OCC_THRESH) {
            continue;
        }

        int seq_no = nspinors - i;
        printf(" [%d] occ number = %.6f\n", seq_no, occ_num_re);
        if (occ_num_im > 1e-13) {
            printf(" imaginary occupation number: %.3e %.3e\n", occ_num_re, occ_num_im);
        }

        print_natural_spinor(stdout, natorb_right + nspinors * i, COEF_THRESH);
        printf("\n");
    }

    /*
     * write natural spinors to formatted files
     */
    double *occ_numbers_re = (double *) cc_calloc(nspinors, sizeof(double));
    get_real_parts(nspinors, occ_numbers, occ_numbers_re);
    write_natural_spinors_formatted(path_to_nat_spinors_file, occ_numbers_re, natorb_right, 1e-16);

    /*
     * clean up
     */
    cc_free(natorb_left);
    cc_free(natorb_right);
    cc_free(occ_numbers);
    cc_free(occ_numbers_re);
    cc_free(dm);
}


/**
 * This subroutine calculates natural transition spinors (NTS) for the pair of electronic states.
 * It also calculates singular values corresponding to each pair (from - to) of transition spinors.
 *
 * Natural transition spinors will be written to the formatted txt files
 * 'natural_spinors_*_from.txt' and 'natural_spinors_*_to.txt'
 */
void calculate_natural_transition_spinors_and_singular_values(
    char *path_to_tran_dm_file,
    char *path_to_nat_spinors_from_file,
    char *path_to_nat_spinors_to_file
)
{
    double const OCC_THRESH = 1e-2;    // threshold for printing natural orbitals
    double const COEF_THRESH = 1e-1;   // threshold for printing model vec-s coeff-s

    /*
     * read density matrix from file
     */
    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);
    read_prop_single_file(nspinors, path_to_tran_dm_file, dm);

    /*
     * calculate natural transition spinors (NTOs)
     * SVD of transition density matrix
     *
     * SVD: T = U*S*VH
     * U = eigenvectors of T*TH
     * V = eigenvectors of TH*T
     * S = sqrt(eigenvalues of U or V)  <<-- singular values
     * pairs (eigenvalue,eigenvector) must be sorted in descending order (by eigenvalue)
     */
    double *lambda = d_zeros(nspinors, 1);      // NTO singular values
    double complex *nat_spinors_from = z_zeros(nspinors, nspinors);
    double complex *nat_spinors_to = z_zeros(nspinors, nspinors);

    // to <-> from swapped to restore the physical meaning
    svd(nspinors, dm, lambda, nat_spinors_to, nat_spinors_from);

    /*
     * print beautiful table with singular values
     */
    printf("\n");
    printf(" squared singular values of natural transition spinors:\n");
    printf(" ------------------------------------------------------\n");
    printf("\n");

    /*
     * print occupation numbers and NTS to stdout
     */
    double sum_lambda_sq = 0.0;  // sum of weights

    for (int i = 0; i < nspinors; i++) {

        double lambda_sq = lambda[i] * lambda[i];
        sum_lambda_sq += lambda_sq;

        // find spinor with max weight => determine NO symmetry
        double max_coef_from = 0.0;
        int max_index_from = 0;
        double max_coef_to = 0.0;
        int max_index_to = 0;

        for (int j = 0; j < nspinors; j++) {
            double complex coef_from = nat_spinors_from[nspinors * i + j];
            double complex coef_to = nat_spinors_to[nspinors * i + j];
            if (cabs(coef_from) > max_coef_from) {
                max_coef_from = cabs(coef_from);
                max_index_from = j;
            }
            if (cabs(coef_to) > max_coef_to) {
                max_coef_to = cabs(coef_to);
                max_index_to = j;
            }
        }

        char *natorb_sym_from = get_irrep_name(get_spinor_irrep(max_index_from));
        char *natorb_sym_to = get_irrep_name(get_spinor_irrep(max_index_to));

        printf("%4d %-6s-> %-6s%10.7f    ", i + 1, natorb_sym_from, natorb_sym_to, lambda_sq);
        if ((i + 1) % 3 == 0) {
            printf("\n");
        }
    }
    if (nspinors % 3 != 0) {
        printf("\n");
    }

    int sect_h = cc_opts->curr_sector_h;
    int sect_p = cc_opts->curr_sector_p;
    printf("\n");
    printf(" sum of weights = %.6f\n", sum_lambda_sq);
    printf(" (should be in range [0,%d])\n", max(sect_h, sect_p));
    printf("\n");

    /*
     * print brief information about expansions of natural spinors to stdout
     */
    printf("\n");
    printf(" natural transition spinors - expansion coefficients\n");
    printf(" ---------------------------------------------------\n");
    printf("\n");
    printf(" threshold for occ numbers : %.2e\n", OCC_THRESH);
    printf(" threshold for coefficients: %.2e\n", COEF_THRESH);
    printf("\n");

    for (int i = 0; i < nspinors; i++) {
        double lambda_sq = lambda[i] * lambda[i];

        if (lambda_sq < OCC_THRESH) {
            continue;
        }

        printf(" [%d] squared singular value (weight) = %.6f\n", i + 1, lambda_sq);
        printf("    from:\n");
        print_natural_spinor(stdout, nat_spinors_from + nspinors * i, COEF_THRESH);
        printf("    to:\n");
        print_natural_spinor(stdout, nat_spinors_to + nspinors * i, COEF_THRESH);
        printf("\n");
    }

    /*
     * write natural spinors to formatted files
     */
    write_natural_spinors_formatted(path_to_nat_spinors_from_file, lambda, nat_spinors_from, 1e-16);
    write_natural_spinors_formatted(path_to_nat_spinors_to_file, lambda, nat_spinors_to, 1e-16);

    /*
     * clean up
     */
    cc_free(dm);
    cc_free(lambda);
    cc_free(nat_spinors_from);
    cc_free(nat_spinors_to);
}


/**
 * prints natural [transition] spinor to an output stream.
 */
void print_natural_spinor(FILE *out, double complex *nat_spinor_coefs, double thresh)
{
    int n_spinors = get_num_spinors();

    fprintf(out, "       coef re     coef im      weight     mol spinor\n");

    for (size_t j = 0; j < n_spinors; j++) {
        double complex coef = nat_spinor_coefs[j];

        if (cabs(coef) < thresh) {
            continue;
        }

        int ispinor = j;
        int irrep = get_spinor_irrep(ispinor);
        char *irrep_name = get_irrep_name(irrep);
        double eps = get_eps(ispinor);
        double weight = cabs(coef) * cabs(coef);
        char *spinor_label = cc_opts->spinor_labels[ispinor];

        fprintf(out, "  %12.6f%12.6f%12.6f  %8s #%4d (%12.6f)  %s\n",
                creal(coef), cimag(coef), weight, irrep_name, ispinor + 1, eps, spinor_label ? spinor_label : "");
    }
}


/**
 * Write natural orbital (spinor) expansions to the formatted file with
 * the standard name natural_spinor_[H]h[P]p_[REP]:[STATE].dat
 */
void write_natural_spinors_formatted(
    char *natorb_file_name,
    double *occ_numbers, double complex *nat_orbs,
    double occ_thresh
)
{
    int n_spinors = get_num_spinors();

    // header contains brief info about spinors with numbering
    FILE *f = fopen(natorb_file_name, "w");
    fprintf(f, "dim %d\n", n_spinors);
    fprintf(f, "spinor info:\n");
    for (int i = 0; i < n_spinors; i++) {
        fprintf(f, "%4d%20.12f%10s\n", i + 1, get_eps(i),
                get_irrep_name(get_spinor_irrep(i)));
    }

    // calculate number of NOs to be printed (with occ_no >= thresh)
    int n_natorb = 0;
    for (int i = 0; i < n_spinors; i++) {
        if (fabs(occ_numbers[i]) >= occ_thresh) {
            n_natorb++;
        }
    }
    fprintf(f, "number of spinors: %d\n", n_natorb);

    // write NOs
    for (int i = n_spinors - 1; i >= 0; i--) {
        double occ = occ_numbers[i];
        if (fabs(occ) < occ_thresh) {
            continue;
        }

        fprintf(f, "occ %25.16e\n", occ);

        // write coefficients (all)
        for (size_t j = 0; j < n_spinors; j++) {
            double complex coef = nat_orbs[n_spinors * i + j];
            double coef_re = creal(coef);
            double coef_im = cimag(coef);
            int ispinor = j;

            if (fabs(coef_re) < 1e-16) {
                coef_re = 0.0;
            }
            if (fabs(coef_im) < 1e-16) {
                coef_im = 0.0;
            }

            if (fabs(coef_re) < 1e-16 && fabs(coef_im) < 1e-16) {
                continue;
            }

            fprintf(f, "%4d  %30.16e%30.16e\n", ispinor + 1, coef_re, coef_im);
        }
    }

    fclose(f);
}

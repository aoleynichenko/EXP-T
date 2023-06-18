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

/*
 * Tools for construction of quasi-natural spinors and quasi-natural
 * transition spinors.
 *
 * NOTE: density matrices computed are not true density matrices since only
 * model vectors are used to compute them. However, this do not seem to be a
 * problem, since density matrices, NOs and NTOs produced are used only for
 * semi-quantitative estimation and visualization purposes.
 *
 * Algorithm notes:
 * (1) Only 'valence' part of the density matrix is constructed =>
 *     dim = nacth * (h > 0) + nactp * (p > 0) for the (h,p) FS sector
 * (2) in all formulas for density matrices we must use global indices of
 *     spinors, however, rows and columns of DMs are enumerated with local
 *     indices; local indices are transformed into global when needed.
 * (3) the 1h1p sector is the special case:
 *     - calculation of approximate natural orbitals for the 0h0p state is
 *       meaningless (NOs = MOs)
 *     - 0h0p - 1h1p transitions ("CIS-like") require separate code for the DM
 *       evaluation (again, two cases: 0h0p->1h1p and 1h1p->0h0p, are not equal
 *       in case of non-hermitian Hamiltonian
 *     - 1h1p - 1h1p transitions can be calculated in the usual way
 *     NOTE: state numbering for the case of the 1h1p sector as a target
 *       must be quite different: state 1 == 0h0p, states from 1h1p: 2,3,4...
 *       (numbering is shifted by 1)
 */

#include "ms_nat_spinors.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "linalg.h"
#include "model_vectors.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "memory.h"
#include "mvcoef.h"


/*
 * functions used in this module
 */
int cmp_natorb_occupations(double complex occ1, double complex occ2);

void write_natural_spinors_formatted(
        char *natorb_file_name, int sect_h, int sect_p,
        int dim, double *occ_numbers, double complex *nat_orbs,
        double occ_thresh
);

int density_matrix_element(int sect_h, int sect_p, int p, int q, slater_det_t *bra, slater_det_t *ket);

void sort_vectors(int n, double complex *ev, double complex *vl, double complex *vr,
                  int (*cmp)(double complex, double complex));

void construct_model_space_natural_spinors(int sect_h, int sect_p, int irrep, int state);

void construct_model_space_natural_transition_spinors(
        int sect_h, int sect_p, int irrep1, int state1, int irrep2, int state2
);

void construct_model_space_density_matrix(int sect_h, int sect_p,
                                          int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
                                          int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
                                          double complex *dm, int *dim_dm);

void print_nat_spinor(FILE *out, int n_active, int *active_spinors, double complex *nat_spinor_coefs, double thresh);


/**
 * Calculates all the required natural spinors (NS) and natural transition spinors (NTS)
 * and writes them to formatted files.
 *
 * (main function of the module).
 */
void calculate_model_space_natural_spinors(int sect_h, int sect_p)
{
    // calculation of density matrices and natural orbitals
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_denmat) {
        for (int ipair = 0; ipair < cc_opts->n_denmat; ipair++) {
            int sect2_h = cc_opts->denmat_query[ipair].sect2[0];
            int sect2_p = cc_opts->denmat_query[ipair].sect2[1];
            if (!(sect2_h == sect_h && sect2_p == sect_p)) {
                continue;
            }
            int rep1 = get_rep_number(cc_opts->denmat_query[ipair].rep1_name);
            int rep2 = get_rep_number(cc_opts->denmat_query[ipair].rep2_name);
            int state1 = cc_opts->denmat_query[ipair].state1;
            int state2 = cc_opts->denmat_query[ipair].state2;

            if (rep1 == rep2 && state1 == state2) {
                construct_model_space_natural_spinors(sect_h, sect_p, rep1, state1);
            }
            else {
                construct_model_space_natural_transition_spinors(sect_h, sect_p, rep1, state1, rep2, state2);
            }
        }
    }
}


/**
 * Calculate natural spinors for the selected electronic state.
 * Write them to a formatted file.
 */
void construct_model_space_natural_spinors(int sect_h, int sect_p, int irrep, int state)
{
    double const OCC_THRESH = 1e-6;    // threshold for printing natural orbitals
    double const COEF_THRESH = 1e-4;   // threshold for printing model vec-s coeff-s
    char *rep_name = get_irrep_name(irrep);

    /*
     * header banner
     */
    printf("\n");
    printf("\n");
    printf(" **\n");
    printf(" ** Natural spinors for the electronic state: irrep %d (%s), state %d ***\n",
           irrep + 1, rep_name, state + 1);
    printf(" **\n");
    printf("\n");

    // only 'valence' parts of DMs are be constructed
    // dim DMs = n_active x n_active
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    get_active_space(sect_h, sect_p, &n_active, active_spinors);

    /*
     * extract target model vector and eigenvalue from the MVCOEF* unformatted file
     */
    int ms_size = 0;
    slater_det_t *det_list = NULL;
    double eigval = 0.0;
    double energy_cm = 0.0;
    double complex *coef_left = NULL;
    double complex *coef_right = NULL;

    char *mvcoef_file_name = NULL;
    if (sect_h == 1 && sect_p == 1 && irrep == get_vacuum_irrep()) {
        printf(" Mixed 0h0p+1h1p vectors will be used\n");
        mvcoef_file_name = cc_strdup("MVCOEF0011");
    }

    mvcoef_read_vectors_unformatted_for_state(
            sect_h, sect_p, mvcoef_file_name, rep_name, state,
            &ms_size, &det_list, &eigval, &energy_cm, &coef_left, &coef_right
    );

    /*
     * density matrix (DM) construction
     */
    printf(" Model density matrix construction...\n");

    double complex *denmat = z_zeros(n_active, n_active);
    construct_model_space_density_matrix(
            sect_h, sect_p,
            ms_size, coef_left, det_list,
            ms_size, coef_right, det_list,
            denmat, &n_active
    );

    /*
     * calculate natural spinors (NS's) by the diagonalization of the density matrix.
     *
     * тут происходит очень неустойчивая диагонализация матрицы с кучей нулевых собственных значений
     * но почему все нормально на неортогонализованных векторах?
     * потому что диагонализуемая матрица плотности неэрмитова!
     */
    double complex *nat_occ = z_zeros(n_active, 1);  // NS occ numbers
    double complex *nat_spinors_left = z_zeros(n_active, n_active);
    double complex *nat_spinors_right = z_zeros(n_active, n_active);

    eigensolver(n_active, denmat, nat_occ, nat_spinors_left, nat_spinors_right);
    sort_vectors(n_active, nat_occ, nat_spinors_left, nat_spinors_right, cmp_natorb_occupations);

    /*
     * print occupation numbers and NOs to stdout
     */
    printf(" Symmetry of this electronic state: %s\n", rep_name);
    printf(" Eigenvalue   = %.8f a.u.\n", eigval);
    printf(" Energy level = %.2f cm^-1\n", energy_cm);
    int count = 1;
    double complex sum_occ = 0.0 + 0.0 * I;  // trace

    for (int i = 0; i < n_active; i++) {
        sum_occ += nat_occ[i];

        if (cimag(nat_occ[i]) > 1e-13) {
            printf(" Imaginary occupation number: %.3e %.3e\n", creal(nat_occ[i]), cimag(nat_occ[i]));
        }
        if (fabs(creal(nat_occ[i])) < OCC_THRESH) {
            continue;
        }

        printf(" [%d] Occ number = %.6f\n", count++, creal(nat_occ[i]));
        print_nat_spinor(stdout, n_active, active_spinors, nat_spinors_right + n_active * i, COEF_THRESH);
    }
    printf(" Sum of occ numbers = %.6f (must be = %d)\n", creal(sum_occ), -sect_h + sect_p);

    /*
     * effective configuration for the state
     */
    printf(" Configuration (weights):\n");
    for (size_t j = 0; j < n_active; j++) {
        double weight = 0.0;
        for (size_t i = 0; i < n_active; i++) {
            double complex coef = nat_spinors_right[n_active * i + j];
            weight += nat_occ[i] * cabs(coef) * cabs(coef);
        }
        printf("  %12.6f", weight);
        int ispinor = active_spinors[j];
        printf("  %4s #%4d (%12.6f)\n", get_irrep_name(spinor_info[ispinor].repno), ispinor + 1,
               spinor_info[ispinor].eps);
    }

    /*
     * flush full information about NOs to the formatted file
     */
    double *nat_occ_re = d_zeros(n_active, 1);
    array_get_real_part(n_active, nat_occ, nat_occ_re);

    char natorb_file_name[256];
    sprintf(natorb_file_name, "NATORB_%dh%dp_%s_%d.dat", sect_h, sect_p, rep_name, state + 1);
    str_replace(natorb_file_name, '/', '|');  // for some irrep names with '/' in order to prevent bugs
    write_natural_spinors_formatted(natorb_file_name, sect_h, sect_p, n_active, nat_occ_re, nat_spinors_right, OCC_THRESH);

    /*
     * cleanup
     */
    cc_free(denmat);
    cc_free(nat_occ);
    cc_free(nat_occ_re);
    cc_free(nat_spinors_left);
    cc_free(nat_spinors_right);
}


/**
 * Calculate natural transition spinors for the pair of electronic states.
 * Write them to a formatted file.
 */
void construct_model_space_natural_transition_spinors(
        int sect_h, int sect_p,
        int irrep1, int state1, int irrep2, int state2
)
{
    double const OCC_THRESH = 1e-6;    // threshold for printing natural orbitals
    double const COEF_THRESH = 1e-4;   // threshold for printing model vec-s coeff-s
    char *rep_name1 = get_irrep_name(irrep1);
    char *rep_name2 = get_irrep_name(irrep2);

    /*
     * header banner
     */
    printf("\n");
    printf("\n");
    printf(" **\n");
    printf(" ** Natural transition spinors for the electronic transition:\n");
    printf(" ** irrep %d (%s), state %d  ->  irrep %d (%s), state %d\n",
           irrep1, rep_name1, state1 + 1, irrep2, rep_name2, state2 + 1);
    printf(" **\n");
    printf("\n");

    /*
     * only 'valence' parts of DMs are be constructed
     * dim DMs = n_active x n_active
     */
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    get_active_space(sect_h, sect_p, &n_active, active_spinors);

    /*
     * extract model vector and eigenvalue from the MVCOEF* unformatted file
     * (bra state)
     */
    int ms_size_1 = 0;
    slater_det_t *det_list_1 = NULL;
    double eigval_1 = 0.0;
    double energy_cm_1 = 0.0;
    double complex *coef_left_1 = NULL;
    double complex *coef_right_1 = NULL;

    char *mvcoef_file_name = NULL;
    if (sect_h == 1 && sect_p == 1 && irrep1 == get_vacuum_irrep()) {
        printf(" Mixed 0h0p+1h1p vectors will be used in <bra|\n");
        mvcoef_file_name = cc_strdup("MVCOEF0011");
    }

    mvcoef_read_vectors_unformatted_for_state(
            sect_h, sect_p, mvcoef_file_name, rep_name1, state1,
            &ms_size_1, &det_list_1, &eigval_1, &energy_cm_1, &coef_left_1, &coef_right_1
    );

    /*
     * extract model vector and eigenvalue from the MVCOEF* unformatted file
     * (ket state)
     */
    int ms_size_2 = 0;
    slater_det_t *det_list_2 = NULL;
    double eigval_2 = 0.0;
    double energy_cm_2 = 0.0;
    double complex *coef_left_2 = NULL;
    double complex *coef_right_2 = NULL;

    mvcoef_file_name = NULL;
    if (sect_h == 1 && sect_p == 1 && irrep2 == get_vacuum_irrep()) {
        printf(" Mixed 0h0p+1h1p vectors will be used in |ket>\n");
        mvcoef_file_name = cc_strdup("MVCOEF0011");
    }

    mvcoef_read_vectors_unformatted_for_state(
            sect_h, sect_p, mvcoef_file_name, rep_name2, state2,
            &ms_size_2, &det_list_2, &eigval_2, &energy_cm_2, &coef_left_2, &coef_right_2
    );

    /*
     * density matrix (DM) construction
     */
    printf(" Model density matrix construction...\n");

    double complex *denmat = z_zeros(n_active, n_active);
    construct_model_space_density_matrix(
            sect_h, sect_p,
            ms_size_1, coef_left_1, det_list_1,
            ms_size_2, coef_right_2, det_list_2,
            denmat, &n_active
    );

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
    double *lambda = d_zeros(n_active, 1);      // NTO singular values
    double complex *nat_spinor_from = z_zeros(n_active, n_active);
    double complex *nat_spinor_to = z_zeros(n_active, n_active);

    svd(n_active, denmat, lambda, nat_spinor_from, nat_spinor_to);

    /*
     * print occupation numbers and NTS to stdout
     */
    printf(" Symmetry: %s --> %s\n", rep_name1, rep_name2);
    printf(" State 1: %10.2f cm^-1 (%s)\n", energy_cm_1, rep_name1);
    printf(" State 2: %10.2f cm^-1 (%s)\n", energy_cm_2, rep_name2);

    double sum_lambda_sq = 0.0;  // sum of weights
    int count = 1;
    for (int i = 0; i < n_active; i++) {

        sum_lambda_sq += lambda[i] * lambda[i];
        if (lambda[i] * lambda[i] < OCC_THRESH) {
            continue;
        }

        printf(" [%d] Squared singular value (weight) = %.6f\n", count++, lambda[i] * lambda[i]);
        printf("    from:\n");
        print_nat_spinor(stdout, n_active, active_spinors, nat_spinor_from + n_active * i, COEF_THRESH);
        printf("    to:\n");
        print_nat_spinor(stdout, n_active, active_spinors, nat_spinor_to + n_active * i, COEF_THRESH);
    }
    printf(" Sum of weights = %.6f (should be in range [0,%d])\n", sum_lambda_sq, max(sect_h, sect_p));

    /*
     * flush full information about NTS's to the formatted file
     */
    char natorb_file_name[256];

    sprintf(natorb_file_name, "NTO_%dh%dp_%s_%d_%s_%d_INITIAL.dat",
            sect_h, sect_p, rep_name1, state1 + 1, rep_name2, state2 + 1);
    str_replace(natorb_file_name, '/', '|');  // for some irrep names with '/' in order to prevent bugs
    write_natural_spinors_formatted(natorb_file_name, sect_h, sect_p, n_active, lambda, nat_spinor_from, OCC_THRESH);

    sprintf(natorb_file_name, "NTO_%dh%dp_%s_%d_%s_%d_FINAL.dat",
            sect_h, sect_p, rep_name1, state1 + 1, rep_name2, state2 + 1);
    str_replace(natorb_file_name, '/', '|');
    write_natural_spinors_formatted(natorb_file_name, sect_h, sect_p, n_active, lambda, nat_spinor_to, OCC_THRESH);

    /*
     * cleanup
     */
    cc_free(denmat);
    cc_free(lambda);
    cc_free(nat_spinor_from);
    cc_free(nat_spinor_to);
}


/**
 * Write natural orbital (spinor) expansions to the formatted file with
 * the standard name NATORB_[H]h[P]p_[REP]:[STATE].dat
 */
void write_natural_spinors_formatted(
        char *natorb_file_name, int sect_h, int sect_p,
        int dim, double *occ_numbers, double complex *nat_orbs,
        double occ_thresh
)
{
    // construct 'local indices -> global indices' mapping
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS];

    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_act_hole(i) && sect_h > 0 || is_act_particle(i) && sect_p > 0) {
            active_spinors[n_active] = i;
            n_active++;
        }
    }

    // header contains brief info about spinors with numbering
    FILE *f = fopen(natorb_file_name, "w");
    fprintf(f, "dim %d\n", n_active);
    fprintf(f, "spinor info:\n");
    for (int i = 0; i < n_active; i++) {
        int ispinor = active_spinors[i];
        fprintf(f, "%4d%20.12f%10s\n", ispinor + 1, spinor_info[ispinor].eps,
                get_irrep_name(spinor_info[ispinor].repno));
    }

    // calculate number of NOs to be printed (with occ_no >= thresh)
    int n_natorb = 0;
    for (int i = 0; i < n_active; i++) {
        if (fabs(creal(occ_numbers[i])) >= occ_thresh) {
            n_natorb++;
        }
    }
    fprintf(f, "NSPINORS %d\n", n_natorb);

    // write NOs
    for (int i = 0; i < n_active; i++) {
        double occ = occ_numbers[i];
        if (fabs(occ) < occ_thresh) {
            continue;
        }

        fprintf(f, "occ %.6f\n", occ);

        // write coefficients (all)
        for (size_t j = 0; j < n_active; j++) {
            double complex coef = nat_orbs[n_active * i + j];
            double coef_re = creal(coef);
            double coef_im = cimag(coef);
            int ispinor = active_spinors[j];
            if (fabs(coef_re) < 1e-16) {
                coef_re = 0.0;
            }
            if (fabs(coef_im) < 1e-16) {
                coef_im = 0.0;
            }
            fprintf(f, "%4d  %24.12e%24.12e\n", ispinor + 1, coef_re, coef_im);
        }
    }
    printf(" written to file: %s\n", natorb_file_name);

    fclose(f);
}


/**
 * Prints natural spinor expanded in the basis of molecular spinors.
 *
 * @param out - output stream.
 * @param n_active - active space dim
 * @param active_spinors - global indices of active spinors
 * @param nat_spinor_coefs - expansion coefficients
 * @param thresh - printing threshold
 */
void print_nat_spinor(FILE *out, int n_active, int *active_spinors, double complex *nat_spinor_coefs, double thresh)
{
    for (size_t j = 0; j < n_active; j++) {
        double complex coef = nat_spinor_coefs[j];

        if (cabs(coef) < thresh) {
            continue;
        }

        int ispinor = active_spinors[j];
        int irrep = get_spinor_irrep(ispinor);
        char *irrep_name = get_irrep_name(irrep);
        double eps = get_eps(ispinor);

        fprintf(out, "  %12.6f%12.6f  %8s #%4d (%12.6f)\n",
                creal(coef), cimag(coef), irrep_name, ispinor + 1, eps);
    }
}


/**
 * Comparator function for rearrangement of natural orbitals by their eigenvalues
 * -- occupation numbers
 *
 * Rules of comparation:
 * (1) negative, then positive
 * (2) negative occ numbers in ascending order: -1 --> 0
 * (3) positive occ numbers in descending order: +1 --> 0
 */
int cmp_natorb_occupations(double complex occ1, double complex occ2)
{
    double occ1_re = creal(occ1);
    double occ2_re = creal(occ2);

    if (occ1_re * occ2_re < 0) {
        return sgn(occ1_re - occ2_re);  // ascending
    }
    else {
        return sgn(fabs(occ2_re) - fabs(occ1_re));  // descending
    }
}

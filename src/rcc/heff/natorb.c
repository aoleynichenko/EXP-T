/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
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
 * 2021 Alexander Oleynichenko
 */

#include "natorb.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "io.h"
#include "linalg.h"
#include "memory.h"
#include "mvcoef.h"
#include "slater_rules.h"
#include "spinors.h"
#include "symmetry.h"
#include "options.h"
#include "utils.h"

// functions used in this module
int cmp_natorb_occupations(double complex occ1, double complex occ2);

void write_NO(char *natorb_file_name, int sect_h, int sect_p,
              int dim, double *occ_numbers, double complex *nat_orbs,
              double occ_thresh);

int density_matrix_element(int sect_h, int sect_p, int p, int q, slater_det_t *bra, slater_det_t *ket);

void sort_vectors(int n, double complex *ev, double complex *vl, double complex *vr,
                  int (*cmp)(double complex, double complex));

void construct_ms_density_matrix(int sect_h, int sect_p,
                                 int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
                                 int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
                                 double complex *dm, int *dim_dm);

void natural_spinor_configuration(int n_active, double complex *natorb_right, double *nat_occ, double *conf);


/*
 * Constructs quasinatural orbitals (NOs) for the given set of model vectors.
 *
 * sect_h, sect_p   Fock space sector
 * ms_size          model space dimension
 * det_list         list of basis Slater determinants (dim = ms_size)
 * coef_left        left model vectors (dim = ms_size)
 * coef_right       right model vectors (dim = ms_size)
 * natorb_left      left NOs (left eigenvectors of the model density matrix),
 *                  matrix of dim = n_active x n_active
 * natorb_right     right NOs, matrix of dim = n_active x n_active
 * nat_occ_numbers  list of occupation numbers of each NO (dim = n_active)
 * config           effective configuration of active spinors (dim = n_active)
 */
void construct_quasi_natural_orbitals(int sect_h, int sect_p, int ms_size, slater_det_t *det_list,
                                      double complex *coef_left, double complex *coef_right,
                                      double complex *natorb_left, double complex *natorb_right,
                                      double *nat_occ_numbers, double *config)
{
    double complex *denmat;
    double complex dm_eig_values[CC_MAX_SPINORS];

    int n_active = 0;
    int act_spinors_indices[CC_MAX_SPINORS];
    get_active_space(sect_h, sect_p, &n_active, act_spinors_indices);

    // density matrix (DM) construction
    denmat = z_zeros(n_active, n_active);
    construct_ms_density_matrix(sect_h, sect_p, ms_size, coef_left, det_list, ms_size, coef_right, det_list, denmat,
                                &n_active);

    // calculate natural spinors (NOs) = diagonalize density matrix
    eigensolver(n_active, denmat, dm_eig_values, natorb_left, natorb_right);
    sort_vectors(n_active, dm_eig_values, natorb_left, natorb_right, cmp_natorb_occupations);

    // check for imaginary parts of occupation numbers
    // and convert occ numbers to reals
    for (size_t i = 0; i < n_active; i++) {
        if (cimag(dm_eig_values[i]) > 1e-13) {
            printf(" Imaginary occupation number: %.3e %.3e\n", creal(dm_eig_values[i]), cimag(dm_eig_values[i]));
        }
    }
    get_real_parts(n_active, dm_eig_values, nat_occ_numbers);

    // effective configuration
    natural_spinor_configuration(n_active, natorb_right, nat_occ_numbers, config);

    cc_free(denmat);
}


void quasi_natural_orbitals_driver(int sect_h, int sect_p, int rep, int state)
{
    double const occ_thresh = 1e-6;   // threshold for printing natural orbitals
    double const coef_thresh = 1e-4;   // threshold for printing model vec-s coeff-s

    int nrep;
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    struct mv_block mv_blocks_0011[CC_MAX_NUM_IRREPS];

    char *rep_name = get_irrep_name(rep);
    printf("\n\n *** DENSITY MATRIX AND NATURAL ORBITALS FOR THE IRREP %d (%s) STATE %d ***\n",
           rep + 1, rep_name, state + 1);

    // only 'valence' parts of DMs are be constructed
    // dim DMs = n_active x n_active
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    get_active_space(sect_h, sect_p, &n_active, active_spinors);

    // allocate working arrays
    double *nat_occ = d_zeros(n_active, 1);  // NO occ numbers
    double *conf = d_zeros(n_active, 1);
    double complex *natorb_left = z_zeros(n_active, n_active);
    double complex *natorb_right = z_zeros(n_active, n_active);

    // extract model vectors and eigenvalues from the MVCOEF* unformatted file
    read_model_vectors_unformatted(sect_h, sect_p, NULL, &nrep, mv_blocks);
    struct mv_block *modvec_block = NULL;

    for (size_t ib = 0; ib < nrep; ib++) {
        if (strcmp(mv_blocks[ib].rep_name, rep_name) == 0) {
            if (sect_h == 1 && sect_p == 1 && rep == get_vacuum_irrep()) {
                // vectors with the intermediate normalization restored (0h0p and 1h1p are mixed)
                // to be used for 00->11, 11->00 and 11->11 matrix elements by default
                int nrep_0011;
                printf(" Mixed 0h0p+1h1p vectors will be used\n");
                read_model_vectors_unformatted(1, 1, "MVCOEF0011", &nrep_0011, mv_blocks_0011);
                modvec_block = &mv_blocks_0011[0];
            }
            else {
                modvec_block = mv_blocks + ib;
            }
        }
    }
    size_t ms_size = modvec_block->ms_size;
    slater_det_t *det_list = modvec_block->dets;
    double eigvalue = creal(modvec_block->eigval[state]);
    double energy_cm = modvec_block->energy_cm[state];
    double complex *coef_left = modvec_block->vl + ms_size * state;
    double complex *coef_right = modvec_block->vr + ms_size * state;
    printf(" Symmetry of this electronic state: %s\n", rep_name);
    printf(" Eigenvalue   = %.8f a.u.\n", eigvalue);
    printf(" Energy level = %.2f cm^-1\n", energy_cm);

    /*
     * Construct (model) density matrix and diagonalize it
     */
    construct_quasi_natural_orbitals(
            sect_h, sect_p, ms_size, det_list, coef_left, coef_right,
            natorb_left, natorb_right, nat_occ, conf
    );

    /*
     * Print natural orbitals to stdout
     */
    int count = 1;
    double sum_occ = sum_doubles(n_active, nat_occ);  // trace
    for (int i = 0; i < n_active; i++) {
        if (fabs(nat_occ[i]) < occ_thresh) {
            continue;
        }

        printf(" [%d] Occ number = %.6f\n", count++, nat_occ[i]);
        for (size_t j = 0; j < n_active; j++) {
            double complex coef = natorb_right[n_active * i + j];
            if (cabs(coef) < coef_thresh) {
                continue;
            }
            printf("  %12.6f%12.6f ", creal(coef), cimag(coef));
            int ispinor = active_spinors[j];
            printf("  %4s #%4d (%12.6f)\n", get_irrep_name(spinor_info[ispinor].repno), ispinor + 1,
                   spinor_info[ispinor].eps);
        }
    }
    printf(" Sum of occ numbers = %.6f (must be = %d)\n", sum_occ, -sect_h + sect_p);
    printf(" Configuration (weights):\n");
    for (size_t j = 0; j < n_active; j++) {
        printf("  %12.6f", conf[j]);
        int ispinor = active_spinors[j];
        printf("  %4s #%4d (%12.6f)\n", get_irrep_name(spinor_info[ispinor].repno), ispinor + 1,
               spinor_info[ispinor].eps);
    }

    /*
     * flush full information about NOs to the formatted file
     */
    char natorb_file_name[256];
    sprintf(natorb_file_name, "NATORB_%dh%dp_%s_%d.dat", sect_h, sect_p, rep_name, state + 1);
    str_replace(natorb_file_name, '/', '|');  // for some irrep names with '/' in order to prevent bugs
    write_NO(natorb_file_name, sect_h, sect_p, n_active, nat_occ, natorb_right, occ_thresh);

    // cleanup
    cc_free(nat_occ);
    cc_free(conf);
    cc_free(natorb_left);
    cc_free(natorb_right);
    for (size_t irep = 0; irep < nrep; irep++) {
        struct mv_block *b = &mv_blocks[irep];
        cc_free(b->dets);
        cc_free(b->eigval);
        cc_free(b->energy_cm);
        cc_free(b->vl);
        cc_free(b->vr);
    }
    printf(" *** end of DENSITY MATRIX AND NATURAL ORBITALS module ***\n");
}


void quasi_natural_transition_orbitals(int sect_h, int sect_p, int rep1, int state1, int rep2, int state2)
{

}


void natural_spinor_configuration(int n_active, double complex *natorb_right, double *nat_occ, double *conf)
{
    for (size_t j = 0; j < n_active; j++) {
        double weight = 0.0;
        for (size_t i = 0; i < n_active; i++) {
            double complex coef = natorb_right[n_active * i + j];
            weight += nat_occ[i] * cabs(coef) * cabs(coef);
        }
        conf[j] = weight;
    }
}


/**
 * Write natural orbital (spinor) expansions to the formatted file with
 * the standard name NATORB_[H]h[P]p_[REP]:[STATE].dat
 */
void write_NO(char *natorb_file_name, int sect_h, int sect_p,
              int dim, double *occ_numbers, double complex *nat_orbs,
              double occ_thresh)
{
    FILE *f;
    int n_natorb;
    int n_active;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping

    // construct 'local indices -> global indices' mapping
    n_active = 0;
    for (int i = 0; i < get_num_spinors(); i++) {
        if (is_act_hole(i) && sect_h > 0 || is_act_particle(i) && sect_p > 0) {
            active_spinors[n_active] = i;
            n_active++;
        }
    }

    // header contains brief info about spinors with numbering
    f = fopen(natorb_file_name, "w");
    fprintf(f, "dim %d\n", n_active);
    fprintf(f, "spinor info:\n");
    for (int i = 0; i < n_active; i++) {
        int ispinor = active_spinors[i];
        fprintf(f, "%4d%20.12f%10s\n", ispinor + 1, spinor_info[ispinor].eps, get_irrep_name(spinor_info[ispinor].repno));
    }

    // calculate number of NOs to be printed (with occ_no >= thresh)
    n_natorb = 0;
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

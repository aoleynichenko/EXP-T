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
 * (1) construction and diagonalization of the effective Hamiltonian matrix;
 * (2) analysis of its eigenvalues and eigenvectors (model vectors);
 * (3) write matrix and vectors to disk.
 * NOTE: ALL THIS CODE IS DESIGNED FOR ABELIAN GROUPS ONLY!
 *
 * 2019-2021 Alexander Oleynichenko
 */

#include "heff.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cc_properies.h"
#include "codata.h"
#include "io.h"
#include "engine.h"
#include "datamodel.h"
#include "linalg.h"
#include "model_space.h"
#include "mvcoef.h"
#include "natorb.h"
#include "options.h"
#include "renorm_omega.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void construct_heff(int sect_h, int sect_p, slater_det_t **det_list, size_t *block_sizes, double complex **heff_blocks,
                    int n_diagrams, char **diagram_names);

void zero_order_heff(int sect_h, int sect_p, size_t dim, slater_det_t *det_list, double complex *heff);

void diagonalize_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff,
                      double complex **eigvalues, double complex **coef_left, double complex **coef_right);

void print_model_vectors_stdout(int sector_h, int sector_p, slater_det_t **det_list, size_t *block_dims,
                                double complex **eigvalues, double complex **coef_left, double complex **coef_right);

void write_model_vectors_unformatted(int sector_h, int sector_p, slater_det_t **det_list, size_t *block_dims,
                                     double complex **eigvalues, double complex **coef_left,
                                     double complex **coef_right);

void print_eigenvalues_table(int sect_h, int sect_p, size_t *block_dims, slater_det_t **det_list,
                             double complex **eigvalues, double degen_thresh, double complex **coef_right);

int get_nroots_for_irrep(char *irrep_name);

void write_formatted_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff);

void write_formatted_heff_0h0p(double total_energy);

size_t first_nonzero_irrep(size_t *block_dims);

void get_nroots(size_t *block_dims, double complex **eigvalues, size_t *nroots);

size_t get_nroots_from_cutoff(size_t irrep_dim, double complex *eigvalues, double emin, double cutoff);

double get_fraction_of_main_space_determinants(
        int sector_h,
        int sector_p,
        size_t dim,
        slater_det_t *det_list,
        double complex *model_vector
);


/**
 * heff_analysis
 *
 * Constructs and diagonalizes effective Hamiltonian.
 * Performs analysis of model vectors.
 * Writes to disk:
 *   HEFF      formatted file with the effective Hamiltonian
 *   MVCOEF**  unformatted files with model vectors and eigenvalues
 *
 * Agruments:
 *   sect_h, sect_p   number of the Fock-space sector
 *   ...              list of names of the diagrams containing 1-particle,
 *                    2-particle etc parts of the effective interaction operator
 */
void heff_analysis(int sect_h, int sect_p, ...)
{
    double complex *heff[CC_MAX_NUM_IRREPS];
    double complex *eigvalues[CC_MAX_NUM_IRREPS];
    double complex *coef_left[CC_MAX_NUM_IRREPS];
    double complex *coef_right[CC_MAX_NUM_IRREPS];
    size_t block_dims[CC_MAX_NUM_IRREPS];
    const int MAX_DIAGRAMS = 64;
    char *diagram_names[MAX_DIAGRAMS];

    /*
     * Construct model space
     */
    slater_det_t **det_basis = construct_model_space(sect_h, sect_p, block_dims);

    if (cc_opts->print_det_list) {
        printf(" List of model space determinants:\n");
        for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
            for (int i = 0; i < block_dims[irrep]; i++) {
                print_slater_det(stdout, sect_h, sect_p, det_basis[irrep] + i);
            }
        }
        printf("\n");
    }

    printf(" Model space dimension:\n");
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] != 0) {
            printf("  %4s [%d]", get_irrep_name(irrep), block_dims[irrep]);
        }
    }
    printf("\n");

    va_list args;
    va_start(args, sect_p);
    int n_diagrams = (sect_h + 1) * (sect_p + 1) - 1; // rectangle - (0,0) sector
    for (int i = 0; i < n_diagrams; i++) {
        char *dg_name = va_arg(args, char *);
        diagram_names[i] = cc_strdup(dg_name);
    }
    va_end(args);
    construct_heff(sect_h, sect_p, det_basis, block_dims, heff, n_diagrams, diagram_names);
    for (int i = 0; i < n_diagrams; i++) {
        cc_free(diagram_names[i]);
    }

    /*
     * Write Heff to the formatted file
     */
    write_formatted_heff(sect_h, sect_p, block_dims, heff);

    /*
     * Diagonalization
     */
    diagonalize_heff(sect_h, sect_p, block_dims, heff, eigvalues, coef_left, coef_right);

    /*
     * For the 1h1p sector only: restoration of the intermediate normalization
     */
    if (sect_h == 1 && sect_p == 1) {
        restore_intermediate_normalization(block_dims, det_basis, heff, eigvalues);
    }

    /*
     * Save model vectors and eigenvalues to the formatted files MVCOEF
     */
    write_model_vectors_unformatted(sect_h, sect_p, det_basis, block_dims,
                                    eigvalues, coef_left, coef_right);

    /*
     * Print model vectors
     */
    print_model_vectors_stdout(sect_h, sect_p, det_basis, block_dims, eigvalues, coef_left, coef_right);

    /*
     * Print table with energy levels
     */
    print_eigenvalues_table(sect_h, sect_p, block_dims, det_basis, eigvalues, cc_opts->degen_thresh, coef_right);

    /*
     * Cleanup
     */
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        cc_free(heff[irrep]);
        cc_free(eigvalues[irrep]);
        cc_free(coef_left[irrep]);
        cc_free(coef_right[irrep]);
        cc_free(det_basis[irrep]);
    }
    cc_free(det_basis);
}


void construct_heff(int sect_h, int sect_p, slater_det_t **det_list, size_t *block_sizes, double complex **heff_blocks,
                    int n_diagrams, char **diagram_names)
{
    /*
     * Allocate memory for the effective Hamiltonian (Heff)
     */
    for (int irep = 0; irep < CC_MAX_NUM_IRREPS; irep++) {
        heff_blocks[irep] = NULL;
    }
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        size_t dim = block_sizes[irep];
        if (dim != 0) {
            heff_blocks[irep] = z_zeros(dim, dim);
        }
    }

    /*
     * Construct symmetry blocks of the Heff matrix
     */
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_sizes[irep] == 0) {
            continue;
        }

        // get pointer to the current block of det-s with symmetry 'irep'
        slater_det_t *rep_dets = det_list[irep];

        // size of the Heff subblock
        size_t ms_size = block_sizes[irep];

        // construct zero-order Hamiltonian
        zero_order_heff(sect_h, sect_p, ms_size, rep_dets, heff_blocks[irep]);

        // loop over n-particle contibutions (n=1,2,...) to the Veff operator
        for (int ioper = 0; ioper < n_diagrams; ioper++) {
            char *dg_name = diagram_names[ioper];
            diagram_t *dg_veff = diagram_stack_find(dg_name);
            if (dg_veff == NULL) {
                errquit("in construct_heff(): diagram '%s' not found", dg_name);
            }
            setup_slater(dg_veff, (matrix_getter_fun) diagram_get, sect_h, sect_p, sect_h, sect_p, rank(dg_name) / 2);

            // construct matrix of the effective interaction (Veff) in the basis
            // of Slater determinants
            for (int i = 0; i < ms_size; i++) {
                for (int j = 0; j < ms_size; j++) {
                    slater_det_t *bra = rep_dets + i;
                    slater_det_t *ket = rep_dets + j;
                    heff_blocks[irep][i * ms_size + j] += slater(bra, ket);
                }
            }
        }
    }
}


/**
 * construct zero-order effective Hamiltonian (H0) matrix:
 *  - only the block belonging to one irrep;
 *  - only one-electron spinor energies are used.
 * It is a diagonal matrix: heff[i,i] = sum{particles} eps - sum{holes} eps
 *
 * @param sect_h, sect_p  Fock space sector signature
 * @param dim dimension of this Heff block (number of model det-s of one symmetry)
 * @param det_list list of model space determinants (with one symmetry)
 * @param heff Heff matrix, of size dim x dim
 */
void zero_order_heff(int sect_h, int sect_p, size_t dim, slater_det_t *det_list, double complex *heff)
{
    for (size_t i = 0; i < dim; i++) {
        slater_det_t *d = det_list + i;

        double H0 = 0.0;
        for (int h = 0; h < sect_h; h++) {
            H0 -= spinor_info[d->indices[h]].eps;
        }
        for (int p = 0; p < sect_p; p++) {
            H0 += spinor_info[d->indices[sect_h + p]].eps;
        }
        heff[i * dim + i] = H0 + 0.0 * I;
    }
}


void diagonalize_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff,
                      double complex **eigvalues, double complex **coef_left, double complex **coef_right)
{
    /*
     * Allocate memory for left, right eigenvectors and eigenvalues
     */
    for (int irep = 0; irep < CC_MAX_NUM_IRREPS; irep++) {
        eigvalues[irep] = NULL;
        coef_left[irep] = NULL;
        coef_right[irep] = NULL;
    }
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        size_t dim = block_dims[irep];
        if (dim != 0) {
            eigvalues[irep] = z_zeros(dim, 1);
            coef_left[irep] = z_zeros(dim, dim);
            coef_right[irep] = z_zeros(dim, dim);
        }
    }

    /*
     * Diagonalize each block of the Heff matrix
     * 1. eigenvalues are sorted in ascending order
     * 2. eigenvectors are resorted too (and biorthonormalized)
     */
    size_t max_dim = size_t_max(get_num_irreps(), block_dims);
    double complex *heff_work = z_zeros(max_dim, max_dim);
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_dims[irep] == 0) {
            continue;
        }
        size_t ms_size = block_dims[irep];
        double complex *heff_block = heff[irep];
        double complex *eigvalues_block = eigvalues[irep];
        double complex *coef_left_block = coef_left[irep];
        double complex *coef_right_block = coef_right[irep];

        memcpy(heff_work, heff_block, sizeof(double complex) * ms_size * ms_size);
        eigensolver(ms_size, heff_work, eigvalues_block, coef_left_block, coef_right_block);
    }
    cc_free(heff_work);

    /*
     * Loewdin orthogonalization of model vectors (if required)
     * NOTE: if no orthogonalization is performed, property matrices will be non-Hermitian
     */
    if (cc_opts->do_hermit == 1) {
        for (int irep = 0; irep < get_num_irreps(); irep++) {
            if (block_dims[irep] == 0) {
                continue;
            }
            size_t ms_size = block_dims[irep];
            double complex *coef_left_block = coef_left[irep];
            double complex *coef_right_block = coef_right[irep];

            loewdin_orth(ms_size, coef_right_block, coef_right_block, cc_opts->print_level == CC_PRINT_DEBUG ? 5 : 0);
            memcpy(coef_left_block, coef_right_block, ms_size * ms_size * sizeof(double complex));
        }
    }
}


double max_energy_of_required_roots(size_t *block_dims, double complex **eigvalues)
{
    double max_energy = -1.0e12;

    size_t nroots[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        if (block_dims[irrep] == 0 || nroots[irrep] == 0) {
            continue;
        }

        // get max number of roots for this irrep
        /*char *irrep_name = get_irrep_name(irrep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : block_dims[irrep];*/

        double complex *ev = eigvalues[irrep];
        if (nroots[irrep] > 0) {
            for (size_t i = 0; i < nroots[irrep]; i++) {
                if (creal(ev[i]) > max_energy) {
                    max_energy = creal(ev[i]);
                }
            }
        }
    }

    return max_energy;
}


/**
 * calculates effective occupation numbers of active spinors
 * for the given electronic state (represented by a model vector)
 */
void get_eff_configuration(
        int sect_h,
        int sect_p,
        int ms_size,
        slater_det_t *det_list,
        double complex *model_vector_coefs,
        double *config
)
{
    // construct list of indices of active spinors
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS]; // local -> global spinor index mapping
    get_active_space(sect_h, sect_p, &n_active, active_spinors);

    memset(config, 0, n_active * sizeof(double));

    // loop over model determinants
    for (int i = 0; i < ms_size; i++) {
        double det_coefft = model_vector_coefs[i];
        double det_weight = cabs(det_coefft) * cabs(det_coefft);
        slater_det_t *det = det_list + i;

        // for each determinant:
        // determine which active spinors are occupied
        for (int j = 0; j < sect_h + sect_p; j++) {
            int act_spinor_index = det->indices[j];
            for (int k = 0; k < n_active; k++) {
                if (active_spinors[k] == act_spinor_index) {
                    config[k] += det_weight;
                    break;
                }
            }
        }
    }
}


/*
 * Prints model vectors to stdout (formatted output).
 * Also prints effective configurations of active spinors calculated via the
 * quasi-natural spinors approach.
 *
 * If the 'nroots' directive is given in the input file, only the
 * vectors specified by this directive will be printed and analyzed.
 */
void print_model_vectors_stdout(int sector_h, int sector_p, slater_det_t **det_list, size_t *block_dims,
                                double complex **eigvalues, double complex **coef_left, double complex **coef_right)
{
    double const COEF_THRESH = 1e-4;   // threshold for printing model vec-s coeff-s
    int first_irrep = first_nonzero_irrep(block_dims);
    double eff_config[CC_MAX_SPINORS];

    double overall_eff_conf[CC_MAX_SPINORS];
    memset(overall_eff_conf, 0, sizeof(overall_eff_conf));

    if (cc_opts->print_model_vectors == 0 &&
        cc_opts->print_eff_config == 0) {
        return;
    }

    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS];
    get_active_space(sector_h, sector_p, &n_active, active_spinors);

    size_t nroots[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots);

    printf("\n Sector (%dh,%dp) -- analysis of model vectors (right vectors)\n", sector_h, sector_p);
    printf(" first line : irrep, state number, total energy, eigenvalue\n");
    printf(" other lines: coefficients of contributing determinants (above a threshold of %.1e)\n",
           COEF_THRESH);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        size_t dim = block_dims[irrep];
        size_t nroots_irrep = nroots[irrep];
        slater_det_t *det_basis = det_list[irrep];

        if (dim == 0 || nroots_irrep == 0) {
            continue;
        }

        // get max number of roots for this irrep
        char *irrep_name = get_irrep_name(irrep);

        for (size_t i = 0; i < nroots_irrep; i++) {
            double complex *left_vector = coef_left[irrep] + dim * i;
            double complex *right_vector = coef_right[irrep] + dim * i;

            print_model_vector(stdout, sector_h, sector_p,
                               irrep - first_irrep + 1, irrep_name, i, eigvalues[irrep][i],
                               dim, coef_right[irrep] + dim * i, det_basis, COEF_THRESH);

            /*
             * calculate effective configurations for the target electronic states
             */
            if (cc_opts->print_eff_config) {

                get_eff_configuration(sector_h, sector_p, dim, det_basis, right_vector, eff_config);

                printf(" effective configuration:\n");
                for (size_t j = 0; j < n_active; j++) {

                    int ispinor = active_spinors[j];

                    // accumulate eff configuration over all electronic states
                    overall_eff_conf[ispinor] += eff_config[j];

                    if (fabs(eff_config[j]) < 1e-2) {
                        continue;
                    }

                    printf("  %12.6f", eff_config[j]);
                    printf("  %4s #%4d (%12.6f)", get_irrep_name(spinor_info[ispinor].repno), ispinor + 1,
                           spinor_info[ispinor].eps);
                    if (cc_opts->spinor_labels[ispinor] != NULL) {
                        printf("  %s", cc_opts->spinor_labels[ispinor]);
                    }
                    printf("\n");
                }
            }

            /*
             * In case of simple Intermediate Hamiltonian:
             * calculate fraction (%) of main space determinants
             */
            if (cc_opts->intham_imms_opts.sectors[sector_h][sector_p]) {
                double fraction_main = get_fraction_of_main_space_determinants(sector_h, sector_p, dim, det_basis,
                                                                               right_vector);
                printf(" fraction of main space determinants: %.1f%%\n", fraction_main * 100);
            }
        }
    }


    if (cc_opts->print_eff_config) {
        printf("\n\n");
        printf(" sum of effective configurations over the electronic states considered\n");
        printf(" (this output can be used to choose a main model space in IH-IMMS calculations)\n");
        printf("\n");

        for (int i = 0; i < get_num_spinors(); i++) {
            if (!is_active(i)) {
                continue;
            }

            printf("  %12.6f", overall_eff_conf[i]);
            printf("  %4s #%4d (%12.6f)", get_irrep_name(spinor_info[i].repno), i + 1, spinor_info[i].eps);
            if (cc_opts->spinor_labels[i] != NULL) {
                printf("  %s", cc_opts->spinor_labels[i]);
            }
            printf("\n");
        }
    }


    printf("\n\n");
}


/**
 * calculates fraction of main space determinants in the given model vector.
 * to be used with the simple Intermediate Hamiltonian technique.
 */
double get_fraction_of_main_space_determinants(
        int sector_h,
        int sector_p,
        size_t dim,
        slater_det_t *det_list,
        double complex *model_vector
)
{
    int intham_imms_is_main_space_det(int sect_h, int sect_p, moindex_t *indices);

    double fraction_main = 0.0;

    for (size_t j = 0; j < dim; j++) {
        if (intham_imms_is_main_space_det(sector_h, sector_p, det_list[j].indices)) {
            double complex coef = model_vector[j];
            fraction_main += cabs(coef) * cabs(coef);
        }
    }

    return fraction_main;
}


double get_lowest_eigenvalue(size_t *block_dims, double complex **eigvalues)
{
    double min_eig_value = 1e9;

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        if (creal(eigvalues[irrep][0]) < min_eig_value) {
            min_eig_value = creal(eigvalues[irrep][0]);
        }
    }

    return min_eig_value;
}


void write_model_vectors_unformatted(int sector_h, int sector_p, slater_det_t **det_list, size_t *block_dims,
                                     double complex **eigvalues, double complex **coef_left,
                                     double complex **coef_right)
{
    // open unformatted file with model vectors
    int f_mvcoef = mvcoef_open(sector_h, sector_p);

    double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);

    size_t nroots[CC_MAX_NUM_IRREPS];
    get_nroots(block_dims, eigvalues, nroots);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        size_t dim = block_dims[irrep];
        size_t nroots_irrep = nroots[irrep];

        if (dim == 0 || nroots_irrep == 0) {
            continue;
        }

        slater_det_t *det_basis = det_list[irrep];
        char *irrep_name = get_irrep_name(irrep);


        // get max number of roots for this irrep
        /*char *irrep_name = get_irrep_name(irrep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : block_dims[irrep];
        if (cc_opts->roots_cutoff_specified) {
            //size_t nroots_from_cutoff = get_nroots_from_cutoff(dim, eigvalues[irrep], lowest_root, cc_opts->roots_cutoff);
        }
        if (nroots == 0) {
            continue;
        }*/

        // print model vectors to the unformatted file MVCOEF**
        mvcoef_write_vectors_unformatted(f_mvcoef, irrep_name, nroots_irrep, dim, det_basis,
                                         eigvalues[irrep], coef_left[irrep], coef_right[irrep]);
    }

    mvcoef_close(f_mvcoef, lowest_root);
}


void print_eigenvalues_table(int sect_h, int sect_p, size_t *block_dims, slater_det_t **det_list,
                             double complex **eigvalues, double degen_thresh, double complex **coef_right)
{
    /*
     * Allocate table for eigenvalues
     */
    size_t total_heff_dim = 0;
    for (int i = 0; i < get_num_irreps(); i++) {
        total_heff_dim += block_dims[i];
    }
    if (sect_h == 1 && sect_p == 1) {
        total_heff_dim += 1;  // for the reference energy
    }
    eigval_t *eigenvalues = (eigval_t *) cc_malloc(total_heff_dim * sizeof(eigval_t));
    size_t n_eigenvalues = 0;

    /*
     * Copy eigenvalues to the table.
     * Vacuum det is always 'main'.
     */
    if (sect_h == 1 && sect_p == 1) {
        eigenvalues[0].eigval = 0.0;
        eigenvalues[0].repno = get_vacuum_irrep();
        eigenvalues[0].percent_main = 100.0;
        n_eigenvalues++;
    }
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t block_dim = block_dims[irrep];
        for (size_t i = 0; i < block_dim; i++) {
            eigenvalues[n_eigenvalues + i].eigval = eigvalues[irrep][i];
            eigenvalues[n_eigenvalues + i].repno = irrep;

            /*
             * In case of simple Intermediate Hamiltonian:
             * calculate fraction (%) of main space determinants
             */
            double complex *right_vector = coef_right[irrep] + block_dim * i;
            if (cc_opts->intham_imms_opts.sectors[sect_h][sect_p]) {
                double fraction_main = get_fraction_of_main_space_determinants(sect_h, sect_p, block_dim, det_list[irrep],
                                                                               right_vector);
                eigenvalues[n_eigenvalues + i].percent_main = fraction_main * 100;
            }
            else {
                eigenvalues[n_eigenvalues + i].percent_main = 0.0;
            }
        }
        n_eigenvalues += block_dim;
    }

    /*
     * Sort table in ascending order
     */
    qsort(eigenvalues, n_eigenvalues, sizeof(eigval_t), eigval_cmp);

    /*
     * Print table with energy levels
     */
    double complex e0 = eigenvalues[0].eigval;
    double reference_energy = cc_opts->eref;   // (see options.h)
    double max_energy = max_energy_of_required_roots(block_dims, eigvalues);
    if (sect_h == 0 && sect_p == 1) {
        cc_opts->ground_energy_0h1p = creal(e0);
    }
    if (sect_h == 0 && sect_p == 2) {
        cc_opts->ground_energy_0h2p = creal(e0);
    }
    if (sect_h == 0 && sect_p == 3) {
        cc_opts->ground_energy_0h3p = creal(e0);
    }
    if (sect_h == 1 && sect_p == 0) {
        cc_opts->ground_energy_1h0p = creal(e0);
    }
    if (sect_h == 2 && sect_p == 0) {
        cc_opts->ground_energy_2h0p = creal(e0);
    }
    if (sect_h == 1 && sect_p == 2) {
        e0 = cc_opts->ground_energy_0h1p + 0.0 * I;
    }

    int intham = 0;
    if (cc_opts->intham_imms_opts.sectors[sect_h][sect_p]) {
        intham = 1;
    }

    printf("\n Heff eigenvalues:\n (degeneracy threshold = %.1e a.u.)\n\n", degen_thresh);

    if (intham == 0) {
        printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  deg  symmetry\n");
        printf(" ------------------------------------------------------------------------------------------------------------------------\n");
    }
    else {
        printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  %% main  deg  symmetry\n");
        printf(" --------------------------------------------------------------------------------------------------------------------------------\n");
    }

    int ilevel = 1;
    for (size_t i = 0; i < n_eigenvalues; i++) {

        // only energy levels with eigenvalue Ei <= max_energy will be printed
        if (creal(eigenvalues[i].eigval) > max_energy + degen_thresh) {
            break;
        }

        // calculate degeneracy
        int deg = 0;
        for (size_t j = i; j < n_eigenvalues; j++) {
            if (cabs(eigenvalues[i].eigval - eigenvalues[j].eigval) > degen_thresh) {
                break;
            }
            else {
                deg++;
            }
        }

        // print info about (maybe degenerate) level
        double complex ei = eigenvalues[i].eigval;
        double abs_energy = reference_energy + creal(ei);
        double rel_energy = creal(ei) - creal(e0);
        double rel_energy_cm = rel_energy * CODATA_AU_TO_CM;
        double rel_energy_ev = rel_energy * CODATA_AU_TO_EV;

        printf("@%5d%16.10f%10.2e%25.17f%16.10f%16.10f%16.6f  ",
               ilevel, creal(ei), cimag(ei), abs_energy, rel_energy,
               rel_energy_ev, rel_energy_cm);

        // percent of main model space determinants (in case of IH-FSCC)
        if (intham) {
            printf("%6.1f  ", eigenvalues[i].percent_main);
        }

        printf("%2d  ", deg);
        for (size_t irrep = 0; irrep < get_num_irreps(); irrep++) {
            int rep_deg = 0;
            for (size_t j = 0; j < deg; j++) {
                if (eigenvalues[i + j].repno == irrep) {
                    rep_deg++;
                }
            }
            if (rep_deg == 0) { continue; }
            printf(" %s", get_irrep_name(irrep));
            if (rep_deg > 1) {
                printf("(%d)", rep_deg);
            }
        }

        printf("\n");

        // go to the next (maybe degenerate) energy level
        i += deg - 1;
        ilevel++;
    }

    printf("\n");

    cc_free(eigenvalues);
}


void model_space_properties_and_natural_orbitals(int sect_h, int sect_p)
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
            // single code for NOs and NTOs:
            // rep1 == rep2 && state1 == state2  =>  single state (NO)
            // else: pair of states (NTO)
            density_matrix(sect_h, sect_p, rep1, state1, rep2, state2);
        }
    }

    // transition dipole moments via the DL-TDM techniques
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->calc_model_space_tdms) {
        dipole_length_tdms(sect_h, sect_p);
    }

    // estimate properties: only for the target sector
    if (cc_opts->sector_h == sect_h &&
        cc_opts->sector_p == sect_p &&
        cc_opts->n_model_space_props > 0) {
        for (int i = 0; i < cc_opts->n_model_space_props; i++) {
            model_space_property(cc_opts->prop_queries + i, CC_PROPERTIES_APPROX_MODEL_SPACE);
            model_space_property(cc_opts->prop_queries + i, CC_PROPERTIES_APPROX_LINEAR);
        }
    }
}


/**
 * Comparator for eigenvalues (for qsort)
 */
int eigval_cmp(const void *aa, const void *bb)
{
    const eigval_t *a = aa, *b = bb;
    return (creal(a->eigval) < creal(b->eigval)) ? -1 : (creal(a->eigval) > creal(b->eigval));
}


/**
 * Prints model vectors.
 * coef_thresh -- threshold for printing model vec-s coeff-s
 */
void print_model_vector(FILE *f_out, int sect_h, int sect_p,
                        int rep_no, char *rep_name, int state_no, double complex energy,
                        size_t len, double complex *coeffs, slater_det_t *det_list, double coef_thresh)
{
    fprintf(f_out, "\n");
    fprintf(f_out, " Irrep %d (%4s) State %d Energy %23.15f Eigenvalue %14.8f%14.4E\n",
            rep_no, rep_name, state_no + 1, cc_opts->eref + creal(energy), creal(energy), cimag(energy));

    if (cc_opts->print_model_vectors == 0) {
        return;
    }

    for (size_t j = 0; j < len; j++) {
        double complex coef = coeffs[j];
        if (cabs(coef) < coef_thresh) {
            continue;
        }
        fprintf(f_out, " %10.5f%10.5f ", creal(coef), cimag(coef));
        print_slater_det(f_out, sect_h, sect_p, det_list + j);
    }
}


/* calculates how many roots in each irrep should be processed.
 * number of roots is given by the two options in the input file:
 * 'nroots' - number of lowest roots for each irrep
 * 'roots_cutoff' - lowest roots are limited by some max energy
 * finally, number of roots = min(nroots,roots_cutoff)
 */
void get_nroots(size_t *block_dims, double complex **eigvalues, size_t *nroots)
{
    double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        size_t nroots_irrep = block_dims[irrep];
        size_t nroots_cutoff = block_dims[irrep];

        // from the 'nroots' option
        if (cc_opts->nroots_specified) {
            char *irrep_name = get_irrep_name(irrep);
            nroots_irrep = get_nroots_for_irrep(irrep_name);
        }
        else if (cc_opts->roots_cutoff_specified) {
            nroots_cutoff = get_nroots_from_cutoff(nroots_irrep, eigvalues[irrep], lowest_root, cc_opts->roots_cutoff);
        }

        nroots[irrep] = MIN(nroots_irrep, nroots_cutoff);
    }
}


size_t get_nroots_from_cutoff(size_t irrep_dim, double complex *eigvalues, double emin, double cutoff)
{
    size_t nroots = 0;

    for (size_t i = 0; i < irrep_dim; i++) {
        double e_i = creal(eigvalues[i]);
        if (e_i - emin <= cutoff) {
            nroots += 1;
        }
    }

    return nroots;
}


int get_nroots_for_irrep(char *irrep_name)
{
    int nroots_irep = 0;

    for (int ii = 0; ii < get_num_irreps(); ii++) {
        if (strcmp(irrep_name, cc_opts->nroots_specs.rep_names[ii]) == 0) {
            nroots_irep = cc_opts->nroots_specs.dim[ii];
            break;
        }
    }

    return nroots_irep;
}


size_t first_nonzero_irrep(size_t *block_dims)
{
    for (size_t i = 0; i < get_num_irreps(); i++) {
        if (block_dims[i] != 0) {
            return i;
        }
    }

    return 0;
}



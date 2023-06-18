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
 * (1) construction and diagonalization of the effective Hamiltonian matrix;
 * (2) analysis of its eigenvalues and eigenvectors (model vectors);
 * (3) write matrix and vectors to disk.
 * NOTE: ALL THIS CODE IS DESIGNED FOR ABELIAN GROUPS ONLY!
 */

#include "heff.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cc_properies.h"
#include "engine.h"
#include "eigenvalues.h"
#include "datamodel.h"
#include "linalg.h"
#include "model_vectors.h"
#include "model_space.h"
#include "options.h"
#include "renorm_omega.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void zero_order_heff(
        int sect_h,
        int sect_p, size_t dim,
        slater_det_t *det_list,
        double complex *heff
);

void construct_heff(
        int sect_h,
        int sect_p,
        slater_det_t **det_list,
        const size_t *block_sizes,
        double complex **heff_blocks,
        int n_diagrams,
        char **diagram_names
);

void diagonalize_heff(
        size_t *block_dims,
        double complex **heff,
        double complex **eigvalues,
        double complex **coef_left,
        double complex **coef_right
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
    print_model_space_info(sect_h, sect_p, block_dims, det_basis);

    /*
     * Construct Heff matrix
     */
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
    diagonalize_heff(block_dims, heff, eigvalues, coef_left, coef_right);

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


/**
 * Constructs blocks of Heff matrix.
 */
void construct_heff(
        int sect_h,
        int sect_p,
        slater_det_t **det_list,
        const size_t *block_sizes,
        double complex **heff_blocks,
        int n_diagrams,
        char **diagram_names
)
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

        // construct zero-order Hamiltonian matrix (H_0)
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
                    heff_blocks[irep][i * ms_size + j] += slater_rule(bra, ket);
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


/**
 * Diagonalizes each block of Heff matrix.
 * Returns model vectors (eigenvectors) and energies wrt vacuum state (eigenvalues).
 */
void diagonalize_heff(
        size_t *block_dims,
        double complex **heff,
        double complex **eigvalues,
        double complex **coef_left,
        double complex **coef_right
)
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





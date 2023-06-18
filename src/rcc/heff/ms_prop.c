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

#include "ms_prop.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cc_properies.h"
#include "mvcoef.h"
#include "options.h"
#include "symmetry.h"
#include "spinors.h"
#include "ms_tdm.h"
#include "engine.h"
#include "linalg.h"


/*
 * functions used in the file
 */
void guess_operator_symmetry(int nspinors, double complex *prop_matrix);

double complex get_matrix_element(double complex *mat, int *indices);

void model_space_property(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
                          cc_ms_prop_query_t *prop_query);

void print_property_matrix(
        int nroots_i, int nroots_f,
        double *energies_i, double *energies_f,
        char *rep_name_i, char *rep_name_f,
        double complex *prop_matrix);

void msprop_transform_slater_to_model(size_t nroots_i, size_t nroots_f,
                                      size_t ms_size_i, size_t ms_size_f,
                                      double complex *coefs_bra, double complex *coefs_ket,
                                      double complex *prop_slater, double complex *prop_model);

void construct_property_matrix_slater(
        int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
        size_t n_det_bra, slater_det_t *bra_dets, size_t n_det_ket, slater_det_t *ket_dets,
        double complex *prp_spinor, double complex *prp_slater
);

int read_prop_two_files(int nspinors, char *file_re, char *file_im, double complex *prop_mat);

int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);


/**
 * Calculates estimates of property operators matrix elements
 * using only model space vectors.
 */
void calculate_model_space_properties(int sect_h, int sect_p)
{
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

            cc_ms_prop_query_t *query = cc_opts->prop_queries + i;

            // here we process only queries without specifications of symmetry
            // and approximation level. Other queries should be processed by the code
            // for direct property matrix element calculations.
            if (!(query->approx_numerator == 0 &&
                query->approx_denominator == 0 &&
                strcmp(query->irrep_name, "") == 0)) {
                continue;
            }

            model_space_property(sect_h, sect_p, sect_h, sect_p, query);

            /*
             * special case of the 1h1p sector: transitions from the ground state
             */
            if (sect_h == 1 && sect_p == 1) {
                model_space_property(0, 0, 1, 1, query);
                model_space_property(1, 1, 0, 0, query);
            }

            /*
             * special case of the 1h2p sector: additional transitions 0h1p - 1h2p
             */
            if (sect_h == 1 && sect_p == 2) {
                model_space_property(0, 1, 1, 2, query);
                model_space_property(1, 2, 0, 1, query);
            }
        }
    }
}


/**
 * Model-space estimation of properties
 */
void model_space_property(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
                          cc_ms_prop_query_t *prop_query)
{
    int nspinors = get_num_spinors();

    // print_header_banner
    printf("\n");
    printf(" **\n");
    printf(" ** Model-space estimates of property matrix elements\n");
    printf(" ** Electronic transitions: %dh%dp - %dh%dp\n", bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);
    if (prop_query->source == CC_PROP_FROM_MDPROP) {
        printf(" ** MDPROP: %s\n", prop_query->prop_name);
    }
    else {
        printf(" ** From text files:\n");
        printf(" ** (real part) %s\n", prop_query->file_real);
        printf(" ** (imag part) %s\n", prop_query->file_imag);
    }
    printf(" **\n");
    printf("\n");

    /*
     * read model vectors from disk
     */
    int nrep_bra = 0;
    struct mv_block mv_blocks_bra[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(bra_sect_h, bra_sect_p, NULL, &nrep_bra, mv_blocks_bra);

    int nrep_ket = 0;
    struct mv_block mv_blocks_ket[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(ket_sect_h, ket_sect_p, NULL, &nrep_ket, mv_blocks_ket);

    /*
     * read property matrix elements (in the basis of molecular spinors)
     */
    double complex *prop_spinor = x_zeros(CC_COMPLEX, nspinors, nspinors);
    if (prop_query->source == CC_PROP_FROM_MDPROP) {
        int code = read_prop_single_file(nspinors, prop_query->prop_name, prop_spinor);
        if (code == EXIT_FAILURE) {
            printf(" Error: file '%s' with property integrals not found\n", prop_query->prop_name);
            printf(" Calculation will be skipped\n");
            goto cleanup;
        }
    }
    else { // from two txt formatted files
        int code = read_prop_two_files(nspinors, prop_query->file_real, prop_query->file_imag, prop_spinor);
        if (code == EXIT_FAILURE) {
            printf(" Error: files '%s' (real) and '%s' (imag) with property integrals not found\n",
                   prop_query->file_real, prop_query->file_imag);
            printf(" Calculation will be skipped\n");
            goto cleanup;
        }
    }
    if (prop_query->do_transpose) {
        xtranspose(CC_COMPLEX, nspinors, nspinors, prop_spinor);
    }

    //xprimat(CC_COMPLEX, prop_spinor, get_num_spinors(), get_num_spinors(), "prop_pp");

    /*
     * print information about non-zero matrix elements of the operator
     */
    guess_operator_symmetry(nspinors, prop_spinor);

    /*
     * loop over pairs of irreducible representations
     */
    for (int irep_bra = 0; irep_bra < nrep_bra; irep_bra++) {
        for (int irep_ket = 0; irep_ket < nrep_ket; irep_ket++) {
            struct mv_block *block_bra = mv_blocks_bra + irep_bra;
            struct mv_block *block_ket = mv_blocks_ket + irep_ket;

            printf(" < %-4s | prop | %-4s >", block_bra->rep_name, block_ket->rep_name);

            double complex *prop_slater = x_zeros(CC_COMPLEX, block_bra->ms_size, block_ket->ms_size);
            double complex *prop = x_zeros(CC_COMPLEX, block_bra->nroots, block_ket->nroots);

            // construct property matrix in the basis of Slater determinants
            construct_property_matrix_slater(
                    bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p,
                    block_bra->ms_size, block_bra->dets, block_ket->ms_size, block_ket->dets,
                    prop_spinor, prop_slater
            );

            // construct property matrix in the basis of model vectors
            msprop_transform_slater_to_model(
                    block_bra->nroots, block_ket->nroots, block_bra->ms_size, block_ket->ms_size,
                    block_bra->vl, block_ket->vr, prop_slater, prop
            );

            // print smart table of matrix elements
            print_property_matrix(
                    block_bra->nroots, block_ket->nroots, block_bra->energy_cm, block_ket->energy_cm,
                    block_bra->rep_name, block_ket->rep_name, prop
            );

            // cleanup
            cc_free(prop_slater);
            cc_free(prop);
        }
    }

    cleanup:
    cc_free(prop_spinor);
    for (size_t irep = 0; irep < nrep_bra; irep++) {
        mvblock_free(mv_blocks_bra + irep);
    }
    for (size_t irep = 0; irep < nrep_ket; irep++) {
        mvblock_free(mv_blocks_ket + irep);
    }
}


/**
 * Transformation of a property matrix from the basis of model Slater determinants
 * to the basis of model vectors.
 *
 * @param nroots_i    number of bra model vectors
 * @param nroots_f    number of ket model vectors
 * @param ms_size_i   dimension of bra vectors
 * @param ms_size_f   dimension of ket vectors
 * @param coefs_bra   bra vectors; row-wise
 * @param coefs_ket   ket vectors; row-wise
 * @param prop_slater property matrix in the basis of model Slater determinants
 * @param prop_model  property matrix in the basis of model vectors
 */
void msprop_transform_slater_to_model(size_t nroots_i, size_t nroots_f,
                                      size_t ms_size_i, size_t ms_size_f,
                                      double complex *coefs_bra, double complex *coefs_ket,
                                      double complex *prop_slater, double complex *prop_model)
{
    for (size_t i = 0; i < nroots_i; i++) {
        for (size_t j = 0; j < nroots_f; j++) {
            double complex prop_if = 0.0 + 0.0 * I;

            for (size_t a = 0; a < ms_size_i; a++) {
                for (size_t b = 0; b < ms_size_f; b++) {
                    double complex coef_bra = coefs_bra[ms_size_i * i + a];
                    double complex coef_ket = coefs_ket[ms_size_f * j + b];
                    double complex prop_ab = prop_slater[ms_size_f * a + b];
                    prop_if += conj(coef_bra) * coef_ket * prop_ab;
                }
            }

            prop_model[i * nroots_f + j] = prop_if;
        }
    }
}


/**
 * Constructs property operator matrix in the basis of model Slater determinants.
 */
void construct_property_matrix_slater(
        int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
        size_t n_det_bra, slater_det_t *bra_dets, size_t n_det_ket, slater_det_t *ket_dets,
        double complex *prp_spinor, double complex *prp_slater
)
{
    memset(prp_slater, 0, sizeof(double complex) * n_det_bra * n_det_ket);

    setup_slater(prp_spinor, (matrix_getter_fun) get_matrix_element, bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p, 1);

    for (size_t i = 0; i < n_det_bra; i++) {
        for (size_t j = 0; j < n_det_ket; j++) {
            slater_det_t *bra = bra_dets + i;
            slater_det_t *ket = ket_dets + j;
            prp_slater[i * n_det_ket + j] += slater_rule(bra, ket);
        }
    }
}


/**
 * Prints property matrix; "i" -- initial (ket) states, "f" -- final (bra) states
 *
 * @param nroots_i     number of bra model vectors
 * @param nroots_f     number of ket model vectors
 * @param energies_i   energies of bra states (cm-1)
 * @param energies_f   energies of ket states (cm-1)
 * @param rep_name_i   symbolic name of the bra irrep
 * @param rep_name_f   symbolic name of the ket irrep
 * @param prop_matrix  property matrix in the basis of model vectors
 */
void print_property_matrix(
        int nroots_i, int nroots_f,
        double *energies_i, double *energies_f,
        char *rep_name_i, char *rep_name_f,
        double complex *prop_matrix
)
{
    int non_zero = 0;

    for (size_t i = 0; i < nroots_i; i++) {
        for (size_t f = 0; f < nroots_f; f++) {
            double complex prop_if = prop_matrix[i * nroots_f + f];
            if (cabs(prop_if) < 1e-6) {
                continue;
            }
            non_zero++;
            if (non_zero == 1) {
                printf("  E1(cm-1)  E2(cm-1)         Re            Im          |prop|\n");
            }
            printf("%2d (%-4s) -> %2d (%-4s) %10.2f %10.2f %14.6f %14.6f %14.6f\n",
                   i + 1, rep_name_i, f + 1, rep_name_f, energies_i[i], energies_f[f],
                   creal(prop_if), cimag(prop_if), cabs(prop_if));
        }
    }

    if (non_zero == 0) {
        printf("     no allowed transitions\n");
    }
    printf("\n");
}


/**
 * "accessor" function
 * returns matrix element mat[i,j], indices=[i,j]
 */
double complex get_matrix_element(double complex *mat, int *indices)
{
    size_t i = indices[0];
    size_t j = indices[1];
    int nspinors = get_num_spinors();

    return mat[i * nspinors + j];
}



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

/*
 * Direct calculation of properties.
 */

#include "finite_order_prop.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cc_properies.h"
#include "engine.h"
#include "../heff/mvcoef.h"
#include "spinors.h"
#include "symmetry.h"
#include "ms_prop.h"
#include "sort.h"
#include "linalg.h"
#include "finite_order_prop_0h3p.h"
#include "finite_order_prop_3h0p.h"


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

void msprop_transform_slater_to_model_vectors(size_t nroots_i, size_t nroots_f,
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

void calc_finite_order_property(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
                                cc_ms_prop_query_t *prop_query);

void construct_cc_wfn_overlap_matrix(int sect_h, int sect_p, int dim, slater_det_t *det_list,
                                     double complex *overlap_slater, double scalar_term,
                                     int n_diagrams, char **list_diagrams);

void construct_finite_order_property_matrix_slater(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p, int bra_dim,
                                                   slater_det_t *bra_dets, int ket_dim, slater_det_t *ket_dets,
                                                   int is_diagonal_block, double complex *prop_slater, double scalar_term,
                                                   int n_diagrams, char **diagram_list);

void print_finite_order_property_settings(cc_ms_prop_query_t *prop_query,
                                          int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p);

void conjugate_cluster_amplitudes_0h0p();

void conjugate_cluster_amplitudes_1h0p();

void conjugate_cluster_amplitudes_0h1p();

void conjugate_cluster_amplitudes_1h1p();

void conjugate_cluster_amplitudes_2h0p();

void conjugate_cluster_amplitudes_0h2p();

void conjugate_cluster_amplitudes_0h3p();

void conjugate_cluster_amplitudes_3h0p();

void construct_effective_property_matrix_irrep_pair(
        int bra_sect_h, int bra_sect_p, struct mv_block *block_bra,
        int ket_sect_h, int ket_sect_p, struct mv_block *block_ket,
        double overlap_scalar_part, int n_overlap_diagrams, char **list_overlap_diagrams,
        double prop_scalar_part, int n_prop_diagrams, char **list_prop_diagrams,
        int scheme, int is_diagonal_block
);


/**
 * Calculation of property matrix elements using the finite-order method.
 */
void calculate_properties_finite_order_method(int sect_h, int sect_p)
{
    // estimate properties: only for the target sector
    if (cc_opts->sector_h == sect_h &&
        cc_opts->sector_p == sect_p &&
        cc_opts->n_model_space_props > 0) {
        for (int i = 0; i < cc_opts->n_model_space_props; i++) {

            cc_ms_prop_query_t *query = cc_opts->prop_queries + i;

            // ignore queries without specifications of symmetry and approximation level
            // these queries should be processed by the specific code for model-space
            // estimates of property operators matrix elements
            if (query->approx_numerator == 0 &&
                query->approx_denominator == 0 &&
                strcmp(query->irrep_name, "") == 0) {
                continue;
            }

            calc_finite_order_property(sect_h, sect_p, sect_h, sect_p, query);
        }
    }
}


/**
 * Direct calculation of property matrix elements
 */
void calc_finite_order_property(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
                                cc_ms_prop_query_t *prop_query)
{
    int nspinors = get_num_spinors();
    int sect_h = bra_sect_h;
    int sect_p = bra_sect_p;

    if ((sect_h == 1 && sect_p == 0) ||
        (sect_h == 0 && sect_p == 1) ||
        (sect_h == 2 && sect_p == 0) ||
        (sect_h == 0 && sect_p == 2) ||
        (sect_h == 1 && sect_p == 1) ||
        (sect_h == 3 && sect_p == 0) ||
        (sect_h == 0 && sect_p == 3)) {
        // OK
    }
    else {
        printf("\n");
        printf(" Warning: property calculations are not available for the %dh%dp sector\n", sect_h, sect_p);
        printf("\n");
        return;
    }

    print_finite_order_property_settings(prop_query, bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);

    /*
     * get complex conjugate for all cluster operators
     */
    conjugate_cluster_amplitudes_0h0p();
    if (sect_p >= 1) {
        conjugate_cluster_amplitudes_0h1p();
    }
    if (sect_h >= 1) {
        conjugate_cluster_amplitudes_1h0p();
    }
    if (sect_h >= 1 && sect_p >= 1) {
        conjugate_cluster_amplitudes_1h1p();
    }
    if (sect_h >= 2) {
        conjugate_cluster_amplitudes_2h0p();
    }
    if (sect_p >= 2) {
        conjugate_cluster_amplitudes_0h2p();
    }
    if (sect_h >= 3) {
        conjugate_cluster_amplitudes_3h0p();
    }
    if (sect_p >= 3) {
        conjugate_cluster_amplitudes_0h3p();
    }

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

    /*
     * print information about non-zero matrix elements of the operator
     */
    guess_operator_symmetry(nspinors, prop_spinor);

    /*
     * integral sorting
     */
    int prop_sym = get_totally_symmetric_irrep();
    if (strcmp(prop_query->irrep_name, "") != 0) {
        prop_sym = get_rep_number(prop_query->irrep_name);
    }

    request_sorting_sym("prop_pp", "pp", "00", "12", prop_sym);
    request_sorting_sym("prop_hp", "hp", "00", "12", prop_sym);
    request_sorting_sym("prop_ph", "ph", "00", "12", prop_sym);
    request_sorting_sym("prop_hh", "hh", "00", "12", prop_sym);
    sort_prop(get_num_spinors(), prop_spinor);

    /*
     * construct effective property operators:
     * <0| {e^T^\dagger} O {e^T} |0>:
     */
    int calc_disconn = prop_query->scheme == CC_DIRECT_PROP_SCHEME_CONNECTED ? 0 : 1;
    int n_prop_diagrams = 0;
    char *prop_diagrams[10];

    double complex prop_0h0p = finite_order_property_0h0p(prop_query->approx_numerator);
    if (sect_h >= 1) {
        printf(" construction of the one-body contribution (1h0p) ...\n");
        finite_order_property_1h0p("prop_1h0p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_1h0p");
    }
    if (sect_p >= 1) {
        printf(" construction of the one-body contribution (0h1p) ...\n");
        finite_order_property_0h1p("prop_0h1p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_0h1p");
    }
    if (sect_h >= 1 && sect_p >= 1) {

        printf(" construction of the one-body contribution (0h0p -> 1h1p) ...\n");
        direct_property_0h0p_1h1p("prop_0h0p_1h1p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_0h0p_1h1p");

        printf(" construction of the one-body contribution (1h1p -> 0h0p) ...\n");
        direct_property_1h1p_0h0p("prop_1h1p_0h0p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_1h1p_0h0p");

        printf(" construction of the two-body contribution (1h1p) ...\n");
        direct_property_1h1p("prop_1h1p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_1h1p");

    }
    if (sect_h >= 2) {
        printf(" construction of the two-body contribution (2h0p) ...\n");
        finite_order_property_2h0p("prop_2h0p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_2h0p");
    }
    if (sect_p >= 2) {
        printf(" construction of the two-body contribution (0h2p) ...\n");
        finite_order_property_0h2p("prop_0h2p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_0h2p");
    }
    if (sect_h >= 3) {
        printf(" construction of the three-body contribution (3h0p) ...\n");
        finite_order_property_3h0p("prop_3h0p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_3h0p");
    }
    if (sect_p >= 3) {
        printf(" construction of the three-body contribution (0h3p) ...\n");
        finite_order_property_0h3p("prop_0h3p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_0h3p");
    }

    printf(" scalar contribution to property = %.8f %.8f\n", creal(prop_0h0p), cimag(prop_0h0p));
    printf(" effective property operator diagrams:\n");
    for (int i = 0; i < n_prop_diagrams; i++) {
        printf(" - %s\n", prop_diagrams[i]);
    }

    /*
     * construct overlap: <0| {e^T^\dagger} {e^T} |0>:
     * scalar part: 1 + fully contracted diagrams without any external lines
     */
    int n_overlap_diagrams = 0;
    char *overlap_diagrams[10];

    double complex overlap_scalar_part = sector_0h0p_overlap(prop_query->approx_denominator);
    if (sect_h >= 1) {
        printf(" construction of the overlap diagram (1h0p) ...\n");
        sector_1h0p_overlap("overlap_1h0p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h0p");
    }
    if (sect_p >= 1) {
        printf(" construction of the overlap diagram (0h1p) ...\n");
        sector_0h1p_overlap("overlap_0h1p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h1p");
    }
    if (sect_h >= 1 && sect_p >= 1) {
        printf(" construction of the overlap diagram (1h1p) ...\n");
        sector_1h1p_overlap("overlap_1h1p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h1p");
    }
    if (sect_h >= 2) {
        printf(" construction of the overlap diagram (2h0p) ...\n");
        sector_2h0p_overlap("overlap_2h0p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_2h0p");
    }
    if (sect_p >= 2) {
        printf(" construction of the overlap diagram (0h2p) ...\n");
        sector_0h2p_overlap("overlap_0h2p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h2p");
    }
    if (sect_h >= 3) {
        printf(" construction of the overlap diagram (3h0p) ...\n");
        sector_3h0p_overlap("overlap_3h0p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_3h0p");
    }
    if (sect_p >= 3) {
        printf(" construction of the overlap diagram (0h3p) ...\n");
        sector_0h3p_overlap("overlap_0h3p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h3p");
    }

    printf(" scalar contribution to overlap = %.8f %.8f\n",
           creal(overlap_scalar_part), cimag(overlap_scalar_part));
    printf(" overlap diagrams:\n");
    for (int i = 0; i < n_overlap_diagrams; i++) {
        printf(" - %s\n", overlap_diagrams[i]);
    }
    printf("\n");

    /*
     * special case: 0h0p -> 1h1p transitions
     */
    if (sect_h == 1 && sect_p == 1) {

        mv_block_t vac_block;

        vac_block.eigval = z_zeros(1, 1);
        vac_block.energy_cm = d_zeros(1, 1);
        vac_block.nroots = 1;
        vac_block.ms_size = 1;

        int vac_irrep = get_vacuum_irrep();
        char *vac_irrep_name = get_irrep_name(vac_irrep);
        strcpy(vac_block.rep_name, vac_irrep_name);

        vac_block.dets = (slater_det_t *) cc_calloc(1, sizeof(slater_det_t));
        vac_block.dets[0].sym = vac_irrep;

        vac_block.vl = z_zeros(1, 1);
        vac_block.vr = z_zeros(1, 1);
        vac_block.vl[0] = 1.0 + 0.0 * I;
        vac_block.vr[0] = 1.0 + 0.0 * I;

        printf("\n 0h0p -> 1h1p transitions\n\n");

        for (int irep_ket = 0; irep_ket < nrep_ket; irep_ket++) {
            struct mv_block *block_ket = mv_blocks_ket + irep_ket;

            construct_effective_property_matrix_irrep_pair(
                    0, 0, &vac_block,
                    1, 1, block_ket,
                    overlap_scalar_part, n_overlap_diagrams, overlap_diagrams,
                    prop_0h0p, n_prop_diagrams, prop_diagrams,
                    prop_query->scheme, 0
            );
        }

        printf("\n 1h1p -> 0h0p transitions\n\n");

        for (int irep_bra = 0; irep_bra < nrep_ket; irep_bra++) {
            struct mv_block *block_bra = mv_blocks_bra + irep_bra;

            construct_effective_property_matrix_irrep_pair(
                    1, 1, block_bra,
                    0, 0, &vac_block,
                    overlap_scalar_part, n_overlap_diagrams, overlap_diagrams,
                    prop_0h0p, n_prop_diagrams, prop_diagrams,
                    prop_query->scheme, 0
            );
        }
    }
    else {
        /*
         * loop over pairs of irreducible representations
         * (except for 1h1p)
         */
        for (int irep_bra = 0; irep_bra < nrep_bra; irep_bra++) {
            for (int irep_ket = 0; irep_ket < nrep_ket; irep_ket++) {
                struct mv_block *block_bra = mv_blocks_bra + irep_bra;
                struct mv_block *block_ket = mv_blocks_ket + irep_ket;

                construct_effective_property_matrix_irrep_pair(
                        sect_h, sect_p, block_bra,
                        sect_h, sect_p, block_ket,
                        overlap_scalar_part, n_overlap_diagrams, overlap_diagrams,
                        prop_0h0p, n_prop_diagrams, prop_diagrams,
                        prop_query->scheme, irep_bra == irep_ket
                );
            }
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
    for (int i = 0; i < n_overlap_diagrams; i++) {
        cc_free(overlap_diagrams[i]);
    }
    for (int i = 0; i < n_prop_diagrams; i++) {
        cc_free(prop_diagrams[i]);
    }
}


void construct_effective_property_matrix_irrep_pair(
        int bra_sect_h, int bra_sect_p, struct mv_block *block_bra,
        int ket_sect_h, int ket_sect_p, struct mv_block *block_ket,
        double overlap_scalar_part, int n_overlap_diagrams, char **list_overlap_diagrams,
        double prop_scalar_part, int n_prop_diagrams, char **list_prop_diagrams,
        int scheme, int is_diagonal_block
)
{
    printf(" < %-4s | prop | %-4s >", block_bra->rep_name, block_ket->rep_name);

    double complex *overlap_bra_slater = x_zeros(CC_COMPLEX, block_bra->ms_size, block_bra->ms_size);
    double complex *overlap_ket_slater = x_zeros(CC_COMPLEX, block_ket->ms_size, block_ket->ms_size);
    double complex *overlap_bra = x_zeros(CC_COMPLEX, block_bra->nroots, block_bra->nroots);
    double complex *overlap_ket = x_zeros(CC_COMPLEX, block_ket->nroots, block_ket->nroots);

    double complex *prop_slater = x_zeros(CC_COMPLEX, block_bra->ms_size, block_ket->ms_size);
    double complex *prop = x_zeros(CC_COMPLEX, block_bra->nroots, block_ket->nroots);

    /*
     * construct overlap integrals: bra vectors
     */
    construct_cc_wfn_overlap_matrix(
            bra_sect_h, bra_sect_p, block_bra->ms_size, block_bra->dets,
            overlap_bra_slater, overlap_scalar_part,
            n_overlap_diagrams, list_overlap_diagrams
    );

    msprop_transform_slater_to_model_vectors(
        block_bra->nroots, block_bra->nroots, block_bra->ms_size, block_bra->ms_size,
        block_bra->vr, block_bra->vr, overlap_bra_slater, overlap_bra
    );

    /*
     * construct overlap integrals: ket vectors
     */
    construct_cc_wfn_overlap_matrix(
            ket_sect_h, ket_sect_p, block_ket->ms_size, block_ket->dets,
            overlap_ket_slater, overlap_scalar_part,
            n_overlap_diagrams, list_overlap_diagrams
    );

    msprop_transform_slater_to_model_vectors(
        block_ket->nroots, block_ket->nroots, block_ket->ms_size, block_ket->ms_size,
        block_ket->vr, block_ket->vr, overlap_ket_slater, overlap_ket
    );

    /*
     * construct property matrix in the basis of Slater determinants
     */
    construct_finite_order_property_matrix_slater(
        bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p,
        block_bra->ms_size, block_bra->dets,
        block_ket->ms_size, block_ket->dets,
        is_diagonal_block, prop_slater,
        prop_scalar_part, n_prop_diagrams, list_prop_diagrams
    );

    /*
     * construct property matrix in the basis of model vectors
     */
    msprop_transform_slater_to_model_vectors(
        block_bra->nroots, block_ket->nroots, block_bra->ms_size, block_ket->ms_size,
        (scheme == CC_DIRECT_PROP_SCHEME_HERMITIAN) ? block_bra->vr : block_bra->vl, block_ket->vr,
        prop_slater, prop
    );

    if (scheme == CC_DIRECT_PROP_SCHEME_HERMITIAN) {
        /*
         * divide by norms of wavefunctions
         */
        for (int i = 0; i < block_bra->nroots; i++) {
            for (int j = 0; j < block_ket->nroots; j++) {

                double complex overlap_ii = overlap_bra[i * block_bra->nroots + i];
                double complex overlap_jj = overlap_ket[j * block_ket->nroots + j];

                double norm_i = sqrt(creal(overlap_ii));
                double norm_j = sqrt(creal(overlap_jj));

                prop[i * block_ket->nroots + j] /= (norm_i * norm_j);
            }
        }
    }
    else if (scheme == CC_DIRECT_PROP_SCHEME_NON_HERMITIAN) {

        double complex *overlap_bra_inv = x_zeros(CC_COMPLEX, block_bra->nroots, block_bra->nroots);
        double complex *prop_tmp = x_zeros(CC_COMPLEX, block_bra->nroots, block_ket->nroots);

        memcpy(prop_tmp, prop, sizeof(double complex) * block_bra->nroots * block_ket->nroots);
        inverse_matrix(block_bra->nroots, overlap_bra, overlap_bra_inv);

        for (int i = 0; i < block_bra->nroots; i++) {
            for (int j = 0; j < block_ket->nroots; j++) {

                double complex prop_ij = 0.0 + 0.0 * I;

                double complex overlap_ii = overlap_bra[i * block_bra->nroots + i];
                double complex overlap_jj = overlap_ket[j * block_ket->nroots + j];
                double norm_i = sqrt(creal(overlap_ii));
                double norm_j = sqrt(creal(overlap_jj));

                for (int m = 0; m < block_bra->nroots; m++) {
                    double complex overlap_inv_im = overlap_bra_inv[i * block_bra->nroots + m];
                    double complex prop_tmp_ij = prop_tmp[m * block_ket->nroots + j];

                    prop_ij += overlap_inv_im * prop_tmp_ij;
                }

                prop[i * block_ket->nroots + j] = prop_ij * norm_i / norm_j;
            }
        }

        cc_free(overlap_bra_inv);
        cc_free(prop_tmp);
    }
    else if (scheme == CC_DIRECT_PROP_SCHEME_CONNECTED) {

        for (int i = 0; i < block_bra->nroots; i++) {
            for (int j = 0; j < block_ket->nroots; j++) {

                double complex overlap_ii = overlap_bra[i * block_bra->nroots + i];
                double complex overlap_jj = overlap_ket[j * block_ket->nroots + j];

                double norm_i = sqrt(creal(overlap_ii));
                double norm_j = sqrt(creal(overlap_jj));

                prop[i * block_ket->nroots + j] = prop[i * block_ket->nroots + j] * norm_i / norm_j;
            }
        }
    }

    /*
     * print smart table of matrix elements
     */
    print_property_matrix(
            block_bra->nroots, block_ket->nroots, block_bra->energy_cm, block_ket->energy_cm,
            block_bra->rep_name, block_ket->rep_name, prop
    );

    /*
     * cleanup
     */
    cc_free(overlap_bra_slater);
    cc_free(overlap_bra);
    cc_free(overlap_ket_slater);
    cc_free(overlap_ket);
    cc_free(prop_slater);
    cc_free(prop);
}


void print_finite_order_property_settings(cc_ms_prop_query_t *prop_query,
                                          int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p)
{
    printf("\n");
    printf(" **\n");
    printf(" ** Finite-order calculation of property matrix elements\n");
    printf(" ** Electronic transitions: %dh%dp - %dh%dp\n", bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);
    if (prop_query->source == CC_PROP_FROM_MDPROP) {
        printf(" ** MDPROP: %s\n", prop_query->prop_name);
    }
    else {
        printf(" ** From text files:\n");
        printf(" ** (real part) %s\n", prop_query->file_real);
        printf(" ** (imag part) %s\n", prop_query->file_imag);
    }
    printf(" ** Operator symmetry: ");
    if (strcmp(prop_query->irrep_name, "") == 0) {
        printf("totally symmetric\n");
    }
    else {
        printf("%s (irrep %d)\n", prop_query->irrep_name, get_rep_number(prop_query->irrep_name));
    }

    printf(" ** Scheme: ");
    if (prop_query->scheme == CC_DIRECT_PROP_SCHEME_HERMITIAN) {
        printf("hermitian");
    }
    else if (prop_query->scheme == CC_DIRECT_PROP_SCHEME_NON_HERMITIAN) {
        printf("non-hermitian");
    }
    else if (prop_query->scheme == CC_DIRECT_PROP_SCHEME_CONNECTED) {
        printf("connected");
    }
    printf("\n");

    printf(" ** Order of numerator: ");
    if (prop_query->approx_numerator == CC_PROPERTIES_APPROX_MODEL_SPACE) {
        printf("zero\n");
    }
    else if (prop_query->approx_numerator == CC_PROPERTIES_APPROX_LINEAR) {
        printf("linear\n");
    }
    else if (prop_query->approx_numerator == CC_PROPERTIES_APPROX_QUADRATIC) {
        printf("quadratic\n");
    }
    else {
        printf("unknown\n");
    }
    printf(" ** Order of denominator: ");
    if (prop_query->approx_denominator == CC_PROPERTIES_APPROX_MODEL_SPACE) {
        printf("zero\n");
    }
    else if (prop_query->approx_denominator == CC_PROPERTIES_APPROX_QUADRATIC) {
        printf("quadratic\n");
    }
    else {
        printf("unknown\n");
    }
    printf(" **\n");
    printf("\n");
}


void conjugate_cluster_amplitudes_0h0p()
{
    diagram_conjugate("t1c", "t1c+");
    diagram_conjugate("t2c", "t2c+");
    if (triples_enabled()) {
        diagram_conjugate("t3c", "t3c+");
    }
}


void conjugate_cluster_amplitudes_1h0p()
{
    diagram_conjugate("h1c", "h1c+");
    diagram_conjugate("h2c", "h2c+");
    if (triples_enabled()) {
        diagram_conjugate("h3c", "h3c+");
    }
}


void conjugate_cluster_amplitudes_0h1p()
{
    diagram_conjugate("s1c", "s1c+");
    diagram_conjugate("s2c", "s2c+");
    if (triples_enabled()) {
        diagram_conjugate("s3c", "s3c+");
    }
}


void conjugate_cluster_amplitudes_2h0p()
{
    diagram_conjugate("g2c", "g2c+");
    if (triples_enabled()) {
        diagram_conjugate("g3c", "g3c+");
    }
}


void conjugate_cluster_amplitudes_0h2p()
{
    diagram_conjugate("x2c", "x2c+");
    if (triples_enabled()) {
        diagram_conjugate("x3c", "x3c+");
    }
}


void conjugate_cluster_amplitudes_1h1p()
{
    diagram_conjugate("e1c", "e1c+");
    diagram_conjugate("e2c", "e2c+");
    if (triples_enabled()) {
        diagram_conjugate("e3c", "e3c+");
    }
}


void conjugate_cluster_amplitudes_0h3p()
{
    if (triples_enabled()) {
        diagram_conjugate("z3c", "z3c+");
    }
}


void conjugate_cluster_amplitudes_3h0p()
{
    if (triples_enabled()) {
        diagram_conjugate("k3c", "k3c+");
    }
}


void construct_finite_order_property_matrix_slater(
        int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
        int bra_dim, slater_det_t *bra_dets, int ket_dim, slater_det_t *ket_dets,
        int is_diagonal_block, double complex *prop_slater,
        double scalar_term, int n_diagrams, char **diagram_list)
{
    memset(prop_slater, 0, sizeof(double complex) * bra_dim * ket_dim);

    if (is_diagonal_block) {
        for (int i = 0; i < bra_dim; i++) {
            prop_slater[i * bra_dim + i] = scalar_term;
        }
    }

    for (int k = 0; k < n_diagrams; k++) {
        char *dg_name = diagram_list[k];
        int npart = rank(dg_name) / 2;
        if (2 * npart > bra_sect_h + bra_sect_p + ket_sect_h + ket_sect_p) {
            continue;
        }

        diagram_t *dg_norm = diagram_stack_find(dg_name);
        setup_slater(dg_norm, (matrix_getter_fun) diagram_get,
                     bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p, npart);

        for (int i = 0; i < bra_dim; i++) {
            for (int j = 0; j < ket_dim; j++) {
                slater_det_t *bra = bra_dets + i;
                slater_det_t *ket = ket_dets + j;

                prop_slater[i * ket_dim + j] += slater_rule(bra, ket);
            }
        }
    }
}


void print_time_mem_usage(char *diagram_name, double t_start, int print_level)
{
    if (print_level < CC_PRINT_HIGH) {
        return;
    }

    const double b2mb = 1.0 / (1024.0 * 1024.0);
    const double b2gb = 1.0 / (1024.0 * 1024.0 * 1024.0);

    size_t max_allocated = cc_get_peak_memory_usage();
    double time_elapsed = abs_time() - t_start;

    printf(" %-20s%14.3f%14.2f%14.2f\n", diagram_name, time_elapsed,
           b2mb * (double) max_allocated, b2gb * (double) max_allocated);
}


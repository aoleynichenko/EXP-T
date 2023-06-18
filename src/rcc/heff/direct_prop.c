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
 * Direct calculation of properties.
 */

#include "direct_prop.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cc_properies.h"
#include "datamodel.h"
#include "engine.h"
#include "mvcoef.h"
#include "spinors.h"
#include "symmetry.h"
#include "ms_prop.h"
#include "sort.h"
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

void direct_property(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
                     cc_ms_prop_query_t *prop_query);

void construct_cc_wfn_overlap_matrix(int sect_h, int sect_p, int dim, slater_det_t *det_list,
                                     double complex *overlap_slater, double scalar_term,
                                     int n_diagrams, char **list_diagrams);

void direct_property_0h1p(char *result_name, int operator_symmetry, int approximation, int disconnected);

void direct_property_1h0p(char *result_name, int operator_symmetry, int approximation, int disconnected);

void direct_property_0h2p(char *result_name, int operator_symmetry, int approximation, int disconnected);

void direct_property_0h0p_1h1p(char *result_name, int operator_symmetry, int approximation, int disconnected);

void direct_property_1h1p_0h0p(char *result_name, int operator_symmetry, int approximation, int disconnected);

void direct_property_1h1p(char *result_name, int operator_symmetry, int approximation, int disconnected);

double complex direct_property_0h0p(int approximation);

void construct_direct_property_matrix(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p, int bra_dim,
                                      slater_det_t *bra_dets, int ket_dim, slater_det_t *ket_dets,
                                      int is_diagonal_block, double complex *overlap_slater, double scalar_term,
                                      int n_diagrams, char **diagram_list);

double complex sector_0h0p_overlap(int approximation);

void sector_1h0p_overlap(char *result_name, int approximation);

void sector_0h1p_overlap(char *result_name, int approximation);

void sector_1h1p_overlap(char *result_name, int approximation);

void sector_0h2p_overlap(char *result_name, int approximation);

void construct_disconnected_rank2_rank2(char *src1_name, char *src2_name, char *dst_name);

void print_direct_property_settings(cc_ms_prop_query_t *prop_query,
                                    int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p);

void conjugate_cluster_amplitudes_0h0p();

void conjugate_cluster_amplitudes_1h0p();

void conjugate_cluster_amplitudes_0h1p();

void conjugate_cluster_amplitudes_1h1p();

void conjugate_cluster_amplitudes_0h2p();

void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence);

void construct_effective_property_matrix_irrep_pair(
        int bra_sect_h, int bra_sect_p, struct mv_block *block_bra,
        int ket_sect_h, int ket_sect_p, struct mv_block *block_ket,
        double overlap_scalar_part, int n_overlap_diagrams, char **list_overlap_diagrams,
        double prop_scalar_part, int n_prop_diagrams, char **list_prop_diagrams,
        int scheme, int is_diagonal_block
);


/**
 * Calculation of property matrix elements using the direct formula.
 */
void calculate_properties_direct(int sect_h, int sect_p)
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

            direct_property(sect_h, sect_p, sect_h, sect_p, query);
        }
    }
}


/**
 * Direct calculation of property matrix elements
 */
void direct_property(int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
                     cc_ms_prop_query_t *prop_query)
{
    int nspinors = get_num_spinors();
    int sect_h = bra_sect_h;
    int sect_p = bra_sect_p;

    if ((sect_h == 1 && sect_p == 0) ||
        (sect_h == 0 && sect_p == 1) ||
        (sect_h == 0 && sect_p == 2) ||
        (sect_h == 1 && sect_p == 1)) {
        // OK
    }
    else {
        printf("\n");
        printf(" Warning: properties calculations are not available for this sector\n");
        printf("\n");
        return;
    }

    print_direct_property_settings(prop_query, bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);

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
    if (sect_p >= 2) {
        conjugate_cluster_amplitudes_0h2p();
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

    double complex prop_0h0p = direct_property_0h0p(prop_query->approx_numerator);
    if (sect_h >= 1) {
        printf(" construction of the one-body contribution (1h0p) ...\n");
        direct_property_1h0p("prop_1h0p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_1h0p");
    }
    if (sect_p >= 1) {
        printf(" construction of the one-body contribution (0h1p) ...\n");
        direct_property_0h1p("prop_0h1p", prop_sym, prop_query->approx_numerator, calc_disconn);
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
    if (sect_p >= 2) {
        printf(" construction of the two-body contribution (0h2p) ...\n");
        direct_property_0h2p("prop_0h2p", prop_sym, prop_query->approx_numerator, calc_disconn);
        prop_diagrams[n_prop_diagrams++] = cc_strdup("prop_0h2p");
    }

    printf(" scalar contribution to property = %.8f %.8f\n", creal(prop_0h0p), cimag(prop_0h0p));
    printf(" effective property operator diagrams: ");
    for (int i = 0; i < n_prop_diagrams; i++) {
        printf("%s ", prop_diagrams[i]);
    }
    printf("\n");

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
    if (sect_p >= 2) {
        printf(" construction of the overlap diagram (0h2p) ...\n");
        sector_0h2p_overlap("overlap_0h2p", prop_query->approx_denominator);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h2p");
    }

    printf(" scalar contribution to overlap = %.8f %.8f\n",
           creal(overlap_scalar_part), cimag(overlap_scalar_part));
    printf(" overlap diagrams: ");
    for (int i = 0; i < n_overlap_diagrams; i++) {
        printf("%s ", overlap_diagrams[i]);
    }
    printf("\n\n");

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

    msprop_transform_slater_to_model(
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

    msprop_transform_slater_to_model(
            block_ket->nroots, block_ket->nroots, block_ket->ms_size, block_ket->ms_size,
            block_ket->vr, block_ket->vr, overlap_ket_slater, overlap_ket
    );

    /*
     * construct property matrix in the basis of Slater determinants
     */
    construct_direct_property_matrix(
            bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p,
            block_bra->ms_size, block_bra->dets,
            block_ket->ms_size, block_ket->dets,
            is_diagonal_block, prop_slater,
            prop_scalar_part, n_prop_diagrams, list_prop_diagrams
    );

    /*
     * construct property matrix in the basis of model vectors
     */
    msprop_transform_slater_to_model(
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


void print_direct_property_settings(cc_ms_prop_query_t *prop_query,
                                    int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p)
{
    printf("\n");
    printf(" **\n");
    printf(" ** Direct calculation of property matrix elements\n");
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


void construct_direct_property_matrix(
        int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
        int bra_dim, slater_det_t *bra_dets, int ket_dim, slater_det_t *ket_dets,
        int is_diagonal_block, double complex *overlap_slater,
        double scalar_term, int n_diagrams, char **diagram_list)
{
    memset(overlap_slater, 0, sizeof(double complex) * bra_dim * ket_dim);

    if (is_diagonal_block) {
        for (int i = 0; i < bra_dim; i++) {
            overlap_slater[i * bra_dim + i] = scalar_term;
        }
    }

    //printf("%d%d -> %d%d\n", bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p);

    for (int k = 0; k < n_diagrams; k++) {
        char *dg_name = diagram_list[k];
        //printf("dg_name = %s\n", dg_name);
        if (strcmp(dg_name, "prop_0h0p_1h1p") == 0) {
            //prt(dg_name);
        }

        int npart = rank(dg_name) / 2;
        if (2 * npart > bra_sect_h + bra_sect_p + ket_sect_h + ket_sect_p) {
            continue;
        }

        diagram_t *dg_norm = diagram_stack_find(dg_name);
        setup_slater(dg_norm, (matrix_getter_fun) diagram_get,
                     bra_sect_h, bra_sect_p, ket_sect_h, ket_sect_p, npart);
        //printf("end setup slater\n");

        //printf("slater_rule = %p\n", slater_rule);

        for (int i = 0; i < bra_dim; i++) {
            for (int j = 0; j < ket_dim; j++) {
                slater_det_t *bra = bra_dets + i;
                slater_det_t *ket = ket_dets + j;
                //printf("i=%d j=%d begin\n", i, j);
                //print_slater_det(stdout, bra_sect_h, bra_sect_p, bra);
                //print_slater_det(stdout, ket_sect_h, ket_sect_p, ket);

                if (strcmp(dg_name, "prop_0h0p_1h1p") == 0) {
                    //print_slater_det(stdout, bra_sect_h, bra_sect_p, bra);
                    //print_slater_det(stdout, ket_sect_h, ket_sect_p, ket);
                    //printf("%.10f  %.10f\n", creal(slater_rule(bra, ket)), cimag(slater_rule(bra, ket)));
                }

                overlap_slater[i * ket_dim + j] += slater_rule(bra, ket);
                //printf("i=%d j=%d end\n", i, j);
            }
        }
    }

    for (int i = 0; i < bra_dim; i++) {
        for (int j = 0; j < ket_dim; j++) {
            slater_det_t *bra = bra_dets + i;
            slater_det_t *ket = ket_dets + j;
            double complex prop_ij = overlap_slater[i * ket_dim + j];
            //printf("[%d,%d] = %.6f %.6f\n", i, j, creal(prop_ij), cimag(prop_ij));
        }
    }
}


void print_time_mem_usage(char *diagram_name, double t_start)
{
    const double b2mb = 1.0 / (1024.0 * 1024.0);
    const double b2gb = 1.0 / (1024.0 * 1024.0 * 1024.0);

    size_t max_allocated = cc_get_peak_memory_usage();
    double time_elapsed = abs_time() - t_start;

    printf(" %-8s%14.3f%14.2f%14.2f\n", diagram_name, time_elapsed,
           b2mb * (double) max_allocated, b2gb * (double) max_allocated);
}


double complex direct_property_0h0p(int approximation)
{
    dg_stack_pos_t pos = get_stack_pos();

    double complex d = 0.0 + 0.0 * I;

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return d;
    }

    // L1a + L2a
    d += 2.0 * scalar_product("C", "N", "prop_hp", "t1c");

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return d;
    }

    // Q1a
    reorder("t1c", "r1", "21");
    mult("prop_hh", "r1", "r2", 1);
    d -= scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    // Q1b
    reorder("prop_pp", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    d += scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    // Q2a + Q2b
    reorder("t2c", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    d += 2.0 * scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    // Q3a
    reorder("t2c", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    reorder("r2", "r3", "3412");
    d += -0.5 * scalar_product("C", "N", "t2c", "r3");
    restore_stack_pos(pos);

    // Q3b
    reorder("prop_pp", "r1", "21");
    mult("t2c", "r1", "r2", 1);
    d += 0.5 * scalar_product("C", "N", "t2c", "r2");
    restore_stack_pos(pos);

    if (triples_enabled()) {

        // Q4a + Q4b
        reorder("t3c", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        d += (0.25 + 0.25) * scalar_product("C", "N", "prop_hp", "r2");
        restore_stack_pos(pos);

        // Q5a
        reorder("t3c", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        reorder("r2", "r3", "456123");
        d -= (1.0 / 12.0) * scalar_product("C", "N", "t3c", "r3");
        restore_stack_pos(pos);

        // Q5b
        reorder("prop_pp", "r1", "21");
        mult("t3c", "r1", "r2", 1);
        d += (1.0 / 12.0) * scalar_product("C", "N", "t3c", "r2");
        restore_stack_pos(pos);
    }

    return d;
}


void direct_property_0h1p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "pp", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c", "t1c_v2", "01", 0);
    restrict_valence("t1c+", "t1c+_v1", "01", 0);
    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);

    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);
    restrict_valence("t2c", "t2c_v2", "0001", 0);
    restrict_valence("t2c+", "t2c+_v2", "0100", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    restrict_valence("s2c", "s2c_v2", "1001", 0);
    restrict_valence("s2c+", "s2c+_v2", "0110", 0);

    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hp", "prop_hv", "01", 0);

    restrict_valence("prop_pp", "prop_pv", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);
    restrict_valence("prop_pp", "prop_vv", "11", 0);

    dg_stack_pos_t pos = get_stack_pos();

    printf(" diagram     time (sec)  max mem (mb)  max mem (gb)\n");

    // O1
    double t_start_O1 = abs_time();
    update(result_name, 1.0, "prop_vv");
    restore_stack_pos(pos);
    print_time_mem_usage("O1", t_start_O1);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    double t_start_L1a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("prop_vh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L1a", t_start_L1a);

    // L1b
    double t_start_L1b = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("t1c+_v", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L1b", t_start_L1b);

    // L2a
    double t_start_L2a = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L2a", t_start_L2a);

    // L2b
    double t_start_L2b = abs_time();
    reorder("s1c+", "r1", "21");
    mult("prop_vp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L2b", t_start_L2b);

    // L3a
    double t_start_L3a = abs_time();
    reorder("s2c_v1", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L3a", t_start_L3a);

    // L3b
    double t_start_L3b = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L3b", t_start_L3b);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a
    double t_start_Q1a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1a", t_start_Q1a);

    // Q1b
    double t_start_Q1b = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1b", t_start_Q1b);

    // Q1c
    double t_start_Q1c = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1c", t_start_Q1c);

    // Q1d - disconnected term
    if (disconnected) {
        double t_start_Q1d = abs_time();
        double t1_t1 = scalar_product("C", "N", "t1c", "t1c");
        update(result_name, t1_t1, "prop_vv");
        restore_stack_pos(pos);
        print_time_mem_usage("Q1d", t_start_Q1d);
    }

    // Q2a
    double t_start_Q2a = abs_time();
    reorder("t2c_v1", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2a", t_start_Q2a);

    // Q2b
    double t_start_Q2b = abs_time();
    reorder("t2c+_v1", "r1", "1342");
    reorder("t1c_v", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2b", t_start_Q2b);

    // Q2c
    double t_start_Q2c = abs_time();
    reorder("t2c_v1", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2c", t_start_Q2c);

    // Q2d
    double t_start_Q2d = abs_time();
    reorder("prop_hv", "r1", "21");
    reorder("t2c+_v1", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2d", t_start_Q2d);

    // Q3a - disconnected term
    if (disconnected) {
        double t_start_Q3a = abs_time();
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        update(result_name, t2_t2, "prop_vv");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3a", t_start_Q3a);
    }

    // Q3b
    double t_start_Q3b = abs_time();
    reorder("t2c+", "r1", "3412");
    reorder("t2c_v2", "r2", "4123");
    mult("prop_vp", "r1", "r3", 1);
    mult("r3", "r2", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3b", t_start_Q3b);

    // Q3c
    double t_start_Q3c = abs_time();
    reorder("prop_pv", "r1", "21");
    reorder("t2c+_v2", "r2", "2341");
    mult("r1", "t2c", "r3", 1);
    mult("r2", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3c", t_start_Q3c);

    // Q3d
    double t_start_Q3d = abs_time();
    reorder("t2c_v1", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("t2c+_v1", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3d", t_start_Q3d);

    // Q3e
    double t_start_Q3e = abs_time();
    reorder("prop_pp", "r1", "21");
    mult("t2c_v1", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v1", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3e", t_start_Q3e);

    // Q4a
    double t_start_Q4a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4a", t_start_Q4a);

    // Q4b
    double t_start_Q4b = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4b", t_start_Q4b);

    // Q4c
    double t_start_Q4c = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4c", t_start_Q4c);

    // Q4d
    double t_start_Q4d = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4d", t_start_Q4d);

    // Q5a
    double t_start_Q5a = abs_time();
    reorder("prop_hh", "r1", "21");
    reorder("s2c_v1", "r2", "1324");
    mult("r1", "t1c+", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5a", t_start_Q5a);

    // Q5b
    double t_start_Q5b = abs_time();
    reorder("t1c", "r1", "21");
    reorder("s2c+_v1", "r2", "1342");
    mult("prop_hh", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5b", t_start_Q5b);

    // Q5c
    double t_start_Q5c = abs_time();
    reorder("t1c+", "r1", "21");
    reorder("s2c_v1", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5c", t_start_Q5c);

    // Q5d
    double t_start_Q5d = abs_time();
    reorder("s2c+_v1", "r1", "1324");
    reorder("prop_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5d", t_start_Q5d);

    // Q6a
    double t_start_Q6a = abs_time();
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("prop_hv", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6a", t_start_Q6a);

    // Q6b
    double t_start_Q6b = abs_time();
    reorder("t2c", "r1", "2341");
    reorder("s2c+", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("prop_vh", "r3", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6b", t_start_Q6b);

    // Q6c
    double t_start_Q6c = abs_time();
    reorder("s2c_v1", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6c", t_start_Q6c);

    // Q6d
    double t_start_Q6d = abs_time();
    reorder("s2c+_v1", "r1", "1324");
    reorder("t2c", "r2", "3142");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6d", t_start_Q6d);

    // Q7
    double t_start_Q7 = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q7", t_start_Q7);

    // Q8a
    double t_start_Q8a = abs_time();
    reorder("s2c", "r1", "1342");
    reorder("s1c+", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8a", t_start_Q8a);

    // Q8b
    double t_start_Q8b = abs_time();
    reorder("s2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8b", t_start_Q8b);

    // Q8c
    double t_start_Q8c = abs_time();
    reorder("s2c", "r1", "1342");
    reorder("prop_pv", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8c", t_start_Q8c);

    // Q8d
    double t_start_Q8d = abs_time();
    reorder("s2c+", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8d", t_start_Q8d);

    // Q9a
    double t_start_Q9a = abs_time();
    reorder("s2c+", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("s2c", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9a", t_start_Q9a);

    // Q9b
    double t_start_Q9b = abs_time();
    reorder("s2c+", "r1", "3124");
    reorder("s2c", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9b", t_start_Q9b);

    if (triples_enabled()) {

        // Q10a
        reorder("t3c", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        reorder("r2", "r3", "21");
        mult("prop_ph", "r3", "r4", 1);
        closed("r4", "r5");
        update(result_name, -0.25, "r5");
        restore_stack_pos(pos);

        // Q10b
        reorder("t3c+", "r1", "145623");
        reorder("prop_hp", "r2", "21");
        mult("r1", "t2c", "r3", 4);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        update(result_name, -0.25, "r5");
        restore_stack_pos(pos);

        // Q10c
        reorder("t3c", "r1", "451263");
        mult("r1", "prop_ph", "r2", 2);
        mult("t2c+", "r2", "r3", 3);
        closed("r3", "r4");
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);

        // Q10d
        reorder("t2c", "r1", "3412");
        reorder("t3c+", "r2", "124563");
        mult("r2", "prop_hp", "r3", 2);
        mult("r3", "r1", "r4", 3);
        closed("r4", "r5");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);

        // Q11a
        reorder("t3c", "r1", "456123");
        reorder("prop_pp", "r2", "21");
        mult("t3c+", "r1", "r3", 5);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        update(result_name, -1.0 / 12.0, "r5");
        restore_stack_pos(pos);

        // Q11b
        reorder("t3c", "r1", "456123");
        mult("r1", "t3c+", "r2", 5);
        mult("prop_pp", "r2", "r3", 1);
        closed("r3", "r4");
        update(result_name, -1.0 / 12.0, "r4");
        restore_stack_pos(pos);

        // Q11c
        reorder("prop_pp", "r1", "21");
        mult("t3c", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("t3c+", "r3", "r4", 5);
        closed("r4", "r5");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);

        // Q11d
        reorder("t3c", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        mult("t3c+", "r2", "r3", 5);
        closed("r3", "r4");
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);

        // Q11e - disconnected term
        if (disconnected) {
            double t3_t3 = (1.0 / 36.0) * scalar_product("C", "N", "t3c", "t3c");
            copy("prop_pp", "r0");
            closed("r0", "r1");
            update(result_name, t3_t3, "r1");
            restore_stack_pos(pos);
        }

        // Q12a
        reorder("s3c", "r1", "145263");
        mult("r1", "prop_ph", "r2", 2);
        mult("r2", "t1c+", "r3", 2);
        closed("r3", "r4");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);

        // Q12b
        reorder("s3c+", "r1", "145263");
        mult("r1", "prop_hp", "r2", 2);
        mult("r2", "t1c", "r3", 2);
        closed("r3", "r4");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);

        // Q13a
        reorder("s3c", "r1", "145623");
        reorder("prop_pp", "r2", "21");
        mult("r1", "t2c+", "r3", 4);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        update(result_name, 0.25, "r5");
        restore_stack_pos(pos);

        // Q13b
        reorder("s3c+", "r1", "415623");
        mult("r1", "t2c", "r2", 4);
        mult("prop_pp", "r2", "r3", 1);
        closed("r3", "r4");
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);

        // Q13c
        reorder("t2c+", "r1", "3412");
        reorder("s3c", "r2", "142356");
        mult("r1", "prop_pp", "r3", 1);
        mult("r2", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q13d
        reorder("s3c+", "r1", "145623");
        reorder("prop_pp", "r2", "21");
        mult("t2c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q13e
        reorder("prop_hh", "r1", "21");
        reorder("s3c", "r2", "145623");
        mult("t2c+", "r1", "r3", 1);
        mult("r2", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);

        // Q13f
        reorder("s3c+", "r1", "142356");
        reorder("t2c", "r2", "3412");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);

        // Q14a
        reorder("t3c+", "r1", "124563");
        reorder("s3c", "r2", "145623");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, 0.25, "r5");
        restore_stack_pos(pos);

        // Q14b
        reorder("s3c+", "r1", "142356");
        reorder("t3c", "r2", "451263");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, 0.25, "r5");
        restore_stack_pos(pos);

        // Q14c
        reorder("s3c+", "r1", "456123");
        mult("r1", "t3c", "r2", 5);
        mult("prop_ph", "r2", "r3", 1);
        closed("r3", "r4");
        update(result_name, -1.0 / 12.0, "r4");
        restore_stack_pos(pos);

        // Q14d
        reorder("t3c+", "r1", "456123");
        reorder("prop_hp", "r2", "21");
        mult("s3c", "r1", "r3", 5);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        update(result_name, -1.0 / 12.0, "r5");
        restore_stack_pos(pos);

        // Q15a
        reorder("s2c+", "r1", "3412");
        reorder("s3c", "r2", "124563");
        mult("r2", "prop_ph", "r3", 2);
        mult("r3", "r1", "r4", 3);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q15b
        reorder("s3c+", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("s2c", "r2", "r3", 3);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);

        // Q16a
        reorder("s3c+", "r1", "456123");
        mult("r1", "prop_pp", "r2", 1);
        mult("s3c", "r2", "r3", 5);
        closed("r3", "r4");
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);

        // Q16b
        reorder("prop_hh", "r1", "21");
        mult("s3c+", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("s3c", "r3", "r4", 5);
        closed("r4", "r5");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_1h0p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "hh", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    dg_stack_pos_t pos = get_stack_pos();

    // O1
    copy("prop_hh", "r1");
    closed("r1", "r2");
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    reorder("prop_ph", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // L1b
    reorder("t1c+", "r1", "21");
    mult("prop_hp", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // L2a
    reorder("h1c", "r1", "21");
    mult("prop_hh", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // L2b
    reorder("prop_hh", "r1", "21");
    mult("h1c+", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // L3a
    reorder("h2c", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // L3b
    reorder("h2c+", "r1", "1342");
    mult("r1", "prop_hp", "r2", 2);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a
    reorder("t1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q1b
    reorder("t1c+", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_hh", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q1c
    reorder("prop_hh", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q1d - disconnected term
    if (disconnected) {
        double t1_t1 = scalar_product("C", "N", "t1c", "t1c");
        copy("prop_hh", "r0");
        closed("r0", "r1");
        update(result_name, t1_t1, "r1");
        restore_stack_pos(pos);
    }

    // Q2a
    reorder("t2c", "r1", "1342");
    reorder("t1c+", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    closed("r4", "r5");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);

    // Q2b
    reorder("t2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("t1c", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q2c
    reorder("t2c", "r1", "1342");
    reorder("prop_ph", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    closed("r4", "r5");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);

    // Q2d
    reorder("t2c+", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_hp", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q3a - disconnected term
    if (disconnected) {
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        copy("prop_hh", "r0");
        closed("r0", "r1");
        update(result_name, t2_t2, "r1");
        restore_stack_pos(pos);
    }

    // Q3b
    reorder("prop_hh", "r1", "21");
    reorder("t2c+", "r2", "4123");
    mult("r1", "r2", "r3", 1);
    mult("t2c", "r3", "r4", 3);
    closed("r4", "r5");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q3c
    reorder("t2c+", "r1", "3412");
    reorder("t2c", "r2", "2341");
    mult("prop_hh", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    closed("r4", "r5");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q3d
    reorder("prop_hh", "r1", "21");
    mult("t2c+", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c", "r3", "r4", 3);
    closed("r4", "r5");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q3e
    reorder("t2c+", "r1", "3412");
    mult("r1", "prop_pp", "r2", 1);
    mult("t2c", "r2", "r3", 3);
    closed("r3", "r4");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4a
    reorder("prop_ph", "r1", "21");
    reorder("h1c", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    closed("r4", "r5");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);

    // Q4b
    reorder("t1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4c
    reorder("prop_ph", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4d
    reorder("h1c", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_hp", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5a
    reorder("h2c", "r1", "1342");
    reorder("prop_hh", "r2", "21");
    mult("t1c+", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    closed("r4", "r5");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);

    // Q5b
    reorder("h2c+", "r1", "1342");
    reorder("t1c", "r2", "21");
    mult("prop_hh", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    closed("r4", "r5");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);

    // Q5c
    reorder("t1c+", "r1", "21");
    reorder("h2c", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    closed("r4", "r5");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);

    // Q5d
    reorder("h2c+", "r1", "1324");
    reorder("prop_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    closed("r4", "r5");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);

    // Q6a
    reorder("h2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("prop_hp", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6b
    reorder("t2c", "r1", "3412");
    reorder("prop_ph", "r2", "21");
    mult("h2c+", "r1", "r3", 3);
    mult("r3", "r2", "r4", 1);
    closed("r4", "r5");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q6c
    reorder("h2c", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    closed("r4", "r5");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);

    // Q6d
    reorder("h2c+", "r1", "1342");
    reorder("t2c", "r2", "1342");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 2);
    closed("r4", "r5");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);

    // Q7
    reorder("h1c", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q8a
    reorder("h2c", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("h1c+", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q8b
    reorder("h2c+", "r1", "1342");
    reorder("h1c", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    closed("r4", "r5");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);

    // Q8c
    reorder("h2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_hh", "r2", "r3", 1);
    closed("r3", "r4");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q8d
    reorder("h2c+", "r1", "1342");
    reorder("prop_hh", "r2", "21");
    mult("r1", "t1c", "r3", 2);
    mult("r3", "r2", "r4", 1);
    closed("r4", "r5");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);

    // Q9a
    reorder("prop_pp", "r1", "21");
    mult("h2c", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("h2c+", "r3", "r4", 3);
    closed("r4", "r5");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q9b
    reorder("h2c", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("h2c+", "r2", "r3", 3);
    closed("r3", "r4");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_0h2p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c", "t1c_v2", "01", 0);
    restrict_valence("t1c+", "t1c+_v1", "01", 0);
    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);

    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);
    restrict_valence("t2c", "t2c_v2", "0001", 0);
    restrict_valence("t2c+", "t2c+_v2", "0100", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    restrict_valence("s2c", "s2c_v2", "1001", 0);
    restrict_valence("s2c+", "s2c+_v2", "0110", 0);

    restrict_valence("x2c", "x2c_v1", "1110", 0);
    restrict_valence("x2c+", "x2c+_v1", "1011", 0);

    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hp", "prop_hv", "01", 0);

    restrict_valence("prop_pp", "prop_pv", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);
    restrict_valence("prop_pp", "prop_vv", "11", 0);


    dg_stack_pos_t pos = get_stack_pos();

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    printf(" diagram     time (sec)  max mem (mb)  max mem (gb)\n");

    // L1a
    double t_start_L1a = abs_time();
    reorder("s2c_v12", "r1", "1342");
    mult("r1", "prop_vh", "r2", 1);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("L1a", t_start_L1a);

    // L1b
    double t_start_L1b = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("s2c+_v12", "r1", "r2", 1);
    perm("r2", "(34)");
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L1b", t_start_L1b);

    // L2a
    double t_start_L2a = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("x2c_v1", "r1", "r2", 1);
    perm("r2", "(34)");
    update(result_name, +1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L2a", t_start_L2a);

    // L2b
    double t_start_L2b = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    mult("r1", "prop_vp", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("L2b", t_start_L2b);

    if (triples_enabled()) {

        // L3a
        reorder("x3c", "r1", "124563");
        mult("r1", "prop_ph", "r2", 2);
        closed("r2", "r3");
        update(result_name, +1.0, "r3");
        restore_stack_pos(pos);

        // L3b
        reorder("x3c+", "r1", "124563");
        mult("r1", "prop_hp", "r2", 2);
        closed("r2", "r3");
        update(result_name, +1.0, "r3");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1 - disconnected term
    if (disconnected) {
        double t_start_Q1 = abs_time();
        reorder("t1c", "r2", "21");
        mult("t1c+", "r2", "r3", 1);
        closed("r3", "r4");
        construct_disconnected_rank2_rank2("r4", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q1", t_start_Q1);
    }

    // Q2a
    double t_start_Q2a = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("prop_vh", "r1", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    perm("r3", "(12)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2a", t_start_Q2a);

    // Q2b
    double t_start_Q2b = abs_time();
    reorder("t1c_v", "r1", "21");
    reorder("prop_hv", "r2", "21");
    mult("r2", "t2c+_v12", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2b", t_start_Q2b);

    // Q3a - disconnected term
    if (disconnected) {
        double t_start_Q3a = abs_time();
        reorder("t2c_v1", "r1", "3412");
        mult("t2c+_v1", "r1", "r2", 3);
        construct_disconnected_rank2_rank2("r2", "prop_vv", "r3");
        perm("r3", "(12|34)");
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3a", t_start_Q3a);
    }

    // Q3b
    double t_start_Q3b = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("r1", "t2c+_v1", "r2", 2);
    mult("r2", "prop_vp", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update(result_name, 0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3b", t_start_Q3b);

    // Q3c
    double t_start_Q3c = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("t2c_v1", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v12", "r3", "r4", 2);
    perm("r4", "(34)");
    update(result_name, 0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3c", t_start_Q3c);

    // Q3d
    double t_start_Q3d = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("t2c+_v12", "r2", "r3", 2);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3d", t_start_Q3d);

    // Q4a - disconnected term
    if (disconnected) {
        double t_start_Q4a = abs_time();
        reorder("s2c_v1", "r1", "1342");
        mult("r1", "t1c+", "r2", 2);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, +1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q4a", t_start_Q4a);
    }

    // Q4b - disconnected term
    if (disconnected) {
        double t_start_Q4b = abs_time();
        reorder("s2c+_v1", "r1", "1342");
        mult("r1", "t1c", "r2", 2);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, +1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q4b", t_start_Q4b);
    }

    // Q4c
    double t_start_Q4c = abs_time();
    reorder("t1c+", "r1", "21");
    reorder("s2c_v12", "r2", "3412");
    mult("prop_vp", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    closed("r5", "r6");
    perm("r6", "(12)");
    update(result_name, -1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4c", t_start_Q4c);

    // Q4d
    double t_start_Q4d = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("s2c+_v12", "r2", "r3", 1);
    closed("r3", "r4");
    perm("r4", "(34)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4d", t_start_Q4d);

    // Q4e
    double t_start_Q4e = abs_time();
    reorder("s2c_v2", "r1", "1423");
    reorder("prop_pv", "r2", "21");
    mult("r2", "r1", "r3", 1);
    mult("r3", "t1c+_v", "r4", 1);
    reorder("r4", "r5", "2413");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4e", t_start_Q4e);

    // Q4f
    double t_start_Q4f = abs_time();
    reorder("s2c+_v2", "r1", "2341");
    reorder("t1c_v", "r2", "21");
    mult("prop_vp", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(12|34)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4f", t_start_Q4f);

    // Q4g
    double t_start_Q4g = abs_time();
    reorder("s2c_v12", "r1", "3412");
    reorder("prop_hh", "r2", "21");
    mult("t1c+_v", "r2", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4g", t_start_Q4g);

    // Q4h
    double t_start_Q4h = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("s2c+_v12", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4h", t_start_Q4h);

    // Q5a
    double t_start_Q5a = abs_time();
    reorder("t2c_v1", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("s2c+_v12", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5a", t_start_Q5a);

    // Q5b
    double t_start_Q5b = abs_time();
    reorder("s2c_v12", "r1", "3412");
    reorder("t2c+_v1", "r2", "1342");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5b", t_start_Q5b);

    // Q5c
    double t_start_Q5c = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    reorder("t2c_v2", "r2", "4213");
    mult("r1", "r2", "r3", 2);
    mult("r3", "prop_vh", "r4", 1);
    reorder("r4", "r5", "1423");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5c", t_start_Q5c);

    // Q5d
    double t_start_Q5d = abs_time();
    reorder("prop_hv", "r1", "21");
    reorder("t2c+_v2", "r2", "2314");
    reorder("s2c_v1", "r3", "1324");
    mult("r1", "r2", "r4", 1);
    mult("r3", "r4", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(12|34)");
    update(result_name, -1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5d", t_start_Q5d);

    // Q6a
    double t_start_Q6a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6a", t_start_Q6a);

    // Q6b
    double t_start_Q6b = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    mult("r1", "prop_hp", "r2", 1);
    mult("r2", "t1c+_v", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6b", t_start_Q6b);

    // Q6c
    double t_start_Q6c = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6c", t_start_Q6c);

    // Q6d
    double t_start_Q6d = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    reorder("t1c", "r2", "21");
    mult("prop_vh", "r2", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6d", t_start_Q6d);

    // Q7 - disconnected term
    if (disconnected) {
        double t_start_Q7 = abs_time();
        reorder("s1c+", "r1", "21");
        mult("s1c", "r1", "r2", 1);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7", t_start_Q7);
    }

    // Q8a
    double t_start_Q8a = abs_time();
    reorder("prop_ph", "r1", "21");
    reorder("s2c_v12", "r2", "3412");
    mult("s1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8a", t_start_Q8a);

    // Q8b
    double t_start_Q8b = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("s2c+_v12", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8b", t_start_Q8b);

    // Q8c
    double t_start_Q8c = abs_time();
    reorder("s1c+", "r1", "21");
    reorder("s2c_v2", "r2", "4123");
    mult("r1", "r2", "r3", 1);
    mult("r3", "prop_vh", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8c", t_start_Q8c);

    // Q8d
    double t_start_Q8d = abs_time();
    reorder("s2c+_v2", "r1", "2341");
    reorder("prop_hv", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(12|34)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8d", t_start_Q8d);

    // Q9a - disconnected term
    if (disconnected) {
        double t_start_Q9a = abs_time();
        reorder("s2c+", "r1", "3412");
        mult("s2c", "r1", "r2", 3);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, +0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9a", t_start_Q9a);
    }

    // Q9b
    double t_start_Q9b = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r2", "2134");
    reorder("prop_pv", "r3", "21");
    mult("r3", "r2", "r4", 1);
    mult("r1", "r4", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(12|34)");
    update(result_name, 1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9b", t_start_Q9b);

    // Q9c
    double t_start_Q9c = abs_time();
    reorder("s2c_v1", "r1", "1324");
    reorder("s2c+", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "4312");
    mult("prop_vp", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9c", t_start_Q9c);

    // Q9d
    double t_start_Q9d = abs_time();
    reorder("s2c_v1", "r1", "1342");
    reorder("s2c+_v1", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "2413");
    mult("r1", "prop_hh", "r3", 1);
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9d", t_start_Q9d);

    // Q9e
    double t_start_Q9e = abs_time();
    reorder("s2c+_v1", "s2c+_21", "2143");
    reorder("s2c+_21", "r1", "2431");
    reorder("s2c_v1", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9e", t_start_Q9e);

    // Q10a
    double t_start_Q10a = abs_time();
    reorder("prop_pv", "r1", "21");
    reorder("s1c+", "r2", "21");
    mult("r1", "x2c", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10a", t_start_Q10a);

    // Q10b
    double t_start_Q10b = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("prop_vp", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    perm("r3", "(12)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10b", t_start_Q10b);

    // Q10c
    double t_start_Q10c = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10c", t_start_Q10c);

    // Q10d
    double t_start_Q10d = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10d", t_start_Q10d);

    // Q11a
    double t_start_Q11a = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    reorder("s2c", "r2", "1342");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11a", t_start_Q11a);

    // Q11b
    double t_start_Q11b = abs_time();
    reorder("s2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11b", t_start_Q11b);

    // Q11c
    double t_start_Q11c = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("s2c", "r1", "r2", 2);
    reorder("r2", "r3", "3412");
    mult("r3", "prop_vh", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11c", t_start_Q11c);

    // Q11d
    double t_start_Q11d = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("s2c+", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("x2c", "r3", "r4", 2);
    perm("r4", "(34)");
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11d", t_start_Q11d);

    // Q12
    double t_start_Q12 = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("r1", "prop_pp", "r2", 1);
    mult("x2c", "r2", "r3", 2);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q12", t_start_Q12);

    if (triples_enabled()) {

        // Q13a
        reorder("t3c", "r1", "451263");
        mult("r1", "prop_ph", "r2", 2);
        mult("t2c+", "r2", "r3", 2);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);

        // Q13b
        reorder("t2c", "r1", "3412");
        reorder("t3c+", "r2", "124563");
        mult("r2", "prop_hp", "r3", 2);
        mult("r3", "r1", "r4", 2);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q13c
        reorder("t3c", "r1", "451623");
        mult("t2c+", "r1", "r2", 3);
        mult("prop_ph", "r2", "r3", 1);
        closed("r3", "r4");
        perm("r4", "(12)");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);

        // Q13d
        reorder("t2c", "r1", "3412");
        reorder("t3c+", "r2", "124356");
        reorder("prop_hp", "r3", "21");
        mult("r1", "r2", "r4", 3);
        mult("r3", "r4", "r5", 1);
        reorder("r5", "r6", "3412");
        closed("r6", "r7");
        perm("r7", "(34)");
        update(result_name, 0.5, "r7");
        restore_stack_pos(pos);

        // Q14a
        reorder("t3c", "r1", "456123");
        reorder("prop_pp", "r2", "21");
        mult("t3c+", "r1", "r3", 4);
        reorder("r3", "r4", "1243");
        mult("r4", "r2", "r5", 1);
        reorder("r5", "r6", "1243");
        closed("r6", "r7");
        perm("r7", "(34)");
        update(result_name, 1.0 / 6.0, "r7");
        restore_stack_pos(pos);

        // Q14b
        reorder("t3c", "r1", "456123");
        mult("t3c+", "r1", "r2", 4);
        reorder("r2", "r3", "2341");
        mult("prop_pp", "r3", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, 1.0 / 6.0, "r5");
        restore_stack_pos(pos);

        // Q14c
        reorder("prop_pp", "r1", "21");
        mult("t3c", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("t3c+", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, 1.0 / 6.0, "r5");
        restore_stack_pos(pos);

        // Q14d
        reorder("t3c", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        mult("t3c+", "r2", "r3", 4);
        closed("r3", "r4");
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);

        // Q14e - disconnected term
        if (disconnected) {
            reorder("t3c", "r1", "456123");
            mult("t3c+", "r1", "r2", 5);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, -1.0 / 12.0, "r5");
            restore_stack_pos(pos);
        }

        // Q15a
        reorder("s3c", "r1", "145263");
        mult("r1", "t1c+", "r2", 2);
        mult("r2", "prop_ph", "r3", 1);
        reorder("r3", "r4", "1423");
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);

        // Q15b
        reorder("s3c+", "r1", "124563");
        reorder("prop_hp", "r2", "21");
        mult("r1", "t1c", "r3", 2);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);

        // Q15c
        reorder("s3c", "r1", "145263");
        mult("r1", "prop_ph", "r2", 2);
        mult("r2", "t1c+", "r3", 1);
        reorder("r3", "r4", "1423");
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);

        // Q15d
        reorder("s3c+", "r1", "124563");
        reorder("t1c", "r2", "21");
        mult("r1", "prop_hp", "r3", 2);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);

        // Q16a
        reorder("s3c", "r1", "145623");
        reorder("t2c+", "r2", "3412");
        mult("r2", "prop_pp", "r3", 1);
        reorder("r3", "r4", "3412");
        mult("r1", "r4", "r5", 3);
        reorder("r5", "r6", "1423");
        closed("r6", "r7");
        perm("r7", "(12)");
        update(result_name, -0.5, "r7");
        restore_stack_pos(pos);

        // Q16b
        reorder("s3c+", "r1", "124356");
        reorder("prop_pp", "r2", "21");
        mult("t2c", "r2", "r3", 1);
        reorder("r3", "r4", "3412");
        mult("r1", "r4", "r5", 3);
        closed("r5", "r6");
        perm("r6", "(34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);

        // Q16c
        reorder("prop_hh", "r1", "21");
        reorder("s3c", "r2", "145623");
        mult("t2c+", "r1", "r3", 1);
        mult("r2", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        closed("r5", "r6");
        perm("r6", "(12)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);

        // Q16d
        reorder("s3c+", "r1", "124356");
        reorder("t2c", "r2", "3412");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 3);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q16e
        reorder("s3c", "r1", "154623");
        reorder("prop_pp", "r2", "21");
        mult("r1", "t2c+", "r3", 3);
        reorder("r3", "r4", "1243");
        mult("r4", "r2", "r5", 1);
        reorder("r5", "r6", "1342");
        closed("r6", "r7");
        perm("r7", "(12|34)");
        update(result_name, -0.5, "r7");
        restore_stack_pos(pos);

        // Q16f
        reorder("s3c+", "r1", "241356");
        reorder("t2c", "r2", "3412");
        mult("r1", "r2", "r3", 3);
        reorder("r3", "r4", "1243");
        mult("prop_pp", "r4", "r5", 1);
        closed("r5", "r6");
        perm("r6", "(12|34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);

        // Q16g - disconnected term
        if (disconnected) {
            reorder("s3c", "r1", "145623");
            mult("r1", "t2c+", "r2", 4);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, 0.25, "r5");
            restore_stack_pos(pos);
        }

        // Q16h - disconnected term
        if (disconnected) {
            reorder("s3c+", "r1", "145623");
            mult("r1", "t2c", "r2", 4);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, 0.25, "r5");
            restore_stack_pos(pos);
        }

        // Q17a
        reorder("t3c+", "r1", "124563");
        reorder("s3c", "r2", "145623");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        closed("r5", "r6");
        perm("r6", "(12)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);

        // Q17b
        reorder("s3c+", "r1", "124356");
        reorder("t3c", "r2", "451263");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 3);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);

        // Q17c
        reorder("t3c+", "r1", "245613");
        reorder("s3c", "r2", "152346");
        reorder("prop_hp", "r3", "21");
        mult("r2", "r1", "r4", 4);
        mult("r4", "r3", "r5", 1);
        reorder("r5", "r6", "1342");
        closed("r6", "r7");
        perm("r7", "(12|34)");
        update(result_name, 0.25, "r7");
        restore_stack_pos(pos);

        // Q17d
        reorder("s3c+", "r1", "245613");
        reorder("t3c", "r2", "512346");
        mult("r1", "r2", "r3", 4);
        mult("prop_ph", "r3", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(12|34)");
        update(result_name, 0.25, "r5");
        restore_stack_pos(pos);

        // Q18a
        reorder("s3c", "r1", "142563");
        reorder("s2c+", "r2", "1342");
        mult("r1", "prop_ph", "r3", 2);
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1324");
        closed("r5", "r6");
        perm("r6", "(12|34)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);

        // Q18b
        reorder("s3c+", "r1", "145263");
        reorder("s2c", "r2", "1324");
        mult("r1", "prop_hp", "r3", 2);
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1324");
        closed("r5", "r6");
        perm("r6", "(12|34)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);

        // Q18c
        reorder("s3c", "r1", "152346");
        reorder("s2c+", "r2", "3412");
        mult("r1", "r2", "r3", 3);
        reorder("r3", "r4", "4213");
        mult("r4", "prop_ph", "r5", 1);
        reorder("r5", "r6", "3412");
        closed("r6", "r7");
        perm("r7", "(12|34)");
        update(result_name, -0.5, "r7");
        restore_stack_pos(pos);

        // Q18d
        reorder("s3c+", "r1", "245613");
        reorder("prop_hp", "r2", "21");
        mult("r1", "s2c", "r3", 3);
        reorder("r3", "r4", "4123");
        mult("r4", "r2", "r5", 1);
        closed("r5", "r6");
        perm("r6", "(12|34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);

        // Q19a
        reorder("prop_pp", "r1", "21");
        reorder("s3c+", "r2", "245613");
        reorder("s3c", "r3", "152346");
        mult("r2", "r3", "r4", 4);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "3124");
        closed("r6", "r7");
        perm("r7", "(12|34)");
        update(result_name, -0.25, "r7");
        restore_stack_pos(pos);

        // Q19b
        reorder("s3c+", "r1", "425613");
        reorder("s3c", "r2", "152346");
        mult("r2", "r1", "r3", 4);
        mult("r3", "prop_pp", "r4", 1);
        reorder("r4", "r5", "1432");
        closed("r5", "r6");
        perm("r6", "(12|34)");
        update(result_name, -0.25, "r6");
        restore_stack_pos(pos);

        // Q19c
        reorder("s3c", "r1", "142356");
        reorder("s3c+", "r2", "145623");
        reorder("prop_pp", "r3", "21");
        mult("r1", "r3", "r4", 1);
        mult("r4", "r2", "r5", 4);
        reorder("r5", "r6", "1324");
        closed("r6", "r7");
        perm("r7", "(12|34)");
        update(result_name, 0.5, "r7");
        restore_stack_pos(pos);

        // Q19d
        reorder("s3c", "r1", "145623");
        reorder("prop_hh", "r2", "21");
        mult("s3c+", "r2", "r3", 1);
        reorder("r3", "r4", "142356");
        mult("r1", "r4", "r5", 4);
        reorder("r5", "r6", "1324");
        closed("r6", "r7");
        perm("r7", "(12|34)");
        update(result_name, -0.5, "r7");
        restore_stack_pos(pos);

        // Q19e - disconnected term
        if (disconnected) {
            reorder("s3c+", "r1", "456123");
            mult("s3c", "r1", "r2", 5);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, +1.0 / 12.0, "r5");
            restore_stack_pos(pos);
        }

        // Q20a
        reorder("x3c", "r1", "124356");
        reorder("s3c+", "r2", "451263");
        mult("r2", "prop_hp", "r3", 2);
        mult("r1", "r3", "r4", 3);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q20b
        reorder("x3c+", "r1", "451623");
        reorder("s3c", "r2", "124563");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 3);
        reorder("r4", "r5", "3412");
        closed("r5", "r6");
        perm("r6", "(12)");
        update(result_name, 0.5, "r6");
        restore_stack_pos(pos);

        // Q20c
        reorder("s3c+", "r1", "456123");
        reorder("prop_hp", "r2", "21");
        mult("x3c", "r1", "r3", 4);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);

        // Q20d
        reorder("x3c+", "r1", "456123");
        mult("r1", "s3c", "r2", 4);
        mult("r2", "prop_ph", "r3", 1);
        reorder("r3", "r4", "3412");
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);

        // Q21a
        reorder("x3c", "r1", "124563");
        reorder("prop_pp", "r2", "21");
        mult("r1", "t1c+", "r3", 2);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q21b
        reorder("x3c+", "r1", "145263");
        mult("r1", "t1c", "r2", 2);
        mult("r2", "prop_pp", "r3", 1);
        reorder("r3", "r4", "1423");
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q21c
        reorder("prop_hh", "r1", "21");
        reorder("x3c", "r2", "124536");
        mult("r1", "t1c+", "r3", 1);
        mult("r2", "r3", "r4", 2);
        closed("r4", "r5");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);

        // Q21d
        reorder("x3c+", "r1", "124536");
        reorder("t1c", "r2", "21");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 2);
        closed("r4", "r5");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);

        // Q21e
        reorder("x3c", "r1", "124536");
        reorder("t1c+", "r2", "21");
        mult("r2", "prop_pp", "r3", 1);
        mult("r1", "r3", "r4", 2);
        closed("r4", "r5");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q21f
        reorder("x3c+", "r1", "124563");
        reorder("prop_pp", "r2", "21");
        mult("t1c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 2);
        closed("r4", "r5");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q22a
        reorder("x3c", "r1", "124536");
        reorder("t2c+", "r2", "3142");
        mult("r2", "prop_hp", "r3", 2);
        mult("r1", "r3", "r4", 2);
        closed("r4", "r5");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q22b
        reorder("x3c+", "r1", "124563");
        reorder("t2c", "r2", "1342");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 2);
        closed("r4", "r5");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q23a
        reorder("x3c", "r1", "124563");
        reorder("s1c+", "r2", "21");
        mult("r1", "prop_ph", "r3", 2);
        mult("r3", "r2", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q23b
        reorder("x3c+", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "3412");
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q24a
        reorder("x3c", "r1", "124356");
        reorder("s2c+", "r2", "3412");
        mult("r2", "prop_pp", "r3", 1);
        mult("r1", "r3", "r4", 3);
        closed("r4", "r5");
        perm("r5", "(34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);

        // Q24b
        reorder("x3c+", "r1", "145623");
        reorder("prop_pp", "r2", "21");
        mult("s2c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        closed("r5", "r6");
        perm("r6", "(12)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);

        // Q24c
        reorder("x3c", "r1", "124356");
        reorder("prop_hh", "r2", "21");
        mult("s2c+", "r2", "r3", 1);
        reorder("r3", "r4", "3412");
        mult("r1", "r4", "r5", 3);
        closed("r5", "r6");
        perm("r6", "(34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);

        // Q24d
        reorder("x3c+", "r1", "145236");
        reorder("s2c", "r2", "1342");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        closed("r5", "r6");
        perm("r6", "(12)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);

        // Q24e
        reorder("x3c", "r1", "124356");
        reorder("s2c+", "r2", "3412");
        reorder("prop_pp", "r3", "21");
        mult("r1", "r2", "r4", 3);
        reorder("r4", "r5", "4123");
        mult("r3", "r5", "r6", 1);
        reorder("r6", "r7", "3412");
        closed("r7", "r8");
        perm("r8", "(34)");
        update(result_name, 0.5, "r8");
        restore_stack_pos(pos);

        // Q24f
        reorder("x3c+", "r1", "145623");
        mult("r1", "s2c", "r2", 3);
        reorder("r2", "r3", "4231");
        mult("prop_pp", "r3", "r4", 1);
        closed("r4", "r5");
        perm("r5", "(12)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q25a
        reorder("x2c+", "r1", "3412");
        reorder("x3c", "r2", "124563");
        mult("r2", "prop_ph", "r3", 2);
        mult("r3", "r1", "r4", 2);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);

        // Q25b
        reorder("x3c+", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("x2c", "r2", "r3", 2);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);

        // Q26a
        reorder("x3c+", "r1", "456123");
        mult("r1", "prop_pp", "r2", 1);
        mult("x3c", "r2", "r3", 4);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);

        // Q26b
        reorder("prop_hh", "r1", "21");
        mult("x3c+", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("x3c", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_0h0p_1h1p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "ph", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c+", "t1c+_v", "10", 0);
    restrict_valence("t1c+", "t1c+_g", "01", 0);
    restrict_valence("t1c+", "t1c+_vg", "11", 0);
    restrict_valence("t1c", "t1c_gv", "11", 0);

    restrict_valence("t2c+", "t2c+_vg", "1010", 0);
    restrict_valence("t2c+", "t2c+_g", "0010", 0);
    restrict_valence("t2c+", "t2c+_v", "1000", 0);

    restrict_valence("prop_ph", "prop_vg", "11", 0);
    restrict_valence("prop_ph", "prop_pg", "01", 0);
    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hh", "prop_hg", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);

    dg_stack_pos_t pos = get_stack_pos();

    // O1
    update(result_name, 1.0, "prop_vg");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    reorder("prop_pg", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L1b
    reorder("h1c", "r1", "21");
    mult("prop_vh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2a
    reorder("prop_hg", "r1", "21");
    mult("t1c+_v", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2b
    reorder("t1c+_g", "r1", "21");
    mult("prop_vp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3a
    reorder("t2c+_vg", "r1", "1342");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3b
    reorder("e2c", "r1", "1432");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a - disconnected term
    if (disconnected) {
        double factor_Q1a = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q1a, "t1c+_vg");
    }

    // Q1b
    reorder("t1c+_g", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q2a - disconnected term
    if (disconnected) {
        double factor_Q2a = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q2a, "e1c");
    }

    // Q2b
    reorder("h1c", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3a - disconnected term
    if (disconnected) {
        double t1_t1 = scalar_product("C", "N", "t1c", "t1c");
        update(result_name, t1_t1, "prop_vg");
    }

    // Q3b - disconnected term
    if (disconnected) {
        double factor_Q3b = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q3b, "t1c+_vg");
    }

    // Q3c - disconnected term
    if (disconnected) {
        double factor_Q3c = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q3c, "e1c");
    }

    // Q3d
    reorder("prop_pg", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3e
    reorder("t1c+_g", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3f
    reorder("prop_hg", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3g
    reorder("h1c", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3h
    reorder("t1c+_g", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q3i
    reorder("h1c", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4a
    reorder("s2c", "r1", "1342");
    reorder("prop_pg", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4b
    reorder("h2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4c
    reorder("s2c", "r1", "1342");
    reorder("t1c+_g", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4d
    reorder("h2c", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4e
    reorder("e2c", "r1", "1432");
    mult("r1", "prop_hh", "r2", 1);
    mult("r2", "t1c+", "r3", 2);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4f
    reorder("e2c", "r1", "1432");
    mult("r1", "t1c+", "r2", 1);
    mult("r2", "prop_pp", "r3", 2);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5a
    reorder("t2c+_g", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q5b
    reorder("prop_hg", "r1", "21");
    reorder("t2c+_v", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5c
    reorder("t2c+_g", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q5d
    reorder("t2c+_v", "r1", "1342");
    reorder("h1c", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5e
    reorder("t1c", "r1", "21");
    reorder("t2c+_vg", "r2", "1342");
    mult("prop_hh", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5f
    reorder("t2c+_vg", "r1", "1324");
    reorder("prop_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6a - disconnected term
    if (disconnected) {
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        update(result_name, t2_t2, "prop_vg");
        restore_stack_pos(pos);
    }

    // Q6b
    reorder("t2c+_vg", "r1", "1324");
    reorder("t2c", "r2", "3142");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6c
    reorder("t2c", "r1", "2341");
    reorder("t2c+_g", "r0", "2143");
    reorder("r0", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("prop_vh", "r3", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6d
    reorder("prop_pg", "r0", "21");
    reorder("t2c", "r1", "4123");
    mult("r0", "r1", "r2", 1);
    mult("t2c+_v", "r2", "r3", 3);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6e
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("prop_hg", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q6f
    reorder("h2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6g
    reorder("t2c+_g", "r1", "3124");
    reorder("s2c", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6h
    reorder("h2c", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("t2c+_v", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q6i
    reorder("t2c+_g", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("s2c", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6k
    reorder("prop_pp", "r1", "21");
    mult("h2c", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6l
    reorder("e2c", "r1", "1423");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    if (triples_enabled()) {

    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_1h1p_0h0p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "hp", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c", "t1c_g", "10", 0);
    restrict_valence("t1c", "t1c_gv", "11", 0);

    restrict_valence("t2c", "t2c_gv", "1010", 0);
    restrict_valence("t2c", "t2c_g", "1000", 0);
    restrict_valence("t2c", "t2c_v", "0010", 0);

    restrict_valence("prop_hh", "prop_gh", "10", 0);
    restrict_valence("prop_hp", "prop_gv", "11", 0);
    restrict_valence("prop_hp", "prop_gp", "10", 0);

    dg_stack_pos_t pos = get_stack_pos();

    // O1
    update(result_name, 1.0, "prop_gv");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    reorder("t1c_v", "r1", "21");
    mult("prop_gh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L1b
    reorder("prop_pv", "r1", "21");
    mult("t1c_g", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L2a
    reorder("prop_hv", "r1", "21");
    mult("h1c+", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2b
    reorder("s1c+", "r1", "21");
    mult("prop_gp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3a
    reorder("t2c_gv", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3b
    reorder("e2c+", "r1", "2341");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a - disconnected term
    if (disconnected) {
        double factor_Q1a = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q1a, "t1c_gv");
    }

    // Q1b
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q2a - disconnected term
    if (disconnected) {
        double factor_Q2a = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q2a, "e1c+");
    }

    // Q2b
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3a - disconnected term
    if (disconnected) {
        double t1_t1_Q3a = scalar_product("C", "N", "t1c", "t1c");
        update(result_name, t1_t1_Q3a, "prop_gv");
        restore_stack_pos(pos);
    }

    // Q3b - disconnected term
    if (disconnected) {
        double factor_Q3b = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q3b, "t1c_gv");
    }

    // Q3c - disconnected term
    if (disconnected) {
        double factor_Q3c = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q3c, "e1c+");
    }

    // Q3d
    reorder("t1c_v", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3e
    reorder("prop_hv", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3f
    reorder("s1c+", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_gh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3g
    reorder("prop_pv", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3h
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q3i
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4a
    reorder("s2c+", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4b
    reorder("prop_hv", "r1", "21");
    reorder("h2c+", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4c
    reorder("s2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4d
    reorder("h2c+", "r1", "1342");
    reorder("t1c_v", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4e
    reorder("t1c", "r1", "21");
    reorder("e2c+", "r2", "2341");
    mult("prop_hh", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4f
    reorder("e2c+", "r1", "2341");
    mult("r1", "prop_pp", "r2", 1);
    mult("r2", "t1c", "r3", 2);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5a
    reorder("t2c_g", "r1", "1342");
    reorder("prop_pv", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q5b
    reorder("t2c_v", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_gh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5c
    reorder("t2c_g", "r1", "1342");
    reorder("s1c+", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q5d
    reorder("t2c_v", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5e
    reorder("prop_hh", "r1", "21");
    reorder("t2c_gv", "r2", "1324");
    mult("r1", "t1c+", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5f
    reorder("t1c+", "r1", "21");
    reorder("t2c_gv", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6a - disconnected term
    if (disconnected) {
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        update(result_name, t2_t2, "prop_gv");
        restore_stack_pos(pos);
    }

    // Q6b
    reorder("t2c_gv", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6c
    reorder("t2c", "r1", "2341");
    reorder("s2c+", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("prop_gh", "r3", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6d
    reorder("t2c", "r1", "3412");
    reorder("prop_pv", "r2", "21");
    mult("h2c+", "r1", "r3", 3);
    mult("r3", "r2", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6e
    reorder("t2c_g", "t2c_21", "2143");
    reorder("t2c_21", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("prop_hv", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q6f
    reorder("t2c_v", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6g
    reorder("s2c+", "r1", "3124");
    reorder("t2c_g", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6h
    reorder("t2c_v", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("h2c+", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q6i
    reorder("s2c+", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("t2c_g", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6k
    reorder("prop_pp", "r1", "21");
    mult("t2c_v", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("h2c+", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6l
    reorder("t2c", "r1", "1342");
    reorder("e2c+", "r2", "2341");
    mult("r1", "prop_ph", "r3", 2);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    if (triples_enabled()) {

    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_1h1p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "phph", "1111", "1234", NOT_PERM_UNIQUE, operator_symmetry);

    /*restrict_valence("t1c", "t1c_v2", "01", 0);
    restrict_valence("t1c+", "t1c+_v1", "01", 0);
    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);

    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);
    restrict_valence("t2c", "t2c_v2", "0001", 0);
    restrict_valence("t2c+", "t2c+_v2", "0100", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    restrict_valence("s2c", "s2c_v2", "1001", 0);
    restrict_valence("s2c+", "s2c+_v2", "0110", 0);

    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hp", "prop_hv", "01", 0);

    restrict_valence("prop_pp", "prop_pv", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);
    restrict_valence("prop_pp", "prop_vv", "11", 0);*/

    dg_stack_pos_t pos = get_stack_pos();


    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }


    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }


    if (triples_enabled()) {

    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}

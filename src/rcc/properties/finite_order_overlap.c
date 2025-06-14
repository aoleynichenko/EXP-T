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
 * Direct calculation of norms and overlaps of correlated wavefunctions
 * in different Fock space sector.
 * The formula <psi_I|psi_J> = < phi_I | {e^T^\dagger} {e^T} | phi_J > is used,
 * where |psi> are exact wavefunctions and |phi> are model vectors.
 * Infinite expansion is truncated at quadratic terms (T^2).
 */

#include "finite_order_overlap.h"

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


void print_time_mem_usage(char *diagram_name, double t_start, int print_level);

void msprop_transform_slater_to_model_vectors(size_t nroots_i, size_t nroots_f,
                                              size_t ms_size_i, size_t ms_size_f,
                                              double complex *coefs_bra, double complex *coefs_ket,
                                              double complex *prop_slater, double complex *prop_model);

void construct_property_matrix_slater(
        int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p,
        size_t n_det_bra, slater_det_t *bra_dets, size_t n_det_ket, slater_det_t *ket_dets,
        double complex *prp_spinor, double complex *prp_slater
);

void print_property_matrix(
        int nroots_i, int nroots_f,
        double *energies_i, double *energies_f,
        char *rep_name_i, char *rep_name_f,
        double complex *prop_matrix);

void conjugate_cluster_amplitudes_0h0p();

void conjugate_cluster_amplitudes_1h0p();

void conjugate_cluster_amplitudes_0h1p();

void conjugate_cluster_amplitudes_2h0p();

void conjugate_cluster_amplitudes_0h2p();

void conjugate_cluster_amplitudes_1h1p();

void construct_cc_wfn_overlap_matrix(int sect_h, int sect_p, int dim, slater_det_t *det_list,
                                     double complex *overlap_slater, double scalar_term,
                                     int n_diagrams, char **list_diagrams);


void calculate_wavefunction_norms_and_overlaps(int sect_h, int sect_p)
{
    int n_overlap_diagrams = 0;
    char *overlap_diagrams[10];

    if (cc_opts->calc_overlap[sect_h][sect_p] == 0) {
        return;
    }

    printf("\n");
    printf(" **\n");
    printf(" ** Overlap of electronic wavefunctions\n");
    printf(" ** Approximation: quadratic\n");
    printf(" **\n");
    printf("\n");

    /*
     * get complex conjugate for all cluster operators
     */
    conjugate_cluster_amplitudes_0h0p();
    if (sect_h >= 1) {
        conjugate_cluster_amplitudes_1h0p();
    }
    if (sect_p >= 1) {
        conjugate_cluster_amplitudes_0h1p();
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

    /*
     * scalar part:
     * 1 + fully contracted diagrams without any external lines
     */
    printf(" construction of the overlap diagram (0h0p) ...\n");
    double complex scalar_part = sector_0h0p_overlap(CC_PROPERTIES_APPROX_QUADRATIC);

    /*
     * one-electron part
     */
    if (sect_h >= 1) {
        printf(" construction of the overlap diagram (1h0p) ...\n");
        sector_1h0p_overlap("overlap_1h0p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h0p");
    }
    if (sect_p >= 1) {
        printf(" construction of the overlap diagram (0h1p) ...\n");
        sector_0h1p_overlap("overlap_0h1p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h1p");
    }

    /*
     * two-electron part
     */
    if (sect_h >= 2) {
        printf(" construction of the overlap diagram (2h0p) ...\n");
        sector_2h0p_overlap("overlap_2h0p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_2h0p");
    }
    if (sect_h >=1 && sect_p >= 1) {
        printf(" construction of the overlap diagram (1h1p) ...\n");
        sector_1h1p_overlap("overlap_1h1p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h1p");
    }
    if (sect_p >= 2) {
        printf(" construction of the overlap diagram (0h2p) ...\n");
        sector_0h2p_overlap("overlap_0h2p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h2p");
    }

    /*
     * three-electron part
     * (only for sectors with three quasiparticles)
     */
    if (sect_h >= 3) {
        printf(" construction of the overlap diagram (3h0p) ...\n");
        sector_3h0p_overlap("overlap_3h0p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_3h0p");
    }
    if (sect_p >= 3) {
        printf(" construction of the overlap diagram (0h3p) ...\n");
        sector_0h3p_overlap("overlap_0h3p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h3p");
    }

    printf(" scalar term: %.6f %.6fi\n", creal(scalar_part), cimag(scalar_part));
    printf(" diagrams:\n");
    for (int i = 0; i < n_overlap_diagrams; i++) {
        printf(" - %s\n", overlap_diagrams[i]);
    }
    printf("\n");

    /*
     * read model vectors from disk
     */
    int n_irreps = 0;
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(sect_h, sect_p, NULL, &n_irreps, mv_blocks);

    for (int irrep = 0; irrep < n_irreps; irrep++) {

        struct mv_block *block = mv_blocks + irrep;

        size_t ms_size = block->ms_size;
        size_t n_roots = block->nroots;
        char *rep_name = block->rep_name;
        slater_det_t *det_list = block->dets;
        double *energies_cm = block->energy_cm;
        double complex *vectors_left = block->vl;
        double complex *vectors_right = block->vr;

        printf(" Irrep %s\n", block->rep_name);

        double complex *overlap_slater = x_zeros(CC_COMPLEX, ms_size, ms_size);
        double complex *overlap = x_zeros(CC_COMPLEX, n_roots, n_roots);

        // obtain overlap matrix in the basis of Slater determinants
        construct_cc_wfn_overlap_matrix(sect_h, sect_p, ms_size, det_list, overlap_slater,
                                        scalar_part, n_overlap_diagrams, overlap_diagrams);

        // construct property matrix in the basis of model vectors
        msprop_transform_slater_to_model_vectors(
            n_roots, n_roots, ms_size, ms_size,
            vectors_right, vectors_right, overlap_slater, overlap
        );

        // print smart table of matrix elements
        print_property_matrix(
                n_roots, n_roots, energies_cm, energies_cm,
                rep_name, rep_name, overlap
        );

        cc_free(overlap_slater);
        cc_free(overlap);
    }

    /*
     * clean up
     */
    for (size_t irep = 0; irep < n_irreps; irep++) {
        mvblock_free(mv_blocks + irep);
    }
}


/**
 * calculate norm of the wavefunction for the given electronic state:
 * N = <psi| \Omega^+ \Omega |psi>^{1/2}
 */
double calculate_wavefunction_norm(int sect_h, int sect_p, char *irrep_name, int state_number)
{
    int n_overlap_diagrams = 0;
    char *overlap_diagrams[10];

    /*
     * load model vector from disk
     */
    int model_space_dim = 0;
    double eigenvalue = 0.0;
    double exc_energy_cm = 0.0;
    slater_det_t *det_list = NULL;
    double complex *coef_left = NULL;
    double complex *coef_right = NULL;
    char mvcoef_file_name[20];

    sprintf(mvcoef_file_name, "MVCOEF%d%d", sect_h, sect_p);
    mvcoef_read_vectors_unformatted_for_state(
        sect_h, sect_p, mvcoef_file_name, irrep_name, state_number,
        &model_space_dim, &det_list, &eigenvalue, &exc_energy_cm, &coef_left, &coef_right
    );

    /*
     * get complex conjugate for all cluster operators
     */
    conjugate_cluster_amplitudes_0h0p();
    if (sect_h >= 1) {
        conjugate_cluster_amplitudes_1h0p();
    }
    if (sect_p >= 1) {
        conjugate_cluster_amplitudes_0h1p();
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

    /*
     * scalar part:
     * 1 + fully contracted diagrams without any external lines
     */
    double complex scalar_part = sector_0h0p_overlap(CC_PROPERTIES_APPROX_QUADRATIC);

    /*
     * one-electron part
     */
    if (sect_h >= 1) {
        sector_1h0p_overlap("overlap_1h0p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h0p");
    }
    if (sect_p >= 1) {
        sector_0h1p_overlap("overlap_0h1p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h1p");
    }

    /*
     * two-electron part
     */
    if (sect_h >=1 && sect_p >= 1) {
        sector_1h1p_overlap("overlap_1h1p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_1h1p");
    }
    if (sect_h >= 2) {
        sector_2h0p_overlap("overlap_2h0p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_2h0p");
    }
    if (sect_p >= 2) {
        sector_0h2p_overlap("overlap_0h2p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h2p");
    }

    /*
     * three-electron part
     * (only for sectors with three quasiparticles)
     */
    if (sect_p >= 3) {
        sector_0h3p_overlap("overlap_0h3p", CC_PROPERTIES_APPROX_QUADRATIC);
        overlap_diagrams[n_overlap_diagrams++] = cc_strdup("overlap_0h3p");
    }

    /*
     * obtain overlap matrix in the basis of Slater determinants
     */
    double complex *overlap_slater = x_zeros(CC_COMPLEX, model_space_dim, model_space_dim);
    construct_cc_wfn_overlap_matrix(sect_h, sect_p, model_space_dim, det_list, overlap_slater,
                                    scalar_part, n_overlap_diagrams, overlap_diagrams);

    /*
     * construct property matrix in the basis of model vectors
     */
    double complex overlap = 0.0 + 0.0 * I;
    msprop_transform_slater_to_model_vectors(
        1, 1, model_space_dim, model_space_dim,
        coef_right, coef_right, overlap_slater, &overlap
    );

    /*
     * clean up
     */
    cc_free(det_list);
    cc_free(coef_left);
    cc_free(coef_right);
    cc_free(overlap_slater);

    return sqrt(creal(overlap));
}


/**
 * constructs overlap (or metric) matrix in the basis of Slater determinants.
 */
void construct_cc_wfn_overlap_matrix(int sect_h, int sect_p, int dim, slater_det_t *det_list,
                                     double complex *overlap_slater, double scalar_term,
                                     int n_diagrams, char **list_diagrams)
{
    memset(overlap_slater, 0, sizeof(double complex) * dim * dim);

    for (int i = 0; i < dim; i++) {
        overlap_slater[i * dim + i] = scalar_term;
    }

    if (sect_h == 0 && sect_p == 0) {
        return;
    }

    for (int k = 0; k < n_diagrams; k++) {
        char *dg_name = list_diagrams[k];
        diagram_t *dg_norm = diagram_stack_find(dg_name);
        setup_slater(dg_norm, (matrix_getter_fun) diagram_get,
                     sect_h, sect_p, sect_h, sect_p, rank(dg_name) / 2);

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                slater_det_t *bra = det_list + i;
                slater_det_t *ket = det_list + j;
                overlap_slater[i * dim + j] += slater_rule(bra, ket);
            }
        }
    }
}


/*
 * overlap of the vacuum-state wave function < 0h0p | 0h0p >.
 * terms up to quadratic are included, S_ij = < i | (1 + (T^+) * T) | j >.
 */
double sector_0h0p_overlap(int approximation)
{
    double overlap = 1.0 + 0.0 * I;

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return overlap;
    }

    double t1t1 = scalar_product("C", "N", "t1c", "t1c");
    double t2t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
    overlap += t1t1 + t2t2;

    /*
     * CCSDT-n and CCSDT models
     */
    if (triples_enabled()) {
        double t3t3 = (1.0 / 36.0) * scalar_product("C", "N", "t3c", "t3c");
        overlap += t3t3;
    }

    return overlap;
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 0h1p sector.
 */
void sector_0h1p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "pp", "11", "12", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    /*
     * prepare tensors for contraction:
     * lines which will be opened in final diagrams are to be turned into active ones
     * (for better performance)
     */
    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);
    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_v1", "000100", 0);
        restrict_valence("t3c+", "t3c+_v1", "100000", 0);
        restrict_valence("s3c", "s3c_v1", "100100", 0);
        restrict_valence("s3c+", "s3c+_v1", "100100", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();

    // Q1
    double t_start_Q1 = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("t1c+_v", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q1", t_start_Q1);

    // Q2
    double t_start_Q2 = abs_time();
    reorder("t2c_v1", "r1", "3412");
    mult("t2c+_v1", "r1", "r2", 3);
    update(result_name, -0.5, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q2", t_start_Q2);

    // Q3
    double t_start_Q3 = abs_time();
    reorder("s1c+", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    update(result_name, +1.0, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q3", t_start_Q3);

    // Q4
    double t_start_Q4 = abs_time();
    reorder("s2c_v1", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    update(result_name, +1.0, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q4", t_start_Q4);

    // Q5
    double t_start_Q5 = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    update(result_name, +1.0, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q5", t_start_Q5);

    // Q6
    double t_start_Q6 = abs_time();
    reorder("s2c+", "r1", "3412");
    mult("s2c", "r1", "r2", 3);
    update(result_name, +0.5, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q6", t_start_Q6);

    /*
     * for CCSDT-n and CCSDT models
     */
    if (triples_enabled()) {

        // Q7
        reorder("t3c_v1", "r1", "456123");
        mult("t3c+_v1", "r1", "r2", 5);
        update(result_name, -1.0 / 12.0, "r2");
        restore_stack_pos(pos);

        // Q8
        reorder("s3c_v1", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        update(result_name, 0.25, "r2");
        restore_stack_pos(pos);

        // Q9
        reorder("s3c+_v1", "r1", "145623");
        mult("r1", "t2c", "r2", 4);
        update(result_name, 0.25, "r2");
        restore_stack_pos(pos);

        // Q10
        reorder("s3c+", "r1", "456123");
        mult("s3c", "r1", "r2", 5);
        update(result_name, +1.0 / 12.0, "r2");
        restore_stack_pos(pos);
    }
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 1h0p sector.
 *
 * TODO: diagrams for the CCSDT-n and CCSDT models
 */
void sector_1h0p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "hh", "11", "12", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_g1", "100000", 0);
        restrict_valence("t3c+", "t3c+_g1", "000100", 0);
        restrict_valence("h3c", "h3c_g1", "100100", 0);
        restrict_valence("h3c+", "h3c+_g1", "100100", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();

    // Q1
    reorder("t1c+", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q2
    reorder("t2c+", "r1", "3412");
    mult("t2c", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, 0.5, "r3");
    restore_stack_pos(pos);

    // Q3
    reorder("h1c", "r1", "21");
    mult("h1c+", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4
    reorder("h2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q5
    reorder("h2c+", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    closed("r2", "r3");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q6
    reorder("h2c", "r1", "3412");
    mult("h2c+", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    /*
     * for CCSDT-n and CCSDT models
     */
    if (triples_enabled()) {

        // Q7
        reorder("t3c+_g1", "r1", "456123");
        mult("t3c_g1", "r1", "r2", 5);
        update(result_name, +1.0 / 12.0, "r2");
        restore_stack_pos(pos);

        // Q8
        reorder("h3c_g1", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        update(result_name, +1.0 / 4.0, "r2");
        restore_stack_pos(pos);

        // Q9
        reorder("h3c+_g1", "r1", "145623");
        mult("r1", "t2c", "r2", 4);
        update(result_name, +1.0 / 4.0, "r2");
        restore_stack_pos(pos);

        // Q10
        reorder("h3c", "r1", "456123");
        mult("h3c+", "r1", "r2", 5);
        update(result_name, -1.0 / 12.0, "r2");
        restore_stack_pos(pos);
    }
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 1h1p sector.
 *
 * TODO: diagrams for the CCSDT-n and CCSDT models
 */
void sector_1h1p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "phph", "1111", "1234", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    restrict_valence("t1c+", "t1c+_g", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);
    restrict_valence("t1c", "t1c_g", "10", 0);
    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t2c", "t2c_vg", "0110", 0);
    restrict_valence("t2c+", "t2c+_vg", "1001", 0);

    restrict_valence("s2c", "s2c_vg", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_vg", "1011", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);

    restrict_valence("h2c", "h2c_gv", "1011", 0);
    restrict_valence("h2c", "h2c_g1", "1010", 0);
    restrict_valence("h2c+", "h2c+_gv", "1110", 0);
    restrict_valence("h2c+", "h2c+_g1", "1010", 0);

    restrict_valence("e2c", "e2c_v", "1011", 0);
    restrict_valence("e2c", "e2c_g", "1101", 0);
    restrict_valence("e2c+", "e2c+_v", "1110", 0);
    restrict_valence("e2c+", "e2c+_g", "0111", 0);

    dg_stack_pos_t pos = get_stack_pos();

    // Q1a
    restrict_valence("t1c+", "t1c+_vg", "11", 0);
    restrict_valence("t1c", "t1c_gv", "11", 0);
    construct_disconnected_rank2_rank2("t1c+_vg", "t1c_gv", "r1");
    reorder("r1", "r2", "1243"); // the sign is changed here
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // Q1b
    construct_disconnected_rank2_rank2("e1c", "e1c+", "r1");
    reorder("r1", "r2", "1243"); // the sign is changed here
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // Q2
    reorder("t2c+_vg", "r1", "1432");
    reorder("t2c_vg", "r2", "2314");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q3a
    reorder("t1c+_g", "r1", "21");
    mult("s2c_vg", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // Q3b
    reorder("h2c_gv", "r1", "3412");
    mult("r1", "t1c+_v", "r2", 1);
    reorder("r2", "r3", "4321"); // + interchange 1 <-> 2
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4a
    reorder("s2c+_vg", "r1", "3412");
    mult("r1", "t1c_g", "r2", 1);
    reorder("r2", "r3", "3412");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4b
    reorder("t1c_v", "r1", "21");
    mult("h2c+_gv", "r1", "r2", 1);
    reorder("r2", "r3", "2143"); // + interchange 1 <-> 2
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5a
    reorder("s2c_v1", "r1", "1324");
    reorder("h2c+_g1", "r2", "1342");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q5b
    reorder("s2c+_v1", "r1", "1342");
    reorder("h2c_g1", "r2", "1324");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6a
    reorder("e2c_v", "r1", "3412");
    mult("r1", "h1c+", "r2", 1);
    reorder("r2", "r3", "3412");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q6b
    reorder("s1c+", "r1", "21");
    reorder("e2c_g", "r2", "1243");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q7a
    reorder("h1c", "r1", "21");
    mult("e2c+_v", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // Q7b
    reorder("e2c+_g", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // Q8
    reorder("e2c+", "r1", "2341");
    reorder("e2c", "r2", "1423");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1342");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    if (triples_enabled()) {

    }
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 2h0p sector.
 */
void sector_2h0p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "hhhh", "1111", "1234", NOT_PERM_UNIQUE);
    return;
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 0h2p sector.
 */
void sector_0h2p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);
    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    restrict_valence("x2c", "x2c_v1", "1110", 0);
    restrict_valence("x2c+", "x2c+_v1", "1011", 0);

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_v12", "000110", 0);
        restrict_valence("t3c+", "t3c+_v12", "110000", 0);

        restrict_valence("s3c", "s3c_v12", "100110", 0);
        restrict_valence("s3c", "s3c_v2", "100010", 0);
        restrict_valence("s3c+", "s3c+_v12", "110100", 0);
        restrict_valence("s3c+", "s3c+_v2", "010100", 0);

        restrict_valence("x3c", "x3c_v1", "110100", 0);
        restrict_valence("x3c", "x3c_v12", "110110", 0);
        restrict_valence("x3c+", "x3c+_v1", "100110", 0);
        restrict_valence("x3c+", "x3c+_v12", "110110", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();

    // Q1
    double t_start_Q1 = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("t2c+_v12", "r1", "r2", 2);
    update(result_name, +0.5, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q1", t_start_Q1);

    // Q2
    double t_start_Q2 = abs_time();
    reorder("s2c_v12", "r1", "1342");
    mult("r1", "t1c+_v", "r2", 1);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q2", t_start_Q2);

    // Q3
    double t_start_Q3 = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("s2c+_v12", "r1", "r2", 1);
    perm("r2", "(34)");
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q3", t_start_Q3);

    // Q4
    double t_start_Q4 = abs_time();
    reorder("s2c+_v1", "s2c+_21", "2143"); // interchange electrons 1 <-> 2
    reorder("s2c+_21", "r1", "2431");
    reorder("s2c_v1", "r2", "1324");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12|34)");
    update(result_name, +1.0, "r4");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q4", t_start_Q4);

    // Q5
    double t_start_Q5 = abs_time();
    reorder("s1c+", "r1", "21");
    mult("x2c_v1", "r1", "r2", 1);
    perm("r2", "(34)");
    update(result_name, +1.0, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q5", t_start_Q5);

    // Q6
    double t_start_Q6 = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    mult("r1", "s1c", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q6", t_start_Q6);

    // Q7
    double t_start_Q7 = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("x2c", "r1", "r2", 2);
    update(result_name, +0.5, "r2");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q7", t_start_Q7);

    /*
     * for CCSDT-n and CCSDT models
     */
    if (triples_enabled()) {

        // Q8
        double t_start_Q8 = abs_time();
        reorder("x3c_v12", "r1", "124563");
        mult("r1", "t1c+", "r2", 2);
        update(result_name, +1.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q8", t_start_Q8, 0);

        // Q9
        double t_start_Q9 = abs_time();
        reorder("x3c+_v12", "r1", "124563");
        mult("r1", "t1c", "r2", 2);
        update(result_name, +1.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9", t_start_Q9, 0);

        // Q10
        double t_start_Q10 = abs_time();
        reorder("s2c+", "r1", "3412");
        reorder("x3c_v1", "r2", "124356");
        mult("r2", "r1", "r3", 3);
        perm("r3", "(34)");
        update(result_name, +0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10", t_start_Q10, 0);

        // Q11
        double t_start_Q11 = abs_time();
        reorder("x3c+_v1", "r1", "145623");
        mult("r1", "s2c", "r2", 3);
        reorder("r2", "r3", "1423");
        perm("r3", "(12)");
        update(result_name, +0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11", t_start_Q11, 0);

        // Q12
        double t_start_Q12 = abs_time();
        reorder("s3c_v12", "r1", "145623");
        mult("r1", "t2c+_v1", "r2", 3);
        reorder("r2", "r3", "1423");
        perm("r3", "(12)");
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12", t_start_Q12, 0);

        // Q13
        double t_start_Q13 = abs_time();
        reorder("s3c+_v12", "r1", "124356");
        reorder("t2c_v1", "r2", "3412");
        mult("r1", "r2", "r3", 3);
        perm("r3", "(34)");
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13", t_start_Q13, 0);

        // Q14
        double t_start_Q14 = abs_time();
        reorder("t3c_v12", "r1", "456123");
        mult("t3c+_v12", "r1", "r2", 4);
        update(result_name, 1.0 / 6.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14", t_start_Q14, 0);

        // Q15
        double t_start_Q15 = abs_time();
        reorder("s3c+_v2", "r1", "245613");
        reorder("s3c_v2", "r2", "152346");
        mult("r2", "r1", "r3", 4);
        reorder("r3", "r4", "1342");
        perm("r4", "(12|34)");
        update(result_name, -0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15", t_start_Q15, 0);

        // Q16
        double t_start_Q16 = abs_time();
        reorder("x3c+", "r1", "456123");
        mult("x3c", "r1", "r2", 4);
        update(result_name, 1.0 / 6.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16", t_start_Q16, 0);
    }
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 3h0p sector.
 */
void sector_3h0p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "hhhhhh", "111111", "123456", NOT_PERM_UNIQUE);
    return;
}


/*
 * calculates diagrams contributing to the overlap (metric) matrix in the 0h3p sector.
 */
void sector_0h3p_overlap(char *result_name, int approximation)
{
    tmplt(result_name, "pppppp", "111111", "123456", NOT_PERM_UNIQUE);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);
    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    restrict_valence("x2c", "x2c_v1", "1110", 0);
    restrict_valence("x2c+", "x2c+_v1", "1011", 0);

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_v12", "000110", 0);
        restrict_valence("t3c+", "t3c+_v12", "110000", 0);
        restrict_valence("t3c", "t3c_v123", "000111", 0);
        restrict_valence("t3c+", "t3c+_v123", "111000", 0);

        restrict_valence("s3c", "s3c_v2", "100010", 0);
        restrict_valence("s3c+", "s3c+_v2", "010100", 0);
        restrict_valence("s3c", "s3c_v12", "100110", 0);
        restrict_valence("s3c+", "s3c+_v12", "110100", 0);
        restrict_valence("s3c", "s3c_v23", "100011", 0);
        restrict_valence("s3c+", "s3c+_v23", "011100", 0);
        restrict_valence("s3c", "s3c_v123", "100111", 0);
        restrict_valence("s3c+", "s3c+_v123", "111100", 0);

        restrict_valence("x3c", "x3c_v1", "110100", 0);
        restrict_valence("x3c+", "x3c+_v1", "100110", 0);
        restrict_valence("x3c", "x3c_v3", "110001", 0);
        restrict_valence("x3c+", "x3c+_v3", "001110", 0);
        restrict_valence("x3c", "x3c_v12", "110110", 0);
        restrict_valence("x3c+", "x3c+_v12", "110110", 0);
        restrict_valence("x3c", "x3c_v123", "110111", 0);
        restrict_valence("x3c+", "x3c+_v123", "111110", 0);

        restrict_valence("z3c", "z3c_v1", "111100", 0);
        restrict_valence("z3c+", "z3c+_v1", "100111", 0);
        restrict_valence("z3c", "z3c_v12", "111110", 0);
        restrict_valence("z3c+", "z3c+_v12", "110111", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();

    /*
     * CCSD terms
     */

    // Q1
    double t_start_Q1 = abs_time();
    reorder("s2c+_v12", "r1", "2134");
    reorder("s2c_v12", "r2", "1342");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "145236");
    perm("r4", "(1/23|6/45)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q1", t_start_Q1);

    // Q2
    double t_start_Q2 = abs_time();
    reorder("x2c+_v1", "x2c+_v2", "2143");
    reorder("x2c+_v2", "r1", "2341");
    mult("x2c_v1", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|4/56)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    //print_time_mem_usage("Q2", t_start_Q2);

    /*
     * for CCSDT-n and CCSDT models
     */
    if (triples_enabled()) {

        // Q3a
        double t_start_Q3a = abs_time();
        reorder("x3c_v123", "r1", "456123");
        mult("r1", "t1c+_v", "r2", 1);
        reorder("r2", "r3", "456123");
        perm("r3", "(3/12)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3a", t_start_Q3a, 0);

        // Q3b
        double t_start_Q3b = abs_time();
        reorder("t1c_v", "r1", "21");
        mult("x3c+_v123", "r1", "r2", 1);
        perm("r2", "(6/45)");
        update(result_name, -1.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3b", t_start_Q3b, 0);

        // Q4a
        double t_start_Q3c = abs_time();
        reorder("s1c+", "r1", "21");
        mult("z3c_v12", "r1", "r2", 1);
        perm("r2", "(6/45)");
        update(result_name, 1.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3c", t_start_Q3c, 0);

        // Q4b
        double t_start_Q3d = abs_time();
        reorder("z3c+_v12", "r1", "456123");
        mult("r1", "s1c", "r2", 1);
        reorder("r2", "r3", "456123");
        perm("r3", "(3/12)");
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3d", t_start_Q3d, 0);

        // Q5a
        double t_start_Q4a = abs_time();
        reorder("s3c_v123", "r1", "456123");
        mult("r1", "t2c+_v12", "r2", 2);
        reorder("r2", "r3", "456123");
        perm("r3", "(1/23)");
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q4a", t_start_Q4a, 0);

        // Q5b
        double t_start_Q4b = abs_time();
        reorder("t2c_v12", "r1", "3412");
        mult("s3c+_v123", "r1", "r2", 2);
        perm("r2", "(4/56)");
        update(result_name, 0.5, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q4b", t_start_Q4b, 0);

        // Q6a
        double t_start_Q5a = abs_time();
        reorder("x3c_v12", "r1", "124536");
        reorder("s2c+_v1", "s2c+_v1_r", "2143");
        reorder("s2c+_v1_r", "r2", "2431");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "125346");
        perm("r4", "(3/12|6/45)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q5a", t_start_Q5a, 0);

        // Q6b
        double t_start_Q5b = abs_time();
        reorder("x3c+_v12", "r1", "124563");
        reorder("s2c_v1", "s2c_v1_r", "2143");
        reorder("s2c_v1_r", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "125346");
        perm("r4", "(3/12|6/45)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q5b", t_start_Q5b, 0);

        // Q7a
        double t_start_Q6a = abs_time();
        reorder("x2c+", "r1", "3412");
        mult("z3c_v1", "r1", "r2", 2);
        perm("r2", "(4/56)");
        update(result_name, 0.5, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q6a", t_start_Q6a, 0);

        // Q7b
        double t_start_Q6b = abs_time();
        reorder("z3c+_v1", "r1", "456123");
        mult("r1", "x2c", "r2", 2);
        reorder("r2", "r3", "456123");
        perm("r3", "(1/23)");
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q6b", t_start_Q6b, 0);

        // Q8
        double t_start_Q7 = abs_time();
        reorder("t3c_v123", "r1", "456123");
        mult("t3c+_v123", "r1", "r2", 3);
        update(result_name, -1.0/6.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7", t_start_Q7, 0);

        // Q9
        double t_start_Q8 = abs_time();
        reorder("s3c+_v23", "r1", "234561");
        reorder("s3c_v23", "r2", "156234");
        mult("r2", "r1", "r3", 3);
        reorder("r3", "r4", "145623");
        perm("r4", "(1/23|4/56)");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q8", t_start_Q8, 0);

        // Q10
        double t_start_Q9 = abs_time();
        reorder("x3c+_v3", "r1", "345612");
        reorder("x3c_v3", "r2", "126345");
        mult("r2", "r1", "r3", 3);
        reorder("r3", "r4", "124563");
        perm("r4", "(3/12|6/45)");
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9", t_start_Q9, 0);

        // Q11
        double t_start_Q10 = abs_time();
        reorder("z3c+", "r1", "456123");
        mult("z3c", "r1", "r2", 3);
        update(result_name, 1.0/6.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10", t_start_Q10, 0);
    }
}









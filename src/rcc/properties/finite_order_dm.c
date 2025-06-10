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

#include "finite_order_dm.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "engine.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/mvcoef.h"
#include "finite_order_overlap.h"
#include "natural_spinors.h"

#include "finite_order_dm_0h0p.h"
#include "finite_order_dm_0h1p.h"
#include "finite_order_dm_1h0p.h"
#include "finite_order_dm_0h2p.h"
#include "finite_order_dm_0h3p.h"


void save_density_matrix(char *path, int nspinors, double complex *dm);

void finite_order_dm_prepare_cluster_operators(int sect_h, int sect_p);

void conjugate_cluster_amplitudes_0h0p();

void conjugate_cluster_amplitudes_0h1p();

void conjugate_cluster_amplitudes_1h0p();

void conjugate_cluster_amplitudes_0h2p();

int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

double complex calculate_vacuum_expectation_value(int nspinors, double complex *prop_matrix);

void finite_order_density_matrix_block_hp(int sect_h, int sect_p, int dm_sym, int transition_dm);

void finite_order_density_matrix_block_ph(int sect_h, int sect_p, int dm_sym, int transition_dm);

void finite_order_density_matrix_block_hh(int sect_h, int sect_p, int dm_sym, int transition_dm);

void finite_order_density_matrix_block_pp(int sect_h, int sect_p, int dm_sym, int transition_dm);

void construct_two_part_active_density_matrix_diagram_0h2p_fast(
    char *diagram_name,
    char *bra_irrep_name, int bra_state_number,
    char *ket_irrep_name, int ket_state_number
);


/**
 * calculates the series of pure density matrices for ALL electronic
 * states defined by the 'roots_cutoff' or 'nroots' directives
 * (these model vectors are stored in the MVCOEF* files).
 */
int finite_order_density_matrix_all_states(int sect_h, int sect_p)
{
    int nrep = 0;
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(sect_h, sect_p, NULL, &nrep, mv_blocks);

    for (int irep = 0; irep < nrep; irep++) {
        for (int iroot = 0; iroot < mv_blocks[irep].nroots; iroot++) {
            finite_order_density_matrix(sect_h, sect_p, mv_blocks[irep].rep_name, iroot,
                                        mv_blocks[irep].rep_name, iroot);
        }
    }

    // deallocate memory
    for (size_t irep = 0; irep < nrep; irep++) {
        mvblock_free(mv_blocks + irep);
    }
}


/**
 * this subroutine constructs approximate density matrix or transition density matrix
 * using the finite-order method.
 * currently this method is implemented only for the CCSD model and the
 * 0h1p, 1h0p, 0h2p, 0h3p sectors.
 *
 * to calculate density matrix for the given pure (ground or excited) state, use:
 * irrep_name_1 == irrep_name_2
 * state_number_1 == state_number_2
 * otherwise transition density matrix will be calculated.
 *
 * density matrices are stored in txt formatted files:
 * density_*.txt and transition_density_*.txt
 *
 * this subroutine also calculates natural spinors or natural transition spinors
 */
int finite_order_density_matrix(
    int sect_h, int sect_p,
    char *irrep_name_1, int state_number_1,
    char *irrep_name_2, int state_number_2
)
{
    if (!((sect_h == 0 && sect_p == 1) ||
          (sect_h == 1 && sect_p == 0) ||
          (sect_h == 0 && sect_p == 2) ||
          (sect_h == 0 && sect_p == 3))) {
        errquit("finite-order for density matrix calculations is not yet implemented for the %dh%dp sector\n",
                sect_h, sect_p);
    }

    /*
     * density matrix for a pure state or a transition density matrix?
     */
    int tran_dm = 0;
    if (strcmp(irrep_name_1, irrep_name_2) != 0 ||
        state_number_1 != state_number_2) {
        tran_dm = 1;
    }

    /*
     * banner for this module
     */
    printf("\n");
    printf(" **\n");
    if (tran_dm == 0) {
        printf(" ** construction of the density matrix for the state %d (irrep %s)\n", state_number_1 + 1,
               irrep_name_1);
    }
    else {
        printf(
            " ** construction of the transition density matrix for the pair of states %d (irrep %s) - %d (irrep %s)\n",
            state_number_1 + 1, irrep_name_1, state_number_2 + 1, irrep_name_2);
    }
    printf(" ** (finite-order method quadratic in cluster amplitudes)\n");
    printf(" **\n");
    printf("\n");

    /*
     * possibly the 1-DM for the 0h0p sector was pre-calculated previously
     * by the exact method (lambda equations). in this case we read it from disk
     */
    double complex *analytic_dm_0h0p = NULL;
    if (cc_opts->calc_density[0][0] == CC_DENSITY_MATRIX_LAMBDA) {
        int nspinors = get_num_spinors();
        analytic_dm_0h0p = z_zeros(nspinors, nspinors);
        printf(" %50s", "reading analytic density matrix from disk...");
        read_prop_single_file(nspinors, "density_0h0p.txt", analytic_dm_0h0p);
        printf("done\n");
    }

    /*
     * only for transition moments: norms of wavefunctions
     */
    double norm_bra = 1.0;
    double norm_ket = 1.0;

    if (tran_dm) {
        printf(" %-50s", "calculation of wavefunction norms...");

        double t_start_norms = abs_time();
        norm_bra = calculate_wavefunction_norm(sect_h, sect_p, irrep_name_1, state_number_1);
        norm_ket = calculate_wavefunction_norm(sect_h, sect_p, irrep_name_2, state_number_2);
        printf("done in %.3f sec\n", abs_time() - t_start_norms);
    }

    /*
     * construct "model space" density matrices:
     * one-particle
     */
    printf(" %-50s", "construction of active-space 1-DM...");
    double t_start_1dm = abs_time();
    if (sect_h == 0 && sect_p == 1) {
        construct_one_part_active_density_matrix_diagram_0h1p("model_01_1dm", irrep_name_1, state_number_1,
                                                              irrep_name_2, state_number_2);
    }
    if (sect_h == 1 && sect_p == 0) {
        construct_one_part_active_density_matrix_diagram_1h0p("model_10_1dm", irrep_name_1, state_number_1,
                                                              irrep_name_2, state_number_2);
    }
    else if (sect_h == 0 && sect_p == 2) {
        construct_n_part_active_density_matrix_diagram_0h2p("model_01_1dm", irrep_name_1, state_number_1, irrep_name_2,
                                                            state_number_2, 1);
    }
    else if (sect_h == 0 && sect_p == 3) {
        construct_n_part_active_density_matrix_diagram_0h3p("model_01_1dm", irrep_name_1, state_number_1, irrep_name_2,
                                                            state_number_2, 1);
    }
    printf("done in %.3f sec\n", abs_time() - t_start_1dm);

    /*
     * construct "model space" density matrices:
     * two-particle
     */
    if ((sect_h == 0 && sect_p == 2) ||
        (sect_h == 0 && sect_p == 3)) {
        printf(" %-50s", "construction of active-space 2-DM...");
        double t_start_2dm = abs_time();
        if (sect_h == 0 && sect_p == 2) {
            construct_two_part_active_density_matrix_diagram_0h2p_fast("model_02_2dm", irrep_name_1, state_number_1, irrep_name_2, state_number_2);
        }
        else if (sect_h == 0 && sect_p == 3) {
            construct_n_part_active_density_matrix_diagram_0h3p("model_02_2dm", irrep_name_1, state_number_1,
                                                                irrep_name_2, state_number_2, 2);
        }
        printf("done in %.3f sec\n", abs_time() - t_start_2dm);
    }

    /*
     * construct "model space" density matrices:
     * three-particle
     */
    if (sect_h == 0 && sect_p == 3) {
        printf(" %-50s", "construction of active-space 3-DM...");
        double t_start_3dm = abs_time();
        construct_n_part_active_density_matrix_diagram_0h3p("model_03_3dm", irrep_name_1, state_number_1, irrep_name_2,
                                                            state_number_2, 3);
        printf("done in %.3f sec\n", abs_time() - t_start_3dm);
    }

    /*
     * prepare cluster operators, T => T^+
     * (complex conjugation)
     */
    printf(" %-50s", "complex conjugation of cluster operators...");
    double t_start_tconj = abs_time();
    conjugate_cluster_amplitudes_0h0p();
    if (sect_h >= 1) {
        conjugate_cluster_amplitudes_1h0p();
    }
    if (sect_p >= 1) {
        conjugate_cluster_amplitudes_0h1p();
    }
    if (sect_p >= 2) {
        conjugate_cluster_amplitudes_0h2p();
    }
    double t_stop_tconj = abs_time();
    printf("done in %.3f sec\n", t_stop_tconj - t_start_tconj);

    /*
     * restrict cluster operators to active lines
     * (in different ways)
     */
    printf(" %-50s", "restrict cluster operators to active lines...");
    double t_start_restrict_active = abs_time();
    finite_order_dm_prepare_cluster_operators(sect_h, sect_p);
    printf("done in %.3f sec\n", abs_time() - t_start_restrict_active);

    /*
     * determine the "symmetry" of the operator representing DM
     */
    int op_sym = detect_operator_symmetry(irrep_name_1, irrep_name_2);

    /*
     * construct diagrams corresponding to different blocks of the (full) one-particle DM
     */
    printf(" %-50s", "construction of the hole-particle block...");
    double t_start_hp = abs_time();
    finite_order_density_matrix_block_hp(sect_h, sect_p, op_sym, tran_dm);   // => "dm_hp"
    printf("done in %.3f sec\n", abs_time() - t_start_hp);

    printf(" %-50s", "construction of the particle-hole block...");
    double t_start_ph = abs_time();
    finite_order_density_matrix_block_ph(sect_h, sect_p, op_sym, tran_dm);   // => "dm_ph"
    printf("done in %.3f sec\n", abs_time() - t_start_ph);

    printf(" %-50s", "construction of the hole-hole block...");
    double t_start_hh = abs_time();
    finite_order_density_matrix_block_hh(sect_h, sect_p, op_sym, tran_dm);   // => "dm_hh"
    printf("done in %.3f sec\n", abs_time() - t_start_hh);

    printf(" %-50s", "construction of the particle-particle block...");
    double t_start_pp = abs_time();
    finite_order_density_matrix_block_pp(sect_h, sect_p, op_sym, tran_dm);   // => "dm_pp"
    printf("done in %.3f sec\n", abs_time() - t_start_pp);

    /*
     * extract matrix elements from diagrams representing 1-DM
     */
    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);

    diagram_t *dg_dm_hp = diagram_stack_find("dm_hp");
    diagram_t *dg_dm_ph = diagram_stack_find("dm_ph");
    diagram_t *dg_dm_hh = diagram_stack_find("dm_hh");
    diagram_t *dg_dm_pp = diagram_stack_find("dm_pp");

    for (int p = 0; p < nspinors; p++) {
        for (int q = 0; q < nspinors; q++) {
            double complex gamma_pq = 0.0 + 0.0 * I;

            /*
             * contribution from the 0h0p sector
             */
            if (op_sym == get_totally_symmetric_irrep() && !tran_dm) {
                if (analytic_dm_0h0p) {
                    gamma_pq = analytic_dm_0h0p[p * nspinors + q];
                }
                else {
                    if (is_hole(p) && is_hole(q)) {
                        gamma_pq = (p == q);
                    }
                }
            }

            if (is_hole(p) && is_hole(q)) { // h-h
                gamma_pq += /*(p == q) +*/ diagram_get_2(dg_dm_hh, p, q);
            }
            else if (is_hole(p) && is_particle(q)) { // h-p
                gamma_pq += diagram_get_2(dg_dm_hp, p, q);
            }
            else if (is_particle(p) && is_hole(q)) { // p-h
                gamma_pq += diagram_get_2(dg_dm_ph, p, q);
            }
            else { // p-p
                gamma_pq += diagram_get_2(dg_dm_pp, p, q);
            }

            dm[p * nspinors + q] = gamma_pq * norm_bra / norm_ket;
        }
    }

    /*
     * save density matrix to disk
     */
    char dm_file_name[256];
    if (tran_dm == 0) {
        sprintf(dm_file_name, "density_%dh%dp_%s_%d.txt", sect_h, sect_p, irrep_name_1, state_number_1 + 1);
    }
    else {
        sprintf(dm_file_name, "transition_density_%dh%dp_%s_%d_%s_%d.txt", sect_h, sect_p,
                irrep_name_1, state_number_1 + 1, irrep_name_2, state_number_2 + 1);
    }
    str_replace(dm_file_name, '/', '.');
    save_density_matrix(dm_file_name, nspinors, dm);
    printf(" %-50s%s\n", "density matrix has written to file...", dm_file_name);

    /*
     * natural spinors and their occupation numbers
     * OR
     * natural transition spinors and their singular values
     */
    if (tran_dm == 0) {
        char nat_spinors_file_name[256];
        sprintf(nat_spinors_file_name, "natural_spinors_%dh%dp_%s_%d.txt", sect_h, sect_p, irrep_name_1,
                state_number_1 + 1);
        str_replace(nat_spinors_file_name, '/', '.');

        calculate_natural_spinors_and_occ_numbers(dm_file_name, nat_spinors_file_name);
    }
    else {
        char nat_spinors_from_file_name[256];
        sprintf(nat_spinors_from_file_name, "natural_spinors_%dh%dp_%s_%d_%s_%d_from.txt",
                sect_h, sect_p, irrep_name_1, state_number_1 + 1, irrep_name_2, state_number_2 + 1);
        str_replace(nat_spinors_from_file_name, '/', '.');

        char nat_spinors_to_file_name[256];
        sprintf(nat_spinors_to_file_name, "natural_spinors_%dh%dp_%s_%d_%s_%d_to.txt",
                sect_h, sect_p, irrep_name_1, state_number_1 + 1, irrep_name_2, state_number_2 + 1);
        str_replace(nat_spinors_to_file_name, '/', '.');

        calculate_natural_transition_spinors_and_singular_values(
            dm_file_name,
            nat_spinors_from_file_name,
            nat_spinors_to_file_name
        );
    }


    // contract DM with property
    /*double complex *prop_matrix = z_zeros(nspinors, nspinors);
    int status = read_prop_single_file(nspinors, "ZDIPLEN", prop_matrix);

    // calculate correlation contribution to the expectation value
    double complex expect_value = 0.0;
    for (int p = 0; p < nspinors; p++) {
        for (int q = 0; q < nspinors; q++) {
            expect_value += dm[p * nspinors + q] * prop_matrix[p * nspinors + q];
        }
    }

    // print expectation value
    double complex scf_contrib = calculate_vacuum_expectation_value(nspinors, prop_matrix);
    double complex ccsd_contrib = expect_value - scf_contrib;

    printf(" SCF contribution          %20.12f%20.12f\n", creal(scf_contrib), cimag(scf_contrib));
    printf(" Correlation contribution  %20.12f%20.12f\n", creal(ccsd_contrib), cimag(ccsd_contrib));
    printf(" Final expectation value   %20.12f%20.12f\n", creal(expect_value), cimag(expect_value));

    cc_free(prop_matrix);*/


    // clean up
    cc_free(dm);
    if (analytic_dm_0h0p) {
        cc_free(analytic_dm_0h0p);
    }

    return EXIT_SUCCESS;
}


/**
 * constructs the 'h-p' block of the one-particle [transition] density matrix.
 *
 * @param dm_sym: in the case of transition density matrix the symmetry of the corresponding
 * diagram (operator) can be non-totally symmetric. 'dm_sym' stands for the number of the
 * corresponding irreducible representation.
 * @param transition_dm is equal to 1 if the transition density matrix is to be calculated;
 * equals to 0 otherwise.
 */
void finite_order_density_matrix_block_hp(int sect_h, int sect_p, int dm_sym, int transition_dm)
{
    tmplt_sym("dm_hp", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);

    if (cc_opts->calc_density[0][0] != CC_DENSITY_MATRIX_LAMBDA && (transition_dm == 0)) {
        finite_order_density_matrix_block_hp_0h0p(dm_sym);
    }

    if (sect_h >= 1) {
        finite_order_density_matrix_block_hp_1h0p(dm_sym);
    }
    if (sect_p >= 1) {
        finite_order_density_matrix_block_hp_0h1p(dm_sym);
    }
    if (sect_p >= 2) {
        finite_order_density_matrix_block_hp_0h2p(dm_sym);
    }
    if (sect_p >= 3) {
        finite_order_density_matrix_block_hp_0h3p(dm_sym);
    }
}


/**
 * constructs the 'p-h' block of the one-particle [transition] density matrix.
 *
 * @param dm_sym: in the case of transition density matrix the symmetry of the corresponding
 * diagram (operator) can be non-totally symmetric. 'dm_sym' stands for the number of the
 * corresponding irreducible representation.
 * @param transition_dm is equal to 1 if the transition density matrix is to be calculated;
 * equals to 0 otherwise.
 */
void finite_order_density_matrix_block_ph(int sect_h, int sect_p, int dm_sym, int transition_dm)
{
    tmplt_sym("dm_ph", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);

    if (cc_opts->calc_density[0][0] != CC_DENSITY_MATRIX_LAMBDA && (transition_dm == 0)) {
        finite_order_density_matrix_block_ph_0h0p(dm_sym);
    }

    if (sect_h >= 1) {
        finite_order_density_matrix_block_ph_1h0p(dm_sym);
    }
    if (sect_p >= 1) {
        finite_order_density_matrix_block_ph_0h1p(dm_sym);
    }
    if (sect_p >= 2) {
        finite_order_density_matrix_block_ph_0h2p(dm_sym);
    }
    if (sect_p >= 3) {
        finite_order_density_matrix_block_ph_0h3p(dm_sym);
    }
}


/**
 * constructs the 'h-h' block of the one-particle [transition] density matrix.
 *
 * @param dm_sym: in the case of transition density matrix the symmetry of the corresponding
 * diagram (operator) can be non-totally symmetric. 'dm_sym' stands for the number of the
 * corresponding irreducible representation.
 * @param transition_dm is equal to 1 if the transition density matrix is to be calculated;
 * equals to 0 otherwise.
 */
void finite_order_density_matrix_block_hh(int sect_h, int sect_p, int dm_sym, int transition_dm)
{
    tmplt_sym("dm_hh", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);

    if (cc_opts->calc_density[0][0] != CC_DENSITY_MATRIX_LAMBDA && (transition_dm == 0)) {
        finite_order_density_matrix_block_hh_0h0p(dm_sym);
    }

    if (sect_h >= 1) {
        finite_order_density_matrix_block_hh_1h0p(dm_sym);
    }
    if (sect_p >= 1) {
        finite_order_density_matrix_block_hh_0h1p(dm_sym);
    }
    if (sect_p >= 2) {
        finite_order_density_matrix_block_hh_0h2p(dm_sym);
    }
    if (sect_p >= 3) {
        finite_order_density_matrix_block_hh_0h3p(dm_sym);
    }
}


/**
 * constructs the 'p-p' block of the one-particle [transition] density matrix.
 *
 * @param dm_sym: in the case of transition density matrix the symmetry of the corresponding
 * diagram (operator) can be non-totally symmetric. 'dm_sym' stands for the number of the
 * corresponding irreducible representation.
 * @param transition_dm is equal to 1 if the transition density matrix is to be calculated;
 * equals to 0 otherwise.
 */
void finite_order_density_matrix_block_pp(int sect_h, int sect_p, int dm_sym, int transition_dm)
{
    tmplt_sym("dm_pp", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);

    if (cc_opts->calc_density[0][0] != CC_DENSITY_MATRIX_LAMBDA && (transition_dm == 0)) {
        finite_order_density_matrix_block_pp_0h0p(dm_sym);
    }

    if (sect_h >= 1) {
        finite_order_density_matrix_block_pp_1h0p(dm_sym);
    }
    if (sect_p >= 1) {
        finite_order_density_matrix_block_pp_0h1p(dm_sym);
    }
    if (sect_p >= 2) {
        finite_order_density_matrix_block_pp_0h2p(dm_sym);
    }
    if (sect_p >= 3) {
        finite_order_density_matrix_block_pp_0h3p(dm_sym);
    }
}


/**
 * prepares operators representing cluster amplitudes for construction of density matrices
 * using the finite-order method. the idea is to turn some 'general' lines of diagrams to
 * 'valence' ones in order to reduce the computational cost.
 */
void finite_order_dm_prepare_cluster_operators(int sect_h, int sect_p)
{
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

    if (sect_h >= 1) {
        restrict_valence("t1c", "t1c_g", "10", 0);
        restrict_valence("t1c+", "t1c+_g", "01", 0);
        restrict_valence("t2c", "t2c_g1", "1000", 0);
        restrict_valence("t2c+", "t2c+_g1", "0010", 0);

        restrict_valence("h2c", "h2c_g1", "1010", 0);
        restrict_valence("h2c+", "h2c+_g1", "1010", 0);
    }

    if (sect_p >= 1) {
        restrict_valence("s2c", "s2c_v12", "1011", 0);
        restrict_valence("s2c+", "s2c+_v12", "1110", 0);
        restrict_valence("s2c", "s2c_v1", "1010", 0);
        restrict_valence("s2c+", "s2c+_v1", "1010", 0);
        restrict_valence("s2c", "s2c_v2", "1001", 0);
        restrict_valence("s2c+", "s2c+_v2", "0110", 0);
    }

    if (sect_p >= 2) {
        restrict_valence("x2c", "x2c_v1", "1110", 0);
        restrict_valence("x2c+", "x2c+_v1", "1011", 0);
    }
}


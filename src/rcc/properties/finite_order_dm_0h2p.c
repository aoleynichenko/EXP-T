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

#include "finite_order_dm_0h2p.h"

#include <assert.h>
#include <string.h>

#include "engine.h"
#include "options.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/mvcoef.h"
#include "finite_order_overlap.h"
#include "cc_properies.h"


complex double one_part_density_matrix_element_0h2p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pq
);

complex double two_part_density_matrix_element_0h2p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pqrs
);


void construct_n_part_active_density_matrix_diagram_0h2p(
    char *diagram_name,
    char *bra_irrep_name, int bra_state_number,
    char *ket_irrep_name, int ket_state_number,
    int n_particles
)
{
    assert(n_particles == 1 || n_particles == 2);

    int tran_dm = 0;
    if (strcmp(bra_irrep_name, ket_irrep_name) != 0 ||
        bra_state_number != ket_state_number) {
        tran_dm = 1;
    }

    /*
     * bra vector
     */
    int bra_dim = 0;
    double bra_eigenvalue = 0.0;
    double bra_exc_energy_cm = 0.0;
    slater_det_t *bra_det_list = NULL;
    double complex *bra_coef_left = NULL;
    double complex *bra_coef_right = NULL;

    mvcoef_read_vectors_unformatted_for_state(
        0, 2, "MVCOEF02", bra_irrep_name, bra_state_number,
        &bra_dim, &bra_det_list, &bra_eigenvalue, &bra_exc_energy_cm, &bra_coef_left, &bra_coef_right
    );

    /*
     * ket vector - if needed
     */
    int ket_dim = 0;
    double ket_eigenvalue = 0.0;
    double ket_exc_energy_cm = 0.0;
    slater_det_t *ket_det_list = NULL;
    double complex *ket_coef_left = NULL;
    double complex *ket_coef_right = NULL;

    if (tran_dm) {
        mvcoef_read_vectors_unformatted_for_state(
            0, 2, "MVCOEF02", ket_irrep_name, ket_state_number,
            &ket_dim, &ket_det_list, &ket_eigenvalue, &ket_exc_energy_cm, &ket_coef_left, &ket_coef_right
        );
    }
    else {
        ket_dim = bra_dim;
        ket_eigenvalue = bra_eigenvalue;
        ket_exc_energy_cm = bra_exc_energy_cm;
        ket_det_list = bra_det_list;
        ket_coef_left = bra_coef_left;
        ket_coef_right = bra_coef_right;
    }

    /*
     * determine the "symmetry" of the operator representing DM
     */
    int op_sym = detect_operator_symmetry(bra_irrep_name, ket_irrep_name);
    //printf("op_sym = %s\n", get_irrep_name(op_sym));

    /*
     * create empty template for the diagram representing density matrix in the active space
     */
    if (n_particles == 1) {
        tmplt_sym(diagram_name, "pp", "11", "12", NOT_PERM_UNIQUE, op_sym);
    }
    else { // n_particles == 2
        tmplt_sym(diagram_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE, op_sym);
    }
    diagram_t *diag_dm = diagram_stack_find(diagram_name);

    for (int iblock = 0; iblock < diag_dm->n_blocks; iblock++) {
        block_t *block = diag_dm->blocks[iblock];
        int *indices = (int *) cc_malloc(block->size * block->rank * sizeof(int));
        block_gen_indices(block, indices);

        for (size_t i = 0; i < block->size; i++) {
            int *idx = indices + i * block->rank;

            // here we interchange pqr... and stu... blocks of indices since in EXP-T
            // the order of indices in the algebraic representation of diagrams is
            // <in|out> instead of <out|in> (which is the standard one)
            int indices_pqrs[4];
            intcpy(indices_pqrs, idx + n_particles, n_particles);
            intcpy(indices_pqrs + n_particles, idx, n_particles);

            double complex dm_elem = 0.0 + 0.0 * I;
            if (n_particles == 1) {
                dm_elem = one_part_density_matrix_element_0h2p(
                    bra_dim, bra_coef_left, bra_det_list,
                    ket_dim, ket_coef_right, ket_det_list,
                    indices_pqrs
                );
            }
            else { // n_particles == 2
                dm_elem = two_part_density_matrix_element_0h2p(
                    bra_dim, bra_coef_left, bra_det_list,
                    ket_dim, ket_coef_right, ket_det_list,
                    indices_pqrs
                );
            }

            if (arith == CC_ARITH_COMPLEX) {
                block->buf[i] = dm_elem;
            }
            else {
                ((double *) block->buf)[i] = creal(dm_elem);
            }
        }

        cc_free(indices);
    }

    /*
     * clean up
     */
    cc_free(bra_det_list);
    cc_free(bra_coef_left);
    cc_free(bra_coef_right);
    if (tran_dm) {
        cc_free(ket_det_list);
        cc_free(ket_coef_left);
        cc_free(ket_coef_right);
    }
}


complex double one_part_density_matrix_element_0h2p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pq
)
{
    double complex matrix_element = 0.0 + 0.0 * I;

    for (size_t i = 0; i < dim_bra; i++) {
        for (size_t j = 0; j < dim_ket; j++) {
            slater_det_t *bra = &dets_bra[i];
            slater_det_t *ket = &dets_ket[j];

            // bra determinant
            int a = bra->indices[0];
            int b = bra->indices[1];

            // ket determinant
            int c = ket->indices[0];
            int d = ket->indices[1];

            // excitation operator (p+ q)
            int p = indices_pq[0];
            int q = indices_pq[1];

            int bra_pq_ket =
                + (b == p) * (a == c) * (d == q)
                - (b == p) * (a == d) * (c == q)
                - (b == c) * (a == p) * (d == q)
                + (b == d) * (a == p) * (c == q);

            matrix_element += conj(coef_bra[i]) * coef_ket[j] * bra_pq_ket;
        }
    }

    return matrix_element;
}


complex double two_part_density_matrix_element_0h2p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pqrs
)
{
    double complex matrix_element = 0.0 + 0.0 * I;

    for (size_t i = 0; i < dim_bra; i++) {
        for (size_t j = 0; j < dim_ket; j++) {
            slater_det_t *bra = &dets_bra[i];
            slater_det_t *ket = &dets_ket[j];

            // bra determinant
            int a = bra->indices[0];
            int b = bra->indices[1];

            // ket determinant
            int c = ket->indices[0];
            int d = ket->indices[1];

            // excitation operator (p+ q)
            int p = indices_pqrs[0];
            int q = indices_pqrs[1];
            int r = indices_pqrs[2];
            int s = indices_pqrs[3];

            int bra_pqrs_ket =
                + (b == p) * (a == q) * (c == s) * (d == r)
                - (b == p) * (a == q) * (d == s) * (c == r)
                - (b == q) * (a == p) * (c == s) * (d == r)
                + (b == q) * (a == p) * (d == s) * (c == r);

            matrix_element += conj(coef_bra[i]) * coef_ket[j] * bra_pqrs_ket;
        }
    }

    return matrix_element;
}


void construct_two_part_active_density_matrix_diagram_0h2p_fast(
    char *diagram_name,
    char *bra_irrep_name, int bra_state_number,
    char *ket_irrep_name, int ket_state_number
)
{
    int tran_dm = 0;
    if (strcmp(bra_irrep_name, ket_irrep_name) != 0 ||
        bra_state_number != ket_state_number) {
        tran_dm = 1;
    }

    /*
     * bra vector
     */
    int bra_dim = 0;
    double bra_eigenvalue = 0.0;
    double bra_exc_energy_cm = 0.0;
    slater_det_t *bra_det_list = NULL;
    double complex *bra_coef_left = NULL;
    double complex *bra_coef_right = NULL;

    mvcoef_read_vectors_unformatted_for_state(
        0, 2, "MVCOEF02", bra_irrep_name, bra_state_number,
        &bra_dim, &bra_det_list, &bra_eigenvalue, &bra_exc_energy_cm, &bra_coef_left, &bra_coef_right
    );

    /*
     * ket vector - if needed
     */
    int ket_dim = 0;
    double ket_eigenvalue = 0.0;
    double ket_exc_energy_cm = 0.0;
    slater_det_t *ket_det_list = NULL;
    double complex *ket_coef_left = NULL;
    double complex *ket_coef_right = NULL;

    if (tran_dm) {
        mvcoef_read_vectors_unformatted_for_state(
            0, 2, "MVCOEF02", ket_irrep_name, ket_state_number,
            &ket_dim, &ket_det_list, &ket_eigenvalue, &ket_exc_energy_cm, &ket_coef_left, &ket_coef_right
        );
    }
    else {
        ket_dim = bra_dim;
        ket_eigenvalue = bra_eigenvalue;
        ket_exc_energy_cm = bra_exc_energy_cm;
        ket_det_list = bra_det_list;
        ket_coef_left = bra_coef_left;
        ket_coef_right = bra_coef_right;
    }

    /*
     * determine the "symmetry" of the operator representing DM
     */
    int op_sym = detect_operator_symmetry(bra_irrep_name, ket_irrep_name);

    tmplt_sym(diagram_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE, op_sym);
    diagram_t *diag_dm = diagram_stack_find(diagram_name);

    for (int i = 0; i < bra_dim; i++) {
        for (int j = 0; j < ket_dim; j++) {
            slater_det_t *det_i = bra_det_list + i;
            slater_det_t *det_j = ket_det_list + j;
            moindex_t r = det_i->indices[0];
            moindex_t s = det_i->indices[1];
            moindex_t p = det_j->indices[0];
            moindex_t q = det_j->indices[1];

            double complex C_ab = bra_coef_left[i];
            double complex C_cd = ket_coef_right[j];
            double complex gamma_rspq = conj(C_ab) * C_cd;

            diagram_set_4(diag_dm, gamma_rspq, p, q, r, s);
            diagram_set_4(diag_dm, gamma_rspq, q, p, s, r);
            diagram_set_4(diag_dm, -gamma_rspq, q, p, r, s);
            diagram_set_4(diag_dm, -gamma_rspq, p, q, s, r);
        }
    }

    /*
     * clean up
     */
    cc_free(bra_det_list);
    cc_free(bra_coef_left);
    cc_free(bra_coef_right);
    if (tran_dm) {
        cc_free(ket_det_list);
        cc_free(ket_coef_left);
        cc_free(ket_coef_right);
    }
}


void finite_order_density_matrix_block_hp_0h2p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_02_L1a
    restrict_valence("s2c", "s2c_v12", "1011", 0);
    reorder("model_02_2dm", "r1", "4123");
    reorder("s2c_v12", "r2", "2341");
    mult("r2", "r1", "r3", 3);
    tmplt_sym("r4", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_hp", -0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q2a
    reorder("t2c_v12", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    mult("r1", "t1c+_v", "r3", 1);
    mult("r3", "r2", "r4", 3);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", 0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q5a
    reorder("s2c+_v12", "r1", "4123");
    reorder("t2c_v1", "r2", "2413");
    reorder("model_02_2dm", "r3", "2341");
    mult("r1", "r3", "r4", 3);
    mult("r2", "r4", "r5", 2);
    tmplt_sym("r6", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_hp", -0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q5c
    reorder("s2c+_v1", "r1", "1342");
    reorder("t2c_v2", "r2", "2413");
    reorder("model_02_2dm", "r3", "4231");
    mult("r2", "r1", "r4", 2);
    mult("r4", "r3", "r5", 3);
    tmplt_sym("r6", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_hp", -1.0, "r6");
    restore_stack_pos(pos);

    // DM_02_Q6a
    reorder("x2c_v1", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    mult("r1", "r2", "r3", 3);
    mult("t1c_v", "r3", "r4", 1);
    update("dm_hp", -0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q6d
    reorder("x2c+_v1", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("t1c", "r3", "r4", 1);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q8a
    reorder("s2c_v12", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    reorder("s1c", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update("dm_hp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q8c
    reorder("s2c_v2", "r1", "2413");
    reorder("model_02_2dm", "r2", "4231");
    mult("r2", "s1c+", "r3", 1);
    mult("r1", "r3", "r4", 3);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", -1.0, "r5");
    restore_stack_pos(pos);

    // DM_02_Q11a
    reorder("s2c", "r1", "2431");
    reorder("model_02_2dm", "r2", "4123");
    reorder("x2c+_v1", "r3", "2341");
    mult("r3", "r2", "r4", 3);
    mult("r1", "r4", "r5", 2);
    tmplt_sym("r6", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_hp", 0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q11c
    reorder("x2c+", "r1", "3412");
    reorder("s2c", "r2", "2134");
    reorder("model_02_2dm", "r3", "4312");
    mult("r2", "r1", "r4", 2);
    mult("r4", "r3", "r5", 3);
    tmplt_sym("r6", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_hp", -0.25, "r6");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_ph_0h2p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_02_L1b
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    reorder("s2c+_v12", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    mult("r2", "r1", "r3", 3);
    tmplt_sym("r4", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_ph", -0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q2b
    reorder("t2c+_v12", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    reorder("t1c_v", "r3", "21");
    mult("r1", "r3", "r4", 1);
    mult("r2", "r4", "r5", 3);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", 0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q5b
    reorder("s2c_v12", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    reorder("t2c+_v1", "r3", "2431");
    mult("r1", "r2", "r4", 3);
    mult("r3", "r4", "r5", 2);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", -0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q5d
    reorder("s2c_v1", "r1", "3124");
    reorder("t2c+_v2", "r2", "4231");
    reorder("model_02_2dm", "r3", "2413");
    mult("r2", "r1", "r4", 2);
    mult("r3", "r4", "r5", 3);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", -1.0, "r6");
    restore_stack_pos(pos);

    // DM_02_Q6b
    reorder("x2c+_v1", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    reorder("t1c+_v", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update("dm_ph", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q6c
    reorder("x2c_v1", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    reorder("t1c+", "r3", "21");
    mult("r2", "r1", "r4", 3);
    mult("r4", "r3", "r5", 1);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", -0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q8b
    reorder("s2c+_v12", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    mult("r1", "r2", "r3", 3);
    mult("s1c+", "r3", "r4", 1);
    update("dm_ph", -0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q8d
    reorder("s2c+_v2", "r1", "4321");
    reorder("model_02_2dm", "r2", "2143");
    reorder("s1c", "r3", "21");
    mult("r2", "r3", "r4", 1);
    mult("r4", "r1", "r5", 3);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", -1.0, "r6");
    restore_stack_pos(pos);

    // DM_02_Q11b
    reorder("x2c_v1", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    reorder("s2c+", "r3", "2431");
    mult("r2", "r1", "r4", 3);
    mult("r3", "r4", "r5", 2);
    update("dm_ph", 0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q11d
    reorder("s2c+", "r1", "4312");
    reorder("x2c", "r2", "3412");
    reorder("model_02_2dm", "r3", "2134");
    mult("r3", "r2", "r4", 2);
    mult("r4", "r1", "r5", 3);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", -0.25, "r6");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_hh_0h2p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_02_Q3d
    reorder("t2c+_v12", "r1", "4312");
    reorder("t2c_v12", "r0", "2134");
    mult("r1", "model_02_2dm", "r2", 2);
    mult("r0", "r2", "r3", 3);
    tmplt_sym("r4", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_hh", -0.25, "r4");
    restore_stack_pos(pos);

    // DM_02_Q4g
    reorder("s2c_v12", "r1", "2341");
    reorder("t1c+_v", "r2", "21");
    mult("r2", "model_02_2dm", "r3", 1);
    mult("r1", "r3", "r4", 3);
    tmplt_sym("r5", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hh", 0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q4h
    reorder("s2c+_v12", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    mult("r1", "r2", "r3", 3);
    mult("t1c_v", "r3", "r4", 1);
    tmplt_sym("r5", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hh", 0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q9d
    reorder("s2c+_v1", "s2c+_v1_r", "2143");
    reorder("s2c+_v1_r", "r1", "3142");
    reorder("s2c_v1", "r2", "2431");
    reorder("model_02_2dm", "r3", "2413");
    mult("r2", "r3", "r4", 2);
    mult("r4", "r1", "r5", 3);
    update("dm_hh", -1.0, "r5");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_pp_0h2p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_02_L2a
    restrict_valence("x2c", "x2c_v1", "1110", 0);
    reorder("x2c_v1", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    mult("r2", "r1", "r3", 3);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", 0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_L2b
    restrict_valence("x2c+", "x2c+_v1", "1011", 0);
    reorder("model_02_2dm", "r1", "4123");
    reorder("x2c+_v1", "r2", "2341");
    mult("r2", "r1", "r3", 3);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", 0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q3b
    reorder("t2c+_v1", "r1", "2134");
    reorder("model_02_2dm", "r2", "4312");
    mult("r2", "t2c_v12", "r3", 2);
    mult("r1", "r3", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 0.25, "r5");
    restore_stack_pos(pos);

    // DM_02_Q3c
    reorder("t2c_v1", "r1", "4312");
    reorder("model_02_2dm", "r2", "2134");
    mult("r1", "t2c+_v12", "r3", 2);
    mult("r2", "r3", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 0.25, "r5");
    restore_stack_pos(pos);

    // DM_02_Q4c
    reorder("s2c_v12", "r1", "2134");
    reorder("model_02_2dm", "r2", "4312");
    mult("r2", "r1", "r3", 3);
    mult("t1c+", "r3", "r4", 1);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q4d
    reorder("s2c+_v12", "r1", "4123");
    reorder("model_02_2dm", "r2", "2341");
    reorder("t1c", "r3", "21");
    mult("r2", "r1", "r4", 3);
    mult("r4", "r3", "r5", 1);
    tmplt_sym("r6", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_pp", -0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q4e
    reorder("s2c_v2", "r1", "3412");
    mult("r1", "t1c+_v", "r2", 1);
    mult("model_02_2dm", "r2", "r3", 3);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_02_Q4f
    reorder("model_02_2dm", "r1", "3412");
    reorder("t1c_v", "r2", "21");
    mult("s2c+_v2", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", -1.0, "r5");
    restore_stack_pos(pos);

    // DM_02_Q9b
    reorder("s2c+_v1", "r1", "1342");
    reorder("s2c", "s2c_r", "2143");
    reorder("s2c_r", "r2", "4213");
    reorder("model_02_2dm", "r3", "2431");
    mult("r2", "r1", "r4", 2);
    mult("r3", "r4", "r5", 3);
    tmplt_sym("r6", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_pp", 1.0, "r6");
    restore_stack_pos(pos);

    // DM_02_Q9c
    reorder("s2c+", "s2c+_r", "2143");
    reorder("s2c+_r", "r1", "2431");
    reorder("s2c_v1", "r2", "3124");
    reorder("model_02_2dm", "r3", "4213");
    mult("r1", "r2", "r4", 2);
    mult("r4", "r3", "r5", 3);
    tmplt_sym("r6", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_pp", 1.0, "r6");
    restore_stack_pos(pos);

    // DM_02_Q9e
    reorder("s2c+_v1", "s2c+_v1_r", "2143");
    reorder("s2c+_v1_r", "r1", "1342");
    reorder("s2c_v1", "r2", "4213");
    reorder("model_02_2dm", "r3", "2431");
    mult("r2", "r3", "r4", 2);
    mult("r1", "r4", "r5", 3);
    update("dm_pp", 1.0, "r5");
    restore_stack_pos(pos);

    // DM_02_Q10a
    reorder("s1c+", "r1", "21");
    reorder("x2c", "r2", "4123");
    reorder("model_02_2dm", "r3", "2341");
    mult("r2", "r1", "r4", 1);
    mult("r3", "r4", "r5", 3);
    tmplt_sym("r6", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_pp", 0.5, "r6");
    restore_stack_pos(pos);

    // DM_02_Q10b
    reorder("x2c+", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    mult("r1", "s1c", "r3", 1);
    mult("r3", "r2", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 0.5, "r5");
    restore_stack_pos(pos);

    // DM_02_Q10c
    reorder("x2c_v1", "r1", "4123");
    reorder("model_02_2dm", "r2", "3412");
    mult("s1c+", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update("dm_pp", 0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q10d
    reorder("x2c+_v1", "r1", "2341");
    reorder("s1c", "r2", "21");
    mult("r2", "model_02_2dm", "r3", 1);
    mult("r1", "r3", "r4", 3);
    update("dm_pp", 0.5, "r4");
    restore_stack_pos(pos);

    // DM_02_Q12
    reorder("x2c+", "r1", "2134");
    reorder("x2c", "r2", "4312");
    mult("r2", "model_02_2dm", "r3", 2);
    mult("r1", "r3", "r4", 3);
    update("dm_pp", 0.25, "r4");
    restore_stack_pos(pos);

    // renormalization term: ((T^+ T)_conn O)_conn
    sector_0h2p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
    reorder("TT_conn", "r1", "2341");
    reorder("model_02_2dm", "r2", "4123");
    mult("r1", "r2", "r3", 3);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", -0.5, "r4");
    restore_stack_pos(pos);
}










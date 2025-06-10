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

#include "finite_order_dm_0h1p.h"

#include "engine.h"
#include "options.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/mvcoef.h"
#include "finite_order_overlap.h"
#include "cc_properies.h"
#include "finite_order_overlap.h"

#include <string.h>
#include <stdlib.h>


complex double one_part_density_matrix_element_0h1p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pq
);


/**
 * constructs diagram representing the active-space one-particle density matrix:
 * D_pq = < bra model vector | p^+ q | ket model vector >
 */
void construct_one_part_active_density_matrix_diagram_0h1p(
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
     * read bra vector
     */
    int bra_dim = 0;
    double bra_eigenvalue = 0.0;
    double bra_exc_energy_cm = 0.0;
    slater_det_t *bra_det_list = NULL;
    double complex *bra_coef_left = NULL;
    double complex *bra_coef_right = NULL;

    mvcoef_read_vectors_unformatted_for_state(
        0, 1, "MVCOEF01", bra_irrep_name, bra_state_number,
        &bra_dim, &bra_det_list, &bra_eigenvalue, &bra_exc_energy_cm, &bra_coef_left, &bra_coef_right
    );

    /*
     * read ket vector - if needed
     */
    int ket_dim = 0;
    double ket_eigenvalue = 0.0;
    double ket_exc_energy_cm = 0.0;
    slater_det_t *ket_det_list = NULL;
    double complex *ket_coef_left = NULL;
    double complex *ket_coef_right = NULL;

    if (tran_dm) {
        mvcoef_read_vectors_unformatted_for_state(
            0, 1, "MVCOEF01", ket_irrep_name, ket_state_number,
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

    tmplt_sym(diagram_name, "pp", "11", "12", NOT_PERM_UNIQUE, op_sym);
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
            int indices_pq[2];
            indices_pq[0] = idx[1];
            indices_pq[1] = idx[0];

            double complex dm_elem = one_part_density_matrix_element_0h1p(
                bra_dim, bra_coef_left, bra_det_list,
                ket_dim, ket_coef_right, ket_det_list,
                indices_pq
            );

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


/**
 * the "D_pq" element of the active-space one-particle density matrix
 * for the given pair of model vectors representing electronic states
 */
complex double one_part_density_matrix_element_0h1p(
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

            // ket determinant
            int b = ket->indices[0];

            // excitation operator (p+ q)
            int p = indices_pq[0];
            int q = indices_pq[1];

            int bra_pq_ket = (a == p) * (b == q);

            matrix_element += conj(coef_bra[i]) * coef_ket[j] * bra_pq_ket;
        }
    }

    return matrix_element;
}


void finite_order_density_matrix_block_hp_0h1p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_01_L1a
    restrict_valence("t1c", "t1c_v", "01", 0);
    reorder("model_01_1dm", "r1", "21");
    mult("t1c_v", "r1", "r2", 1);
    tmplt_sym("r3", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r2", "r3");
    update("dm_hp", -1.0, "r3");
    restore_stack_pos(pos);

    // DM_01_L3a
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    reorder("s2c_v1", "r1", "2431");
    mult("r1", "model_01_1dm", "r2", 2);
    update("dm_hp", 1.0, "r2");
    restore_stack_pos(pos);

    // DM_01_Q2a
    reorder("t2c_v2", "r1", "1324");
    reorder("t1c+_v", "r2", "21");
    mult("r2", "model_01_1dm", "r3", 1);
    mult("r1", "r3", "r4", 2);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", -1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q2c
    reorder("t2c_v1", "r1", "1342");
    reorder("model_01_1dm", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", -1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q4a
    reorder("s1c", "r1", "21");
    mult("r1", "model_01_1dm", "r2", 1);
    mult("t1c_v", "r2", "r3", 1);
    tmplt_sym("r4", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_hp", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q4c
    reorder("model_01_1dm", "r1", "21");
    mult("r1", "s1c+", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    tmplt_sym("r4", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_hp", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q6b
    reorder("s2c+", "r1", "3412");
    reorder("model_01_1dm", "r2", "21");
    mult("t2c", "r1", "r3", 3);
    mult("r3", "r2", "r4", 1);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_01_Q6d
    reorder("s2c+_v1", "r1", "4231");
    reorder("t2c", "r2", "2413");
    mult("r1", "model_01_1dm", "r3", 2);
    mult("r2", "r3", "r4", 2);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", 1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_8a
    reorder("model_01_1dm", "r1", "21");
    reorder("s2c", "r2", "2431");
    mult("s1c+", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", 1.0, "r5");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_ph_0h1p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_01_L1b
    restrict_valence("t1c+", "t1c+_v", "10", 0);
    reorder("t1c+_v", "r1", "21");
    mult("model_01_1dm", "r1", "r2", 1);
    tmplt_sym("r3", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r2", "r3");
    update("dm_ph", -1.0, "r3");
    restore_stack_pos(pos);

    // DM_01_L3b
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    reorder("s2c+_v1", "r1", "2431");
    mult("r1", "model_01_1dm", "r2", 2);
    update("dm_ph", 1.0, "r2");
    restore_stack_pos(pos);

    // DM_01_Q2b
    reorder("model_01_1dm", "r1", "21");
    reorder("t2c+_v2", "r2", "1342");
    mult("t1c_v", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    tmplt_sym("r5", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_ph", -1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q2d
    reorder("t2c+_v1", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("model_01_1dm", "r2", "r3", 1);
    tmplt_sym("r4", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_ph", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q4b
    reorder("t1c+_v", "r1", "21");
    mult("r1", "model_01_1dm", "r2", 1);
    mult("s1c+", "r2", "r3", 1);
    tmplt_sym("r4", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_ph", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q4d
    reorder("t1c+", "r1", "21");
    mult("r1", "s1c", "r2", 1);
    mult("model_01_1dm", "r2", "r3", 1);
    tmplt_sym("r4", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_ph", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q6a
    reorder("t2c+", "r1", "3412");
    mult("r1", "s2c", "r2", 3);
    mult("model_01_1dm", "r2", "r3", 1);
    tmplt_sym("r4", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_ph", -0.5, "r4");
    restore_stack_pos(pos);

    // DM_01_Q6c
    reorder("t2c+", "r1", "2431");
    reorder("s2c_v1", "r2", "2431");
    mult("r2", "model_01_1dm", "r3", 2);
    mult("r1", "r3", "r4", 2);
    tmplt_sym("r5", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_ph", 1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q8b
    reorder("s2c+", "r1", "2413");
    reorder("s1c", "r2", "21");
    mult("r2", "model_01_1dm", "r3", 1);
    mult("r1", "r3", "r4", 2);
    tmplt_sym("r5", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_ph", 1.0, "r5");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_hh_0h1p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_01_Q1a
    reorder("t1c+_v", "r1", "21");
    mult("r1", "model_01_1dm", "r2", 1);
    mult("t1c_v", "r2", "r3", 1);
    update("dm_hh", 1.0, "r3");
    restore_stack_pos(pos);

    // DM_01_Q3d
    reorder("t2c+_v2", "r1", "3412");
    mult("r1", "model_01_1dm", "r2", 1);
    mult("t2c_v2", "r2", "r3", 3);
    tmplt_sym("r4", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_hh", 1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q5a
    reorder("s2c_v1", "r1", "2431");
    reorder("t1c+", "r2", "21");
    mult("r1", "model_01_1dm", "r3", 2);
    mult("r3", "r2", "r4", 1);
    tmplt_sym("r5", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hh", -1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q5b
    reorder("s2c+_v1", "r1", "4231");
    mult("r1", "model_01_1dm", "r2", 2);
    mult("t1c", "r2", "r3", 1);
    tmplt_sym("r4", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_hh", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q9b
    reorder("s2c+", "r1", "4123");
    reorder("s2c", "r2", "2341");
    mult("r2", "model_01_1dm", "r3", 1);
    mult("r3", "r1", "r4", 3);
    tmplt_sym("r5", "hh", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hh", -0.5, "r5");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_pp_0h1p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM_01_O1
    tmplt_sym("r1", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("model_01_1dm", "r1");
    update("dm_pp", 1.0, "r1");
    restore_stack_pos(pos);

    // DM_01_L2a
    reorder("s1c", "r1", "21");
    mult("model_01_1dm", "r1", "r2", 1);
    tmplt_sym("r3", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r2", "r3");
    update("dm_pp", 1.0, "r3");
    restore_stack_pos(pos);

    // DM_01_L2b
    reorder("model_01_1dm", "r1", "21");
    mult("s1c+", "r1", "r2", 1);
    tmplt_sym("r3", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r2", "r3");
    update("dm_pp", 1.0, "r3");
    restore_stack_pos(pos);

    // DM_01_Q1b
    reorder("t1c", "r1", "21");
    mult("r1", "t1c+_v", "r2", 1);
    mult("model_01_1dm", "r2", "r3", 1);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q1c
    reorder("model_01_1dm", "r1", "21");
    mult("r1", "t1c_v", "r2", 1);
    mult("t1c+", "r2", "r3", 1);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", -1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q3b
    reorder("t2c_v1", "r1", "3412");
    reorder("model_01_1dm", "r2", "21");
    mult("t2c+", "r1", "r3", 3);
    mult("r3", "r2", "r4", 1);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_01_Q3c
    reorder("t2c", "r1", "3412");
    mult("r1", "t2c+_v1", "r2", 3);
    mult("model_01_1dm", "r2", "r3", 1);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", -0.5, "r4");
    restore_stack_pos(pos);

    // DM_01_Q3e
    reorder("model_01_1dm", "r1", "21");
    mult("t2c_v2", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v2", "r3", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM_01_Q5c
    reorder("s2c_v1", "r1", "4231");
    mult("r1", "model_01_1dm", "r2", 2);
    mult("t1c+", "r2", "r3", 1);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", 1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q5d
    reorder("s2c+_v1", "r1", "2431");
    reorder("t1c", "r2", "21");
    mult("r1", "model_01_1dm", "r3", 2);
    mult("r3", "r2", "r4", 1);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q7
    reorder("s1c", "r1", "21");
    mult("r1", "model_01_1dm", "r2", 1);
    mult("s1c+", "r2", "r3", 1);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", 1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q8c
    reorder("s2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("model_01_1dm", "r2", "r3", 1);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", 1.0, "r4");
    restore_stack_pos(pos);

    // DM_01_Q8d
    reorder("s2c+", "r1", "1342");
    reorder("model_01_1dm", "r2", "21");
    mult("r1", "t1c", "r3", 2);
    mult("r3", "r2", "r4", 1);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 1.0, "r5");
    restore_stack_pos(pos);

    // DM_01_Q9a
    reorder("s2c+", "r1", "2413");
    reorder("s2c", "r2", "4231");
    mult("r2", "model_01_1dm", "r3", 1);
    mult("r1", "r3", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 1.0, "r5");
    restore_stack_pos(pos);

    // renormalization connected term
    sector_0h1p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
    reorder("model_01_1dm", "r1", "21");
    mult("TT_conn", "r1", "r2", 1);
    tmplt_sym("r3", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r2", "r3");
    update("dm_pp", -1.0, "r3");
    restore_stack_pos(pos);
}



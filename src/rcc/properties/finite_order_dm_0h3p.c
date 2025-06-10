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

#include "finite_order_dm_0h3p.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../models/ccutils.h"
#include "engine.h"
#include "../models/diis.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/mvcoef.h"
#include "finite_order_overlap.h"
#include "cc_properies.h"


complex double one_part_density_matrix_element_0h3p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pq
);


complex double two_part_density_matrix_element_0h3p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pqrs
);


complex double three_part_density_matrix_element_0h3p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pqrstu
);


void construct_n_part_active_density_matrix_diagram_0h3p(
    char *diagram_name,
    char *bra_irrep_name, int bra_state_number,
    char *ket_irrep_name, int ket_state_number,
    int n_particles
)
{
    assert(n_particles == 1 || n_particles == 2 || n_particles == 3);

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
        0, 3, "MVCOEF03", bra_irrep_name, bra_state_number,
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
            0, 3, "MVCOEF03", ket_irrep_name, ket_state_number,
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

    /*
     * create empty template for the diagram representing density matrix in the active space
     */
    if (n_particles == 1) {
        tmplt_sym(diagram_name, "pp", "11", "12", NOT_PERM_UNIQUE, op_sym);
    }
    else if (n_particles == 2) {
        tmplt_sym(diagram_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE, op_sym);
    }
    else {
        tmplt_sym(diagram_name, "pppppp", "111111", "123456", NOT_PERM_UNIQUE, op_sym);
    }
    diagram_t *diag_dm = diagram_stack_find(diagram_name);

    /*
     * loop over matrix elements of a diagram
     */
    for (int iblock = 0; iblock < diag_dm->n_blocks; iblock++) {
        block_t *block = diag_dm->blocks[iblock];
        int *indices = (int *) cc_malloc(block->size * block->rank * sizeof(int));
        block_gen_indices(block, indices);

        for (size_t i = 0; i < block->size; i++) {
            int *idx = indices + i * block->rank;

            // here we interchange pqr... and stu... blocks of indices since in EXP-T
            // the order of indices in the algebraic representation of diagrams is
            // <in|out> instead of <out|in> (which is the standard one)
            int indices_pqrstu[6];
            intcpy(indices_pqrstu, idx + n_particles, n_particles);
            intcpy(indices_pqrstu + n_particles, idx, n_particles);

            double complex dm_elem = 0.0;
            if (n_particles == 1) {
                dm_elem = one_part_density_matrix_element_0h3p(
                    bra_dim, bra_coef_left, bra_det_list,
                    ket_dim, ket_coef_right, ket_det_list,
                    indices_pqrstu
                );
            }
            else if (n_particles == 2) {
                dm_elem = two_part_density_matrix_element_0h3p(
                    bra_dim, bra_coef_left, bra_det_list,
                    ket_dim, ket_coef_right, ket_det_list,
                    indices_pqrstu
                );
            }
            else {
                dm_elem = three_part_density_matrix_element_0h3p(
                    bra_dim, bra_coef_left, bra_det_list,
                    ket_dim, ket_coef_right, ket_det_list,
                    indices_pqrstu
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


complex double one_part_density_matrix_element_0h3p(
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
            int c = bra->indices[2];

            // ket determinant
            int d = ket->indices[0];
            int e = ket->indices[1];
            int f = ket->indices[2];

            // excitation operator (p+ q)
            int p = indices_pq[0];
            int q = indices_pq[1];

            int bra_pq_ket =
                -(c == p) * (b == d) * (a == e) * (f == q)
                + (c == p) * (b == d) * (a == f) * (e == q)
                + (c == p) * (b == e) * (a == d) * (f == q)
                - (c == p) * (b == e) * (a == f) * (d == q)
                - (c == p) * (b == f) * (a == d) * (e == q)
                + (c == p) * (b == f) * (a == e) * (d == q)
                + (c == d) * (b == p) * (a == e) * (f == q)
                - (c == d) * (b == p) * (a == f) * (e == q)
                - (c == d) * (b == e) * (a == p) * (f == q)
                + (c == d) * (b == f) * (a == p) * (e == q)
                - (c == e) * (b == p) * (a == d) * (f == q)
                + (c == e) * (b == p) * (a == f) * (d == q)
                + (c == e) * (b == d) * (a == p) * (f == q)
                - (c == e) * (b == f) * (a == p) * (d == q)
                + (c == f) * (b == p) * (a == d) * (e == q)
                - (c == f) * (b == p) * (a == e) * (d == q)
                - (c == f) * (b == d) * (a == p) * (e == q)
                + (c == f) * (b == e) * (a == p) * (d == q);

            matrix_element += conj(coef_bra[i]) * coef_ket[j] * bra_pq_ket;
        }
    }

    return matrix_element;
}


complex double two_part_density_matrix_element_0h3p(
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
            int c = bra->indices[2];

            // ket determinant
            int d = ket->indices[0];
            int e = ket->indices[1];
            int f = ket->indices[2];

            // excitation operator (p+ q+ s r)
            int p = indices_pqrs[0];
            int q = indices_pqrs[1];
            int r = indices_pqrs[2];
            int s = indices_pqrs[3];

            int bra_pqrs_ket =
                +(c == p) * (b == q) * (a == d) * (e == s) * (f == r)
                - (c == p) * (b == q) * (a == d) * (f == s) * (e == r)
                - (c == p) * (b == q) * (a == e) * (d == s) * (f == r)
                + (c == p) * (b == q) * (a == e) * (f == s) * (d == r)
                + (c == p) * (b == q) * (a == f) * (d == s) * (e == r)
                - (c == p) * (b == q) * (a == f) * (e == s) * (d == r)
                - (c == p) * (b == d) * (a == q) * (e == s) * (f == r)
                + (c == p) * (b == d) * (a == q) * (f == s) * (e == r)
                + (c == p) * (b == e) * (a == q) * (d == s) * (f == r)
                - (c == p) * (b == e) * (a == q) * (f == s) * (d == r)
                - (c == p) * (b == f) * (a == q) * (d == s) * (e == r)
                + (c == p) * (b == f) * (a == q) * (e == s) * (d == r)
                - (c == q) * (b == p) * (a == d) * (e == s) * (f == r)
                + (c == q) * (b == p) * (a == d) * (f == s) * (e == r)
                + (c == q) * (b == p) * (a == e) * (d == s) * (f == r)
                - (c == q) * (b == p) * (a == e) * (f == s) * (d == r)
                - (c == q) * (b == p) * (a == f) * (d == s) * (e == r)
                + (c == q) * (b == p) * (a == f) * (e == s) * (d == r)
                + (c == q) * (b == d) * (a == p) * (e == s) * (f == r)
                - (c == q) * (b == d) * (a == p) * (f == s) * (e == r)
                - (c == q) * (b == e) * (a == p) * (d == s) * (f == r)
                + (c == q) * (b == e) * (a == p) * (f == s) * (d == r)
                + (c == q) * (b == f) * (a == p) * (d == s) * (e == r)
                - (c == q) * (b == f) * (a == p) * (e == s) * (d == r)
                + (c == d) * (b == p) * (a == q) * (e == s) * (f == r)
                - (c == d) * (b == p) * (a == q) * (f == s) * (e == r)
                - (c == d) * (b == q) * (a == p) * (e == s) * (f == r)
                + (c == d) * (b == q) * (a == p) * (f == s) * (e == r)
                - (c == e) * (b == p) * (a == q) * (d == s) * (f == r)
                + (c == e) * (b == p) * (a == q) * (f == s) * (d == r)
                + (c == e) * (b == q) * (a == p) * (d == s) * (f == r)
                - (c == e) * (b == q) * (a == p) * (f == s) * (d == r)
                + (c == f) * (b == p) * (a == q) * (d == s) * (e == r)
                - (c == f) * (b == p) * (a == q) * (e == s) * (d == r)
                - (c == f) * (b == q) * (a == p) * (d == s) * (e == r)
                + (c == f) * (b == q) * (a == p) * (e == s) * (d == r);

            matrix_element += conj(coef_bra[i]) * coef_ket[j] * bra_pqrs_ket;
        }
    }

    return matrix_element;
}


complex double three_part_density_matrix_element_0h3p(
    int dim_bra, double complex *coef_bra, slater_det_t *dets_bra,
    int dim_ket, double complex *coef_ket, slater_det_t *dets_ket,
    int *indices_pqrstu
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
            int c = bra->indices[2];

            // ket determinant
            int d = ket->indices[0];
            int e = ket->indices[1];
            int f = ket->indices[2];

            // excitation operator (p+ q+ r+ u t s)
            int p = indices_pqrstu[0];
            int q = indices_pqrstu[1];
            int r = indices_pqrstu[2];
            int s = indices_pqrstu[3];
            int t = indices_pqrstu[4];
            int u = indices_pqrstu[5];

            int bra_pqrstu_ket =
                +(c == p) * (b == q) * (a == r) * (d == u) * (e == t) * (f == s)
                - (c == p) * (b == q) * (a == r) * (d == u) * (f == t) * (e == s)
                - (c == p) * (b == q) * (a == r) * (e == u) * (d == t) * (f == s)
                + (c == p) * (b == q) * (a == r) * (e == u) * (f == t) * (d == s)
                + (c == p) * (b == q) * (a == r) * (f == u) * (d == t) * (e == s)
                - (c == p) * (b == q) * (a == r) * (f == u) * (e == t) * (d == s)
                - (c == p) * (b == r) * (a == q) * (d == u) * (e == t) * (f == s)
                + (c == p) * (b == r) * (a == q) * (d == u) * (f == t) * (e == s)
                + (c == p) * (b == r) * (a == q) * (e == u) * (d == t) * (f == s)
                - (c == p) * (b == r) * (a == q) * (e == u) * (f == t) * (d == s)
                - (c == p) * (b == r) * (a == q) * (f == u) * (d == t) * (e == s)
                + (c == p) * (b == r) * (a == q) * (f == u) * (e == t) * (d == s)
                - (c == q) * (b == p) * (a == r) * (d == u) * (e == t) * (f == s)
                + (c == q) * (b == p) * (a == r) * (d == u) * (f == t) * (e == s)
                + (c == q) * (b == p) * (a == r) * (e == u) * (d == t) * (f == s)
                - (c == q) * (b == p) * (a == r) * (e == u) * (f == t) * (d == s)
                - (c == q) * (b == p) * (a == r) * (f == u) * (d == t) * (e == s)
                + (c == q) * (b == p) * (a == r) * (f == u) * (e == t) * (d == s)
                + (c == q) * (b == r) * (a == p) * (d == u) * (e == t) * (f == s)
                - (c == q) * (b == r) * (a == p) * (d == u) * (f == t) * (e == s)
                - (c == q) * (b == r) * (a == p) * (e == u) * (d == t) * (f == s)
                + (c == q) * (b == r) * (a == p) * (e == u) * (f == t) * (d == s)
                + (c == q) * (b == r) * (a == p) * (f == u) * (d == t) * (e == s)
                - (c == q) * (b == r) * (a == p) * (f == u) * (e == t) * (d == s)
                + (c == r) * (b == p) * (a == q) * (d == u) * (e == t) * (f == s)
                - (c == r) * (b == p) * (a == q) * (d == u) * (f == t) * (e == s)
                - (c == r) * (b == p) * (a == q) * (e == u) * (d == t) * (f == s)
                + (c == r) * (b == p) * (a == q) * (e == u) * (f == t) * (d == s)
                + (c == r) * (b == p) * (a == q) * (f == u) * (d == t) * (e == s)
                - (c == r) * (b == p) * (a == q) * (f == u) * (e == t) * (d == s)
                - (c == r) * (b == q) * (a == p) * (d == u) * (e == t) * (f == s)
                + (c == r) * (b == q) * (a == p) * (d == u) * (f == t) * (e == s)
                + (c == r) * (b == q) * (a == p) * (e == u) * (d == t) * (f == s)
                - (c == r) * (b == q) * (a == p) * (e == u) * (f == t) * (d == s)
                - (c == r) * (b == q) * (a == p) * (f == u) * (d == t) * (e == s)
                + (c == r) * (b == q) * (a == p) * (f == u) * (e == t) * (d == s);

            matrix_element += conj(coef_bra[i]) * coef_ket[j] * bra_pqrstu_ket;
        }
    }

    return matrix_element;
}


void finite_order_density_matrix_block_hp_0h3p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // Q1a
    reorder("model_03_3dm", "r1", "451623");
    reorder("s2c_v12", "s2c_v12_r", "2143");
    reorder("x2c_v1", "r2", "4123");
    mult("s2c_v12_r", "r1", "r3", 3);
    mult("r3", "r2", "r4", 3);
    update("dm_hp", -0.25, "r4");
    restore_stack_pos(pos);

    // Q2b
    reorder("model_03_3dm", "r1", "412356");
    reorder("s2c+_v12", "r2", "4321");
    reorder("t2c_v12", "r3", "1342");
    mult("r1", "r2", "r4", 3);
    mult("r3", "r4", "r5", 3);
    tmplt_sym("r6", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_hp", 0.25, "r6");
    restore_stack_pos(pos);

    // Q3a
    reorder("x2c+_v1", "x2c+_v1_r", "2143");
    reorder("model_03_3dm", "r1", "541623");
    reorder("s2c_v1", "r2", "2134");
    mult("r1", "x2c+_v1_r", "r3", 3);
    mult("r2", "r3", "r4", 3);
    tmplt_sym("r5", "hp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_hp", -0.5, "r5");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_ph_0h3p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // Q1b
    reorder("model_03_3dm", "r1", "124356");
    reorder("x2c+_v1", "r2", "2341");
    reorder("s2c+_v12", "r3", "4321");
    mult("r3", "r1", "r4", 3);
    mult("r2", "r4", "r5", 3);
    update("dm_ph", -0.25, "r5");
    restore_stack_pos(pos);

    // Q2c
    reorder("model_03_3dm", "r1", "356124");
    reorder("s2c_v12", "r2", "2341");
    reorder("t2c+_v12", "r3", "4123");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 3);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", 0.25, "r6");
    restore_stack_pos(pos);

    // Q3b
    reorder("model_03_3dm", "r1", "263451");
    reorder("x2c_v1", "r2", "4123");
    reorder("s2c+_v1", "s2c+_v1_r", "2143");
    reorder("s2c+_v1_r", "r3", "3241");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 3);
    tmplt_sym("r6", "ph", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_ph", -0.5, "r6");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_hh_0h3p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // Q2a
    reorder("model_03_3dm", "r1", "124356");
    reorder("s2c+_v12", "r2", "4321");
    reorder("s2c_v12", "r3", "2341");
    mult("r2", "r1", "r4", 3);
    mult("r3", "r4", "r5", 3);
    update("dm_hh", 0.25, "r5");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_pp_0h3p(int dm_sym)
{
    dg_stack_pos_t pos = get_stack_pos();

    // Q2d
    reorder("model_03_3dm", "r1", "124356");
    reorder("s2c_v2", "r2", "3412");
    reorder("s2c+_v12", "r3", "4321");
    mult("r1", "r3", "r4", 3);
    mult("r4", "r2", "r5", 3);
    tmplt_sym("r6", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_pp", -0.5, "r6");
    restore_stack_pos(pos);

    // Q2e
    reorder("model_03_3dm", "r1", "653124");
    reorder("s2c_v12", "r2", "2341");
    mult("r1", "r2", "r3", 3);
    mult("s2c+_v2", "r3", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", -0.5, "r5");
    restore_stack_pos(pos);

    // Q4a
    reorder("model_03_3dm", "r1", "145623");
    reorder("x2c", "r2", "3124");
    reorder("x2c+_v1", "x2c+_v1_r", "2143");
    mult("r1", "x2c+_v1_r", "r3", 3);
    mult("r3", "r2", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 0.25, "r5");
    restore_stack_pos(pos);

    // Q4b
    reorder("model_03_3dm", "r1", "623451");
    reorder("x2c_v1", "r2", "4123");
    reorder("x2c+", "r3", "2341");
    mult("r1", "r2", "r4", 3);
    mult("r3", "r4", "r5", 3);
    tmplt_sym("r6", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r5", "r6");
    update("dm_pp", 0.25, "r6");
    restore_stack_pos(pos);

    // Q4c
    reorder("model_03_3dm", "r1", "623451");
    reorder("x2c_v1", "r2", "4123");
    reorder("x2c+_v1", "x2c+_v1_r", "2143");
    mult("r2", "r1", "r3", 3);
    mult("x2c+_v1_r", "r3", "r4", 3);
    tmplt_sym("r5", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r4", "r5");
    update("dm_pp", 0.25, "r5");
    restore_stack_pos(pos);

    // renormalization term
    sector_0h3p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
    reorder("TT_conn", "r1", "345612");
    reorder("model_03_3dm", "r2", "612345");
    mult("r1", "r2", "r3", 5);
    tmplt_sym("r4", "pp", "00", "12", NOT_PERM_UNIQUE, dm_sym);
    expand_diagram("r3", "r4");
    update("dm_pp", -1.0/12.0, "r4");
    restore_stack_pos(pos);
}



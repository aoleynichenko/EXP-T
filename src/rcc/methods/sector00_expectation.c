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

/**
 * Implementation of the density matrix construction in the 0h0p sector
 * via the expectation value formula:
 * D_pq = <0|(e^{-T^+} {p^+ q} e^T)_conn|0>
 */

#include "cc_properies.h"
#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ccutils.h"
#include "engine.h"
#include "datamodel.h"
#include "diis.h"
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void sector_0h0p_construct_density_matrix_block_hp(int pt_order);

void sector_0h0p_construct_density_matrix_block_hh(int pt_order);

void sector_0h0p_construct_density_matrix_block_pp(int pt_order);

void sector_0h0p_construct_density_matrix_block_hp_T1();

void sector_0h0p_construct_density_matrix_block_hp_T1x_T2();

void sector_0h0p_construct_density_matrix_block_hp_T1x_T1_T1();

void sector_0h0p_construct_density_matrix_block_hp_T2x_T1_T2();

void sector_0h0p_construct_density_matrix_block_hp_T2x_T3();

void sector_0h0p_construct_density_matrix_block_hp_T2x_T1_T1_T1();

void sector_0h0p_construct_density_matrix_block_hp_T1x_T1x_T3();

void sector_0h0p_construct_density_matrix_block_hp_T1x_T1x_T1_T2();

void sector_0h0p_construct_density_matrix_block_hp_T3x_T2_T2();

void sector_0h0p_construct_density_matrix_block_hp_T1x_T2x_T2_T2();

void sector_0h0p_construct_density_matrix_block_hh_T1x_T1();

void sector_0h0p_construct_density_matrix_block_hh_T2x_T2();

void sector_0h0p_construct_density_matrix_block_hh_T2x_T1_T1();

void sector_0h0p_construct_density_matrix_block_hh_T1x_T1x_T2();

void sector_0h0p_construct_density_matrix_block_hh_T1x_T1x_T1_T1();

void sector_0h0p_construct_density_matrix_block_hh_T3x_T3();

void sector_0h0p_construct_density_matrix_block_hh_T3x_T1_T2();

void sector_0h0p_construct_density_matrix_block_hh_T1x_T2x_T3();

void sector_0h0p_construct_density_matrix_block_hh_T1x_T2x_T1_T2();

void sector_0h0p_construct_density_matrix_block_hh_T2x_T2x_T2_T2();

void sector_0h0p_construct_density_matrix_block_pp_T1x_T1();

void sector_0h0p_construct_density_matrix_block_pp_T2x_T2();

void sector_0h0p_construct_density_matrix_block_pp_T2x_T1_T1();

void sector_0h0p_construct_density_matrix_block_pp_T1x_T1x_T2();

void sector_0h0p_construct_density_matrix_block_pp_T1x_T1x_T1_T1();

void sector_0h0p_construct_density_matrix_block_pp_T3x_T3();

void sector_0h0p_construct_density_matrix_block_pp_T3x_T1_T2();

void sector_0h0p_construct_density_matrix_block_pp_T1x_T2x_T3();

void sector_0h0p_construct_density_matrix_block_pp_T1x_T2x_T1_T2();

void sector_0h0p_construct_density_matrix_block_pp_T2x_T2x_T2_T2();

double complex sector_0h0p_prop_linear();

double complex norm_0h0p_T1x_T1();

double complex norm_0h0p_T2x_T2();

double complex property_0h0p_D1();

double complex property_0h0p_D2();

double complex property_0h0p_D3();

double complex property_0h0p_D4();

double complex property_0h0p_D5();

double complex property_0h0p_D6();

double complex property_0h0p_D7();

double complex property_0h0p_D8();

double complex property_0h0p_D9();

double complex property_0h0p_D10();

double complex property_0h0p_D11();

double complex property_0h0p_D12();

double complex property_0h0p_D13();

double complex property_0h0p_D14();

double complex property_0h0p_D15();

double complex property_0h0p_D16();

double complex property_0h0p_D17();

double complex property_0h0p_D18();

double complex property_0h0p_D19();

double complex property_0h0p_D20();

void save_density_matrix(char *path, int nspinors, double complex *dm);

int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

void guess_operator_symmetry(int nspinors, double complex *prop_matrix);

double complex calculate_vacuum_expectation_value(int nspinors, double complex *prop_matrix);

double complex sector_0h0p_calculate_overlap(int nspinors, double complex *dm);

void diagram_conjugate(char *source_name, char *target_name);

void sector_0h0p_analytic_density_matrix_expectation(int pt_order);

void sector_0h0p_calculate_natural_spinors();

void conjugate_t2();


/**
 * Constructs matrix elements of the 1-particle density matrix in the basis of
 * one-electron functions. Terms are included up to the given order of perturbation
 * theory ('pt_order'). Currently max PT order is 4.
 */
void sector_0h0p_analytic_density_matrix_expectation(int pt_order)
{
    double t_start = abs_time();
    printf(" begin construction of density matrix\n");
    printf(" max PT order: %d\n", pt_order);

    int nspinors = get_num_spinors();
    double complex *dm = z_zeros(nspinors, nspinors);
    memset(dm, 0, sizeof(double complex) * nspinors * nspinors);

    /*
     * preparation step
     */
    diagram_conjugate("t1c", "t1c+");
    //diagram_conjugate("t2c", "t2c+");
    conjugate_t2();
    if (triples_enabled()) {
        diagram_conjugate("t3c", "t3c+");
    }

    /*
     * construct matrices for the hp, hh and pp blocks of the DM
     */
    sector_0h0p_construct_density_matrix_block_hp(pt_order);   // => "dm00_hp"
    sector_0h0p_construct_density_matrix_block_hh(pt_order);   // => "dm00_hh"
    sector_0h0p_construct_density_matrix_block_pp(pt_order);   // => "dm00_pp"

    /*
     * extract matrix elements from diagrams and save them into matrix
     */
    diagram_t *dg_dm00_hp = diagram_stack_find("dm00_hp");
    diagram_t *dg_dm00_hh = diagram_stack_find("dm00_hh");
    diagram_t *dg_dm00_pp = diagram_stack_find("dm00_pp");

    for (int p = 0; p < nspinors; p++) {
        for (int q = 0; q < nspinors; q++) {
            double complex gamma_pq = 0.0 + 0.0 * I;

            if (is_hole(p) && is_hole(q)) { // h-h
                gamma_pq = (p == q) + diagram_get_2(dg_dm00_hh, p, q);
            }
            else if (is_hole(p) && is_particle(q)) { // h-p
                gamma_pq = diagram_get_2(dg_dm00_hp, p, q);
            }
            else if (is_particle(p) && is_hole(q)) { // p-h
                double complex gamma_qp = diagram_get_2(dg_dm00_hp, q, p);
                gamma_pq = conj(gamma_qp);
            }
            else { // p-p
                gamma_pq = diagram_get_2(dg_dm00_pp, p, q);
            }

            dm[p * nspinors + q] = gamma_pq;
        }
    }

    /*
     * save density matrix to disk
     */
    save_density_matrix("density_0h0p.txt", nspinors, dm);

    /*
     * calculate overlap integral <psi|psi> and wavefunction norm
     */
    double complex trace = ztrace(nspinors, dm);
    double complex norm2 = sector_0h0p_calculate_overlap(nspinors, dm);
    double norm = sqrt(creal(norm2));

    printf(" density matrix constructed in %.2f seconds\n", abs_time() - t_start);
    printf(" sum of diagonal elements   = %20.12f %20.12f\n", creal(trace), cimag(trace));
    printf(" overlap integral <psi|psi> = %20.12f %20.12f\n", creal(norm2), cimag(norm2));
    printf(" norm                       = %20.12f\n", norm);
    printf("\n");

    sector_0h0p_calculate_natural_spinors();

    // clean up
    cc_free(dm);
}


/**
 * hole-particle block of the 1-electron density matrix
 * (terms up to 4th order of PT).
 * the particle-hole block can be obtained by complex conjugation of the hole-particle block.
 */
void sector_0h0p_construct_density_matrix_block_hp(int pt_order)
{
    tmplt("dm00_hp", "hp", "00", "12", NOT_PERM_UNIQUE);

    /*
     * 1st PT order
     */
    if (pt_order >= PT_1) {
        sector_0h0p_construct_density_matrix_block_hp_T1();
    }

    /*
     * 2nd PT order
     */
    if (pt_order >= PT_2) {
        sector_0h0p_construct_density_matrix_block_hp_T1x_T2();
    }

    /*
     * 3rd PT order
     */
    if (pt_order >= PT_3) {
        sector_0h0p_construct_density_matrix_block_hp_T1x_T1_T1();
        sector_0h0p_construct_density_matrix_block_hp_T2x_T1_T2();

        if (triples_enabled()) {
            sector_0h0p_construct_density_matrix_block_hp_T2x_T3();
        }
    }

    /*
     * 4th PT order
     */
    if (pt_order >= PT_4) {
        sector_0h0p_construct_density_matrix_block_hp_T2x_T1_T1_T1();
        sector_0h0p_construct_density_matrix_block_hp_T1x_T1x_T1_T2();
        sector_0h0p_construct_density_matrix_block_hp_T1x_T2x_T2_T2();

        if (triples_enabled()) {
            sector_0h0p_construct_density_matrix_block_hp_T1x_T1x_T3();
            sector_0h0p_construct_density_matrix_block_hp_T3x_T2_T2();
        }
    }
}


void sector_0h0p_construct_density_matrix_block_hp_T1()
{
    // DM-1
    update("dm00_hp", 1.0, "t1c");
}


void sector_0h0p_construct_density_matrix_block_hp_T1x_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-2
    reorder("t2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    update("dm00_hp", 1.0, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T1x_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-3
    reorder("t1c", "r1", "21");
    reorder("t1c+", "r0", "21");
    mult("t1c", "r0", "r2", 1);
    mult("r2", "r1", "r3", 1);
    update("dm00_hp", -1.0, "r3");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T2x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-4
    reorder("t2c", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "t1c", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update("dm00_hp", 1.0, "r4");
    restore_stack_pos(pos);

    // DM-5
    reorder("t1c", "t1cr", "21");
    reorder("t2c+", "r0", "3412");
    mult("t2c", "r0", "r1", 3);
    mult("r1", "t1cr", "r2", 1);
    update("dm00_hp", -0.5, "r2");
    restore_stack_pos(pos);

    // DM-6
    reorder("t2c", "r0", "3412");
    mult("r0", "t2c+", "r1", 3);
    mult("t1c", "r1", "r2", 1);
    update("dm00_hp", -0.5, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T2x_T3()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-7
    reorder("t3c", "r1", "145623");
    mult("r1", "t2c+", "r2", 4);
    update("dm00_hp", 0.25, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T2x_T1_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-8
    reorder("t1c", "t1cr", "21");
    reorder("t2c+", "r1", "2413");
    mult("r1", "t1cr", "r2", 2);
    reorder("r2", "r3", "21");
    mult("t1c", "r3", "r4", 1);
    mult("r4", "t1cr", "r5", 1);
    update("dm00_hp", -1.0, "r5");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T1x_T1x_T3()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-9
    reorder("t3c", "r1", "145263");
    mult("r1", "t1c+", "r2", 2);
    mult("r2", "t1c+", "r3", 2);
    update("dm00_hp", 0.5, "r3");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T1x_T1x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-10
    reorder("t2c", "r1", "1324");
    reorder("t1c+", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r3", "t1c+", "r4", 1);
    mult("r1", "r4", "r5", 2);
    update("dm00_hp", -1.0, "r5");
    restore_stack_pos(pos);

    // DM-11
    reorder("t2c", "r2", "1342");
    reorder("t1c", "r1", "21");
    reorder("t1c+", "r3", "21");
    mult("r2", "t1c+", "r4", 2);
    mult("r4", "r3", "r5", 1);
    mult("r5", "r1", "r6", 1);
    update("dm00_hp", -0.5, "r6");
    restore_stack_pos(pos);

    // DM-12
    reorder("t2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("r2", "t1c+", "r3", 1);
    mult("t1c", "r3", "r4", 1);
    update("dm00_hp", -0.5, "r4");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T3x_T2_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-13 - CCSDT
    reorder("t2c", "r1", "1324");
    reorder("t3c+", "r2", "415623");
    mult("r2", "t2c", "r3", 4);
    mult("r1", "r3", "r4", 2);
    update("dm00_hp", 0.25, "r4");
    restore_stack_pos(pos);

    // DM-14 - CCSDT
    reorder("t2c", "r1", "2134");
    reorder("t2c", "r2", "3412");
    reorder("t3c+", "r3", "412356");
    mult("r3", "r2", "r4", 3);
    reorder("r4", "r5", "4123");
    mult("r1", "r5", "r6", 3);
    update("dm00_hp", -0.25, "r6");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hp_T1x_T2x_T2_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-15
    reorder("t2c", "r1", "1324");
    reorder("t2c", "r2", "1342");
    reorder("t2c+", "r3", "3142");
    mult("r2", "t1c+", "r4", 2);
    mult("r3", "r4", "r5", 2);
    mult("r1", "r5", "r6", 2);
    update("dm00_hp", 1.0, "r6");
    restore_stack_pos(pos);

    // DM-16
    reorder("t2c+", "r1", "3412");
    reorder("t2c", "r2", "1324");
    mult("r1", "t2c", "r3", 3);
    mult("r3", "t1c+", "r4", 1);
    mult("r2", "r4", "r5", 2);
    update("dm00_hp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM-17
    reorder("t2c", "r1", "3412");
    reorder("t2c", "r2", "1324");
    reorder("t1c+", "r3", "21");
    mult("r1", "t2c+", "r4", 3);
    mult("r2", "r4", "r5", 1);
    mult("r5", "r3", "r6", 2);
    update("dm00_hp", -0.5, "r6");
    restore_stack_pos(pos);

    // DM-18
    reorder("t2c+", "r1", "3412");
    reorder("t2c", "r2", "4231");
    mult("t2c", "r1", "r3", 3);
    mult("r2", "t1c+", "r4", 2);
    mult("r3", "r4", "r5", 1);
    update("dm00_hp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM-19
    reorder("t2c", "r1", "3412");
    reorder("t2c", "r2", "2431");
    mult("r1", "t2c+", "r3", 3);
    mult("r2", "t1c+", "r4", 2);
    mult("r4", "r3", "r5", 1);
    update("dm00_hp", -0.5, "r5");
    restore_stack_pos(pos);
}


/**
 * hole-hole block of the 1-electron density matrix
 * (terms up to 4th order of PT)
 */
void sector_0h0p_construct_density_matrix_block_hh(int pt_order)
{
    tmplt("dm00_hh", "hh", "00", "12", NOT_PERM_UNIQUE);

    /*
     * 2nd PT order
     */
    if (pt_order >= PT_2) {
        sector_0h0p_construct_density_matrix_block_hh_T1x_T1();
        sector_0h0p_construct_density_matrix_block_hh_T2x_T2();
    }

    /*
     * 3rd PT order
     */
    if (pt_order >= PT_3) {
        sector_0h0p_construct_density_matrix_block_hh_T2x_T1_T1();
        sector_0h0p_construct_density_matrix_block_hh_T1x_T1x_T2();
    }

    /*
     * 4th PT order
     */
    if (pt_order >= PT_4) {
        sector_0h0p_construct_density_matrix_block_hh_T1x_T1x_T1_T1();
        sector_0h0p_construct_density_matrix_block_hh_T1x_T2x_T1_T2();
        sector_0h0p_construct_density_matrix_block_hh_T2x_T2x_T2_T2();

        if (triples_enabled()) {
            sector_0h0p_construct_density_matrix_block_hh_T3x_T3();
            sector_0h0p_construct_density_matrix_block_hh_T3x_T1_T2();
            sector_0h0p_construct_density_matrix_block_hh_T1x_T2x_T3();
        }
    }
}


void sector_0h0p_construct_density_matrix_block_hh_T1x_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-20
    reorder("t1c+", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    update("dm00_hh", -1.0, "r2");
    restore_stack_pos(pos);
}


void conjugate_t2()
{
    tmplt("t2c+", "pphh", "0000", "1234", IS_PERM_UNIQUE);

    diagram_t *diag_t2c = diagram_stack_find("t2c");
    diagram_t *diag_t2conj = diagram_stack_find("t2c+");



    int rank = 4;

    for (size_t iblock = 0; iblock < diag_t2conj->n_blocks; iblock++) {
        block_t *block = diag_t2conj->blocks[iblock];
        block_load(block);
        if (block->is_unique == 0) {
            continue;
        }

        int *indices = (int *) cc_malloc(block->size * rank * sizeof(int));
        block_gen_indices(block, indices);

        for (size_t i = 0; i < block->size; i++) {
            int idx_a = indices[4 * i + 0];
            int idx_b = indices[4 * i + 1];
            int idx_i = indices[4 * i + 2];
            int idx_j = indices[4 * i + 3];

            double complex t_ijab = diagram_get_4(diag_t2c, idx_i, idx_j, idx_a, idx_b);
            double complex t_abij = conj(t_ijab);

            if (arith == CC_ARITH_COMPLEX) {
                block->buf[i] = t_abij;
            }
            else {
                ((double *) block->buf)[i] = creal(t_abij);
            }
        }

        cc_free(indices);

        block_store(block);
    }
}


void sector_0h0p_construct_density_matrix_block_hh_T2x_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-21
    reorder("t2c+", "r1", "3412");
    mult("t2c", "r1", "r2", 3);
    update("dm00_hh", -0.5, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hh_T2x_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-22
    reorder("t2c+", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("t1c", "r2", "r3", 1);
    update("dm00_hh", -1.0, "r3");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hh_T1x_T1x_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-23
    reorder("t1c+", "r1", "21");
    reorder("t2c", "r2", "1342");
    mult("r2", "t1c+", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update("dm00_hh", -1.0, "r4");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hh_T1x_T1x_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-24
    reorder("t1c", "r1", "21");
    reorder("t1c+", "r2", "21");
    reorder("t1c+", "r3", "21");
    mult("t1c", "r3", "r4", 1);
    mult("r4", "r1", "r5", 1);
    mult("r5", "r2", "r6", 1);
    update("dm00_hh", 1.0, "r6");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hh_T3x_T3()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-25 - CCSDT
    reorder("t3c+", "r1", "456123");
    mult("t3c", "r1", "r2", 5);
    update("dm00_hh", -1.0 / 12, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hh_T3x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-26
    reorder("t3c+", "r1", "451263");
    mult("r1", "t1c", "r2", 2);
    mult("t2c", "r2", "r3", 3);
    update("dm00_hh", -0.5, "r3");
    restore_stack_pos(pos);

    // DM-27
    reorder("t3c+", "r1", "415623");
    mult("r1", "t2c", "r2", 4);
    mult("t1c", "r2", "r3", 1);
    update("dm00_hh", -0.25, "r3");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_hh_T1x_T2x_T3()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-28
    reorder("t2c+", "r1", "3412");
    reorder("t3c", "r2", "124563");
    mult("r2", "t1c+", "r3", 2);
    mult("r3", "r1", "r4", 3);
    update("dm00_hh", -0.5, "r4");
    restore_stack_pos(pos);

    // DM-29
    reorder("t3c", "r1", "145623");
    reorder("t1c+", "r2", "21");
    mult("r1", "t2c+", "r3", 4);
    mult("r3", "r2", "r4", 1);
    update("dm00_hh", -0.25, "r4");
    restore_stack_pos(pos);
}


/**
 * Note that the factors for the DM-38 and DM-39 in the Bartlett's book are incorrect:
 *         Bartlett      EXP-T
 * DM-33      +1         +1/2
 * DM-34    +1/2           +1
 */
void sector_0h0p_construct_density_matrix_block_hh_T1x_T2x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-30
    reorder("t2c+", "r1", "3142");
    reorder("t2c", "r2", "1342");
    mult("r2", "t1c+", "r3", 2);
    mult("r1", "r3", "r4", 2);
    mult("t1c", "r4", "r5", 1);
    update("dm00_hh", -1.0, "r5");
    restore_stack_pos(pos);

    // DM-31
    reorder("t1c+", "r1", "21");
    reorder("t2c", "r2", "1324");
    reorder("t2c+", "r3", "3142");
    mult("r3", "t1c", "r4", 2);
    mult("r2", "r4", "r5", 2);
    mult("r5", "r1", "r6", 1);
    update("dm00_hh", -1.0, "r6");
    restore_stack_pos(pos);

    // DM-32
    reorder("t2c+", "r1", "3142");
    reorder("t2c", "r2", "2431");
    mult("r1", "t1c", "r3", 2);
    mult("r2", "t1c+", "r4", 2);
    mult("r4", "r3", "r5", 1);
    update("dm00_hh", -1.0, "r5");
    restore_stack_pos(pos);

    // DM-33
    reorder("t2c", "r1", "1342");
    reorder("t2c+", "r2", "3124");
    mult("r1", "t1c+", "r3", 1);
    mult("r3", "t1c", "r4", 1);
    mult("r4", "r2", "r5", 3);
    update("dm00_hh", 0.5, "r5");
    restore_stack_pos(pos);

    // DM-34
    reorder("t2c+", "r1", "3412");
    mult("r1", "t1c", "r2", 1);
    mult("r2", "t1c+", "r3", 1);
    mult("t2c", "r3", "r4", 3);
    update("dm00_hh", 1.0, "r4");
    restore_stack_pos(pos);

    // DM-35
    reorder("t2c", "r1", "3412");
    reorder("t1c+", "r2", "21");
    mult("t2c+", "r1", "r3", 3);
    mult("r2", "r3", "r4", 1);
    mult("t1c", "r4", "r5", 1);
    update("dm00_hh", 0.5, "r5");
    restore_stack_pos(pos);

    // DM-36
    reorder("t2c+", "r1", "4123");
    reorder("t1c", "r2", "21");
    reorder("t1c+", "r3", "21");
    mult("r2", "r1", "r4", 1);
    mult("t2c", "r4", "r5", 3);
    mult("r5", "r3", "r6", 1);
    update("dm00_hh", 0.5, "r6");
    restore_stack_pos(pos);

    // DM-37
    reorder("t2c+", "r1", "3412");
    mult("r1", "t2c", "r2", 3);
    mult("r2", "t1c+", "r3", 1);
    mult("t1c", "r3", "r4", 1);
    update("dm00_hh", 0.5, "r4");
    restore_stack_pos(pos);
}


/**
 * Note that the factors for the DM-38 and DM-39 in the Bartlett's book are incorrect:
 *         Bartlett      EXP-T
 * DM-38    +1/2          +1/4
 * DM-39    +1/4          +1/2
 */
void sector_0h0p_construct_density_matrix_block_hh_T2x_T2x_T2_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-38
    reorder("t2c+", "r0", "3412");
    reorder("t2c", "r1", "1342");
    reorder("t2c+", "r2", "3124");
    mult("t2c", "r0", "r3", 3);
    mult("r1", "r3", "r4", 1);
    mult("r4", "r2", "r5", 3);
    update("dm00_hh", 0.25, "r5");
    restore_stack_pos(pos);

    // DM-39
    reorder("t2c+", "r1", "3412");
    reorder("t2c", "r2", "3412");
    mult("r2", "t2c+", "r3", 3);
    mult("t2c", "r3", "r4", 1);
    mult("r4", "r1", "r5", 3);
    update("dm00_hh", 0.5, "r5");
    restore_stack_pos(pos);

    // DM-40
    reorder("t2c", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("t2c+", "r3", "3412");
    mult("r1", "r2", "r4", 3);
    mult("r3", "t2c", "r5", 3);
    mult("r4", "r5", "r6", 1);
    update("dm00_hh", 0.25, "r6");
    restore_stack_pos(pos);

    // DM-41
    reorder("t2c+", "r1", "4231");
    reorder("t2c", "r2", "1324");
    reorder("t2c+", "r3", "3142");
    mult("r1", "r2", "r4", 2);
    mult("r4", "r3", "r5", 2);
    reorder("r5", "r6", "3142");
    mult("t2c", "r6", "r7", 3);
    update("dm00_hh", -1.0, "r7");
    restore_stack_pos(pos);
}


/**
 * particle-particle block of the 1-electron density matrix
 * (terms up to 4th order of PT)
 */
void sector_0h0p_construct_density_matrix_block_pp(int pt_order)
{
    tmplt("dm00_pp", "pp", "00", "12", NOT_PERM_UNIQUE);

    /*
     * 2nd PT order
     */
    if (pt_order >= PT_2) {
        sector_0h0p_construct_density_matrix_block_pp_T1x_T1();
        sector_0h0p_construct_density_matrix_block_pp_T2x_T2();
    }

    /*
     * 3rd PT order
     */
    if (pt_order >= PT_3) {
        sector_0h0p_construct_density_matrix_block_pp_T2x_T1_T1();
        sector_0h0p_construct_density_matrix_block_pp_T1x_T1x_T2();
    }

    /*
     * 4th PT order
     */
    if (pt_order >= PT_4) {
        sector_0h0p_construct_density_matrix_block_pp_T1x_T1x_T1_T1();
        sector_0h0p_construct_density_matrix_block_pp_T1x_T2x_T1_T2();
        sector_0h0p_construct_density_matrix_block_pp_T2x_T2x_T2_T2();

        if (triples_enabled()) {
            sector_0h0p_construct_density_matrix_block_pp_T3x_T3();
            sector_0h0p_construct_density_matrix_block_pp_T3x_T1_T2();
            sector_0h0p_construct_density_matrix_block_pp_T1x_T2x_T3();
        }
    }
}


void sector_0h0p_construct_density_matrix_block_pp_T1x_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-20
    reorder("t1c", "r1", "21");
    mult("t1c+", "r1", "r2", 1);
    update("dm00_pp", 1.0, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T2x_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-21
    reorder("t2c", "r1", "3412");
    mult("t2c+", "r1", "r2", 3);
    update("dm00_pp", 0.5, "r2");

    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T2x_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-22
    reorder("t2c+", "r1", "2431");
    reorder("t1c", "r2", "21");
    mult("r1", "t1c", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update("dm00_pp", 1.0, "r4");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T1x_T1x_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-23
    reorder("t2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("t1c+", "r2", "r3", 1);
    update("dm00_pp", 1.0, "r3");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T1x_T1x_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-24
    reorder("t1c", "r1", "21");
    reorder("t1c+", "r2", "21");
    reorder("t1c", "r3", "21");
    mult("t1c+", "r1", "r4", 1);
    mult("r4", "r2", "r5", 1);
    mult("r5", "r3", "r6", 1);
    update("dm00_pp", -1.0, "r6");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T3x_T3()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-25 - CCSDT
    reorder("t3c", "r1", "456123");
    mult("t3c+", "r1", "r2", 5);
    update("dm00_pp", +1.0 / 12, "r2");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T3x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-26
    reorder("t3c+", "r1", "124563");
    reorder("t2c", "r2", "3412");
    mult("r1", "t1c", "r3", 2);
    mult("r3", "r2", "r4", 3);
    update("dm00_pp", 0.5, "r4");
    restore_stack_pos(pos);

    // DM-27
    reorder("t3c+", "r1", "145623");
    reorder("t1c", "r2", "21");
    mult("r1", "t2c", "r3", 4);
    mult("r3", "r2", "r4", 1);
    update("dm00_pp", 0.25, "r4");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T1x_T2x_T3()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-28
    reorder("t3c", "r1", "451263");
    mult("r1", "t1c+", "r2", 2);
    mult("t2c+", "r2", "r3", 3);
    update("dm00_pp", 0.5, "r3");
    restore_stack_pos(pos);

    // DM-29
    reorder("t3c", "r1", "415623");
    mult("r1", "t2c+", "r2", 4);
    mult("t1c+", "r2", "r3", 1);
    update("dm00_pp", 0.25, "r3");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T1x_T2x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-30
    reorder("t2c+", "r1", "1324");
    reorder("t2c", "r2", "3142");
    reorder("t1c", "r3", "21");
    mult("r2", "t1c+", "r4", 2);
    mult("r1", "r4", "r5", 2);
    mult("r5", "r3", "r6", 1);
    update("dm00_pp", 1.0, "r6");
    restore_stack_pos(pos);

    // DM-31
    reorder("t2c", "r1", "3142");
    reorder("t2c+", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r1", "r3", "r4", 2);
    mult("t1c+", "r4", "r5", 1);
    update("dm00_pp", 1.0, "r5");
    restore_stack_pos(pos);

    // DM-32
    reorder("t2c+", "r1", "1342");
    reorder("t2c", "r2", "4231");
    mult("r1", "t1c", "r3", 2);
    mult("r2", "t1c+", "r4", 2);
    mult("r3", "r4", "r5", 1);
    update("dm00_pp", 1.0, "r5");
    restore_stack_pos(pos);

    // DM-33
    reorder("t1c+", "r1", "21");
    reorder("t2c", "r2", "3412");
    mult("t1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    mult("t2c+", "r4", "r5", 3);
    update("dm00_pp", -1.0, "r5");
    restore_stack_pos(pos);

    // DM-34
    reorder("t1c", "r1", "21");
    reorder("t2c", "r2", "3124");
    reorder("t2c+", "r3", "1342");
    mult("t1c+", "r1", "r4", 1);
    mult("r3", "r4", "r5", 1);
    mult("r5", "r2", "r6", 3);
    update("dm00_pp", -0.5, "r6");
    restore_stack_pos(pos);

    // DM-35
    reorder("t1c", "r1", "21");
    reorder("t2c+", "r2", "3412");
    mult("t2c", "r2", "r3", 3);
    mult("r1", "r3", "r4", 1);
    mult("t1c+", "r4", "r5", 1);
    update("dm00_pp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM-36
    reorder("t2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("r2", "t1c", "r3", 1);
    mult("r3", "t1c+", "r4", 1);
    update("dm00_pp", -0.5, "r4");
    restore_stack_pos(pos);

    // DM-37
    reorder("t2c", "r1", "3412");
    reorder("t1c", "r2", "21");
    reorder("t1c+", "r3", "21");
    mult("t2c+", "r1", "r4", 3);
    mult("r4", "r3", "r5", 1);
    mult("r5", "r2", "r6", 1);
    update("dm00_pp", -0.5, "r6");
    restore_stack_pos(pos);
}


void sector_0h0p_construct_density_matrix_block_pp_T2x_T2x_T2_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // DM-38
    reorder("t2c+", "r1", "3412");
    reorder("t2c", "r2", "3412");
    mult("t2c", "r1", "r3", 3);
    mult("r2", "r3", "r4", 1);
    mult("t2c+", "r4", "r5", 3);
    update("dm00_pp", -0.5, "r5");
    restore_stack_pos(pos);

    // DM-39
    reorder("t2c+", "r1", "1342");
    reorder("t2c", "r2", "3124");
    reorder("t2c", "r3", "3412");
    mult("r3", "t2c+", "r4", 3);
    mult("r2", "r4", "r5", 1);
    mult("r1", "r5", "r6", 3);
    update("dm00_pp", -0.25, "r6");
    restore_stack_pos(pos);

    // DM-40
    reorder("t2c", "r1", "3412");
    reorder("t2c", "r2", "4123");
    reorder("t2c+", "r3", "2341");
    mult("t2c+", "r1", "r4", 3);
    mult("r2", "r3", "r5", 3);
    mult("r4", "r5", "r6", 1);
    update("dm00_pp", -0.25, "r6");
    restore_stack_pos(pos);

    // DM-41
    reorder("t2c", "r1", "3412");
    reorder("t2c+", "r2", "1342");
    reorder("t2c", "r3", "2413");
    reorder("t2c+", "r4", "2431");
    mult("r2", "r3", "r5", 2);
    mult("r5", "r4", "r6", 2);
    reorder("r6", "r7", "1324");
    mult("r7", "r1", "r8", 3);
    update("dm00_pp", 1.0, "r8");
    restore_stack_pos(pos);
}


double complex sector_0h0p_norm(int approximation)
{
    if (approximation == CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return 1.0 + 0.0 * I;
    }

    return 1.0 + norm_0h0p_T1x_T1() + norm_0h0p_T2x_T2();
}


double complex sector_0h0p_prop(int approximation)
{
    double complex prop_value = 0.0;

    if (approximation >= CC_PROPERTIES_APPROX_LINEAR) {
        prop_value += sector_0h0p_prop_linear();
    }

    return prop_value;
}


double complex sector_0h0p_prop_linear()
{
    double complex d1 = property_0h0p_D1();
    double complex d2 = property_0h0p_D2();

    return d1 + d2;
}


double complex property_0h0p_D1()
{
    double complex d1 = scalar_product("C", "N", "op_hp", "t1c");

    return d1;
}


double complex property_0h0p_D2()
{
    double complex d2 = scalar_product("C", "N", "t1c", "op_hp");

    return d2;
}


double complex property_0h0p_D3()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("op_pp", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    double complex d3 = +1.0 * scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    return d3;
}


double complex property_0h0p_D4()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r1", "21");
    mult("op_hh", "r1", "r2", 1);
    double complex d4 = -1.0 * scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    return d4;
}


double complex property_0h0p_D5()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("op_pp", "r1", "21");
    mult("t2c", "r1", "r2", 1);
    double complex d5 = +0.5 * scalar_product("C", "N", "t2c", "r2");
    restore_stack_pos(pos);

    return d5;
}


double complex property_0h0p_D6()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c", "r1", "3412");
    mult("r1", "op_hh", "r2", 1);
    reorder("r2", "r3", "3412");
    double complex d6 = -0.5 * scalar_product("C", "N", "t2c", "r3");
    restore_stack_pos(pos);

    return d6;
}


double complex property_0h0p_D7()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c", "r1", "1324");
    reorder("op_ph", "r2", "21");
    mult("r1", "r2", "r3", 2);
    double complex d7 = scalar_product("C", "N", "t1c", "r3");
    restore_stack_pos(pos);

    return d7;
}


double complex property_0h0p_D8()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c+", "r1", "4231");
    mult("r1", "t1c", "r2", 2);
    double complex d8 = scalar_product("N", "N", "op_hp", "r2");
    restore_stack_pos(pos);

    return d8;
}


double complex property_0h0p_D9()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r1", "21");
    reorder("t2c+", "r2", "3124");
    mult("r1", "op_hh", "r3", 1);
    mult("r2", "r3", "r4", 2);
    double complex d9 = -1.0 * scalar_product("N", "N", "t1c", "r4");
    restore_stack_pos(pos);

    return d9;
}


double complex property_0h0p_D10()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c+", "r1", "3124");
    reorder("op_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    double complex d10 = 1.0 * scalar_product("N", "N", "t1c", "r4");
    restore_stack_pos(pos);

    return d10;
}


double complex property_0h0p_D11()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("op_hh", "r1", "21");
    reorder("t2c", "r2", "1324");
    mult("r1", "t1c+", "r3", 1);
    mult("r2", "r3", "r4", 2);
    double complex d11 = -1.0 * scalar_product("C", "N", "t1c", "r4");
    restore_stack_pos(pos);

    return d11;
}


double complex property_0h0p_D12()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c+", "r1", "21");
    reorder("t2c", "r2", "1324");
    mult("r1", "op_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    double complex d12 = +1.0 * scalar_product("C", "N", "t1c", "r4");
    restore_stack_pos(pos);

    return d12;
}


double complex property_0h0p_D13()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r0", "21");
    reorder("t1c+", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    mult("r2", "r0", "r3", 1);
    double complex d13 = -1.0 * scalar_product("C", "N", "op_hp", "r3");
    restore_stack_pos(pos);

    return d13;
}


double complex property_0h0p_D14()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c+", "r1", "21");
    reorder("op_hp", "r2", "21");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    double complex d14 = -1.0 * scalar_product("C", "N", "t1c", "r4");
    restore_stack_pos(pos);

    return d14;
}


double complex property_0h0p_D15()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r0", "21");
    reorder("t2c+", "r1", "2413");
    mult("r1", "r0", "r2", 2);
    reorder("r2", "r3", "21");
    reorder("t2c", "r4", "2413");
    mult("r4", "r3", "r5", 2);
    double complex d15 = 1.0 * scalar_product("C", "N", "op_hp", "r5");
    restore_stack_pos(pos);

    return d15;
}


double complex property_0h0p_D16()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("op_hp", "r0", "21");
    reorder("t2c+", "r1", "2413");
    mult("r1", "r0", "r2", 2);
    reorder("r2", "r3", "21");
    reorder("t2c", "r4", "2413");
    mult("r4", "r3", "r5", 2);
    double complex d16 = 1.0 * scalar_product("C", "N", "t1c", "r5");
    restore_stack_pos(pos);

    return d16;
}


double complex property_0h0p_D17()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r1", "21");
    reorder("t2c+", "r2", "3412");
    mult("t2c", "r2", "r3", 3);
    mult("r3", "r1", "r4", 1);
    double complex d17 = -0.5 * scalar_product("C", "N", "op_hp", "r4");
    restore_stack_pos(pos);

    return d17;
}


double complex property_0h0p_D18()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("t1c", "r2", "r3", 1);
    double complex d18 = -0.5 * scalar_product("C", "N", "op_hp", "r3");
    restore_stack_pos(pos);

    return d18;
}


double complex property_0h0p_D19()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("op_hp", "r1", "21");
    reorder("t2c+", "r2", "3412");
    mult("t2c", "r2", "r3", 3);
    mult("r3", "r1", "r4", 1);
    double complex d19 = -0.5 * scalar_product("C", "N", "t1c", "r4");
    restore_stack_pos(pos);

    return d19;
}


double complex property_0h0p_D20()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("op_hp", "r2", "r3", 1);
    double complex d20 = -0.5 * scalar_product("C", "N", "t1c", "r3");
    restore_stack_pos(pos);

    return d20;
}


void guess_operator_symmetry(int nspinors, double complex *prop_matrix)
{
    int nsym = get_num_irreps();

    int *nonzero_blocks = cc_malloc(sizeof(int) * nsym * nsym);
    memset(nonzero_blocks, 0, sizeof(int) * nsym * nsym);

/*
 * detect pairs of irreps for which matrix elements are non-zero
 */
    for (int i = 0; i < nspinors; i++) {
        for (int j = 0; j < nspinors; j++) {
            double complex prop_ij = prop_matrix[i * nspinors + j];
            if (cabs(prop_ij) > 1e-10) {
                int irrep_i = get_spinor_irrep(i);
                int irrep_j = get_spinor_irrep(j);
                nonzero_blocks[irrep_i * nsym + irrep_j] = 1;
            }
        }
    }

/*
 * try to guess irrep of the operator
 */
    printf(" list of non-zero symmetry blocks:\n");
    for (int irrep_i = 0; irrep_i < nsym; irrep_i++) {
        for (int irrep_j = 0; irrep_j < nsym; irrep_j++) {
            if (nonzero_blocks[irrep_i * nsym + irrep_j]) {
                printf(" %s - %s\n", get_irrep_name(irrep_i), get_irrep_name(irrep_j));

                /*for (int op_rep = 0; op_rep < nsym; op_rep++) {
                    if (irrep_i == mulrep2_abelian(op_rep, irrep_j)) {
                        printf("%s ", get_irrep_name(op_rep));
                    }
                }*/
            }
        }
    }
    printf("\n");

    cc_free(nonzero_blocks);
}


double complex calculate_vacuum_expectation_value(int nspinors, double complex *prop_matrix)
{
    double complex expect_value = 0.0 + 0.0 * I;

    for (int ispinor = 0; ispinor < nspinors; ispinor++) {
        if (is_hole(ispinor)) {
            expect_value += prop_matrix[nspinors * ispinor + ispinor];
        }
    }

    return expect_value;
}


double complex norm_0h0p_T1x_T1()
{
    return scalar_product("C", "N", "t1c", "t1c");
}


double complex norm_0h0p_T2x_T2()
{
    return 0.25 * scalar_product("C", "N", "t2c", "t2c");
}


double complex norm_0h0p_T1x_T1x_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("t2c", "r1", "1324");
    reorder("t1c+", "r2", "21");
    mult("r1", "r2", "r3", 2);
    double complex s = 0.5 * scalar_product("C", "N", "t1c", "r3");

    restore_stack_pos(pos);

    return s;
}


double complex norm_0h0p_T1x_T1x_T1_T1()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D1
    double complex one_loop = norm_0h0p_T1x_T1();
    double complex d1 = 0.25 * one_loop * one_loop;

    // D2
    reorder("t1c+", "r0", "21");
    reorder("t1c", "r1", "21");
    mult("t1c", "r0", "r2", 1);
    mult("r2", "r1", "r3", 1);
    double complex d2 = -1.0 * scalar_product("C", "N", "t1c", "r3");
    restore_stack_pos(pos);

    // D3
    /*reorder("t1c", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    double complex d3 = - 0.25 * scalar_product("C", "N", "t1c", "r3");
    restore_stack_pos(pos);*/

    printf("\nT1x T1x - T1 T1\n");
    printf("D1 = %12.8f %12.8f\n", creal(d1), cimag(d1));
    printf("D2 = %12.8f %12.8f\n", creal(d2), cimag(d2));
    //printf("D3 = %12.8f %12.8f\n", creal(d3), cimag(d3));

    return d1 + d2;
}


double complex norm_0h0p_T1x_T2x_T1_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D1
    double complex one_loop = norm_0h0p_T1x_T1();
    double complex two_loop = norm_0h0p_T2x_T2();
    double complex d1 = two_loop * one_loop;

    // D2
    reorder("t2c+", "r1", "3412");
    reorder("t1c", "r3", "21");
    mult("t2c", "r1", "r2", 3);
    mult("r2", "r3", "r4", 1);
    double complex d2 = -0.5 * scalar_product("C", "N", "t1c", "r4");
    restore_stack_pos(pos);

    // D3
    reorder("t2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("t1c", "r2", "r3", 1);
    double complex d3 = -0.5 * scalar_product("C", "N", "t1c", "r3");
    restore_stack_pos(pos);

    // D4
    reorder("t2c", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "t1c", "r3", 2);
    mult("r1", "r3", "r4", 2);
    double complex d4 = +1.0 * scalar_product("C", "N", "t1c", "r4");
    restore_stack_pos(pos);

    printf("\nT1x T2x - T1 T2\n");
    printf("D1 = %12.8f %12.8f\n", creal(d1), cimag(d1));
    printf("D2 = %12.8f %12.8f\n", creal(d2), cimag(d2));
    printf("D3 = %12.8f %12.8f\n", creal(d3), cimag(d3));
    printf("D4 = %12.8f %12.8f\n", creal(d4), cimag(d4));

    return d1 + d2 + d3 + d4;
}


double complex norm_0h0p_T2x_T2x_T2_T2()
{
    dg_stack_pos_t pos = get_stack_pos();

    // D1
    double complex two_loop = norm_0h0p_T2x_T2();
    double complex d1 = two_loop * two_loop;

    // D2
    reorder("t2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("t2c", "r2", "r3", 1);
    double complex d2 = -0.25 * scalar_product("C", "N", "t2c", "r3");
    restore_stack_pos(pos);

    // D3
    reorder("t2c", "r1", "3412");
    reorder("t2c+", "r2", "3412");
    mult("t2c", "r2", "r3", 3);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    double complex d3 = -0.25 * scalar_product("C", "N", "t2c", "r5");
    restore_stack_pos(pos);

    // D4
    reorder("t2c", "r3", "3412");
    reorder("t2c+", "r4", "3412");
    mult("t2c", "r4", "r1", 2);
    mult("r1", "r3", "r2", 2);
    double complex d4 = (0.25 * 0.25) * scalar_product("C", "N", "t2c", "r2");
    restore_stack_pos(pos);

    // D5
    reorder("t2c", "r1", "1324");
    reorder("t2c+", "r2", "4231");
    mult("r1", "r2", "r3", 2);
    mult("r3", "r1", "r4", 2);
    reorder("r4", "r5", "1324");
    double complex d5 = 0.25 * scalar_product("C", "N", "t2c", "r5");
    restore_stack_pos(pos);

    // D6
    /*reorder("t2c+", "r0", "3412");
    mult("t2c", "r0", "r1", 3);
    reorder("t2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    double complex d6 = - 0.25 * scalar_product("C", "N", "t2c", "r3");
    restore_stack_pos(pos);

    // D7
    reorder("t2c", "r0", "3412");
    mult("r0", "pphh", "r1", 3);
    mult("t2c", "r1", "r2", 1);
    double complex d7 = - 0.25 * scalar_product("C", "N", "t2c", "r2");
    restore_stack_pos(pos);*/

    printf("\nT2x T2x - T2 T2\n");
    printf("D1 = %12.8f %12.8f\n", creal(d1), cimag(d1));
    printf("D2 = %12.8f %12.8f\n", creal(d2), cimag(d2));
    printf("D3 = %12.8f %12.8f\n", creal(d3), cimag(d3));
    printf("D4 = %12.8f %12.8f\n", creal(d4), cimag(d4));
    printf("D5 = %12.8f %12.8f\n", creal(d5), cimag(d5));
    //printf("D6 = %12.8f %12.8f\n", creal(d6), cimag(d6));
    //printf("D7 = %12.8f %12.8f\n", creal(d7), cimag(d7));

    return d1 + d2 + d3 + d4 + d5;
}


void sector00_calc_wf_norm()
{
    diagram_conjugate("t1c", "t1c+");
    diagram_conjugate("t2c", "t2c+");

    double complex norm2 = 1.0;

    // < T1 | T1 >
    double complex norm2_T1x_T1 = norm_0h0p_T1x_T1();
    norm2 += norm2_T1x_T1;

    // < T2 | T2 >
    double complex norm2_T2x_T2 = norm_0h0p_T2x_T2();
    norm2 += norm2_T2x_T2;

    // < T1^2 | T2 >
    double complex norm2_T1x_T1x_T2 = norm_0h0p_T1x_T1x_T2();
    norm2 += norm2_T1x_T1x_T2;

    // < T2 | T1^2 >
    double complex norm2_T2x_T1_T1 = norm2_T1x_T1x_T2;
    norm2 += norm2_T2x_T1_T1;

    // < T1 T1 | T1 T1 >
    double complex norm2_T1x_T1x_T1_T1 = norm_0h0p_T1x_T1x_T1_T1();
    norm2 += norm2_T1x_T1x_T1_T1;

    // < T1 T2 | T1 T2 >
    double complex norm2_T1x_T2x_T1_T2 = norm_0h0p_T1x_T2x_T1_T2();
    norm2 += norm2_T1x_T2x_T1_T2;

    // < T2 T2 | T2 T2 >
    double complex norm2_T2x_T2x_T2_T2 = norm_0h0p_T2x_T2x_T2_T2();
    norm2 += norm2_T2x_T2x_T2_T2;

    printf("\n");
    printf("  %20s%12.8f\n", "1 : 1", 1.0);
    printf("  %20s%12.8f%12.8f\n", "T1 : T1", creal(norm2_T1x_T1), cimag(norm2_T1x_T1));
    printf("  %20s%12.8f%12.8f\n", "T2 : T2", creal(norm2_T2x_T2), cimag(norm2_T2x_T2));
    printf("  %20s%12.8f%12.8f\n", "T1 T1 : T2", creal(norm2_T1x_T1x_T2), cimag(norm2_T1x_T1x_T2));
    printf("  %20s%12.8f%12.8f\n", "T2 : T1 T1", creal(norm2_T2x_T1_T1), cimag(norm2_T2x_T1_T1));
    printf("  %20s%12.8f%12.8f\n", "T1 T1 : T1 T1", creal(norm2_T1x_T1x_T1_T1), cimag(norm2_T1x_T1x_T1_T1));
    printf("  %20s%12.8f%12.8f\n", "T1 T2 : T1 T2", creal(norm2_T1x_T2x_T1_T2), cimag(norm2_T1x_T2x_T1_T2));
    printf("  %20s%12.8f%12.8f\n", "T2 T2 : T2 T2", creal(norm2_T2x_T2x_T2_T2), cimag(norm2_T2x_T2x_T2_T2));
    printf("  ------------------------------------------------\n");
    printf("  %20s%12.8f%12.8f\n", "< psi | psi >", creal(norm2), cimag(norm2));
    printf("  %20s%12.8f\n", "norm ||psi||", sqrt(creal(norm2)));
    printf("\n");
}


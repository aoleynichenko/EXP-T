/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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

/*******************************************************************************
 * t3corr.c
 * ==========
 *
 * Perturbative corrections to account for triple excitations: +T(CCSD) and (T)
 * For a detailed description of the theoretical part, see:
 *   (1) I. Shavitt, R. J. Bartlett, "Many-Body Methods in Chemistry and
 *       Physics: MBPT and Coupled-Cluster Theory", 1st edition
 *   (2) M. J. O. Deegan, P. J. Knowles, Chem. Phys. Lett. V. 227, P. 321 (1994)
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "engine.h"
#include "platform.h"
#include "datamodel.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"

void calc_T3(int pt_order);

void t3_0h0p_contrib_to_singles(int pt_order);

void t3_0h0p_contrib_to_doubles(int pt_order);


/*******************************************************************************
 * t3corr
 *
 * Calculates the (T) and [T] perturbative corrections for the ground-state CCSD
 * (FS sector 0h0p) and prints them out.
 * NOTE: this is an experimental code; permutational symmetry of T3 amplitudes
 * is not employed and T3 amplitudes by default are stored on disk;
 * so the performance of this code is still very low
 ******************************************************************************/
void t3corr()
{
    printf("\n");
    printf(" Perturbative corrections to account for triple excitations\n");
    printf(" This is an experimental code (is not optimized at all)\n");

    // extract one-electron energies
    double *eps = (double *) cc_malloc(NSPINORS * sizeof(double));
    for (size_t i = 0; i < NSPINORS; i++) {
        eps[i] = spinor_info[i].eps;
    }

    // construct approximate triples amplitudes
    dg_stack_pos_t pos;

    timer_new_entry("00-(T)", "0h0p -- CCSD(T) perturbative triples");
    timer_start("00-(T)");

    tmplt("t3c", "hhhppp", "000000", "123456", NOT_PERM_UNIQUE);
    pos = get_stack_pos();

    // T1a
    printf(" construction of T3 amplitudes:\n");
    double t0 = abs_time();

    reorder("phpp", "r1", "2341");
    double t1 = abs_time();
    printf("  %-25s%8.3f sec\n", "reorder <ph||pp>", t1 - t0);
    mult("t2c", "r1", "r2", 1);
    double t2 = abs_time();
    printf("  %-25s%8.3f sec\n", "contraction T2 - phpp", t2 - t1);
    reorder("r2", "r3", "124356");
    double t3 = abs_time();
    printf("  %-25s%8.3f sec\n", "reorder intermediate", t3 - t2);
    perm("r3", "(3/12|4/56)");
    double t4 = abs_time();
    printf("  %-25s%8.3f sec\n", "perm P(k/ij|a/bc)", t4 - t3);
    update("t3c", 1.0, "r3");
    double t5 = abs_time();
    printf("  %-25s%8.3f sec\n", "update T3 amplitudes", t5 - t4);
    restore_stack_pos(pos);

    // T1b
    reorder("t2c", "r1", "3412");
    double t6 = abs_time();
    printf("  %-25s%8.3f sec\n", "reorder T2 amplitudes", t6 - t5);
    reorder("hhhp", "r2", "1243");
    double t7 = abs_time();
    printf("  %-25s%8.3f sec\n", "reorder <hh||hp>", t7 - t6);
    mult("r1", "r2", "r3", 1);
    double t8 = abs_time();
    printf("  %-25s%8.3f sec\n", "contraction T2 - hhhp", t8 - t7);
    reorder("r3", "r4", "345126");
    double t9 = abs_time();
    printf("  %-25s%8.3f sec\n", "reorder intermediate", t9 - t8);
    perm("r4", "(1/23|6/45)");
    double t10 = abs_time();
    printf("  %-25s%8.3f sec\n", "perm P(i/jk|c/ab)", t10 - t9);
    update("t3c", -1.0, "r4");
    double t11 = abs_time();
    printf("  %-25s%8.3f sec\n", "update T3 amplitudes", t11 - t10);
    restore_stack_pos(pos);

    diveps("t3c");
    double t12 = abs_time();
    printf("  %-25s%8.3f sec\n", "divide T3 by denom-s", t12 - t11);

    double de_t = 0.0;
    double de_st = 0.0;
    double de_dt = 0.0;

    // calculate 4-order [T] correction (blockwise)
    // NOTE: permutational symmetry is not employed
    diagram_t *dg_t3 = diagram_stack_find("t3c");
    for (size_t isb = 0; isb < dg_t3->n_blocks; isb++) {
        block_t *block = dg_t3->blocks[isb];
        symblock_load(block);
        double complex *buf = block->buf;
        double denom = 0;
        int dims_0 = block->indices[0][0];
        int dims_1 = block->indices[1][0];
        int dims_2 = block->indices[2][0];
        int dims_3 = block->indices[3][0];
        int dims_4 = block->indices[4][0];
        int dims_5 = block->indices[5][0];
        int *block_indices_0 = block->indices[0] + 1; // +1 to skip first element (length)
        int *block_indices_1 = block->indices[1] + 1;
        int *block_indices_2 = block->indices[2] + 1;
        int *block_indices_3 = block->indices[3] + 1;
        int *block_indices_4 = block->indices[4] + 1;
        int *block_indices_5 = block->indices[5] + 1;

        // energy denominators are constructed stepwise
        // to reduce the number of memory access operations
        // denom = ei + ej + ek - ea - eb - ec
        size_t index = 0;
        for (int i0 = 0; i0 < dims_0; i0++) {
            int idx_0 = block_indices_0[i0]; // relative index to spinor indices
            double denom_0 = eps[idx_0];
            for (int i1 = 0; i1 < dims_1; i1++) {
                int idx_1 = block_indices_1[i1];
                double denom_1 = denom_0 + eps[idx_1];
                for (int i2 = 0; i2 < dims_2; i2++) {
                    int idx_2 = block_indices_2[i2];
                    double denom_2 = denom_1 + eps[idx_2];
                    for (int i3 = 0; i3 < dims_3; i3++) {
                        int idx_3 = block_indices_3[i3];
                        double denom_3 = denom_2 - eps[idx_3];
                        for (int i4 = 0; i4 < dims_4; i4++) {
                            int idx_4 = block_indices_4[i4];
                            double denom_4 = denom_3 - eps[idx_4];
                            for (int i5 = 0; i5 < dims_5; i5++) {
                                int idx_5 = block_indices_5[i5];
                                double denom = denom_4 - eps[idx_5];

                                // [T] += |t_ijkabc|^2 * e_ijkabc
                                double t3abs;
                                if (carith) {
                                    t3abs = cabs(buf[index]);
                                }
                                else {
                                    t3abs = fabs(((double *) buf)[index]);
                                }
                                de_t += t3abs * t3abs * denom;

                                index++;
                            }
                        }
                    }
                }
            }
        }

        symblock_unload(block);
    }
    de_t /= 36.0;
    double t13 = abs_time();
    printf("  %-25s%8.3f sec\n", "dE_t contribution", t13 - t12);

    // E_st correction
    reorder("t3c", "r1", "145623");
    mult("r1", "pphh", "r2", 4);
    de_st = 0.25 * creal(scalar_product("C", "N", "t1c", "r2"));
    double t14 = abs_time();
    printf("  %-25s%8.3f sec\n", "dE_st contribution", t14 - t13);
    restore_stack_pos(pos);

    // E_dt correction
    reorder("t3c", "r1", "124563");
    mult("r1", "ph", "r2", 2);
    de_dt = 0.25 * creal(scalar_product("C", "N", "t2c", "r2"));
    double t15 = abs_time();
    printf("  %-25s%8.3f sec\n", "dE_dt contribution", t15 - t14);
    restore_stack_pos(pos);

    // prepare diagram T3 for subsequent use in other sectors: rename it
    //renamedg("t3c", "t3");

    printf("\n");
    printf("          SCF reference energy = %20.12f\n", cc_opts->escf);
    printf("       CCSD correlation energy = %20.12f\n", cc_opts->eref - cc_opts->escf);
    printf("             Total CCSD energy = %20.12f\n", cc_opts->eref);
    printf("      4th order triples corr-n = %20.12f\n", de_t);
    printf("  5th order triples (T) corr-n = %20.12f\n", de_st + de_dt);
    printf("     Total CCSD+T(CCSD) energy = %20.12f\n", cc_opts->eref + de_t);
    printf("          Total CCSD(T) energy = %20.12f\n", cc_opts->eref + de_t + de_st + de_dt);
    printf("\n");

    timer_stop("00-(T)");

    cc_free(eps);
}


/**
 * Constructs triples amplitudes which are correct up to 2nd PT order
 * (energy is correct up to 4th PT order)
 */
void sector_0h0p_ccsd_t3()
{
    timer_new_entry("00-T(3)", "0h0p -- CCSD+T(3) perturbative triples correction");
    timer_start("00-T(3)");

    printf(" 0h0p: 2nd order triples\n");

    tmplt("t3nw", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);
    calc_T3(PT_2);
    rename_diagram("t3nw", "t3c");
    diveps("t3c");

    timer_stop("00-T(3)");
}


/**
 * Constructs triples amplitudes which are correct up to 3rd PT order
 * (energy is correct up to 5th PT order)
 */
void sector_0h0p_ccsd_t4()
{
    dg_stack_pos_t pos;
    timer_new_entry("00-T(4)", "0h0p -- CCSD+T(4) perturbative triples correction");
    timer_start("00-T(4)");

    copy("t1c", "t1c_ccsd");
    copy("t2c", "t2c_ccsd");

    // calculate PT2 perturbative triples
    sector_0h0p_ccsd_t3();

    // calculate PT3 corrections
    tmplt("t1nw", "hp", "00", "12", NOT_PERM_UNIQUE);
    tmplt("t2nw", "hhpp", "0000", "1234", IS_PERM_UNIQUE);
    printf(" 0h0p: 3rd order correction to singles\n");
    t3_0h0p_contrib_to_singles(PT_3);
    printf(" 0h0p: 3rd order correction to doubles\n");
    t3_0h0p_contrib_to_doubles(PT_3);
    rename_diagram("t1nw", "t1_pt3_corr");
    rename_diagram("t2nw", "t2_pt3_corr");
    diveps("t1_pt3_corr");
    diveps("t2_pt3_corr");

    // update T1 and T2 with PT3 corrections
    update("t1c", 1.0, "t1_pt3_corr");
    update("t2c", 1.0, "t2_pt3_corr");
    // at this point 't1c' and 't2c' contain amplitudes which are correct up to PT3

    // calculate PT3 perturbative triples
    printf(" 0h0p: 3rd order triples\n");
    tmplt("t3nw", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);
    calc_T3(PT_3);
    rename_diagram("t3nw", "t3c");
    diveps("t3c");
    // at this point 't3c' contain amplitudes which are correct up to PT3

    timer_stop("00-T(4)");
}

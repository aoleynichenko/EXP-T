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
 * sector01.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with one particle (Fock space sector 0h1p)
 *
 * Cluster Operator: T = T1 + T2 (+ T3) + S1 + S2 (+ S3)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip hole creation lines to turn them to valence particle
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef VERSION_DEVEL
#error This file is a part of the development version only!
#endif

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ccutils.h"
#include "diis.h"
#include "engine.h"
#include "datamodel.h"
#include "options.h"
#include "methods.h"
#include "sort.h"
#include "heff.h"
#include "spinors.h"
#include "utils.h"

void calc_T3(int pt_order);

void const_terms_0h1p();

void calc_S1();

void calc_S3(int pt_order);

void folded_0h1p();

void t3_0h0p_contrib_to_singles(int pt_order);

void t3_0h0p_contrib_to_doubles(int pt_order);

void t3_0h1p_contrib_to_doubles(int pt_order);

void diag_T3_0h1p_i_jk_c_ab_inv1(int pt_order);
void diag_T3_0h1p_i_jk_c_ab_inv3(int pt_order);

void diag_T3_0h1p_ijk_c_ab_inv1(int pt_order);
void diag_T3_0h1p_ijk_c_ab_inv2(int pt_order);
void diag_T3_0h1p_ijk_c_ab_inv3(int pt_order);

void diag_T3_0h1p_k_ij_a_bc_inv1(int pt_order);
void diag_T3_0h1p_k_ij_a_bc_inv3(int pt_order);

void diag_T3_0h1p_k_ij_abc_inv1(int pt_order);
void diag_T3_0h1p_k_ij_abc_inv3(int pt_order);

void diag_T3_0h1p_c_ab(int pt_order);

void diag_T3_0h1p_a_bc(int pt_order);

void diag_T3_0h1p_abc(int pt_order);

void diag_T3_0h1p_k_ij_c_ab_inv1(int pt_order);
void diag_T3_0h1p_k_ij_c_ab_inv3(int pt_order);

void diag_T3_0h1p_ijk_inv1(int pt_order);
void diag_T3_0h1p_ijk_inv2(int pt_order);
void diag_T3_0h1p_ijk_inv3(int pt_order);

void diag_T3_0h1p_i_jk_inv1(int pt_order);
void diag_T3_0h1p_i_jk_inv3(int pt_order);

void diag_T3_0h1p_k_ij_inv1(int pt_order);
void diag_T3_0h1p_k_ij_inv3(int pt_order);


void t3_0h1p_const_contrib_to_triples(int pt_order)
{
    dg_stack_pos_t pos = get_stack_pos();

    // T1a-2
    reorder("vppp", "r1", "3412");
    reorder("t2c", "r2", "1243");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "612453");
    diagram_stack_erase("r3");
    perm("r4", "(6/45)");
    update("s3_0", 1.0, "r4");
    restore_stack_pos(pos);

    // T1b-2
    reorder("t2c", "r1", "2341");
    mult("vhph", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(23|4/56)");
    update("s3_0", -1.0, "r3");
    restore_stack_pos(pos);
}


void t3_0h1p_contrib_to_singles(int pt_order)
{
    if (pt_order <= PT_2) {
        return;
    }

    timer_new_entry("01-T1-S7", "0h1p -- Triples contribution to singles");
    timer_start("01-T1-S7");

    // S7
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s3c", "r1", "145623");
    mult("r1", "pphh", "r2", 4);
    update("s1nw", 0.25, "r2");
    restore_stack_pos(pos);

    timer_stop("01-T1-S7");
}


void t3_0h1p_contrib_to_doubles(int pt_order)
{
    if (pt_order <= PT_2) {
        return;
    }

    timer_new_entry("01-T2-Triples", "0h1p -- Triples contribution to doubles");
    timer_start("01-T2-Triples");

    dg_stack_pos_t pos;
    pos = get_stack_pos();

    // D10a, D11a
    // T3 reordering: 124536
    reorder("s3c", "r1", "124536");
    tmplt("i1", "hp", "00", "21", NOT_PERM_UNIQUE);
    {
        dg_stack_pos_t pos2 = get_stack_pos();

        // D10a
        reorder("ph", "phr", "21");
        update("i1", 1.0, "phr");
        restore_stack_pos(pos2);

        if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_1A) ||
            (pt_order != PT_INF && pt_order >= PT_4)) {
            // D11a
            reorder("pphh", "_r1", "3142");
            mult("_r1", "t1c", "_r2", 2);
            update("i1", 1.0, "_r2");
            restore_stack_pos(pos2);
        }
    }
    mult("r1", "i1", "r2", 2);
    update("s2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10b, D11b
    // T3 reordering: 124356
    // permutation: P(34)
    reorder("s3c", "r1", "124356");
    tmplt("i1", "phpp", "0000", "3412", NOT_PERM_UNIQUE);
    {
        dg_stack_pos_t pos2 = get_stack_pos();

        // D10b
        reorder("ppph", "_r1", "3412");
        update("i1", 0.5, "_r1");
        restore_stack_pos(pos2);

        if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_1A) ||
            (pt_order != PT_INF && pt_order >= PT_4)) {
            // D11b
            reorder("pphh", "_r1", "4123");
            reorder("t1c", "_r2", "21");
            mult("_r2", "_r1", "_r3", 1);
            update("i1", -0.5, "_r3");
            restore_stack_pos(pos2);
        }
    }
    mult("r1", "i1", "r2", 3);
    perm("r2", "(34)");
    update("s2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10c, D11c
    // T3 reordering: 145236
    // permutation: (12)
    // valence line in the T3 diagram
    reorder("s3c", "r1", "145236");
    tmplt("i1", "hhhp", "0000", "1342", NOT_PERM_UNIQUE);
    {
        dg_stack_pos_t pos2 = get_stack_pos();

        // D10c
        reorder("hphh", "_r1", "1342");
        update("i1", -0.5, "_r1");
        restore_stack_pos(pos2);

        if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_1A) ||
            (pt_order != PT_INF && pt_order >= PT_4)) {
            // D11c
            reorder("pphh", "_r1", "3421");
            mult("t1c", "_r1", "_r2", 1);
            update("i1", -0.5, "_r2");
            restore_stack_pos(pos2);
        }
    }
    mult("r1", "i1", "r2", 3);
    reorder("r2", "r3", "1423");
    update("s2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // D10c, D11с
    // T3 reordering: 145236
    // valence line in the intermediate diagram
    reorder("t3c", "r1", "145236");
    tmplt("i1", "phhp", "1000", "1342", NOT_PERM_UNIQUE);
    {
        dg_stack_pos_t pos2 = get_stack_pos();

        // D10c
        reorder("vphh", "_r1", "1342");
        update("i1", -0.5, "_r1");
        restore_stack_pos(pos2);

        if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_1A) ||
            (pt_order != PT_INF && pt_order >= PT_4)) {
            // D11c
            reorder("pphh", "_r1", "3421");
            mult("s1c", "_r1", "_r2", 1);
            update("i1", -0.5, "_r2");
            restore_stack_pos(pos2);
        }
    }
    mult("r1", "i1", "r2", 3);
    reorder("r2", "r3", "4132");
    update("s2nw", 1.0, "r3");
    restore_stack_pos(pos);

    timer_stop("01-T2-Triples");
}


void calc_S3(int pt_order)
{
    timer_new_entry("01-S3", "0h1p -- Triples equations (T3)");
    timer_start("01-S3");

    copy("s3_0", "s3nw");

    diag_T3_0h1p_i_jk_c_ab_inv1(pt_order);
    diag_T3_0h1p_i_jk_c_ab_inv3(pt_order);

    diag_T3_0h1p_ijk_c_ab_inv1(pt_order);
    diag_T3_0h1p_ijk_c_ab_inv2(pt_order);
    diag_T3_0h1p_ijk_c_ab_inv3(pt_order);

    diag_T3_0h1p_k_ij_a_bc_inv1(pt_order);
    diag_T3_0h1p_k_ij_a_bc_inv3(pt_order);

    diag_T3_0h1p_k_ij_abc_inv1(pt_order);
    diag_T3_0h1p_k_ij_abc_inv3(pt_order);

    diag_T3_0h1p_k_ij_c_ab_inv1(pt_order);
    diag_T3_0h1p_k_ij_c_ab_inv3(pt_order);

    diag_T3_0h1p_c_ab(pt_order);

    diag_T3_0h1p_a_bc(pt_order);

    diag_T3_0h1p_abc(pt_order);

    diag_T3_0h1p_ijk_inv1(pt_order);
    diag_T3_0h1p_ijk_inv2(pt_order);
    diag_T3_0h1p_ijk_inv3(pt_order);

    diag_T3_0h1p_i_jk_inv1(pt_order);
    diag_T3_0h1p_i_jk_inv3(pt_order);

    diag_T3_0h1p_k_ij_inv1(pt_order);
    diag_T3_0h1p_k_ij_inv3(pt_order);

    timer_stop("01-S3");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_0h1p_i_jk_c_ab_inv1(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhph", "0000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hhhp", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3a
        reorder("ph", "phr_", "21");
        reorder("t2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphp", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("t1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphp", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("t2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("t2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("s3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("t2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        diagram_stack_erase("t5d_r1");
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        diagram_stack_erase("t5d_r4");
        reorder("t5d_r5", "t5d_r6", "156234");
        diagram_stack_erase("t5d_r5");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(6/45)");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_i_jk_c_ab_inv3(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3a
        reorder("ph", "phr_", "21");
        reorder("s2c", "r1", "2134");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("s2c", "s2c_21", "2143");
        reorder("pphp", "r2", "4312");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("t1c", "r2", "21");
        mult("hvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        reorder("s2c", "s2c_21", "2143");
        mult("s2c_21", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("t3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        diagram_stack_erase("t5d_r1");
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        diagram_stack_erase("t5d_r4");
        reorder("t5d_r5", "t5d_r6", "156234");
        diagram_stack_erase("t5d_r5");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(12|6/45)");
    reorder("r4", "r5", "321654");
    diagram_stack_erase("r4");
    update("s3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_0h1p_ijk_c_ab_inv1(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhph", "0000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        mult("t1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("t1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("t2c", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    perm("r4", "(23|6/45)");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_ijk_c_ab_inv2(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "phph", "1000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        mult("s1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("t1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("t2c", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    perm("r4", "(13|6/45)");
    reorder("r4", "r5", "213546");
    diagram_stack_erase("r4");
    update("s3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_ijk_c_ab_inv3(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("pvhp", "r1", "2341");
        mult("t1c", "r1", "r2", 1);
        reorder("r2", "r3", "1243");
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("t1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    perm("r4", "(12|6/45)");
    reorder("r4", "r5", "321654");
    update("s3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_0h1p_k_ij_a_bc_inv1(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hppp", "0000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phpp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3e
        reorder("t2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppppr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("t1c", "r1", "21");
        mult("r1", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("s3c", "t5e_r1", "134562");
        reorder("pphh", "t5e_r2", "4123");
        reorder("t2c", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r4");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(23|4/56)");
    update("s3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_k_ij_a_bc_inv3(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3e
        reorder("t2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppppr", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("t1c", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d --- ?
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b --- ?
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    diagram_stack_erase("i1");
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("t3c", "t5e_r1", "134562");
        reorder("pphh", "t5e_r2", "4123");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r4");
        reorder("t5e_r5", "t5e_r6", "124356");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(4/56)");
    reorder("r3", "r4", "321654");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_0h1p_k_ij_abc_inv1(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hppp", "0000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("t2c", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("phhp", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("t2c", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(23|456)");
    update("s3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_k_ij_abc_inv3(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("pvhp", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    diagram_stack_erase("i1");
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    perm("r3", "(456)");
    reorder("r3", "r4", "321654");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_0h1p_k_ij_c_ab_inv1(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hphp", "0000", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        update("i1", 1.0, "phhp");

        // T5a
        reorder("pphh", "r2", "3142");
        reorder("t2c", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("pphp", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T6d
        reorder("phhh", "r1", "2314");
        reorder("t1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "2134");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9e
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3124");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(23|6/45)");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_k_ij_c_ab_inv3(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pphp", "1000", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        reorder("pvhp", "r1", "2431");
        update("i1", 1.0, "r1");

        // T5a --- ?
        reorder("pphh", "r2", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T6d --- ?
        reorder("pvhh", "r1", "2314");
        reorder("t1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "2134");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9e
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3124");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(6/45)");
    reorder("r4", "r5", "321654");
    update("s3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_0h1p_c_ab(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pp", "00", "21", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2a
        reorder("pp", "ppr_", "21");
        update("i1", 1.0, "ppr_");
        restore_stack_pos(pos2);

        // T5c
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 3);
        update("i1", -0.5, "r2");
        restore_stack_pos(pos2);

        // T6b
        reorder("t1c", "r1", "21");
        mult("r1", "ph", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T6e
        reorder("ppph", "r1", "3142");
        mult("r1", "t1c", "r2", 2);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9b
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        mult("r1", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s3c", "i1", "r1", 1);
    perm("r1", "(6/45)");
    update("s3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_0h1p_a_bc(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pppp", "0000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2c
        update("i1", 0.5, "ppppr");

        // T5g
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        update("i1", 0.25, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9d
        reorder("t1c", "r1", "21");
        reorder("t1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s3c", "i1", "r1", 2);
    diagram_stack_erase("i1");
    perm("r1", "(4/56)");
    update("s3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_0h1p_abc(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_3) {
        return;
    }

    tmplt("i1", "pppp", "0000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6g
        reorder("t1c", "r1", "21");
        mult("ppph", "r1", "r2", 1);
        reorder("r2", "r3", "3412");
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);
    }
    mult("s3c", "i1", "r1", 2);
    diagram_stack_erase("i1");
    perm("r1", "(456)");
    update("s3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_0h1p_ijk_inv1(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_3) {
        return;
    }

    tmplt("i1", "hhhh", "0000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6h
        reorder("hphh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(23)");
    update("s3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_ijk_inv2(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_3) {
        return;
    }

    tmplt("i1", "phhh", "1000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6h
        reorder("vphh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(13)");
    reorder("r3", "r4", "213546");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_ijk_inv3(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_3) {
        return;
    }

    tmplt("i1", "hphh", "0100", "1234", NOT_PERM_UNIQUE);
    {
        pos2 = get_stack_pos();

        // T6h
        reorder("hphh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(12)");
    reorder("r3", "r4", "321654");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_0h1p_i_jk_inv1(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhh", "0000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2d
        update("i1", 0.5, "hhhh");

        // T5f
        reorder("pphh", "r2", "3412");
        mult("t2c", "r2", "r3", 2);
        update("i1", 0.25, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9c
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    update("s3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_i_jk_inv3(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "phhh", "1000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2d
        update("i1", 0.5, "vhhh");

        // T5f
        reorder("pphh", "r2", "3412");
        mult("s2c", "r2", "r3", 2);
        update("i1", 0.25, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9c
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(13)");
    reorder("r3", "r4", "213546");
    update("s3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_0h1p_k_ij_inv1(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hh", "00", "12", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2b
        update("i1", -1.0, "hh");

        // T5b
        reorder("pphh", "r2", "3412");
        mult("t2c", "r2", "r3", 3);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T6a
        reorder("ph", "r2", "21");
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T6f
        reorder("hphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9a
        reorder("pphh", "r2", "3142");
        mult("r2", "t1c", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(23)");
    update("s3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h1p_k_ij_inv3(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "ph", "10", "12", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2b
        update("i1", -1.0, "vh");

        // T5b
        reorder("pphh", "r2", "3412");
        mult("s2c", "r2", "r3", 3);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T6a
        reorder("ph", "r2", "21");
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T6f
        reorder("vphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9a
        reorder("pphh", "r2", "3142");
        mult("r2", "t1c", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "654321"); // 456123, then 321654
    update("s3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void t3_0h1p_contrib_to_folded(int pt_order)
{
    if (pt_order <= PT_2) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();
    reorder("s3c", "r1", "234561");
    mult("r1", "veff01", "r2", 1);
    reorder("r2", "r3", "612345");
    update("s3nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void sector_0h1p_ccsd_t3()
{
    // estimate triples amplitudes
    dg_stack_pos_t pos;

    timer_new_entry("01-S3-corr", "0h1p -- Triples pert correction +T(3)");
    timer_start("01-S3-corr");

    printf("\n 3th order perturbative correction\n");
    printf(" ---------------------------------\n");

    copy("s1c", "s1c_ccsd");
    copy("s2c", "s2c_ccsd");
    copy("veff01", "veff01_ccsd");

    tmplt("s3_0", "phhppp", "100000", "123456", IS_PERM_UNIQUE);
    tmplt("s3nw", "phhppp", "100000", "123456", IS_PERM_UNIQUE);
    t3_0h1p_const_contrib_to_triples(PT_2);
    calc_S3(PT_2);
    // folded are PT3+ => no folded here
    rename_diagram("s3nw", "s3c");
    diveps("s3c");

    // contribution to effective interaction
    copy("veff01", "veff01_pt3");
    pos = get_stack_pos();
    reorder("s3c", "r1", "145623");
    mult("r1", "pphh", "r2", 4);
    closed("r2", "r3");
    update("veff01_pt3", 0.25, "r3");
    restore_stack_pos(pos);

    // construction and diagonalization of the effective Hamiltonian
    diag_heff(0, 1, "veff01_pt3");

    timer_stop("01-S3-corr");
}


void sector_0h1p_ccsd_t4()
{
    dg_stack_pos_t pos;

    timer_new_entry("01-T(4)", "0h1p -- CCSD+T(4) triples correction");
    timer_start("01-T(4)");
    double t0 = abs_time();

    printf("\n 4th order perturbative correction\n");
    printf(" ---------------------------------\n");

    copy("t1c", "t1c_ccsd");
    copy("t2c", "t2c_ccsd");
    copy("s1c", "s1c_ccsd");
    copy("s2c", "s2c_ccsd");
    copy("veff01", "veff01_ccsd");

    /*
     * 1. Calculate 2nd order T3 in the 0h0p and 0h1p sectors
     */
    printf(" [1] calculation of 2nd order triples\n");
    // 0h0p
    tmplt("t3nw", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);
    calc_T3(PT_2);
    rename("t3nw", "t3c");
    diveps("t3c");
    copy("t3c", "t3c_pt2");
    double t1 = abs_time();
    printf("     0h0p %.3f sec\n", t1 - t0);
    // 0h1p
    tmplt("s3_0", "phhppp", "100000", "123456", IS_PERM_UNIQUE);
    tmplt("s3nw", "phhppp", "100000", "123456", IS_PERM_UNIQUE);
    t3_0h1p_const_contrib_to_triples(PT_2);
    calc_S3(PT_2);
    t3_0h1p_contrib_to_folded(PT_2);  // however, folded are PT3+
    rename_diagram("s3nw", "s3c");
    diveps("s3c");
    copy("s3c", "s3c_pt2");
    double t2 = abs_time();
    printf("     0h1p %.3f sec\n", t2 - t1);
    // at this point 't3c' and 's3c' contain triples which are correct up to PT2


    /*
     * 2. Calculate 3rd order T1 and T2 in the 0h0p and 0h1p sectors
     *    and the corresponding Veff
     */
    printf(" [2] calculation of PT3 corrections to singles and doubles\n");
    // 0h0p
    tmplt("t1nw", "hp", "00", "12", NOT_PERM_UNIQUE);
    tmplt("t2nw", "hhpp", "0000", "1234", IS_PERM_UNIQUE);
    t3_0h0p_contrib_to_singles(PT_3);
    t3_0h0p_contrib_to_doubles(PT_3);
    rename_diagram("t1nw", "t1_pt3_corr");
    rename_diagram("t2nw", "t2_pt3_corr");
    diveps("t1_pt3_corr");
    diveps("t2_pt3_corr");
    update("t1c", 1.0, "t1_pt3_corr");
    update("t2c", 1.0, "t2_pt3_corr");
    double t3 = abs_time();
    printf("     0h0p %.3f sec\n", t3 - t2);
    // 0h1p
    tmplt("s1nw", "pp", "10", "12", NOT_PERM_UNIQUE);
    tmplt("s2nw", "phpp", "1000", "1234", IS_PERM_UNIQUE);
    //const_terms_0h1p();
    //calc_S1();
    //calc_S2();
    //folded_0h1p();
    t3_0h1p_contrib_to_singles(PT_3);
    t3_0h1p_contrib_to_doubles(PT_3);
    pos = get_stack_pos();
    closed("s1nw", "x");
    restore_stack_pos(pos);
    rename_diagram("s1nw", "s1_pt3_corr");
    rename_diagram("s2nw", "s2_pt3_corr");
    diveps("s1_pt3_corr");
    diveps("s2_pt3_corr");
    update("s1c", 1.0, "s1_pt3_corr");
    update("s2c", 1.0, "s2_pt3_corr");
    double t4 = abs_time();
    printf("     0h1p %.3f sec\n", t4 - t3);
    // at this point 't1c', 't2c', 's1c', 's2c' are correct up to PT3


    /*
     * 3. Calculate 3nd order T3 in the 0h1p sector
     *    (0h0p 3rd order triples are not required in the 0h1p sector)
     */
    printf(" [3] calculation of 3rd order triples\n");
    // 0h1p
    clear("s3_0");
    t3_0h1p_const_contrib_to_triples(PT_3);
    calc_S3(PT_3);
    t3_0h1p_contrib_to_folded(PT_3);
    copy("s3nw", "s3c");
    diveps("s3c");
    double t5 = abs_time();
    printf("     0h1p %.3f sec\n", t5 - t4);
    // at this point 't3c' and 's3c' contain triples which are correct up to PT3

    /*
     * 4. Calculate 4th order Veff the 0h1p sector
     */
    printf(" [4] calculation of 4th order Heff\n");
    tmplt("s1nw", "pp", "10", "12", NOT_PERM_UNIQUE);
    const_terms_0h1p();
    calc_S1();
    t3_0h1p_contrib_to_singles(PT_3);
    closed("s1nw", "veff01");
    double t6 = abs_time();
    printf("     0h1p %.3f sec\n", t6 - t5);

    /*
     * construct and diagonalize Heff
     */
    diag_heff(0, 1, "veff01");

    /*
     * save PT-corrected diagrams for future uses in the higher sectors
     */
    copy("t1c", "t1c_pt3");
    copy("t2c", "t2c_pt3");
    rename_diagram("t3c", "t3c_pt3");
    copy("s1c", "s1c_pt3");
    copy("s2c", "s2c_pt3");
    rename_diagram("s3c", "s3c_pt3");
    copy("veff01", "veff01_pt4");
    // and restore CCSD amplitudes and Heff
    copy("t1c_ccsd", "t1c");
    copy("t2c_ccsd", "t2c");
    copy("s1c_ccsd", "s1c");
    copy("s2c_ccsd", "s2c");
    copy("veff01_ccsd", "veff01");

    timer_stop("01-T(4)");
}




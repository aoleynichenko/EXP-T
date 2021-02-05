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
 * sector02.c
 * ==========
 *
 * Multireference Fock-space coupled cluster method for model spaces
 * with two particles (Fock space sector 0h2p)
 *
 * Cluster Operator: T = T1 + T2 + S1 + S2 + X2
 * X2 =    |     |
 *         ^     ^
 *         |     |
 *         =======
 *         ^     ^
 *         ^     ^
 *         |     |
 * (excluding amplitudes with four valence indices)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip hole creation lines to turn them to valence particle
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 * Only two-particle diagrams are required.
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
#include "datamodel.h"
#include "diis.h"
#include "engine.h"
#include "heff.h"
#include "methods.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"

void init_amplitudes_0h2p();

void const_terms_0h2p();

void t3corr_0h2p();

void calc_X2();

void calc_X3();

void folded_0h2p();

void t3_0h2p_contrib_to_doubles(int pt_order);

void sector_0h2p_ccsd_t3();

void diag_T3_0h2p_i_jk_c_ab_inv12(int pt_order);
void diag_T3_0h2p_i_jk_c_ab_inv23(int pt_order);

void diag_T3_0h2p_k_ij_a_bc_inv12(int pt_order);
void diag_T3_0h2p_k_ij_a_bc_inv13(int pt_order);

void diag_T3_0h2p_k_ij_c_ab_inv12(int pt_order);
void diag_T3_0h2p_k_ij_c_ab_inv13(int pt_order);

void diag_T3_0h2p_a_bc_inv12(int pt_order);

void diag_T3_0h2p_c_ab_inv12(int pt_order);

void diag_T3_0h2p_abc_inv12(int pt_order);

void diag_T3_0h2p_k_ij_abc_inv12(int pt_order);
void diag_T3_0h2p_k_ij_abc_inv13(int pt_order);

void diag_T3_0h2p_ijk_c_ab_inv12(int pt_order);
void diag_T3_0h2p_ijk_c_ab_inv13(int pt_order);
void diag_T3_0h2p_ijk_c_ab_inv23(int pt_order);

void diag_T3_0h2p_ijk_inv12(int pt_order);
void diag_T3_0h2p_ijk_inv13(int pt_order);
void diag_T3_0h2p_ijk_inv23(int pt_order);

void diag_T3_0h2p_k_ij_inv12(int pt_order);
void diag_T3_0h2p_k_ij_inv13(int pt_order);

void diag_T3_0h2p_i_jk_inv12(int pt_order);
void diag_T3_0h2p_i_jk_inv23(int pt_order);


void t3_0h2p_const_contrib_to_triples(int pt_order)
{
    diag_T3_0h2p_i_jk_c_ab_inv12(pt_order);
    diag_T3_0h2p_k_ij_a_bc_inv13(pt_order);
    diag_T3_0h2p_k_ij_c_ab_inv13(pt_order);
    diag_T3_0h2p_k_ij_abc_inv13(pt_order);
    diag_T3_0h2p_ijk_c_ab_inv12(pt_order);
    diag_T3_0h2p_ijk_c_ab_inv13(pt_order);
    diag_T3_0h2p_ijk_c_ab_inv23(pt_order);
    diag_T3_0h2p_ijk_inv12(pt_order);
    diag_T3_0h2p_ijk_inv13(pt_order);
    diag_T3_0h2p_ijk_inv23(pt_order);
    diag_T3_0h2p_k_ij_inv13(pt_order);
    diag_T3_0h2p_i_jk_inv12(pt_order);
}


void t3_0h2p_contrib_to_doubles(int pt_order)
{
    if (pt_order <= PT_2) {
        return;
    }

    timer_new_entry("02-T2-Triples", "0h2p -- Triples contribution to doubles");
    timer_start("02-T2-Triples");

    dg_stack_pos_t pos;
    pos = get_stack_pos();

    // D10a, D11a
    // T3 reordering: 124536
    reorder("x3c", "r1", "124536");
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
    update("x2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10b, D11b
    // T3 reordering: 124356
    // permutation: P(34)
    reorder("x3c", "r1", "124356");
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
    update("x2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10c, D11c
    // T3 reordering: 145236
    // permutation: (12)
    reorder("s3c", "r1", "145236");
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
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update("x2nw", 1.0, "r3");
    restore_stack_pos(pos);

    timer_stop("02-T2-Triples");
}


void calc_X3(int pt_order)
{
    // estimate triples amplitudes
    dg_stack_pos_t pos;

    timer_new_entry("02-X3", "0h2p -- Triples equations (X3)");
    timer_start("02-X3");

    copy("x3_0", "x3nw");
    // tmplt("x3nw", "pphppp", "110000", "123456");

    //diag_T3_0h2p_i_jk_c_ab_inv12();
    diag_T3_0h2p_i_jk_c_ab_inv23(pt_order);

    diag_T3_0h2p_k_ij_a_bc_inv12(pt_order);
    //diag_T3_0h2p_k_ij_a_bc_inv13();

    diag_T3_0h2p_k_ij_c_ab_inv12(pt_order);
    //diag_T3_0h2p_k_ij_c_ab_inv13();

    diag_T3_0h2p_a_bc_inv12(pt_order);

    diag_T3_0h2p_c_ab_inv12(pt_order);

    diag_T3_0h2p_abc_inv12(pt_order);

    diag_T3_0h2p_k_ij_abc_inv12(pt_order);
    //diag_T3_0h2p_k_ij_abc_inv13();

    //diag_T3_0h2p_ijk_c_ab_inv12();
    //diag_T3_0h2p_ijk_c_ab_inv13();
    //diag_T3_0h2p_ijk_c_ab_inv23();

    //diag_T3_0h2p_ijk_inv12();
    //diag_T3_0h2p_ijk_inv13();
    //diag_T3_0h2p_ijk_inv23();

    diag_T3_0h2p_k_ij_inv12(pt_order);
    //diag_T3_0h2p_k_ij_inv13();

    //diag_T3_0h2p_i_jk_inv12();
    diag_T3_0h2p_i_jk_inv23(pt_order);

    timer_stop("02-X3");
}


void t3_0h2p_contrib_to_folded(int pt_order)
{
    timer_new_entry("folded_X3", "Folded contribution to S{02}_3");
    timer_start("folded_X3");

    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model == CC_MODEL_CCSDT_1B_PRIME) {
        // CCSDT-1'
        // here: undressed Veff
        copy("vvpp", "r1");
        closed("r1", "r2");
        //copy("veff02", "r2");
        reorder("s2c", "r3", "2341");
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "124356");
        tmplt("r6", "pphppp", "110000", "123456", NOT_PERM_UNIQUE);
        expand_diagram("r5", "r6");
        perm("r6", "(4/56)");
        update("x3nw", -1.0, "r6");
        restore_stack_pos(pos);
    }
    else {
        // F5
        // the only folded diagram to be included at the CCSDT-1 and CCSDT-2 levels of theory
        // dressed Veff
        reorder("s2c", "r1", "2341");
        mult("veff02", "r1", "r2", 1);
        reorder("r2", "r3", "124356");
        tmplt("r4", "pphppp", "110000", "123456", NOT_PERM_UNIQUE);
        expand_diagram("r3", "r4");
        perm("r4", "(4/56)");
        update("x3nw", -1.0, "r4");
        restore_stack_pos(pos);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto end_folded;
        }
        if (pt_order <= PT_2) {
            goto end_folded;
        }

        // F8
        reorder("s2c", "r1", "2341");
        reorder("veff02", "r2", "1243");
        reorder("s1c", "r3", "21");
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "124356");
        perm("r7", "(4/56)");
        update("x3nw", -1.0, "r7");
        restore_stack_pos(pos);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_3) {
            goto end_folded;
        }

        // only in the CCSDT model:
        // F6
        reorder("x3c", "r1", "234561");
        mult("veff01", "r1", "r2", 1);
        perm("r2", "(12)");
        update("x3nw", -1.0, "r2");
        restore_stack_pos(pos);

        // F7
        reorder("x3c", "r1", "345612");
        mult("veff02", "r1", "r2", 2);
        update("x3nw", -0.5, "r2");
        restore_stack_pos(pos);
    }

    end_folded:
    timer_stop("folded_X3");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_0h2p_i_jk_c_ab_inv12(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phph", "1000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vhhp", "r2", "1243");
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
        reorder("s2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        //reorder("t2c", "r1", "3412");
        reorder("pphp", "r2", "4312");
        mult("s2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("t1c", "r2", "21");
        mult("vhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphp", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("s2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("s2c", "r1", "r2", 2);
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
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("s3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("s2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(12|6/45)");
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_i_jk_c_ab_inv23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhp", "r2", "1243");
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
        reorder("x2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphp", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("t1c", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("x2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("x2c", "r1", "r2", 2);
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
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("t3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("x2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(6/45)");
    reorder("r4", "r5", "321654"); // 1 <-> 3
    update("x3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_0h2p_k_ij_a_bc_inv12(int pt_order)
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
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        // 123456 -> 132465 ---[2-3]-->
        reorder("x3c", "t5e_r1", "124653");  // interchange 2 <-> 3, then reorder 134562
        reorder("pphh", "t5e_r2", "4123");
        reorder("t2c", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        reorder("t5e_r5", "t5e_r6", "124356");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(4/56)");
    update("x3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_k_ij_a_bc_inv13(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a -- constant term
        reorder("pvpp", "r1", "2341");
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
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("s3c", "t5e_r1", "134562");
        reorder("pphh", "t5e_r2", "4123");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        reorder("t5e_r5", "t5e_r6", "124356");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(13|4/56)");
    reorder("r3", "r4", "132465"); // interchange 2 <-> 3
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_0h2p_k_ij_c_ab_inv12(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    reorder("x3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    reorder("r3", "r4", "125346");
    perm("r4", "(6/45)");
    update("x3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_k_ij_c_ab_inv13(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    reorder("s3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    reorder("r3", "r4", "125346");
    perm("r4", "(13|6/45)");
    reorder("r4", "r5", "132465");  // interchange 2<->3
    update("x3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_0h2p_a_bc_inv12(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    mult("x3c", "i1", "r1", 2);
    perm("r1", "(4/56)");
    update("x3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_0h2p_c_ab_inv12(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    mult("x3c", "i1", "r1", 1);
    perm("r1", "(6/45)");
    update("x3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_0h2p_abc_inv12(int pt_order)
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
    mult("x3c", "i1", "r1", 2);
    perm("r1", "(456)");
    update("x3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_0h2p_k_ij_abc_inv12(int pt_order)
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
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(456)");
    update("x3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_k_ij_abc_inv13(int pt_order)
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
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|456)");
    reorder("r3", "r4", "132465");  // interchange 2 <-> 3
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_0h2p_ijk_c_ab_inv12(int pt_order)
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
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12|6/45)");
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_ijk_c_ab_inv13(int pt_order)
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
        reorder("pvhp", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
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
        reorder("s2c", "s2c_21", "2143");
        reorder("pphh", "r1", "3142");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13|6/45)");
    reorder("r4", "r5", "132465");  // interchange 2 <-> 3
    update("x3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_ijk_c_ab_inv23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "ppph", "1100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("pvhp", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("t1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("s2c", "s2c_21", "2143");
        reorder("pphh", "r1", "3142");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23|6/45)");
    reorder("r4", "r5", "321654");  // interchange 1 <-> 3
    update("x3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_0h2p_ijk_inv12(int pt_order)
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

    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    reorder("r2", "r3", "456123");
    perm("r3", "(12)");
    update("x3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_ijk_inv13(int pt_order)
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

    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    reorder("r2", "r3", "456123");
    perm("r3", "(13)");
    reorder("r3", "r4", "132465");  // interchange 2 <-> 3
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_ijk_inv23(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_3) {
        return;
    }

    tmplt("i1", "pphh", "1100", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6h
        reorder("vphh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    reorder("r2", "r3", "456123");
    perm("r3", "(23)");
    reorder("r3", "r4", "321654");  // interchange 1 <-> 3
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_0h2p_k_ij_inv12(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    reorder("x3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    reorder("r2", "r3", "456123");
    update("x3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_k_ij_inv13(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    reorder("r2", "r3", "456123");
    perm("r3", "(13)");
    reorder("r3", "r4", "132465");  // interchange 2 <-> 3
    update("x3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_0h2p_i_jk_inv12(int pt_order)
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    reorder("s3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    reorder("r2", "r3", "456123");
    perm("r3", "(12)");
    update("x3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_0h2p_i_jk_inv23(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pphh", "1100", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2d
        update("i1", 0.5, "vvhh");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5f
        reorder("pphh", "r2", "3412");
        mult("x2c", "r2", "r3", 2);
        update("i1", 0.25, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9c
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    reorder("r2", "r3", "654321");  // 456123, then interchange 1 <-> 3 = 654321
    update("x3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void sector_0h2p_ccsd_t3()
{
    // estimate triples amplitudes
    dg_stack_pos_t pos;

    timer_new_entry("02-T(3)", "0h1p -- CCSD+T(3) triples correction");
    timer_start("02-T(3)");

    printf("\n 3th order perturbative correction\n");
    printf(" ---------------------------------\n");

    copy("x2c", "x2c_ccsd");
    copy("veff02", "veff02_ccsd");

    tmplt("x3_0", "pphppp", "110000", "123456", IS_PERM_UNIQUE);
    tmplt("x3nw", "pphppp", "110000", "123456", IS_PERM_UNIQUE);
    t3_0h2p_const_contrib_to_triples(PT_2);
    calc_X3(PT_2);
    t3_0h2p_contrib_to_folded(PT_2);
    rename_diagram("x3nw", "x3c");
    diveps("x3c");

    // contribution to effective interaction
    tmplt("x2nw", "pppp", "1111", "1234", IS_PERM_UNIQUE);
    const_terms_0h2p();
    calc_X2();
    t3_0h2p_contrib_to_doubles(PT_2);  // here: T3 contrib to Veff
    closed("x2nw", "veff02_pt3");

    // construction and diagonalization of the effective Hamiltonian
    diag_heff(0, 2, "veff01_pt3", "veff02_pt3");

    timer_stop("02-T(3)");
}


void sector_0h2p_ccsd_t4()
{
    dg_stack_pos_t pos;

    timer_new_entry("02-T(4)", "0h2p -- CCSD+T(4) triples correction");
    timer_start("02-T(4)");
    double t0 = abs_time();

    printf("\n 4th order perturbative correction\n");
    printf(" ---------------------------------\n");

    copy("x2c", "x2c_ccsd");
    copy("veff02", "veff02_ccsd");

    /*
     * 1. Calculate 2nd order T3 in the 0h2p sector.
     *    Triples for 0h0p and 0h1p are borrowed from previous calculation
     *    (see subroutine sector_0h1p_ccsd_t4()).
     */
    printf(" [1] calculation of 2nd order triples\n");
    // 0h2p
    tmplt("x3_0", "pphppp", "110000", "123456", IS_PERM_UNIQUE);
    tmplt("x3nw", "pphppp", "110000", "123456", IS_PERM_UNIQUE);
    t3_0h2p_const_contrib_to_triples(PT_2);
    calc_X3(PT_2);
    t3_0h2p_contrib_to_folded(PT_2);
    rename_diagram("x3nw", "x3c");
    diveps("x3c");
    copy("x3c", "x3c_pt2");
    copy("t3c_pt2", "t3c");
    copy("s3c_pt2", "s3c");
    double t1 = abs_time();
    printf("     0h2p %.3f sec\n", t1 - t0);
    // at this point 't3c', 's3c', 'x3c' contain triples which are correct up to PT2

    /*
     * 2. Calculate 3rd order T2 in the 0h2p sector
     *    and the corresponding Veff
     */
    printf(" [2] calculation of PT3 corrections to doubles\n");
    // 0h2p
    tmplt("x2nw", "pppp", "1100", "1234", IS_PERM_UNIQUE);
    //const_terms_0h1p();
    //calc_S1();
    //calc_S2();
    //folded_0h1p();
    t3_0h2p_contrib_to_doubles(PT_3);
    pos = get_stack_pos();
    closed("x2nw", "x");
    restore_stack_pos(pos);
    rename_diagram("x2nw", "x2_pt3_corr");
    diveps("x2_pt3_corr");
    update("x2c", 1.0, "x2_pt3_corr");
    double t2 = abs_time();
    printf("     0h2p %.3f sec\n", t2 - t1);
    // at this point 't1c', 't2c', 's1c', 's2c', 'x2c' are correct up to PT3

    /*
     * 3. Calculate 3nd order T3 in the 0h1p sector
     *    (0h0p 3rd order triples are not required in the 0h1p sector)
     */
    printf(" [3] calculation of 3rd order triples\n");
    // 0h2p
    clear("x3_0");
    t3_0h2p_const_contrib_to_triples(PT_3);
    calc_X3(PT_3);
    t3_0h2p_contrib_to_folded(PT_3);
    copy("x3nw", "x3c");
    diveps("x3c");
    copy("t3c_pt3", "t3c");
    copy("s3c_pt3", "s3c");
    double t3 = abs_time();
    printf("     0h2p %.3f sec\n", t3 - t2);
    // at this point 't3c', 's3c', 'x3c' contain triples which are correct up to PT3

    /*
     * 4. Calculate 4th order Veff the 0h1p sector
     */
    // Effective Hamiltonian from PT3 singles and doubles
    printf(" [4] calculation of 4th order Heff\n");
    tmplt("x2nw", "pppp", "1100", "1234", IS_PERM_UNIQUE);
    const_terms_0h2p();
    calc_X2();
    t3_0h2p_contrib_to_doubles(PT_3);
    closed("x2nw", "veff02");
    double t4 = abs_time();
    printf("     0h2p %.3f sec\n", t4 - t3);

    /*
     * construct and diagonalize Heff
     */
    diag_heff(0, 2, "veff01", "veff02");

    /*
     * save PT-corrected diagrams for future uses in the higher sectors
     */
    copy("x2c", "x2c_pt3");
    rename_diagram("x3c", "x3c_pt3");
    copy("veff02", "veff02_pt4");
    // and restore CCSD amplitudes and Heff
    copy("t1c_ccsd", "t1c");
    copy("t2c_ccsd", "t2c");
    copy("s1c_ccsd", "s1c");
    copy("s2c_ccsd", "s2c");
    copy("x2c_ccsd", "x2c");
    copy("veff01_ccsd", "veff01");
    copy("veff02_ccsd", "veff02");

    timer_stop("02-T(4)");
}







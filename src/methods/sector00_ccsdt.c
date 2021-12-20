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
 * sector00_devel.c
 *
 * Single-reference coupled cluster method (Fock space sector 0h0p):
 * triple excitations.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

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
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"

void calc_T3(int pt_order);

void t3_0h0p_contrib_to_singles(int pt_order);
void t3_0h0p_contrib_to_doubles(int pt_order);

void diag_T3_0h0p_c_ab(int pt_order);
void diag_T3_0h0p_a_bc(int pt_order);
void diag_T3_0h0p_abc(int pt_order);
void diag_T3_0h0p_k_ij_c_ab(int pt_order);
void diag_T3_0h0p_i_jk(int pt_order);
void diag_T3_0h0p_ijk(int pt_order);
void diag_T3_0h0p_k_ij(int pt_order);
void diag_T3_0h0p_k_ij_a_bc(int pt_order);
void diag_T3_0h0p_k_ij_abc(int pt_order);
void diag_T3_0h0p_i_jk_c_ab(int pt_order);
void diag_T3_0h0p_ijk_c_ab(int pt_order);


/**
 * T3 contributions to the Singles (T2) equations in the 0h0p sector.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void t3_0h0p_contrib_to_singles(int pt_order)
{
    timer_new_entry("00-T1-S7", "0h0p -- Triples contibution to singles");
    timer_start("00-T1-S7");

    if (pt_order <= PT_2) {
        return;
    }

    // S7
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t3c", "r1", "145623");
    mult("r1", "pphh", "r2", 4);
    update("t1nw", 0.25, "r2");
    restore_stack_pos(pos);

    timer_stop("00-T1-S7");
}


/**
 * T3 contributions to the Doubles (T2) equations in the 0h0p sector.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void t3_0h0p_contrib_to_doubles(int pt_order)
{
    timer_new_entry("00-T2-Triples", "0h0p -- Triples contribution to doubles");
    timer_start("00-T2-Triples");

    if (pt_order <= PT_2) {
        return;
    }

    dg_stack_pos_t pos;
    pos = get_stack_pos();

    // D10a, D11a
    // T3 reordering: 124536
    reorder("t3c", "r1", "124536");
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
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10b, D11b
    // T3 reordering: 124356
    // permutation: P(34)
    reorder("t3c", "r1", "124356");
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
    update("t2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10c, D11c
    // T3 reordering: 145236
    // permutation: (12)
    reorder("t3c", "r1", "145236");
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
    perm("r3", "(12)");
    update("t2nw", 1.0, "r3");
    restore_stack_pos(pos);

    timer_stop("00-T2-Triples");
}


/**
 * calc_T3
 *
 * Full triples (CCSDT) amplitude equations in the 0h0p sector.
 * To construct any other models by the neglection of some diagrams, the CC model
 * of the min PT order must be given.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void calc_T3(int pt_order)
{
    timer_new_entry("00-T3", "0h0p -- Triples equations (T3)");
    timer_start("00-T3");

    tmplt("t3nw", "hhhppp", "000000", "123456", IS_PERM_UNIQUE);

    diag_T3_0h0p_k_ij_a_bc(pt_order);
    diag_T3_0h0p_k_ij_abc(pt_order);
    diag_T3_0h0p_i_jk_c_ab(pt_order);
    diag_T3_0h0p_ijk_c_ab(pt_order);

    diag_T3_0h0p_c_ab(pt_order);
    diag_T3_0h0p_a_bc(pt_order);
    diag_T3_0h0p_abc(pt_order);

    diag_T3_0h0p_k_ij_c_ab(pt_order);
    diag_T3_0h0p_i_jk(pt_order);
    diag_T3_0h0p_k_ij(pt_order);
    diag_T3_0h0p_ijk(pt_order);

    timer_stop("00-T3");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_0h0p_i_jk_c_ab(int pt_order)
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME && pt_order == PT_INF) {
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
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
        reorder("t2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(1/23|6/45)");
    update("t3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_0h0p_ijk_c_ab(int pt_order)
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
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(123|6/45)");
    update("t3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_0h0p_k_ij_a_bc(int pt_order)
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
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("t3c", "t5e_r1", "134562");
        reorder("pphh", "t5e_r2", "4123");
        reorder("t2c", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(3/12|4/56)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_0h0p_k_ij_abc(int pt_order)
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
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|456)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_0h0p_c_ab(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    mult("t3c", "i1", "r1", 1);
    perm("r1", "(6/45)");
    update("t3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_0h0p_a_bc(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    mult("t3c", "i1", "r1", 2);
    perm("r1", "(4/56)");
    update("t3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_0h0p_abc(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
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
    mult("t3c", "i1", "r1", 2);
    perm("r1", "(456)");
    update("t3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_0h0p_k_ij_c_ab(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    reorder("t3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(3/12|6/45)");
    update("t3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_0h0p_ijk(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
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

    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(123)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_0h0p_i_jk(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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

        if (pt_order <= PT_3) {
            goto finish;
        }

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
    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(1/23)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_0h0p_k_ij(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    reorder("t3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(3/12)");
    update("t3nw", 1.0, "r3");
    restore_stack_pos(pos);
}



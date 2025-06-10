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
 * Fock space coupled cluster method for excitation energies: triples
 * (Fock space sector 1h1p).
 *
 * Symbolic names used for diagrams:
 * e1c   current approximation to T^{1h1p}_1 amplitudes
 * e2c   current approximation to T^{1h1p}_2 amplitudes
 * e3c   current approximation to T^{1h1p}_3 amplitudes
 * e1nw  next (new) approximation to T^{1h1p}_1 amplitudes
 * e2nw  next (new) approximation to T^{1h1p}_2 amplitudes
 * e3nw  next (new) approximation to T^{1h1p}_3 amplitudes
 */

#include "methods.h"

#include "engine.h"
#include "options.h"

void construct_triples_1h1p(int pt_order);

void diag_T3_1h1p_i_jk_c_ab_inv_12(int pt_order);
void diag_T3_1h1p_i_jk_c_ab_inv_13_const(int pt_order);
void diag_T3_1h1p_i_jk_c_ab_inv_23(int pt_order);
void diag_T3_1h1p_i_jk_c_ab_inv_21_const(int pt_order);

void diag_T3_1h1p_ijk_c_ab_inv_12(int pt_order);
void diag_T3_1h1p_ijk_c_ab_inv_21_const(int pt_order);
void diag_T3_1h1p_ijk_c_ab_inv_31_const(int pt_order);
void diag_T3_1h1p_ijk_c_ab_inv_13_const(int pt_order);
void diag_T3_1h1p_ijk_c_ab_inv_23_const(int pt_order);
void diag_T3_1h1p_ijk_c_ab_inv_33(int pt_order);

void diag_T3_1h1p_k_ij_a_bc_inv_21(int pt_order);
void diag_T3_1h1p_k_ij_a_bc_inv_12_const(int pt_order);
void diag_T3_1h1p_k_ij_a_bc_inv_31_const(int pt_order);
void diag_T3_1h1p_k_ij_a_bc_inv_32(int pt_order);

void diag_T3_1h1p_k_ij_abc_inv_12_const(int pt_order);
void diag_T3_1h1p_k_ij_abc_inv_21(int pt_order);
void diag_T3_1h1p_k_ij_abc_inv_13_const(int pt_order);
void diag_T3_1h1p_k_ij_abc_inv_32_const(int pt_order);
void diag_T3_1h1p_k_ij_abc_inv_33(int pt_order);
void diag_T3_1h1p_k_ij_abc_inv_31_const(int pt_order);

void diag_T3_1h1p_c_ab_inv_12(int pt_order);
void diag_T3_1h1p_c_ab_inv_13_const(int pt_order);

void diag_T3_1h1p_a_bc_inv_12_const(int pt_order);
void diag_T3_1h1p_a_bc_inv_21(int pt_order);

void diag_T3_1h1p_abc_inv_21(int pt_order);
void diag_T3_1h1p_abc_inv_12_const(int pt_order);
void diag_T3_1h1p_abc_inv_13_const(int pt_order);

void diag_T3_1h1p_k_ij_c_ab_inv_12(int pt_order);
void diag_T3_1h1p_k_ij_c_ab_inv_13_const(int pt_order);
void diag_T3_1h1p_k_ij_c_ab_inv_31_const(int pt_order);
void diag_T3_1h1p_k_ij_c_ab_inv_33(int pt_order);

void diag_T3_1h1p_ijk_inv_12(int pt_order);
void diag_T3_1h1p_ijk_inv_21_const(int pt_order);
void diag_T3_1h1p_ijk_inv_31_const(int pt_order);

void diag_T3_1h1p_i_jk_inv_12(int pt_order);
void diag_T3_1h1p_i_jk_inv_21_const(int pt_order);

void diag_T3_1h1p_k_ij_inv_12(int pt_order);
void diag_T3_1h1p_k_ij_inv_31_const(int pt_order);

void t3_contrib_to_singles_1h1p(int pt_order);
void t3_contrib_to_doubles_1h1p(int pt_order);


/**
 * T3 contributions to the Singles (T1) equations in the 1h1p sector.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void t3_contrib_to_singles_1h1p(int pt_order)
{
    timer_new_entry("11-T1-S7", "1h1p -- Triples contibution to singles");
    timer_start("11-T1-S7");

    if (pt_order <= PT_2) {
        return;
    }

    // S7
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e3c", "e3c_m", "123546");
    reorder("e3c_m", "r1", "145623");
    mult("r1", "pphh", "r2", 4);
    update("e1nw", -0.25, "r2"); // sign is changed!
    restore_stack_pos(pos);

    timer_stop("11-T1-S7");
}


/**
 * T3 contributions to the Doubles (T2) equations in the 1h1p sector.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void t3_contrib_to_doubles_1h1p(int pt_order)
{
    timer_new_entry("11-T2-Triples", "1h1p -- Triples contribution to doubles");
    timer_start("11-T2-Triples");

    if (pt_order <= PT_2) {
        return;
    }

    dg_stack_pos_t pos;
    pos = get_stack_pos();

    // D10a, D11a
    // T3 reordering: 124536
    reorder("e3c", "r1", "124536");
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
    update("e2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10b - inv_11
    reorder("e3c", "r1", "215346");  // 124356, then interchange electrons 1 <-> 2
    reorder("ppph", "i1", "3412");
    mult("r1", "i1", "r2", 3);
    reorder("r2", "r3", "2143");
    update("e2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D10b - inv_12
    reorder("s3c", "r1", "124356");
    reorder("ppgh", "i1", "3412");
    mult("r1", "i1", "r2", 3);
    update("e2nw", 0.5, "r2");
    restore_stack_pos(pos);

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_1A) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {

        // D11b - inv_11
        reorder("s3c", "r1", "124356");
        reorder("pphh", "_r1", "4123");
        reorder("h1c", "_r2", "21");
        mult("_r2", "_r1", "i1", 1);
        mult("r1", "i1", "r2", 3);
        update("e2nw", -0.5, "r2");
        restore_stack_pos(pos);

        // D11b - inv_12
        reorder("e3c", "r1", "125346");
        reorder("pphh", "r2", "4123");
        reorder("t1c", "r3", "21");
        mult("r3", "r2", "r4", 1);
        mult("r1", "r4", "r5", 3);
        reorder("r5", "r6", "1243");
        update("e2nw", -0.5, "r6");
        restore_stack_pos(pos);
    }

    // D10c - inv_12
    reorder("e3c", "r1", "145236");
    reorder("hphh", "i1", "1342");
    mult("r1", "i1", "r2", 3);
    reorder("r2", "r3", "1423");
    update("e2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D10c - inv_21
    reorder("h3c", "r1", "145236");
    reorder("vphh", "i1", "1342");
    mult("r1", "i1", "r2", 3);
    reorder("r2", "r3", "1423");
    reorder("r3", "r4", "2143");  // interchange electrons: 1 <-> 2
    update("e2nw", -0.5, "r4");
    restore_stack_pos(pos);

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_1A) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {

        // D11c - inv_12
        reorder("h3c", "r1", "145236");
        reorder("pphh", "_r1", "3421");
        mult("s1c", "_r1", "i1", 1);
        mult("r1", "i1", "r2", 3);
        reorder("r2", "r3", "1423");
        reorder("r3", "r4", "2143");
        update("e2nw", -0.5, "r4");
        restore_stack_pos(pos);

        // D11c - inv_21
        reorder("e3c", "r1", "145236");
        reorder("pphh", "_r1", "3421");
        mult("t1c", "_r1", "i1", 1);
        mult("r1", "i1", "r2", 3);
        reorder("r2", "r3", "1423");
        update("e2nw", -0.5, "r3");
        restore_stack_pos(pos);
    }

    timer_stop("11-T2-Triples");
}


/**
 * construct_triples_1h1p
 *
 * Full triples (CCSDT) amplitude equations in the 1h1p sector.
 * To construct any other models by the neglection of some diagrams, the CC model
 * of the min PT order must be given.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void construct_triples_1h1p(int pt_order)
{
    timer_new_entry("11-T3", "1h1p -- Triples equations (T3)");
    timer_start("11-T3");

    copy("e3_0", "e3nw");
    //tmplt("e3nw", "phhphp", "100010", "123456", IS_PERM_UNIQUE);

    diag_T3_1h1p_i_jk_c_ab_inv_12(pt_order);
    //diag_T3_1h1p_i_jk_c_ab_inv_13_const(pt_order);
    diag_T3_1h1p_i_jk_c_ab_inv_23(pt_order);
    //diag_T3_1h1p_i_jk_c_ab_inv_21_const(pt_order);

    diag_T3_1h1p_ijk_c_ab_inv_12(pt_order);
    //diag_T3_1h1p_ijk_c_ab_inv_21_const(pt_order);
    //diag_T3_1h1p_ijk_c_ab_inv_31_const(pt_order);
    //diag_T3_1h1p_ijk_c_ab_inv_13_const(pt_order);
    //diag_T3_1h1p_ijk_c_ab_inv_23_const(pt_order);
    diag_T3_1h1p_ijk_c_ab_inv_33(pt_order);

    diag_T3_1h1p_k_ij_a_bc_inv_21(pt_order);
    //diag_T3_1h1p_k_ij_a_bc_inv_12_const(pt_order);
    //diag_T3_1h1p_k_ij_a_bc_inv_31_const(pt_order);
    diag_T3_1h1p_k_ij_a_bc_inv_32(pt_order);

    //diag_T3_1h1p_k_ij_abc_inv_12_const(pt_order);
    diag_T3_1h1p_k_ij_abc_inv_21(pt_order);
    //diag_T3_1h1p_k_ij_abc_inv_13_const(pt_order);
    //diag_T3_1h1p_k_ij_abc_inv_32_const(pt_order);
    diag_T3_1h1p_k_ij_abc_inv_33(pt_order);
    //diag_T3_1h1p_k_ij_abc_inv_31_const(pt_order);

    diag_T3_1h1p_c_ab_inv_12(pt_order);
    //diag_T3_1h1p_c_ab_inv_13_const(pt_order);

    //diag_T3_1h1p_a_bc_inv_12_const(pt_order);
    diag_T3_1h1p_a_bc_inv_21(pt_order);

    diag_T3_1h1p_abc_inv_21(pt_order);
    //diag_T3_1h1p_abc_inv_12_const(pt_order);
    //diag_T3_1h1p_abc_inv_13_const(pt_order);

    diag_T3_1h1p_k_ij_c_ab_inv_12(pt_order);
    //diag_T3_1h1p_k_ij_c_ab_inv_13_const(pt_order);
    //diag_T3_1h1p_k_ij_c_ab_inv_31_const(pt_order);
    diag_T3_1h1p_k_ij_c_ab_inv_33(pt_order);

    diag_T3_1h1p_ijk_inv_12(pt_order);
    //diag_T3_1h1p_ijk_inv_21_const(pt_order);
    //diag_T3_1h1p_ijk_inv_31_const(pt_order);

    diag_T3_1h1p_i_jk_inv_12(pt_order);
    //diag_T3_1h1p_i_jk_inv_21_const(pt_order);

    diag_T3_1h1p_k_ij_inv_12(pt_order);
    //diag_T3_1h1p_k_ij_inv_31_const(pt_order);

    timer_stop("11-T3");
}


void const_terms_1h1p_triples(int pt_order)
{
    timer_new_entry("11-T3-const", "1h1p -- Const contributions to triples");
    timer_start("11-T3-const");

    //diag_T3_1h1p_i_jk_c_ab_inv_12(pt_order);
    diag_T3_1h1p_i_jk_c_ab_inv_13_const(pt_order);
    //diag_T3_1h1p_i_jk_c_ab_inv_23(pt_order);
    diag_T3_1h1p_i_jk_c_ab_inv_21_const(pt_order);

    //diag_T3_1h1p_ijk_c_ab_inv_12(pt_order);
    diag_T3_1h1p_ijk_c_ab_inv_21_const(pt_order);
    diag_T3_1h1p_ijk_c_ab_inv_31_const(pt_order);
    diag_T3_1h1p_ijk_c_ab_inv_13_const(pt_order);
    diag_T3_1h1p_ijk_c_ab_inv_23_const(pt_order);
    //diag_T3_1h1p_ijk_c_ab_inv_33(pt_order);

    //diag_T3_1h1p_k_ij_a_bc_inv_21(pt_order);
    diag_T3_1h1p_k_ij_a_bc_inv_12_const(pt_order);
    diag_T3_1h1p_k_ij_a_bc_inv_31_const(pt_order);
    //diag_T3_1h1p_k_ij_a_bc_inv_32(pt_order);

    diag_T3_1h1p_k_ij_abc_inv_12_const(pt_order);
    //diag_T3_1h1p_k_ij_abc_inv_21(pt_order);
    diag_T3_1h1p_k_ij_abc_inv_13_const(pt_order);
    diag_T3_1h1p_k_ij_abc_inv_32_const(pt_order);
    //diag_T3_1h1p_k_ij_abc_inv_33(pt_order);
    diag_T3_1h1p_k_ij_abc_inv_31_const(pt_order);

    //diag_T3_1h1p_c_ab_inv_12(pt_order);
    diag_T3_1h1p_c_ab_inv_13_const(pt_order);

    diag_T3_1h1p_a_bc_inv_12_const(pt_order);
    //diag_T3_1h1p_a_bc_inv_21(pt_order);

    //diag_T3_1h1p_abc_inv_21(pt_order);
    diag_T3_1h1p_abc_inv_12_const(pt_order);
    diag_T3_1h1p_abc_inv_13_const(pt_order);

    //diag_T3_1h1p_k_ij_c_ab_inv_12(pt_order);
    diag_T3_1h1p_k_ij_c_ab_inv_13_const(pt_order);
    diag_T3_1h1p_k_ij_c_ab_inv_31_const(pt_order);
    //diag_T3_1h1p_k_ij_c_ab_inv_33(pt_order);

    //diag_T3_1h1p_ijk_inv_12(pt_order);
    diag_T3_1h1p_ijk_inv_21_const(pt_order);
    diag_T3_1h1p_ijk_inv_31_const(pt_order);

    //diag_T3_1h1p_i_jk_inv_12(pt_order);
    diag_T3_1h1p_i_jk_inv_21_const(pt_order);

    //diag_T3_1h1p_k_ij_inv_12(pt_order);
    diag_T3_1h1p_k_ij_inv_31_const(pt_order);

    timer_stop("11-T3-const");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_1h1p_i_jk_c_ab_inv_12(int pt_order)
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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        // TODO: merge two reorderings of "e3c"
        reorder("e3c", "e3c_23", "132465");
        reorder("e3c_23", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("t2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(46)");
    update("e3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_i_jk_c_ab_inv_13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hhhg", "r2", "1243");
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
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphg", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
            goto finish;
        }

        // T4b
        reorder("h1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("t2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10a
        reorder("h1c", "r1", "21");
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
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    reorder("r4", "r4_23", "132465");
    update("e3_0", 1.0, "r4_23");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_i_jk_c_ab_inv_23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vhhg", "r2", "1243");
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
        reorder("e2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        // TODO: "r1" - for what?
        //reorder("t2c", "r1", "3412");
        reorder("pphg", "r2", "4312");
        mult("s2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
            goto finish;
        }

        // T4b
        reorder("h1c", "r2", "21");
        mult("vhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("e2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("s2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
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
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("t3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("e2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(13)");
    reorder("r4", "r4_231", "231564");
    update("e3nw", 1.0, "r4_231");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_i_jk_c_ab_inv_21_const(int pt_order)
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME && pt_order == PT_INF) {
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("h3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("s2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(13|56)");
    reorder("r4", "r4_12","213546"); // interchange electrons: 1 <-> 2
    update("e3_0", 1.0, "r4_12");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_1h1p_ijk_c_ab_inv_12(int pt_order)
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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23|46)");
    update("e3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_c_ab_inv_21_const(int pt_order)
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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13|56)");
    reorder("r4", "r4_21", "213546");
    update("e3_0", 1.0, "r4_21");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_c_ab_inv_31_const(int pt_order)
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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12|56)");
    reorder("r4", "r4_312", "312645");
    update("e3_0", 1.0, "r4_312");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_c_ab_inv_13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("phhg", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("h1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23)");
    reorder("r4", "r4_23", "132465");
    update("e3_0", 1.0, "r4_23");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_c_ab_inv_23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("phhg", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13)");
    reorder("r4", "r4_231", "231564");
    update("e3_0", 1.0, "r4_231");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_c_ab_inv_33(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hphh", "0110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("e2c", "e2c_m", "2134"); // minus sign here!
        reorder("e2c_m", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        //reorder("pvhg", "r0", "2431");
        mult("t1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("h1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("e2c", "e2c_m", "2134"); // minus sign here!
        reorder("e2c_m", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12)");
    reorder("r4", "r5", "312465"); // minus sign here!
    update("e3nw", -1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_1h1p_k_ij_a_bc_inv_21(int pt_order)
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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        // TODO: merge two reorderings of e3c
        reorder("e3c", "e3c_231", "231564");  // interchange electrons: [231]
        reorder("e3c_231", "t5e_r1", "134562");
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

    perm("r3", "(13)");
    reorder("r3", "r4", "213546"); // interchange electrons 1 <-> 2
    update("e3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_a_bc_inv_12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhpp", "0100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phgp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3e
        reorder("h2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppgp", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("t1c", "t1c_21", "21");
        reorder("h1c", "h1c_21", "21");
        mult("t1c_21", "phhh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("t1c", "t1c_21", "21");
        reorder("h1c", "h1c_21", "21");
        mult("t1c_21", "pphh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
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
        reorder("h2c", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(23|46)");
    update("e3_0", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_a_bc_inv_31_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
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
        reorder("t1c", "t1c_21", "21");
        mult("t1c_21", "pvhh", "r2", 1);
        mult("t1c_21", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("h3c", "t5e_r1", "134562");
        reorder("pphh", "t5e_r2", "4123");
        reorder("s2c", "s2c_21", "2143");  // interchange electrons 1 <-> 2
        reorder("s2c_21", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    reorder("r3", "r4", "312645");  // move 3rd electron to the 1st position
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_a_bc_inv_32(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvgp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3e
        reorder("h2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppgp", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("t1c", "t1c_21", "21");
        reorder("h1c", "h1c_21", "21");
        mult("t1c_21", "pvhh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("t1c", "t1c_21", "21");
        reorder("h1c", "h1c_21", "21");
        mult("t1c_21", "pphh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
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
        reorder("e2c", "e2c_21", "2143");  // interchange electrons 1 <-> 2
        reorder("e2c_21", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(46)");
    reorder("r3", "r4", "321654");  // interchange electrons: 321
    update("e3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_1h1p_k_ij_abc_inv_12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhpp", "0100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppgh", "r1", "1342");
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
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
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
    perm("r3", "(23|46)");
    update("e3_0", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_abc_inv_21(int pt_order)
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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|56)");
    reorder("r3", "r4", "213546"); // interchange electrons 1 <-> 2
    update("e3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_abc_inv_13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppph", "r1", "1342");
	    reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
	    reorder("phhg", "r0", "2431");
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
        reorder("h1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
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
    perm("r3", "(23|45)");
    reorder("r3", "r4", "132465"); // interchange electrons: 2 <-> 3
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_abc_inv_32_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppgh", "r1", "1342");
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
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
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
    reorder("r2", "r3", "124356");
    perm("r3", "(46)");
    reorder("r3", "r4", "321654");  // interchange electrons 1 <-> 3
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_abc_inv_33(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pphp", "1010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("e2c", "e2c_m_21", "2134");  // sign is changed
        reorder("e2c_m_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("e2c", "e2c_m_21", "2134");  // sign is changed
        reorder("e2c_m_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(45)");
    reorder("r3", "r4", "312465"); // sign is changed again
    update("e3nw", -1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_abc_inv_31_const(int pt_order)
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
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(56)");
    // interchange electrons: 3rd electron to the 1st position
    reorder("r3", "r4", "312645");
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_1h1p_c_ab_inv_12(int pt_order)
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
    mult("e3c", "i1", "r1", 1);
    perm("r1", "(46)");
    update("e3nw", 1.0, "r1");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_c_ab_inv_13_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hp", "10", "21", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2a
        reorder("pg", "pgr", "21");
        update("i1", 1.0, "pgr");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5c
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 3);
        update("i1", -0.5, "r2");
        restore_stack_pos(pos2);

        // T6b
        reorder("h1c", "r1", "21");
        mult("r1", "ph", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T6e
        reorder("ppgh", "r1", "3142");
        mult("r1", "t1c", "r2", 2);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9b
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        mult("r1", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s3c", "i1", "r1", 1);
    reorder("r1", "r2", "132465");
    update("e3_0", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_1h1p_a_bc_inv_12_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hppp", "1000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2c
        update("i1", 0.5, "ppgp");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5g
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        update("i1", 0.25, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9d
        reorder("t1c", "r1", "21");
        reorder("h1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s3c", "i1", "r1", 2);
    perm("r1", "(46)");
    update("e3_0", 1.0, "r1");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_a_bc_inv_21(int pt_order)
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
    reorder("e3c", "e3c_21", "213546");
    mult("e3c_21", "i1", "r1", 2);
    reorder("r1", "r2", "213546");
    update("e3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_1h1p_abc_inv_21(int pt_order)
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
    reorder("e3c", "e3c_21", "213546");
    mult("e3c_21", "i1", "r1", 2);
    perm("r1", "(56)");
    reorder("r1", "r2", "213546");  // interchange electrons 1 <-> 2
    update("e3nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_abc_inv_12_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
        return;
    }

    tmplt("i1", "hppp", "1000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6g
        reorder("t1c", "r1", "21");
        mult("ppgh", "r1", "r2", 1);
        reorder("r2", "r3", "3412");
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);
    }
    mult("s3c", "i1", "r1", 2);
    perm("r1", "(46)");
    update("e3_0", 1.0, "r1");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_abc_inv_13_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
        return;
    }

    tmplt("i1", "phpp", "0100", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6g
        reorder("h1c", "r1", "21");
        mult("ppph", "r1", "r2", 1);
        reorder("r2", "r3", "3412");
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);
    }
    mult("s3c", "i1", "r1", 2);
    perm("r1", "(45)");
    reorder("r1", "r2", "132465");
    update("e3_0", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_1h1p_k_ij_c_ab_inv_12(int pt_order)
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
    reorder("e3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(23|46)");
    update("e3nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_c_ab_inv_13_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhp", "0100", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        reorder("phhg", "r1", "2431");
        update("i1", 1.0, "r1");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5a
        reorder("pphh", "r2", "3142");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T6d
        reorder("phhh", "r1", "2314");
        reorder("h1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "2134");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9e
        reorder("h1c", "r1", "21");
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
    perm("r4", "(23)");
    reorder("r4", "r5", "132465");
    update("e3_0", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_c_ab_inv_31_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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

        // T5a
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

        // T6d
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
    reorder("h3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(56)");
    reorder("r4", "r5", "312645");
    update("e3_0", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_c_ab_inv_33(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "phhp", "1100", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        update("i1", 1.0, "pvhg");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5a
        reorder("pphh", "r2", "3142");
        reorder("e2c", "e2c_m_21", "2134"); // sign is changed
        reorder("e2c_m_21", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("pphg", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T6d
        reorder("pvhh", "r1", "2314");
        reorder("h1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "2134");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9e
        reorder("h1c", "r1", "21");
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
    reorder("r4", "r5", "312465"); // sign is changed
    update("e3nw", -1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_1h1p_ijk_inv_12(int pt_order)
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

    reorder("e3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(23)");
    update("e3nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_inv_21_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
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

    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(13)");
    reorder("r3", "r4", "213546");
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_ijk_inv_31_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
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

    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(12)");
    // place 3rd electron to the 1st position
    reorder("r3", "r4", "312645");
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_1h1p_i_jk_inv_12(int pt_order)
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
    reorder("e3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    update("e3nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_i_jk_inv_21_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(13)");
    // interchange electrons 1 <-> 2
    reorder("r3", "r4", "213546");
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_1h1p_k_ij_inv_12(int pt_order)
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
    reorder("e3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(23)");
    update("e3nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_T3_1h1p_k_ij_inv_31_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    reorder("r3", "r4", "312645");
    update("e3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Seven folded diagrams.
 * You can find them in (for example):
 * M. Musial, R. J. Bartlett, J. Chem. Phys. 121, 1670 (2004); doi: 10.1063/1.1765096
 */
void folded_diagrams_triples_1h1p(int pt_order)
{
    timer_new_entry("folded_E3", "Folded contribution to S{11}_3");
    timer_start("folded_E3");

    dg_stack_pos_t pos = get_stack_pos();

    if (pt_order == PT_INF && cc_opts->cc_model == CC_MODEL_CCSDT_1B_PRIME) {
        // CCSDT-1'
        // here: undressed Veff = V_closed is used in contractions instead of Veff{1h1p}

        // F4
        // T{1h0p}_2 * V_closed
        copy("vhpg", "r0");
        closed("r0", "V_closed");
        reorder("h2c", "r1", "1243");
        reorder("V_closed", "r2", "1342");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "145236");
        tmplt("r5", "phhphp", "100010", "123456", NOT_PERM_UNIQUE);
        expand_diagram("r4", "r5");
        perm("r5", "(46)");
        update("e3nw", 1.0, "r5");
        restore_stack_pos(pos);

        // F5
        // T{0h1p}_2 * V_closed
        copy("vhpg", "r0");
        closed("r0", "V_closed");
        reorder("s2c", "r1", "2341");
        reorder("V_closed", "V_closed_21", "2143");
        mult("V_closed_21", "r1", "r2", 1);
        reorder("r2", "r3", "124356");
        tmplt("r4", "hphhpp", "010100", "213546", NOT_PERM_UNIQUE);
        expand_diagram("r3", "r4");
        perm("r4", "(13)");
        reorder("r4", "r5", "213546");
        update("e3nw", -1.0, "r5");
        restore_stack_pos(pos);

    }
    else {
        // the only two folded diagrams to be included at the CCSDT-1 and CCSDT-2 levels of theory
        // Veff is dressed

        // F4
        // T{1h0p}_2 * veff{1h1p}
        reorder("h2c", "r1", "1243");
        reorder("veff11", "r2", "1342");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "145236");
        tmplt("r5", "phhphp", "100010", "123456", NOT_PERM_UNIQUE);
        expand_diagram("r4", "r5");
        perm("r5", "(46)");
        update("e3nw", 1.0, "r5");
        restore_stack_pos(pos);

        // F5
        // T{0h1p}_2 * veff{1h1p}
        reorder("s2c", "r1", "2341");
        reorder("veff11", "veff11_21", "2143");
        mult("veff11_21", "r1", "r2", 1);
        reorder("r2", "r3", "124356");
        tmplt("r4", "hphhpp", "010100", "213546", NOT_PERM_UNIQUE);
        expand_diagram("r3", "r4");
        perm("r4", "(13)");
        reorder("r4", "r5", "213546");
        update("e3nw", -1.0, "r5");
        restore_stack_pos(pos);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto end_folded;
        }
        if (pt_order <= PT_2) {
            goto end_folded;
        }

        // F6
        // T{0h1p}_1 * T{1h0p}_2 * veff{1h1p}
        reorder("s1c", "r1", "21");
        reorder("h2c", "r2", "1243");
        reorder("veff11", "r3", "1423");
        mult("r3", "r1", "r4", 1);
        reorder("r4", "r5", "1423");
        mult("r5", "r2", "r6", 1);
        reorder("r6", "r7", "145236");
        perm("r7", "(46)");
        update("e3nw", 1.0, "r7");
        restore_stack_pos(pos);

        // F7
        // T{1h0p}_1 * T{0h1p}_2 * veff{1h1p}
        reorder("veff11", "veff11_21", "2143");
        reorder("veff11_21", "r1", "2341");
        reorder("s2c", "r2", "2341");
        mult("h1c", "r1", "r3", 1);
        mult("r3", "r2", "r4", 1);
        reorder("r4", "r5", "124356");
        perm("r5", "(13)");
        reorder("r5", "r6", "213546");
        update("e3nw", 1.0, "r6");
        restore_stack_pos(pos);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_3) {
            goto end_folded;
        }

        // F1
        // T{1h1p}_3 * veff{0h1p}
        reorder("e3c", "r1", "234561");
        mult("veff01", "r1", "r2", 1);
        update("e3nw", -1.0, "r2");
        restore_stack_pos(pos);

        // F2
        // T{1h1p}_3 * veff{1h0p}
        reorder("e3c", "r1", "123465");
        reorder("veff10", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "123465");
        update("e3nw", 1.0, "r4");
        restore_stack_pos(pos);

        // F3
        // T{1h1p}_3 * veff{1h1p}
        reorder("e3c", "r1", "234615");
        reorder("veff11", "r2", "1432");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "512364");
        update("e3nw", 1.0, "r4");
        restore_stack_pos(pos);
    }

    end_folded:
    timer_stop("folded_E3");
}



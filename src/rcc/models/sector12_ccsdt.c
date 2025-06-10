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
 * Fock space sector 1h2p. CCSDT model
 */


#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ccutils.h"
#include "engine.h"
#include "diis.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"


/*
 * declarations of functions used in this file
 */

void t3_contrib_to_doubles_1h2p(int pt_order);

void t3_contrib_to_folded_1h2p(int pt_order);

void diag_T3_1h2p_k_ij_a_bc__inv_h1_p23_const(int pt_order);
void diag_T3_1h2p_k_ij_a_bc__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_k_ij_a_bc__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_k_ij_a_bc__inv_h1_p12(int pt_order);

void diag_T3_1h2p_k_ij_abc__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_k_ij_abc__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_k_ij_abc__inv_h2_p12_const(int pt_order);
void diag_T3_1h2p_k_ij_abc__inv_h1_p23_const(int pt_order);
void diag_T3_1h2p_k_ij_abc__inv_h1_p12(int pt_order);
void diag_T3_1h2p_k_ij_abc__inv_h3_p23_const(int pt_order);

void diag_T3_1h2p_i_jk_c_ab__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_i_jk_c_ab__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_i_jk_c_ab__inv_h1_p23_const(int pt_order);
void diag_T3_1h2p_i_jk_c_ab__inv_h3_p23(int pt_order);

void diag_T3_1h2p_ijk_c_ab__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_ijk_c_ab__inv_h3_p23_const(int pt_order);
void diag_T3_1h2p_ijk_c_ab__inv_h1_p23_const(int pt_order);
void diag_T3_1h2p_ijk_c_ab__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_ijk_c_ab__inv_h3_p13_const(int pt_order);
void diag_T3_1h2p_ijk_c_ab__inv_h2_p12_const(int pt_order);

void diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_k_ij_c_ab__inv_h3_p13_const(int pt_order);
void diag_T3_1h2p_k_ij_c_ab__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_m3c(int pt_order);

void diag_T3_1h2p_k_ij__inv_h3_p12(int pt_order);
void diag_T3_1h2p_k_ij__inv_h2_p13_const(int pt_order);

void diag_T3_1h2p_i_jk__inv_h1_p23_const(int pt_order);
void diag_T3_1h2p_i_jk__inv_h2_p13_const(int pt_order);

void diag_T3_1h2p_ijk__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_ijk__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_ijk__inv_h1_p23_const(int pt_order);

void diag_T3_1h2p_c_ab__inv_h3_p12_const(int pt_order);
void diag_T3_1h2p_c_ab__inv_h1_p23(int pt_order);

void diag_T3_1h2p_a_bc__inv_h1_p23(int pt_order);
void diag_T3_1h2p_a_bc__inv_h3_p12_const(int pt_order);

void diag_T3_1h2p_abc__inv_h1_p23(int pt_order);
void diag_T3_1h2p_abc__inv_h2_p13_const(int pt_order);
void diag_T3_1h2p_abc__inv_h3_p12_const(int pt_order);


void ampl_equations_1h2p_folded_triples_F1a();
void ampl_equations_1h2p_folded_triples_F1b();
void ampl_equations_1h2p_folded_triples_F2a();
void ampl_equations_1h2p_folded_triples_F2b();
void ampl_equations_1h2p_folded_triples_F2c();
void ampl_equations_1h2p_folded_triples_F3a();
void ampl_equations_1h2p_folded_triples_F3b();
void ampl_equations_1h2p_folded_triples_F4a();
void ampl_equations_1h2p_folded_triples_F4b();
void ampl_equations_1h2p_folded_triples_F4c();
void ampl_equations_1h2p_folded_triples_F5a();
void ampl_equations_1h2p_folded_triples_F5b();
void ampl_equations_1h2p_folded_triples_F6a();
void ampl_equations_1h2p_folded_triples_F6b();
void ampl_equations_1h2p_folded_triples_F7a();
void ampl_equations_1h2p_folded_triples_F7b();
void ampl_equations_1h2p_folded_triples_F8();
void ampl_equations_1h2p_folded_triples_F9a();
void ampl_equations_1h2p_folded_triples_F9b();
void ampl_equations_1h2p_folded_triples_F10a();
void ampl_equations_1h2p_folded_triples_F10b();
void ampl_equations_1h2p_folded_triples_F11();
void ampl_equations_1h2p_folded_triples_F12();


/**
 * T3 contributions to the Doubles (T2) equations in the 1h2p sector.
 *
 * pt_order == PT_INF:
 *   hierarchy of CC models is used for the selection of diagrams.
 * pt_order != PT_INF:
 *   only PT arguments are used for the selection of diagrams. A diagram
 *   contributes if it arises at the 'pt_order'th order of PT.
 */
void t3_contrib_to_doubles_1h2p(int pt_order)
{
    timer_new_entry("12-T2-Triples", "1h2p -- Triples contribution to doubles");
    timer_start("12-T2-Triples");

    if (pt_order <= PT_2) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // D10a
    reorder("m3c", "r1", "124653");
    mult("r1", "ph", "r2", 2);
    update("m2c", -1.0, "r2");
    restore_stack_pos(pos);

    // D10b-1
    reorder("x3c", "r1", "124356");
    reorder("ppgh", "r2", "3412");
    mult("r1", "r2", "r3", 3);
    update("m2c", 0.5, "r3");
    restore_stack_pos(pos);

    // D10b-2
    reorder("m3c", "r1", "126345");
    reorder("ppph", "r2", "3412");
    mult("r1", "r2", "r3", 3);
    reorder("r3", "r4", "1243");
    update("m2c", -0.5, "r4");
    restore_stack_pos(pos);

    // D10c
    reorder("e3c", "r1", "145623");
    mult("r1", "vphh", "r2", 3);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update("m2c", -0.5, "r3");
    restore_stack_pos(pos);

    // D11a
    reorder("m3c", "r1", "124635");
    reorder("pphh", "r2", "3142");
    mult("r2", "t1c", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update("m2c", -1.0, "r4");
    restore_stack_pos(pos);

    // D11b-1
    reorder("m3c", "r1", "126345");
    reorder("pphh", "r2", "3412");
    reorder("t1c", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "1243");
    update("m2c", 0.5, "r6");
    restore_stack_pos(pos);

    // D11b-2
    reorder("x3c", "r1", "124356");
    reorder("pphh", "r2", "4123");
    reorder("h1c", "r3", "21");
    mult("r3", "r2", "r4", 1);
    mult("r1", "r4", "r5", 3);
    update("m2c", -0.5, "r5");
    restore_stack_pos(pos);

    // D11c
    reorder("e3c", "r1", "145623");
    reorder("pphh", "r2", "2341");
    mult("s1c", "r2", "r3", 1);
    mult("r1", "r3", "r4", 3);
    reorder("r4", "r5", "1423");
    perm("r5", "(12)");
    update("m2c", -0.5, "r5");
    restore_stack_pos(pos);

    timer_stop("12-T2-Triples");
}


/**
 * T3 contributions to folded diagrams in the 1h2p sector.
 */
void t3_contrib_to_folded_1h2p(int pt_order)
{
    timer_new_entry("folded_M3", "Folded contribution to T{1h2p}_3");
    timer_start("folded_M3");

    ampl_equations_1h2p_folded_triples_F1a();
    ampl_equations_1h2p_folded_triples_F1b();
    ampl_equations_1h2p_folded_triples_F2a();
    ampl_equations_1h2p_folded_triples_F2b();
    ampl_equations_1h2p_folded_triples_F2c();
    ampl_equations_1h2p_folded_triples_F3a();
    ampl_equations_1h2p_folded_triples_F3b();
    ampl_equations_1h2p_folded_triples_F4a();
    ampl_equations_1h2p_folded_triples_F4b();
    ampl_equations_1h2p_folded_triples_F4c();
    ampl_equations_1h2p_folded_triples_F5a();
    ampl_equations_1h2p_folded_triples_F5b();
    ampl_equations_1h2p_folded_triples_F6a();
    ampl_equations_1h2p_folded_triples_F6b();
    ampl_equations_1h2p_folded_triples_F7a();
    ampl_equations_1h2p_folded_triples_F7b();
    ampl_equations_1h2p_folded_triples_F8();
    ampl_equations_1h2p_folded_triples_F9a();
    ampl_equations_1h2p_folded_triples_F9b();
    ampl_equations_1h2p_folded_triples_F10a();
    ampl_equations_1h2p_folded_triples_F10b();
    ampl_equations_1h2p_folded_triples_F11();
    ampl_equations_1h2p_folded_triples_F12();

    timer_stop("folded_M3");
}


void construct_triples_1h2p_const_terms(int pt_order)
{
    diag_T3_1h2p_k_ij_a_bc__inv_h1_p23_const(PT_INF);
    diag_T3_1h2p_k_ij_a_bc__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_k_ij_a_bc__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_k_ij_a_bc__inv_h1_p12(PT_INF); // depends on m2c and m3c

    diag_T3_1h2p_k_ij_abc__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_k_ij_abc__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_k_ij_abc__inv_h2_p12_const(PT_INF);
    diag_T3_1h2p_k_ij_abc__inv_h1_p23_const(PT_INF);
    //diag_T3_1h2p_k_ij_abc__inv_h1_p12(PT_INF); // depends on m2c
    diag_T3_1h2p_k_ij_abc__inv_h3_p23_const(PT_INF);

    diag_T3_1h2p_i_jk_c_ab__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_i_jk_c_ab__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_i_jk_c_ab__inv_h1_p23_const(PT_INF);
    //diag_T3_1h2p_i_jk_c_ab__inv_h3_p23(PT_INF); // depends on m2c

    diag_T3_1h2p_ijk_c_ab__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_ijk_c_ab__inv_h3_p23_const(PT_INF);
    diag_T3_1h2p_ijk_c_ab__inv_h1_p23_const(PT_INF);
    diag_T3_1h2p_ijk_c_ab__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_ijk_c_ab__inv_h3_p13_const(PT_INF);
    diag_T3_1h2p_ijk_c_ab__inv_h2_p12_const(PT_INF);

    diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_k_ij_c_ab__inv_h3_p13_const(PT_INF);
    diag_T3_1h2p_k_ij_c_ab__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_m3c(PT_INF); // depends on m3c

    //diag_T3_1h2p_k_ij__inv_h3_p12(PT_INF); // depends on m3c
    diag_T3_1h2p_k_ij__inv_h2_p13_const(PT_INF);

    diag_T3_1h2p_i_jk__inv_h1_p23_const(PT_INF);
    diag_T3_1h2p_i_jk__inv_h2_p13_const(PT_INF);

    diag_T3_1h2p_ijk__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_ijk__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_ijk__inv_h1_p23_const(PT_INF);

    diag_T3_1h2p_c_ab__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_c_ab__inv_h1_p23(PT_INF); // depends on m3c

    //diag_T3_1h2p_a_bc__inv_h1_p23(PT_INF); // depends on m3c
    diag_T3_1h2p_a_bc__inv_h3_p12_const(PT_INF);

    //diag_T3_1h2p_abc__inv_h1_p23(PT_INF); // depends on m3c
    diag_T3_1h2p_abc__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_abc__inv_h3_p12_const(PT_INF);
}


/**
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
void construct_triples_1h2p(int pt_order)
{
    timer_new_entry("12-T3", "1h2p -- Triples equations (T3)");
    timer_start("12-T3");

    //tmplt("m3nw", "pphpph", "110001", "123456", IS_PERM_UNIQUE);
    copy("m3_0", "m3nw");

    //diag_T3_1h2p_k_ij_a_bc__inv_h1_p23_const(PT_INF);
    //diag_T3_1h2p_k_ij_a_bc__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_k_ij_a_bc__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_k_ij_a_bc__inv_h1_p12(PT_INF); // depends on m2c and m3c

    //diag_T3_1h2p_k_ij_abc__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_k_ij_abc__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_k_ij_abc__inv_h2_p12_const(PT_INF);
    //diag_T3_1h2p_k_ij_abc__inv_h1_p23_const(PT_INF);
    diag_T3_1h2p_k_ij_abc__inv_h1_p12(PT_INF); // depends on m2c
    //diag_T3_1h2p_k_ij_abc__inv_h3_p23_const(PT_INF);

    //diag_T3_1h2p_i_jk_c_ab__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_i_jk_c_ab__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_i_jk_c_ab__inv_h1_p23_const(PT_INF);
    diag_T3_1h2p_i_jk_c_ab__inv_h3_p23(PT_INF); // depends on m2c

    //diag_T3_1h2p_ijk_c_ab__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_ijk_c_ab__inv_h3_p23_const(PT_INF);
    //diag_T3_1h2p_ijk_c_ab__inv_h1_p23_const(PT_INF);
    //diag_T3_1h2p_ijk_c_ab__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_ijk_c_ab__inv_h3_p13_const(PT_INF);
    //diag_T3_1h2p_ijk_c_ab__inv_h2_p12_const(PT_INF);

    //diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_k_ij_c_ab__inv_h3_p13_const(PT_INF);
    //diag_T3_1h2p_k_ij_c_ab__inv_h2_p13_const(PT_INF);
    diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_m3c(PT_INF); // depends on m3c

    diag_T3_1h2p_k_ij__inv_h3_p12(PT_INF); // depends on m3c
    //diag_T3_1h2p_k_ij__inv_h2_p13_const(PT_INF);

    //diag_T3_1h2p_i_jk__inv_h1_p23_const(PT_INF);
    //diag_T3_1h2p_i_jk__inv_h2_p13_const(PT_INF);

    //diag_T3_1h2p_ijk__inv_h3_p12_const(PT_INF);
    //diag_T3_1h2p_ijk__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_ijk__inv_h1_p23_const(PT_INF);

    //diag_T3_1h2p_c_ab__inv_h3_p12_const(PT_INF);
    diag_T3_1h2p_c_ab__inv_h1_p23(PT_INF); // depends on m3c

    diag_T3_1h2p_a_bc__inv_h1_p23(PT_INF); // depends on m3c
    //diag_T3_1h2p_a_bc__inv_h3_p12_const(PT_INF);

    diag_T3_1h2p_abc__inv_h1_p23(PT_INF); // depends on m3c
    //diag_T3_1h2p_abc__inv_h2_p13_const(PT_INF);
    //diag_T3_1h2p_abc__inv_h3_p12_const(PT_INF);

    timer_stop("12-T3");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_1h2p_i_jk_c_ab__inv_h2_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hvhp", "r2", "1243");
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
        reorder("s2c", "r1", "2134"); // [1243] + interchange electrons
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphp", "r2", "4312");
        reorder("s2c", "s2c_21", "2143");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
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
        reorder("s2c", "r3", "2134"); // [1243] + interchange electrons
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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("e3c", "t5d_r1", "145326");
        reorder("pphh", "t5d_r2", "2341");
        reorder("s2c", "t5d_r3", "3214");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons: 2 <-> 3
    update("m3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_i_jk_c_ab__inv_h3_p12_const(int pt_order)
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
        reorder("e2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(12)");
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_i_jk_c_ab__inv_h1_p23_const(int pt_order)
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME && pt_order == PT_INF) {
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
        reorder("pphp", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
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
        reorder("x2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(56)");
    reorder("r4", "r5", "321654"); // interchange electrons 1 <-> 3
    update("m3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_i_jk_c_ab__inv_h3_p23(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhg", "r2", "1243");
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
        reorder("m2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphg", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
            goto finish;
        }

        // T4b
        reorder("h1c", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphg", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("m2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("x2c", "r1", "r2", 2);
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
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("t3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("m2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    reorder("r4", "r5", "321456"); // the sign will be changed
    update("m3nw", -1.0, "r5"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_1h2p_ijk_c_ab__inv_h3_p12_const(int pt_order)
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
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12)");
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk_c_ab__inv_h3_p23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "pphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("e2c", "e2c_2134", "2134"); // the sign will be changed
        reorder("e2c_2134", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        mult("s1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("e2c", "e2c_2134", "2134"); // the sign will be changed
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23)");
    reorder("r4", "r5", "321456"); // the sign will be changed
    update("m3_0", -1.0, "r5"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk_c_ab__inv_h1_p23_const(int pt_order)
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
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23|56)");
    reorder("r4", "r5", "321654"); // interchange electrons 1 <-> 3
    update("m3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk_c_ab__inv_h2_p13_const(int pt_order)
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
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465");
    update("m3_0", 1.0, "r5"); // interchange electrons 2 <-> 3
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk_c_ab__inv_h3_p13_const(int pt_order)
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
        reorder("e2c", "e2c_2134", "2134"); // the sign will be changed
        reorder("e2c_2134", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
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
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13)");
    reorder("r4", "r5", "132456"); // the sign will be changed
    update("m3_0", -1.0, "r5"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk_c_ab__inv_h2_p12_const(int pt_order)
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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12|46)");
    reorder("r4", "r5", "123465"); // the sign will be changed
    update("m3_0", -1.0, "r5"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_1h2p_k_ij_a_bc__inv_h1_p23_const(int pt_order)
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
        reorder("t1c", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("e3c", "t5e_r1", "215643"); // interchange electrons [231] + reorder [134562]
        reorder("pphh", "t5e_r2", "4123");
        reorder("s2c", "s2c_21", "2143");
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

    perm("r3", "(23)");
    reorder("r3", "r4", "321654"); // interchange electrons 1 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_a_bc__inv_h3_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phpg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3e
        reorder("h2c", "r1", "4321"); // [3412] + interchange electrons
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("pppgr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "r1", "4321"); // [3412] + interchange electrons
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("x3c", "t5e_r1", "124653"); // interchange electrons 2 <-> 3 + reorder [134562]
        reorder("pphh", "t5e_r2", "4123");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(45)");
    update("m3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_a_bc__inv_h2_p13_const(int pt_order)
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
        reorder("t1c", "r0", "21");
        reorder("h1c", "r1", "21");
        mult("r0", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
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
        reorder("t1c", "r0", "21");
        reorder("h1c", "r1", "21");
        mult("r0", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
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
        reorder("e2c", "e2c_21", "2143");
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

    perm("r3", "(13|46)");
    reorder("r3", "r4", "132465"); // interchange electrons: 2 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_a_bc__inv_h1_p12(int pt_order)
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
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("m3c", "t5e_r1", "126453"); // the sign will be changed
        reorder("pphh", "t5e_r2", "4123");
        reorder("t2c", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", 0.5, "t5e_r6"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    reorder("r3", "r4", "123654"); // the sign will be changed
    update("m3nw", -1.0, "r4"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_1h2p_k_ij_abc__inv_h3_p12_const(int pt_order)
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
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(45)");
    update("m3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_abc__inv_h2_p13_const(int pt_order)
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
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|46)");
    reorder("r3", "r4", "132465"); // interchange electrons 2 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_abc__inv_h2_p12_const(int pt_order)
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
    mult("x2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(46)");
    reorder("r3", "r4", "123465"); // the sign will be changed
    update("m3_0", -1.0, "r4"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_abc__inv_h1_p23_const(int pt_order)
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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(23|56)");
    reorder("r3", "r4", "321654"); // interchange electrons 1 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_abc__inv_h1_p12(int pt_order)
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
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(56)");
    reorder("r3", "r4", "123654"); // the sign will be changed
    update("m3nw", -1.0, "r4"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_abc__inv_h3_p23_const(int pt_order)
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
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5"); // the sign was changed, +1 -> -1
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
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c", "s2c_21", "2143");
    mult("s2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(23|45)");
    reorder("r3", "r4", "321456"); // the sign will be changed, +1 -> -1
    update("m3_0", -1.0, "r4"); // the sign was changed
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_1h2p_c_ab__inv_h3_p12_const(int pt_order)
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
        reorder("pg", "pgr_", "21");
        update("i1", 1.0, "pgr_");
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
    mult("x3c", "i1", "r1", 1);
    update("m3_0", 1.0, "r1");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_c_ab__inv_h1_p23(int pt_order)
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
    reorder("m3c", "r0", "321654"); // 1 <-> 3
    mult("r0", "i1", "r1", 1);
    perm("r1", "(56)");
    reorder("r1", "r2", "321654"); // interchange electrons 1 <-> 3
    update("m3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_1h2p_a_bc__inv_h1_p23(int pt_order)
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
    reorder("m3c", "r0", "321654"); // 1 <-> 3
    mult("r0", "i1", "r1", 2);
    reorder("r1", "r2", "321654"); // interchange electrons 1 <-> 3
    update("m3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_a_bc__inv_h3_p12_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "phpp", "0100", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2c
        update("i1", 0.5, "pppgr");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5g
        reorder("h2c", "r1", "4321"); // [3412] + interchange electrons
        mult("r1", "pphh", "r2", 2);
        update("i1", 0.25, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9d
        reorder("h1c", "r1", "21");
        reorder("t1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x3c", "i1", "r1", 2);
    perm("r1", "(45)");
    update("m3_0", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_1h2p_abc__inv_h1_p23(int pt_order)
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
    reorder("m3c", "r0", "321654");
    mult("r0", "i1", "r1", 2);
    perm("r1", "(56)");
    reorder("r1", "r2", "321654"); // interchange electrons: 1 <-> 3
    update("m3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_abc__inv_h2_p13_const(int pt_order)
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
    reorder("x3c", "r0", "132465");
    mult("r0", "i1", "r1", 2);
    diagram_stack_erase("r0");
    perm("r1", "(46)");
    reorder("r1", "r2", "132465"); // interchange electrons 2 <-> 3
    update("m3_0", 1.0, "r2");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_abc__inv_h3_p12_const(int pt_order)
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
    mult("x3c", "i1", "r1", 2);
    perm("r1", "(45)");
    update("m3_0", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_const(int pt_order)
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
        reorder("phhg", "r0", "2431");
        update("i1", 1.0, "r0");

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
    reorder("x3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_c_ab__inv_h3_p13_const(int pt_order)
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
        reorder("e2c", "e2c_2134", "2134"); // the sign will be changed
        reorder("e2c_2134", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", -1.0, "r4"); // the sign was changed, +1 -> -1
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
    reorder("s3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(13)");
    reorder("r4", "r5", "132456"); // the sign will be changed
    update("m3_0", -1.0, "r5"); // the sign was changed, +1 -> -1
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_c_ab__inv_h2_p13_const(int pt_order)
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
        reorder("pvhp", "r0", "2431");
        update("i1", 1.0, "r0");

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
    reorder("e3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons 2 <-> 3
    update("m3_0", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij_c_ab__inv_h3_p12_m3c(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i0", "phph", "0000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        update("i0", -1.0, "phph");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5a
        reorder("t2c", "r1", "2314");
        reorder("pphh", "r2", "4132");
        mult("r2", "r1", "r3", 2);
        reorder("r3", "r4", "2341");
        update("i0", 1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("ppph", "r1", "3412");
        mult("r1", "t1c", "r2", 1);
        reorder("r2", "r3", "3412");
        update("i0", -1.0, "r3");
        restore_stack_pos(pos2);

        // T6d
        reorder("phhh", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "3412");
        update("i0", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9e
        reorder("t1c", "r0", "21");
        reorder("pphh", "r1", "1432");
        mult("t1c", "r1", "r2", 1);
        mult("r0", "r2", "r3", 1);
        reorder("r3", "r4", "3214");
        update("i0", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("i0", "i1", "2341");
    reorder("m3c", "i2", "124635");
    mult("i2", "i1", "i3", 2);
    reorder("i3", "i4", "125364");
    perm("i4", "(45)");
    update("m3nw", 1.0, "i4");

    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_1h2p_ijk__inv_h3_p12_const(int pt_order)
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

    reorder("e3c", "r1", "465132");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(12)");
    update("m3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk__inv_h2_p13_const(int pt_order)
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

    reorder("e3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(13)");
    reorder("r3", "r4", "132465"); // interchange electrons 2 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_ijk__inv_h1_p23_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
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

    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(23)");
    reorder("r3", "r4", "321654"); // interchange electrons 1 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_1h2p_i_jk__inv_h1_p23_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
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
    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    reorder("r3", "r4", "321654"); // interchange electrons 1 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_i_jk__inv_h2_p13_const(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hphh", "0100", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2d
        update("i1", 0.5, "hvhh");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5f
        reorder("pphh", "r2", "3412");
        reorder("s2c", "s2c_21", "2143");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", 0.25, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9c
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
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
    perm("r3", "(13)");
    reorder("r3", "r4", "132465"); // interchange electrons 2 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_1h2p_k_ij__inv_h3_p12(int pt_order)
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
    reorder("m3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    update("m3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_1h2p_k_ij__inv_h2_p13_const(int pt_order)
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
    reorder("e3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(13)");
    reorder("r3", "r4", "132465"); // interchange electrons 2 <-> 3
    update("m3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Folded diagrams for triples equations
 */
void ampl_equations_1h2p_folded_triples_F1a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("m3c", "r1", "234561");
    mult("veff01", "r1", "r3", 1);
    diagram_stack_erase("r1");
    perm("r3", "(12)");
    update("m3nw", -1.0, "r3");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F1b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff10", "r1", "21");
    mult("m3c", "r1", "r2", 1);
    update("m3nw", 1.0, "r2");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F2a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "r1", "2341");
    mult("veff02", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    tmplt("r4", "pphpph", "110001", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r4");
    diagram_stack_erase("r3");

    perm("r4", "(45)");
    update("m3nw", -1.0, "r4");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F2b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "e2c_21", "2143");
    reorder("e2c_21", "r1", "1243");
    reorder("veff11", "r2", "3412");

    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    tmplt("r5", "phpphp", "101010", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r4", "r5");
    diagram_stack_erase("r4");

    perm("r5", "(13|46)");
    reorder("r5", "r6", "132465"); // interchange electrons 2 <-> 3
    update("m3nw", 1.0, "r6");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F2c()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("x2c", "r1", "2341");
    reorder("veff11", "veff11_21", "2143");
    mult("veff11_21", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    tmplt("r4", "hpphpp", "011100", "213546", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r4");
    diagram_stack_erase("r3");

    perm("r4", "(23)");
    reorder("r4", "r5", "321654"); // interchange electrons 1 <-> 3
    update("m3nw", -1.0, "r5");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F3a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("m3c", "r1", "345612");
    mult("veff02", "r1", "r2", 2);
    update("m3nw", -0.5, "r2");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F3b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("m3c", "r1", "134562");
    reorder("veff11", "r2", "1423");
    mult("r1", "r2", "r3", 2);
    diagram_stack_erase("r1");

    reorder("r3", "r4", "152346");
    diagram_stack_erase("r3");

    perm("r4", "(12)");
    update("m3nw", 1.0, "r4");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F4a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff02", "r1", "1243");
    reorder("e2c", "r2", "2341");
    reorder("s1c", "r3", "21");

    mult("r3", "r1", "r4", 1);
    mult("r4", "r2", "r5", 1);
    reorder("r5", "r6", "234156");
    diagram_stack_erase("r5");

    perm("r6", "(45)");
    update("m3nw", -1.0, "r6");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F4b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("e2c", "r1", "2134"); // interchange electrons + [1243]
    reorder("veff11", "r2", "1423");
    reorder("s1c", "r3", "21");

    mult("r3", "r2", "r4", 1);
    mult("r4", "r1", "r5", 1);
    reorder("r5", "r6", "245136");
    diagram_stack_erase("r5");

    perm("r6", "(13|46)");
    reorder("r6", "r7", "132465"); // interchange electrons 2 <-> 3
    update("m3nw", 1.0, "r7");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F4c()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("x2c", "r1", "2341");
    reorder("veff11", "veff11_21", "2143");
    reorder("veff11_21", "r2", "2341");

    mult("h1c", "r2", "r3", 1);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "124356");
    diagram_stack_erase("r4");

    perm("r5", "(23)");
    reorder("r5", "r6", "321654");
    update("m3nw", 1.0, "r6");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F5a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff02", "r1", "1243");
    reorder("e1c", "r2", "21");
    reorder("s2c", "r3", "2341");

    mult("r2", "r1", "r4", 1);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "234156");
    diagram_stack_erase("r5");

    reorder("r6", "r7", "123654"); // the sign will be changed
    update("m3nw", 1.0, "r7"); // the sign was changed, -1 -> +1

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F5b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("s2c", "r1", "2341");
    reorder("veff11", "veff11_21", "2143");
    reorder("veff11_21", "r2", "2341");

    mult("e1c", "r2", "r3", 1);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "124356");
    diagram_stack_erase("r4");

    perm("r5", "(12)");
    reorder("r5", "r6", "123654"); // the sign will be changed
    update("m3nw", -1.0, "r6"); // the sign was changed, +1 -> -1

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F6a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "561234");
    reorder("s1c", "r2", "21");
    
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "456123");
    diagram_stack_erase("r3");

    tmplt("r5", "pphpph", "110001", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r4", "r5");
    diagram_stack_erase("r4");
    perm("r5", "(45)");
    update("m3nw", -1.0, "r5");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F6b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "456123");
    mult("r1", "h1c", "r2", 1);
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");

    tmplt("r4", "pphpph", "110001", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r4");
    update("m3nw", 1.0, "r4");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F7a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "612345");
    reorder("x2c", "r2", "3412");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "456123");
    diagram_stack_erase("r3");

    tmplt("r5", "pphpph", "110001", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r4", "r5");
    update("m3nw", -0.5, "r5");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F7b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "124635");
    reorder("e2c", "r2", "2341");

    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "125364");
    diagram_stack_erase("r3");

    tmplt("r5", "pphpph", "110001", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r4", "r5");
    perm("r5", "(45)");
    update("m3nw", 1.0, "r5");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F8()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "126345");
    reorder("m3c", "r2", "345612");
    mult("r1", "r2", "r3", 3);
    diagram_stack_erase("r1");
    diagram_stack_erase("r2");
    reorder("r3", "r4", "124563");
    update("m3nw", 0.5, "r4");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F9a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "612345");
    reorder("s1c", "r2", "21");

    mult("r2", "r1", "r3", 1);
    diagram_stack_erase("r1");
    mult("r2", "r3", "r4", 1);
    diagram_stack_erase("r3");

    reorder("r4", "r5", "456123");
    diagram_stack_erase("r4");
    tmplt("r6", "pphpph", "110001", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r5", "r6");
    update("m3nw", -1.0, "r6");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F9b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "124653");
    reorder("s1c", "r2", "21");

    mult("h1c", "r1", "r3", 1);
    diagram_stack_erase("r1");
    mult("r3", "r2", "r4", 1);
    diagram_stack_erase("r3");
    reorder("r4", "r5", "231465");
    diagram_stack_erase("r4");

    tmplt("r6", "pphpph", "110001", "231564", NOT_PERM_UNIQUE);
    expand_diagram("r5", "r6");
    diagram_stack_erase("r5");
    perm("r6", "(45)");
    update("m3nw", 1.0, "r6");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F10a()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "126453");
    reorder("x2c", "r2", "3412");

    mult("h1c", "r1", "r3", 1);
    diagram_stack_erase("r1");
    mult("r3", "r2", "r4", 2);
    diagram_stack_erase("r3");

    reorder("r4", "r5", "231564");
    update("m3nw", 0.5, "r5");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F10b()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "126435");
    reorder("e2c", "r2", "2341");
    reorder("s1c", "r5", "21");

    mult("r1", "r2", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125634");
    diagram_stack_erase("r3");
    mult("r4", "r5", "r6", 1);
    diagram_stack_erase("r4");
    reorder("r6", "r7", "123645");
    diagram_stack_erase("r6");

    perm("r7", "(45)");
    update("m3nw", 1.0, "r7");

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F11()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "126453");
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r2", "1342");

    mult("r1", "e1c", "r3", 2);
    diagram_stack_erase("r1");
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "145236");
    diagram_stack_erase("r4");

    reorder("r5", "r6", "321456"); // the sign will be changed
    update("m3nw", -1.0, "r6"); // the sign was changed, +1 -> -1

    restore_stack_pos(pos);
}


void ampl_equations_1h2p_folded_triples_F12()
{
    dg_stack_pos_t pos = get_stack_pos();

    reorder("veff12", "r1", "456123");
    reorder("s1c", "r2", "21");

    mult("r1", "h1c", "r3", 1);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "456312");
    diagram_stack_erase("r3");
    mult("r2", "r4", "r5", 1);
    diagram_stack_erase("r4");
    mult("r2", "r5", "r6", 1);
    diagram_stack_erase("r5");

    reorder("r6", "r7", "345126");
    update("m3nw", 1.0, "r7");

    restore_stack_pos(pos);
}




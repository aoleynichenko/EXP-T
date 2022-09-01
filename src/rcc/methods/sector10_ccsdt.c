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

/*
 * Multireference Fock-space coupled cluster method for model spaces
 * with one hole (Fock space sector 1h0p)
 *
 * Symbolic names used for diagrams:
 * h1c   current approximation to T^{1h0p}_1 amplitudes
 * h2c   current approximation to T^{1h0p}_2 amplitudes
 * h3c   current approximation to T^{1h0p}_3 amplitudes
 * h1nw  next (new) approximation to T^{1h0p}_1 amplitudes
 * h2nw  next (new) approximation to T^{1h0p}_2 amplitudes
 * h3nw  next (new) approximation to T^{1h0p}_3 amplitudes
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip particle creation lines to turn them to valence hole
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 *
 * 2019-2022 Alexander Oleynichenko
 */

#include "methods.h"

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
#include "heff.h"
#include "linalg.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void init_amplitudes_1h0p();

void const_terms_1h0p();

void construct_singles_1h0p();

void construct_doubles_1h0p();

void construct_triples_1h0p();

void construct_folded_1h0p();


void t3_contrib_to_singles_1h0p()
{
    timer_new_entry("10-T1-S7", "1h0p -- Triples contibution to singles");
    timer_start("10-T1-S7");

    // S7
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h3c", "r1", "145623");
    mult("r1", "pphh", "r2", 4);
    update("h1nw", 0.25, "r2");
    restore_stack_pos(pos);

    timer_stop("10-T1-S7");
}


void t3_contrib_to_doubles_1h0p()
{
    timer_new_entry("10-T2-Triples", "1h0p -- Triples contibution to doubles");
    timer_start("10-T2-Triples");

    dg_stack_pos_t pos = get_stack_pos();

    // D10a
    reorder("h3c", "r1", "124536");
    reorder("ph", "phr", "21");
    mult("r1", "phr", "r2", 2);
    update("h2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // D10b-1
    reorder("ppph", "r1", "3412");
    reorder("h3c", "r2", "124356");
    mult("r2", "r1", "r3", 3);
    update("h2nw", 0.5, "r3");
    restore_stack_pos(pos);

    // D10b-2 -- constant term
    reorder("ppgh", "r1", "3412");
    reorder("t3c", "r2", "124356");
    mult("r2", "r1", "r3", 3);
    reorder("r3", "r4", "2143");
    update("h2nw", 0.5, "r4");
    restore_stack_pos(pos);

    // D10c
    reorder("h3c", "r1", "145623");
    mult("r1", "hphh", "r2", 3);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update("h2nw", -0.5, "r3");
    restore_stack_pos(pos);

    if (cc_opts->cc_model <= CC_MODEL_CCSDT_1A) {
        timer_stop("10-T2-Triples");
        timer_stop("10-T2");
        return;
    }

    // D11a
    reorder("pphh", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    reorder("h3c", "r3", "124536");
    mult("r3", "r2", "r4", 2);
    update("h2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // D11b-1
    reorder("t1c", "t1r_", "21");
    mult("t1r_", "pphh", "r1", 1);
    reorder("h3c", "r2", "134562");
    mult("r2", "r1", "r3", 3);
    update("h2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // D11b-2
    reorder("h1c", "h1cr_", "21");
    mult("h1cr_", "pphh", "r1", 1);
    reorder("t3c", "r2", "134562");
    mult("r2", "r1", "r3", 3);
    reorder("r3", "r4", "2143");  // interchange electrons 1 <-> 2
    update("h2nw", -0.5, "r4");
    restore_stack_pos(pos);

    // D11c
    reorder("pphh", "r1", "3412");
    mult("t1c", "r1", "r2", 1);
    reorder("h3c", "r3", "146235");
    mult("r3", "r2", "r4", 3);
    reorder("r4", "r5", "1423");
    perm("r5", "(12)");
    update("h2nw", -0.5, "r5");
    restore_stack_pos(pos);

    timer_stop("10-T2-Triples");
}


void t3_contrib_to_folded_1h0p()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h3c", "r1", "123564");
    reorder("veff10", "r2", "21");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "123645");
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_k_ij_a_bc_inv1();
void diag_T3_1h0p_k_ij_a_bc_inv2();
void diag_T3_1h0p_i_jk_c_ab_inv1();
void diag_T3_1h0p_i_jk_c_ab_inv3();
void diag_T3_1h0p_k_ij_c_ab_inv1();
void diag_T3_1h0p_k_ij_c_ab_inv3();
void diag_T3_1h0p_k_ij_abc_inv1();
void diag_T3_1h0p_k_ij_abc_inv2();
void diag_T3_1h0p_k_ij_abc_inv3();
void diag_T3_1h0p_ijk_c_ab_inv1();
void diag_T3_1h0p_ijk_c_ab_inv3();

void diag_T3_1h0p_a_bc_inv1();
void diag_T3_1h0p_a_bc_inv2();
void diag_T3_1h0p_c_ab_inv1();
void diag_T3_1h0p_c_ab_inv3();
void diag_T3_1h0p_abc_inv1();
void diag_T3_1h0p_abc_inv2();
void diag_T3_1h0p_abc_inv3();

void diag_T3_1h0p_i_jk();
void diag_T3_1h0p_k_ij();
void diag_T3_1h0p_ijk();


/**
 * Triples amplitude equations.
 * This code is invoked only if cc_model >= CCSDT-1a
 */
void construct_triples_1h0p()
{
    timer_new_entry("10-T3", "1h0p -- Triples equations (T3)");
    timer_start("10-T3");

    tmplt("h3nw", "hhhhpp", "000100", "123456", IS_PERM_UNIQUE);

    diag_T3_1h0p_k_ij_a_bc_inv1();
    diag_T3_1h0p_k_ij_a_bc_inv2();
    diag_T3_1h0p_i_jk_c_ab_inv1();
    diag_T3_1h0p_i_jk_c_ab_inv3();
    diag_T3_1h0p_k_ij_c_ab_inv1();
    diag_T3_1h0p_k_ij_c_ab_inv3();
    diag_T3_1h0p_k_ij_abc_inv1();
    diag_T3_1h0p_k_ij_abc_inv2();
    diag_T3_1h0p_k_ij_abc_inv3();
    diag_T3_1h0p_ijk_c_ab_inv1();
    diag_T3_1h0p_ijk_c_ab_inv3();

    diag_T3_1h0p_a_bc_inv1();
    diag_T3_1h0p_a_bc_inv2();
    diag_T3_1h0p_c_ab_inv1();
    diag_T3_1h0p_c_ab_inv3();
    diag_T3_1h0p_abc_inv1();
    diag_T3_1h0p_abc_inv2();
    diag_T3_1h0p_abc_inv3();

    diag_T3_1h0p_i_jk();
    diag_T3_1h0p_k_ij();
    diag_T3_1h0p_ijk();

    timer_stop("10-T3");
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_1h0p_k_ij_a_bc_inv1()
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b

        // T3e
        reorder("t2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppppr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T4d

        // T7a

        // T7d
        reorder("t1c", "r1", "21");
        mult("r1", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8c

        // T8d
        reorder("t2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if (cc_opts->cc_model > CC_MODEL_CCSDT_3) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("h3c", "t5e_r1", "134562");
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

    perm("r3", "(3/12)");
    update("h3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_k_ij_a_bc_inv2()
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b

        // T3e
        reorder("h2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        mult("ppgp", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T4d

        // T7a

        // T7d
        reorder("t1c", "r1", "21");
        reorder("h1c", "h1cr", "21");
        mult("r1", "phhh", "r2", 1);
        mult("h1cr", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8c

        // T8d
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c", "r1", "21");
        reorder("h1c", "h1cr", "21");
        mult("r1", "pphh", "r2", 1);
        mult("h1cr", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if (cc_opts->cc_model > CC_MODEL_CCSDT_3) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("t3c", "t5e_r1", "134562");
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

    perm("r3", "(3/12|46)");
    reorder("r3", "r4", "213546"); // interchange 1 <-> 2
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_1h0p_i_jk_c_ab_inv1()
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3a
        reorder("ph", "phr_", "21");
        reorder("t2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3c

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphp", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("t1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T4c

        // T7b

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

        // T8b

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("t2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if (cc_opts->cc_model > CC_MODEL_CCSDT_3) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("h3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("t2c", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(1/23|56)");
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_i_jk_c_ab_inv3()
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3a
        reorder("ph", "phr_", "21");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3c

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphg", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4b
        reorder("h1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T4c

        // T7b

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

        // T8b

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("t2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

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
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if (cc_opts->cc_model > CC_MODEL_CCSDT_3) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("t3c", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "t5d_r3", "4123");
        mult("t5d_r1", "t5d_r2", "t5d_r4", 3);
        mult("t5d_r4", "t5d_r3", "t5d_r5", 1);
        reorder("t5d_r5", "t5d_r6", "156234");
        update("r4", -0.5, "t5d_r6");
        restore_stack_pos(pos2);
    }

    perm("r4", "(1/23)");
    reorder("r4", "r5", "321654");   // interchange 1 <-> 3
    update("h3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_1h0p_k_ij_c_ab_inv1()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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

        // T9e
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3124");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }
    reorder("h3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(3/12|56)");
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_k_ij_c_ab_inv3()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "hhhp", "0100", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        reorder("phhg", "r0", "2431");
        update("i1", 1.0, "r0");

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

        // T9e
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3124");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }
    reorder("t3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(3/12)");
    reorder("r4", "r5", "321654");  // interchange 1 <-> 3
    update("h3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_1h0p_k_ij_abc_inv1()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hppp", "0000", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("phhp", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|56)");
    update("h3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_k_ij_abc_inv2()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhpp", "0100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("phhp", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|46)");   // interchange 1 <-> 2
    reorder("r3", "r4", "213546");
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_k_ij_abc_inv3()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
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
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|45)");
    reorder("r3", "r4", "321654");  // interchange 1 <-> 3
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_1h0p_ijk_c_ab_inv1()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhph", "0000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        mult("t1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(123|56)");
    update("h3nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_ijk_c_ab_inv3()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
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

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("phhg", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
    reorder("t2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(123)");
    reorder("r4", "r5", "321654");  // interchange 1 <-> 3
    update("h3nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_1h0p_a_bc_inv1()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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

        // T9d
        reorder("t1c", "r1", "21");
        reorder("t1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);

        // T6g
    }
    mult("h3c", "i1", "r1", 2);
    update("h3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_a_bc_inv2()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "hppp", "1000", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2c
        update("i1", 0.5, "ppgp");

        // T5g
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        update("i1", 0.25, "r2");
        restore_stack_pos(pos2);

        // T9d
        reorder("t1c", "r1", "21");
        reorder("h1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);

        // T6g
    }
    mult("t3c", "i1", "r1", 2);
    perm("r1", "(46)");
    reorder("r1", "r2", "213546");  // interchange 1 <-> 2
    update("h3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_1h0p_c_ab_inv1()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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

        // T9b
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        mult("r1", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    mult("h3c", "i1", "r1", 1);
    perm("r1", "(56)");
    update("h3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_c_ab_inv3()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
        return;
    }

    tmplt("i1", "hp", "10", "21", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2a
        reorder("pg", "pgr_", "21");
        update("i1", 1.0, "pgr_");
        restore_stack_pos(pos2);

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

        // T9b
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        mult("r1", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    mult("t3c", "i1", "r1", 1);
    reorder("r1", "r2", "321654");  // interchange 1 <-> 3
    update("h3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_1h0p_abc_inv1()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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
    mult("h3c", "i1", "r1", 2);
    perm("r1", "(56)");
    update("h3nw", 1.0, "r1");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_abc_inv2()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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
    mult("t3c", "i1", "r1", 2);
    perm("r1", "(46)");
    reorder("r1", "r2", "213546");   // interchange 1 <-> 2
    update("h3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


void diag_T3_1h0p_abc_inv3()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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
    mult("t3c", "i1", "r1", 2);
    perm("r1", "(45)");
    reorder("r1", "r2", "321654");  // interchange 1 <-> 3
    update("h3nw", 1.0, "r2");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_1h0p_ijk()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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

    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(123)");
    update("h3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_1h0p_i_jk()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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

        // T9c
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);

        // T6h
    }

    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(1/23)");
    update("h3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_1h0p_k_ij()
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT) {
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

        // T9a
        reorder("pphh", "r2", "3142");
        mult("r2", "t1c", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    reorder("h3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(3/12)");
    update("h3nw", 1.0, "r3");
    restore_stack_pos(pos);
}



/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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
 * Fock space sector 2h1p. CCSDT model
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


void diag_T3_2h1p_i_jk_c_ab__inv_h12_p3_const();
void diag_T3_2h1p_i_jk_c_ab__inv_h23_p1_const();
void diag_T3_2h1p_i_jk_c_ab__inv_h13_p2_const();
void diag_T3_2h1p_i_jk_c_ab__inv_h12_p1_w2c();

void diag_T3_2h1p_ijk_c_ab__inv_h23_p1_const();
void diag_T3_2h1p_ijk_c_ab__inv_h13_p2_const();
void diag_T3_2h1p_ijk_c_ab__inv_h12_p2_const();
void diag_T3_2h1p_ijk_c_ab__inv_h13_p3_const();
void diag_T3_2h1p_ijk_c_ab__inv_h12_p3_const();
void diag_T3_2h1p_ijk_c_ab__inv_h12_p1_w2c();

void diag_T3_2h1p_k_ij_a_bc__inv_h12_p3_const();
void diag_T3_2h1p_k_ij_a_bc__inv_h23_p1_const();
void diag_T3_2h1p_k_ij_a_bc__inv_h13_p2_const();
void diag_T3_2h1p_k_ij_a_bc__inv_h23_p3_const();

void diag_T3_2h1p_k_ij_abc__inv_h12_p2_const();
void diag_T3_2h1p_k_ij_abc__inv_h23_p1_const();
void diag_T3_2h1p_k_ij_abc__inv_h13_p2_const();
void diag_T3_2h1p_k_ij_abc__inv_h12_p3_const();
void diag_T3_2h1p_k_ij_abc__inv_h23_p3_const();
void diag_T3_2h1p_k_ij_abc__inv_h13_p3_const();


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_2h1p_i_jk_c_ab__inv_h12_p3_const()
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

        // T3a
        reorder("ph", "phr_", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("s2c", "s2c_21", "2143");
        reorder("pphp", "r2", "4312");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c", "r2", "21");
        mult("hvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

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
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    perm("r4", "(12)");
    closed("r4", "r5");
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_i_jk_c_ab__inv_h23_p1_const()
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

        // T4b
        reorder("h1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    perm("r4", "(56)");
    closed("r4", "r5");
    reorder("r5", "r6", "321654");
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_i_jk_c_ab__inv_h13_p2_const()
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

        // T4b
        reorder("h1c", "r2", "21");
        mult("vhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    closed("r4", "r5");
    perm("r5", "(13|46)");
    reorder("r5", "r6", "132465"); // interchange electrons 2 <-> 3
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


// contribution to T3 dependent on w2c (ampl-s in 2h1p)
void diag_T3_2h1p_i_jk_c_ab__inv_h12_p1_w2c()
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

        // T4b
        reorder("t1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

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
    reorder("w2c", "r1", "4321"); // [3412] + interchange electrons 1 <-> 2
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    closed("r4", "r5");
    reorder("r5", "r6", "321456"); // and the sign must be changed, +1 -> -1
    update("veff21", -1.0, "r6");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_2h1p_ijk_c_ab__inv_h23_p1_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
        reorder("h2c", "h2c_21", "2143");
        reorder("pphh", "r1", "3142");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(23|56)");
    reorder("r5", "r6", "321654"); // interchange electrons 1 <-> 3
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_ijk_c_ab__inv_h13_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("phhg", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(13|46)");
    reorder("r5", "r6", "132465"); // interchange electrons 2 <-> 3
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_ijk_c_ab__inv_h12_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phph", "1000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("s1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(13)");
    reorder("r5", "r6", "132456"); // + sign will be changed!
    update("veff21_const", -1.0, "r6");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_ijk_c_ab__inv_h13_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphh", "0110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("t1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(12|46)");
    reorder("r5", "r6", "123465"); // the sign will be changed!
    update("veff21_const", -1.0, "r6");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_ijk_c_ab__inv_h12_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pvhp", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(12)");
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


// TODO: contribution to T3 dependent on w2c (ampl-s in 2h1p)
void diag_T3_2h1p_ijk_c_ab__inv_h12_p1_w2c()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhph", "0000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
    reorder("w2c", "r1", "4321"); // [3412] + interchange electrons 1 <-> 2
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(23)");
    reorder("r5", "r6", "321456"); // the sign will be changed, +1 -> -1
    update("veff21", -1.0, "r6");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_2h1p_k_ij_a_bc__inv_h12_p3_const()
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

        // T3e
        reorder("h2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        mult("ppgp", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

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
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    perm("r4", "(45)");
    update("veff21_const", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_a_bc__inv_h23_p1_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhp", "0110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phgg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("g2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppgg", "ppggr", "3412");
        mult("ppggr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r1", "21");
        mult("r1", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("g2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "r1", "21");
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
    closed("r3", "r4");
    perm("r4", "(23)");
    reorder("r4", "r5", "321654"); // interchange electrons: 1 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_a_bc__inv_h13_p2_const()
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

        // T3e
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("pppg", "pppgr", "3412");
        mult("pppgr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons: 2 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_a_bc__inv_h23_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvgg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("g2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppgg", "ppggr", "3412");
        mult("ppggr", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("g2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    reorder("r4", "r5", "123654"); // the sign of the prefactor will be changed, +1 -> -1
    update("veff21_const", -1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_2h1p_k_ij_abc__inv_h12_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhpp", "0100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("t2c", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(13|45)");
    reorder("r4", "r5", "132456"); // the sign will be changed
    update("veff21_const", -1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_abc__inv_h23_p1_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhp", "0110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("phhg", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
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
    closed("r3", "r4");
    perm("r4", "(23|56)");
    reorder("r4", "r5", "321654"); // interchange electrons: 1 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_abc__inv_h13_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons: 2 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_abc__inv_h12_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhp", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(45)");
    update("veff21_const", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_abc__inv_h23_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5");  // the sign was changed, +1 -> -1
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7");  // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(56)");
    reorder("r4", "r5", "123654"); // the sign will be changed, +1 -> -1
    update("veff21_const", -1.0, "r5");
    restore_stack_pos(pos);
}


void diag_T3_2h1p_k_ij_abc__inv_h13_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphp", "1010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5"); // the sign was changed, +1 -> -1
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(46)");
    reorder("r4", "r5", "123465"); // the sign will be changed, +1 -> -1
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


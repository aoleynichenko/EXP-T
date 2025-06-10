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


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cc_properies.h"
#include "engine.h"
#include "../heff/mvcoef.h"
#include "spinors.h"
#include "symmetry.h"
#include "ms_prop.h"
#include "sort.h"
#include "linalg.h"


void direct_property_0h0p_1h1p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "ph", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c+", "t1c+_v", "10", 0);
    restrict_valence("t1c+", "t1c+_g", "01", 0);
    restrict_valence("t1c+", "t1c+_vg", "11", 0);
    restrict_valence("t1c", "t1c_gv", "11", 0);

    restrict_valence("t2c+", "t2c+_vg", "1010", 0);
    restrict_valence("t2c+", "t2c+_g", "0010", 0);
    restrict_valence("t2c+", "t2c+_v", "1000", 0);

    restrict_valence("prop_ph", "prop_vg", "11", 0);
    restrict_valence("prop_ph", "prop_pg", "01", 0);
    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hh", "prop_hg", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);

    dg_stack_pos_t pos = get_stack_pos();

    // O1
    update(result_name, 1.0, "prop_vg");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    reorder("prop_pg", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L1b
    reorder("h1c", "r1", "21");
    mult("prop_vh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2a
    reorder("prop_hg", "r1", "21");
    mult("t1c+_v", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2b
    reorder("t1c+_g", "r1", "21");
    mult("prop_vp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3a
    reorder("t2c+_vg", "r1", "1342");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3b
    reorder("e2c", "r1", "1432");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a - disconnected term
    if (disconnected) {
        double factor_Q1a = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q1a, "t1c+_vg");
    }

    // Q1b
    reorder("t1c+_g", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q2a - disconnected term
    if (disconnected) {
        double factor_Q2a = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q2a, "e1c");
    }

    // Q2b
    reorder("h1c", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3a - disconnected term
    if (disconnected) {
        double t1_t1 = scalar_product("C", "N", "t1c", "t1c");
        update(result_name, t1_t1, "prop_vg");
    }

    // Q3b - disconnected term
    if (disconnected) {
        double factor_Q3b = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q3b, "t1c+_vg");
    }

    // Q3c - disconnected term
    if (disconnected) {
        double factor_Q3c = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q3c, "e1c");
    }

    // Q3d
    reorder("prop_pg", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3e
    reorder("t1c+_g", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3f
    reorder("prop_hg", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3g
    reorder("h1c", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3h
    reorder("t1c+_g", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q3i
    reorder("h1c", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4a
    reorder("s2c", "r1", "1342");
    reorder("prop_pg", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4b
    reorder("h2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4c
    reorder("s2c", "r1", "1342");
    reorder("t1c+_g", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4d
    reorder("h2c", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4e
    reorder("e2c", "r1", "1432");
    mult("r1", "prop_hh", "r2", 1);
    mult("r2", "t1c+", "r3", 2);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4f
    reorder("e2c", "r1", "1432");
    mult("r1", "t1c+", "r2", 1);
    mult("r2", "prop_pp", "r3", 2);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5a
    reorder("t2c+_g", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q5b
    reorder("prop_hg", "r1", "21");
    reorder("t2c+_v", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5c
    reorder("t2c+_g", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q5d
    reorder("t2c+_v", "r1", "1342");
    reorder("h1c", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5e
    reorder("t1c", "r1", "21");
    reorder("t2c+_vg", "r2", "1342");
    mult("prop_hh", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5f
    reorder("t2c+_vg", "r1", "1324");
    reorder("prop_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6a - disconnected term
    if (disconnected) {
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        update(result_name, t2_t2, "prop_vg");
        restore_stack_pos(pos);
    }

    // Q6b
    reorder("t2c+_vg", "r1", "1324");
    reorder("t2c", "r2", "3142");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6c
    reorder("t2c", "r1", "2341");
    reorder("t2c+_g", "r0", "2143");
    reorder("r0", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("prop_vh", "r3", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6d
    reorder("prop_pg", "r0", "21");
    reorder("t2c", "r1", "4123");
    mult("r0", "r1", "r2", 1);
    mult("t2c+_v", "r2", "r3", 3);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6e
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("prop_hg", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q6f
    reorder("h2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6g
    reorder("t2c+_g", "r1", "3124");
    reorder("s2c", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6h
    reorder("h2c", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("t2c+_v", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q6i
    reorder("t2c+_g", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("s2c", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6k
    reorder("prop_pp", "r1", "21");
    mult("h2c", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6l
    reorder("e2c", "r1", "1423");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    if (triples_enabled()) {

    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_1h1p_0h0p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "hp", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c", "t1c_g", "10", 0);
    restrict_valence("t1c", "t1c_gv", "11", 0);

    restrict_valence("t2c", "t2c_gv", "1010", 0);
    restrict_valence("t2c", "t2c_g", "1000", 0);
    restrict_valence("t2c", "t2c_v", "0010", 0);

    restrict_valence("prop_hh", "prop_gh", "10", 0);
    restrict_valence("prop_hp", "prop_gv", "11", 0);
    restrict_valence("prop_hp", "prop_gp", "10", 0);

    dg_stack_pos_t pos = get_stack_pos();

    // O1
    update(result_name, 1.0, "prop_gv");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    reorder("t1c_v", "r1", "21");
    mult("prop_gh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L1b
    reorder("prop_pv", "r1", "21");
    mult("t1c_g", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L2a
    reorder("prop_hv", "r1", "21");
    mult("h1c+", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2b
    reorder("s1c+", "r1", "21");
    mult("prop_gp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3a
    reorder("t2c_gv", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3b
    reorder("e2c+", "r1", "2341");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a - disconnected term
    if (disconnected) {
        double factor_Q1a = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q1a, "t1c_gv");
    }

    // Q1b
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q2a - disconnected term
    if (disconnected) {
        double factor_Q2a = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q2a, "e1c+");
    }

    // Q2b
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3a - disconnected term
    if (disconnected) {
        double t1_t1_Q3a = scalar_product("C", "N", "t1c", "t1c");
        update(result_name, t1_t1_Q3a, "prop_gv");
        restore_stack_pos(pos);
    }

    // Q3b - disconnected term
    if (disconnected) {
        double factor_Q3b = scalar_product("C", "N", "t1c", "prop_hp");
        update(result_name, factor_Q3b, "t1c_gv");
    }

    // Q3c - disconnected term
    if (disconnected) {
        double factor_Q3c = scalar_product("C", "N", "prop_hp", "t1c");
        update(result_name, factor_Q3c, "e1c+");
    }

    // Q3d
    reorder("t1c_v", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3e
    reorder("prop_hv", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3f
    reorder("s1c+", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_gh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3g
    reorder("prop_pv", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q3h
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q3i
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4a
    reorder("s2c+", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4b
    reorder("prop_hv", "r1", "21");
    reorder("h2c+", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4c
    reorder("s2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4d
    reorder("h2c+", "r1", "1342");
    reorder("t1c_v", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4e
    reorder("t1c", "r1", "21");
    reorder("e2c+", "r2", "2341");
    mult("prop_hh", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q4f
    reorder("e2c+", "r1", "2341");
    mult("r1", "prop_pp", "r2", 1);
    mult("r2", "t1c", "r3", 2);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5a
    reorder("t2c_g", "r1", "1342");
    reorder("prop_pv", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q5b
    reorder("t2c_v", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_gh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5c
    reorder("t2c_g", "r1", "1342");
    reorder("s1c+", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q5d
    reorder("t2c_v", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5e
    reorder("prop_hh", "r1", "21");
    reorder("t2c_gv", "r2", "1324");
    mult("r1", "t1c+", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5f
    reorder("t1c+", "r1", "21");
    reorder("t2c_gv", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6a - disconnected term
    if (disconnected) {
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        update(result_name, t2_t2, "prop_gv");
        restore_stack_pos(pos);
    }

    // Q6b
    reorder("t2c_gv", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6c
    reorder("t2c", "r1", "2341");
    reorder("s2c+", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("prop_gh", "r3", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6d
    reorder("t2c", "r1", "3412");
    reorder("prop_pv", "r2", "21");
    mult("h2c+", "r1", "r3", 3);
    mult("r3", "r2", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6e
    reorder("t2c_g", "t2c_21", "2143");
    reorder("t2c_21", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("prop_hv", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);

    // Q6f
    reorder("t2c_v", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6g
    reorder("s2c+", "r1", "3124");
    reorder("t2c_g", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6h
    reorder("t2c_v", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("h2c+", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q6i
    reorder("s2c+", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("t2c_g", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6k
    reorder("prop_pp", "r1", "21");
    mult("t2c_v", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("h2c+", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6l
    reorder("t2c", "r1", "1342");
    reorder("e2c+", "r2", "2341");
    mult("r1", "prop_ph", "r3", 2);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    if (triples_enabled()) {

    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


void direct_property_1h1p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "phph", "1111", "1234", NOT_PERM_UNIQUE, operator_symmetry);

    /*restrict_valence("t1c", "t1c_v2", "01", 0);
    restrict_valence("t1c+", "t1c+_v1", "01", 0);
    restrict_valence("t1c", "t1c_v", "01", 0);
    restrict_valence("t1c+", "t1c+_v", "10", 0);

    restrict_valence("t2c", "t2c_v12", "0011", 0);
    restrict_valence("t2c+", "t2c+_v12", "1100", 0);
    restrict_valence("t2c", "t2c_v1", "0010", 0);
    restrict_valence("t2c+", "t2c+_v1", "1000", 0);
    restrict_valence("t2c", "t2c_v2", "0001", 0);
    restrict_valence("t2c+", "t2c+_v2", "0100", 0);

    restrict_valence("s2c", "s2c_v12", "1011", 0);
    restrict_valence("s2c+", "s2c+_v12", "1110", 0);
    restrict_valence("s2c", "s2c_v1", "1010", 0);
    restrict_valence("s2c+", "s2c+_v1", "1010", 0);
    restrict_valence("s2c", "s2c_v2", "1001", 0);
    restrict_valence("s2c+", "s2c+_v2", "0110", 0);

    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hp", "prop_hv", "01", 0);

    restrict_valence("prop_pp", "prop_pv", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);
    restrict_valence("prop_pp", "prop_vv", "11", 0);*/

    dg_stack_pos_t pos = get_stack_pos();


    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }


    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }


    if (triples_enabled()) {

    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}

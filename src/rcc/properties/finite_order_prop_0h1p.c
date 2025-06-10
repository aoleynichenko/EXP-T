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

#include "finite_order_prop_0h1p.h"

#include "cc_properies.h"
#include "engine.h"
#include "spinors.h"

#include "finite_order_prop.h"
#include "finite_order_overlap.h"


void finite_order_property_0h1p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "pp", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    restrict_valence("t1c", "t1c_v2", "01", 0);
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
    restrict_valence("prop_pp", "prop_vv", "11", 0);

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_v1", "000100", 0);
        restrict_valence("t3c+", "t3c+_v1", "100000", 0);

        restrict_valence("s3c", "s3c_v1", "100100", 0);
        restrict_valence("s3c+", "s3c+_v1", "100100", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();

    int print_level = cc_opts->print_level;
    if (print_level >= CC_PRINT_HIGH) {
        printf(" diagram                 time (sec)  max mem (mb)  max mem (gb)\n");
    }

    // O1
    double t_start_O1 = abs_time();
    update(result_name, 1.0, "prop_vv");
    restore_stack_pos(pos);
    print_time_mem_usage("O1", t_start_O1, print_level);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    double t_start_L1a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("prop_vh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L1a", t_start_L1a, print_level);

    // L1b
    double t_start_L1b = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("t1c+_v", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L1b", t_start_L1b, print_level);

    // L2a
    double t_start_L2a = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L2a", t_start_L2a, print_level);

    // L2b
    double t_start_L2b = abs_time();
    reorder("s1c+", "r1", "21");
    mult("prop_vp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L2b", t_start_L2b, print_level);

    // L3a
    double t_start_L3a = abs_time();
    reorder("s2c_v1", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L3a", t_start_L3a, print_level);

    // L3b
    double t_start_L3b = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L3b", t_start_L3b, print_level);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a
    double t_start_Q1a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1a", t_start_Q1a, print_level);

    // Q1b
    double t_start_Q1b = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1b", t_start_Q1b, print_level);

    // Q1c
    double t_start_Q1c = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1c", t_start_Q1c, print_level);

    // Q1d - disconnected term
    if (disconnected) {
        double t_start_Q1d = abs_time();
        double t1_t1 = scalar_product("C", "N", "t1c", "t1c");
        update(result_name, t1_t1, "prop_vv");
        restore_stack_pos(pos);
        print_time_mem_usage("Q1d", t_start_Q1d, print_level);
    }

    // Q2a
    double t_start_Q2a = abs_time();
    reorder("t2c_v1", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2a", t_start_Q2a, print_level);

    // Q2b
    double t_start_Q2b = abs_time();
    reorder("t2c+_v1", "r1", "1342");
    reorder("t1c_v", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2b", t_start_Q2b, print_level);

    // Q2c
    double t_start_Q2c = abs_time();
    reorder("t2c_v1", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2c", t_start_Q2c, print_level);

    // Q2d
    double t_start_Q2d = abs_time();
    reorder("prop_hv", "r1", "21");
    reorder("t2c+_v1", "r2", "1342");
    mult("r2", "t1c", "r3", 2);
    mult("r3", "r1", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2d", t_start_Q2d, print_level);

    // Q3a - disconnected term
    if (disconnected) {
        double t_start_Q3a = abs_time();
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        update(result_name, t2_t2, "prop_vv");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3a", t_start_Q3a, print_level);
    }

    // Q3b
    double t_start_Q3b = abs_time();
    reorder("t2c+", "r1", "3412");
    reorder("t2c_v2", "r2", "4123");
    mult("prop_vp", "r1", "r3", 1);
    mult("r3", "r2", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3b", t_start_Q3b, print_level);

    // Q3c
    double t_start_Q3c = abs_time();
    reorder("prop_pv", "r1", "21");
    reorder("t2c+_v2", "r2", "2341");
    mult("r1", "t2c", "r3", 1);
    mult("r2", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3c", t_start_Q3c, print_level);

    // Q3d
    double t_start_Q3d = abs_time();
    reorder("t2c_v1", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("t2c+_v1", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3d", t_start_Q3d, print_level);

    // Q3e
    double t_start_Q3e = abs_time();
    reorder("prop_pp", "r1", "21");
    mult("t2c_v1", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v1", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3e", t_start_Q3e, print_level);

    // Q4a
    double t_start_Q4a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4a", t_start_Q4a, print_level);

    // Q4b
    double t_start_Q4b = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4b", t_start_Q4b, print_level);

    // Q4c
    double t_start_Q4c = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_vh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4c", t_start_Q4c, print_level);

    // Q4d
    double t_start_Q4d = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4d", t_start_Q4d, print_level);

    // Q5a
    double t_start_Q5a = abs_time();
    reorder("prop_hh", "r1", "21");
    reorder("s2c_v1", "r2", "1324");
    mult("r1", "t1c+", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5a", t_start_Q5a, print_level);

    // Q5b
    double t_start_Q5b = abs_time();
    reorder("t1c", "r1", "21");
    reorder("s2c+_v1", "r2", "1342");
    mult("prop_hh", "r1", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5b", t_start_Q5b, print_level);

    // Q5c
    double t_start_Q5c = abs_time();
    reorder("t1c+", "r1", "21");
    reorder("s2c_v1", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5c", t_start_Q5c, print_level);

    // Q5d
    double t_start_Q5d = abs_time();
    reorder("s2c+_v1", "r1", "1324");
    reorder("prop_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5d", t_start_Q5d, print_level);

    // Q6a
    double t_start_Q6a = abs_time();
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r1", "2341");
    reorder("t2c+", "r2", "4123");
    reorder("prop_hv", "r3", "21");
    mult("r1", "r2", "r4", 3);
    mult("r4", "r3", "r5", 1);
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6a", t_start_Q6a, print_level);

    // Q6b
    double t_start_Q6b = abs_time();
    reorder("t2c", "r1", "2341");
    reorder("s2c+", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "4123");
    mult("r2", "r1", "r3", 3);
    mult("prop_vh", "r3", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6b", t_start_Q6b, print_level);

    // Q6c
    double t_start_Q6c = abs_time();
    reorder("s2c_v1", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6c", t_start_Q6c, print_level);

    // Q6d
    double t_start_Q6d = abs_time();
    reorder("s2c+_v1", "r1", "1324");
    reorder("t2c", "r2", "3142");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6d", t_start_Q6d, print_level);

    // Q7
    double t_start_Q7 = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q7", t_start_Q7, print_level);

    // Q8a
    double t_start_Q8a = abs_time();
    reorder("s2c", "r1", "1342");
    reorder("s1c+", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8a", t_start_Q8a, print_level);

    // Q8b
    double t_start_Q8b = abs_time();
    reorder("s2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("s1c", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8b", t_start_Q8b, print_level);

    // Q8c
    double t_start_Q8c = abs_time();
    reorder("s2c", "r1", "1342");
    reorder("prop_pv", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8c", t_start_Q8c, print_level);

    // Q8d
    double t_start_Q8d = abs_time();
    reorder("s2c+", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_vp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8d", t_start_Q8d, print_level);

    // Q9a
    double t_start_Q9a = abs_time();
    reorder("s2c+", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("s2c", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9a", t_start_Q9a, print_level);

    // Q9b
    double t_start_Q9b = abs_time();
    reorder("s2c+", "r1", "3124");
    reorder("s2c", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9b", t_start_Q9b, print_level);

    if (triples_enabled()) {

        // Q10a
        double t_start_Q10a = abs_time();
        reorder("t3c_v1", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        reorder("r2", "r3", "21");
        mult("prop_vh", "r3", "r4", 1);
        update(result_name, -0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10a", t_start_Q10a, print_level);

        // Q10b
        double t_start_Q10b = abs_time();
        reorder("t3c+_v1", "r1", "145623");
        reorder("prop_hv", "r2", "21");
        mult("r1", "t2c", "r3", 4);
        mult("r3", "r2", "r4", 1);
        update(result_name, -0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10b", t_start_Q10b, print_level);

        // Q10c
        double t_start_Q10c = abs_time();
        reorder("t3c_v1", "r1", "451263");
        mult("r1", "prop_ph", "r2", 2);
        mult("t2c+_v1", "r2", "r3", 3);
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10c", t_start_Q10c, print_level);

        // Q10d
        double t_start_Q10d = abs_time();
        reorder("t2c_v1", "r1", "3412");
        reorder("t3c+_v1", "r2", "124563");
        mult("r2", "prop_hp", "r3", 2);
        mult("r3", "r1", "r4", 3);
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10d", t_start_Q10d, print_level);

        // Q11a
        double t_start_Q11a = abs_time();
        reorder("t3c", "r1", "456123");
        reorder("prop_pv", "r2", "21");
        mult("t3c+_v1", "r1", "r3", 5);
        mult("r3", "r2", "r4", 1);
        update(result_name, -1.0 / 12.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11a", t_start_Q11a, print_level);

        // Q11b
        double t_start_Q11b = abs_time();
        reorder("t3c_v1", "r1", "456123");
        mult("r1", "t3c+", "r2", 5);
        mult("prop_vp", "r2", "r3", 1);
        update(result_name, -1.0 / 12.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11b", t_start_Q11b, print_level);

        // Q11c
        double t_start_Q11c = abs_time();
        reorder("prop_pp", "r1", "21");
        mult("t3c_v1", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("t3c+_v1", "r3", "r4", 5);
        update(result_name, -1.0 / 6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11c", t_start_Q11c, print_level);

        // Q11d
        double t_start_Q11d = abs_time();
        reorder("t3c_v1", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        mult("t3c+_v1", "r2", "r3", 5);
        update(result_name, 0.25, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11d", t_start_Q11d, print_level);

        // Q11e - disconnected term
        if (disconnected) {
            double t3_t3 = (1.0 / 36.0) * scalar_product("C", "N", "t3c", "t3c");
            copy("prop_pp", "r0");
            closed("r0", "r1");
            update(result_name, t3_t3, "r1");
            restore_stack_pos(pos);
        }

        // Q12a
        double t_start_Q12a = abs_time();
        reorder("s3c_v1", "r1", "145263");
        mult("r1", "prop_ph", "r2", 2);
        mult("r2", "t1c+", "r3", 2);
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12a", t_start_Q12a, print_level);

        // Q12b
        double t_start_Q12b = abs_time();
        reorder("s3c+_v1", "r1", "145263");
        mult("r1", "prop_hp", "r2", 2);
        mult("r2", "t1c", "r3", 2);
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12b", t_start_Q12b, print_level);

        // Q13a
        double t_start_Q13a = abs_time();
        reorder("s3c", "r1", "145623");
        reorder("prop_pv", "r2", "21");
        mult("r1", "t2c+", "r3", 4);
        mult("r3", "r2", "r4", 1);
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13a", t_start_Q13a, print_level);

        // Q13b
        double t_start_Q13b = abs_time();
        reorder("s3c+", "r1", "415623");
        mult("r1", "t2c", "r2", 4);
        mult("prop_vp", "r2", "r3", 1);
        update(result_name, 0.25, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13b", t_start_Q13b, print_level);

        // Q13c
        double t_start_Q13c = abs_time();
        reorder("t2c+", "r1", "3412");
        reorder("s3c_v1", "r2", "142356");
        mult("r1", "prop_pp", "r3", 1);
        mult("r2", "r3", "r4", 4);
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13c", t_start_Q13c, print_level);

        // Q13d
        double t_start_Q13d = abs_time();
        reorder("s3c+_v1", "r1", "145623");
        reorder("prop_pp", "r2", "21");
        mult("t2c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 4);
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13d", t_start_Q13d, print_level);

        // Q13e
        double t_start_Q13e = abs_time();
        reorder("prop_hh", "r1", "21");
        reorder("s3c_v1", "r2", "145623");
        mult("t2c+", "r1", "r3", 1);
        mult("r2", "r3", "r4", 4);
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13e", t_start_Q13e, print_level);

        // Q13f
        double t_start_Q13f = abs_time();
        reorder("s3c+_v1", "r1", "142356");
        reorder("t2c", "r2", "3412");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 4);
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13f", t_start_Q13f, print_level);

        // Q14a
        double t_start_Q14a = abs_time();
        reorder("t3c+", "r1", "124563");
        reorder("s3c_v1", "r2", "145623");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 4);
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14a", t_start_Q14a, print_level);

        // Q14b
        double t_start_Q14b = abs_time();
        reorder("s3c+_v1", "r1", "142356");
        reorder("t3c", "r2", "451263");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 4);
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14b", t_start_Q14b, print_level);

        // Q14c
        double t_start_Q14c = abs_time();
        reorder("s3c+", "r1", "456123");
        mult("r1", "t3c", "r2", 5);
        mult("prop_vh", "r2", "r3", 1);
        update(result_name, -1.0 / 12.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14c", t_start_Q14c, print_level);

        // Q14d
        double t_start_Q14d = abs_time();
        reorder("t3c+", "r1", "456123");
        reorder("prop_hv", "r2", "21");
        mult("s3c", "r1", "r3", 5);
        mult("r3", "r2", "r4", 1);
        update(result_name, -1.0 / 12.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14d", t_start_Q14d, print_level);

        // Q15a
        double t_start_Q15a = abs_time();
        reorder("s2c+", "r1", "3412");
        reorder("s3c", "r2", "124563");
        mult("r2", "prop_ph", "r3", 2);
        mult("r3", "r1", "r4", 3);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15a", t_start_Q15a, print_level);

        // Q15b
        double t_start_Q15b = abs_time();
        reorder("s3c+", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("s2c", "r2", "r3", 3);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15b", t_start_Q15b, print_level);

        // Q16a
        double t_start_Q16a = abs_time();
        reorder("s3c+", "r1", "456123");
        mult("r1", "prop_pp", "r2", 1);
        mult("s3c", "r2", "r3", 5);
        closed("r3", "r4");
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16a", t_start_Q16a, print_level);

        // Q16b
        double t_start_Q16b = abs_time();
        reorder("prop_hh", "r1", "21");
        mult("s3c+", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("s3c", "r3", "r4", 5);
        closed("r4", "r5");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16b", t_start_Q16b, print_level);
    }

    /*
     * for the "connected" scheme only: connected diagrams ((T^+T)_conn O)_conn.
     * In fact, "overlap" (or metric) operator (T^+ T)_conn is contracted with
     * the property operator. The sign standing before this term is "minus".
     */
    if (disconnected == 0) {
        pos = get_stack_pos();
        sector_0h1p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
        reorder("TT_conn", "r1", "21");
        mult("prop_vv", "r1", "r2", 1);
        update(result_name, -1.0, "r2");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}

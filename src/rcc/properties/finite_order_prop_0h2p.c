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

#include "finite_order_prop_0h2p.h"

#include "cc_properies.h"
#include "engine.h"
#include "spinors.h"

#include "finite_order_prop.h"
#include "finite_order_overlap.h"


void finite_order_property_0h2p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE, operator_symmetry);

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

    restrict_valence("x2c", "x2c_v1", "1110", 0);
    restrict_valence("x2c+", "x2c+_v1", "1011", 0);

    restrict_valence("prop_ph", "prop_vh", "10", 0);
    restrict_valence("prop_hp", "prop_hv", "01", 0);

    restrict_valence("prop_pp", "prop_pv", "01", 0);
    restrict_valence("prop_pp", "prop_vp", "10", 0);
    restrict_valence("prop_pp", "prop_vv", "11", 0);

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_v1", "000100", 0);
        restrict_valence("t3c+", "t3c+_v1", "100000", 0);
        restrict_valence("t3c", "t3c_v2", "000010", 0);
        restrict_valence("t3c+", "t3c+_v2", "010000", 0);
        restrict_valence("t3c", "t3c_v12", "000110", 0);
        restrict_valence("t3c+", "t3c+_v12", "110000", 0);

        restrict_valence("s3c", "s3c_v1", "100100", 0);
        restrict_valence("s3c+", "s3c+_v1", "100100", 0);
        restrict_valence("s3c", "s3c_v2", "100010", 0);
        restrict_valence("s3c+", "s3c+_v2", "010100", 0);
        restrict_valence("s3c", "s3c_v12", "100110", 0);
        restrict_valence("s3c+", "s3c+_v12", "110100", 0);

        restrict_valence("x3c", "x3c_v1", "110100", 0);
        restrict_valence("x3c+", "x3c+_v1", "100110", 0);
        restrict_valence("x3c", "x3c_v12", "110110", 0);
        restrict_valence("x3c+", "x3c+_v12", "110110", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();
    int print_level = cc_opts->print_level;

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    if (print_level >= CC_PRINT_HIGH) {
        printf(" diagram                 time (sec)  max mem (mb)  max mem (gb)\n");
    }

    // L1a
    double t_start_L1a = abs_time();
    reorder("s2c_v12", "r1", "1342");
    mult("r1", "prop_vh", "r2", 1);
    reorder("r2", "r3", "1423");
    perm("r3", "(12)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("L1a", t_start_L1a, print_level);

    // L1b
    double t_start_L1b = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("s2c+_v12", "r1", "r2", 1);
    perm("r2", "(34)");
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L1b", t_start_L1b, print_level);

    // L2a
    double t_start_L2a = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("x2c_v1", "r1", "r2", 1);
    perm("r2", "(34)");
    update(result_name, +1.0, "r2");
    restore_stack_pos(pos);
    print_time_mem_usage("L2a", t_start_L2a, print_level);

    // L2b
    double t_start_L2b = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    mult("r1", "prop_vp", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("L2b", t_start_L2b, print_level);

    if (triples_enabled()) {

        // L3a
        reorder("x3c", "r1", "124563");
        mult("r1", "prop_ph", "r2", 2);
        closed("r2", "r3");
        update(result_name, +1.0, "r3");
        restore_stack_pos(pos);

        // L3b
        reorder("x3c+", "r1", "124563");
        mult("r1", "prop_hp", "r2", 2);
        closed("r2", "r3");
        update(result_name, +1.0, "r3");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1 - disconnected term
    if (disconnected) {
        double t_start_Q1 = abs_time();
        reorder("t1c", "r2", "21");
        mult("t1c+", "r2", "r3", 1);
        closed("r3", "r4");
        construct_disconnected_rank2_rank2("r4", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q1", t_start_Q1, print_level);
    }

    // Q2a
    double t_start_Q2a = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("prop_vh", "r1", "r2", 1);
    mult("t1c+_v", "r2", "r3", 1);
    perm("r3", "(12)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2a", t_start_Q2a, print_level);

    // Q2b
    double t_start_Q2b = abs_time();
    reorder("t1c_v", "r1", "21");
    reorder("prop_hv", "r2", "21");
    mult("r2", "t2c+_v12", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2b", t_start_Q2b, print_level);

    // Q3a - disconnected term
    if (disconnected) {
        double t_start_Q3a = abs_time();
        reorder("t2c_v1", "r1", "3412");
        mult("t2c+_v1", "r1", "r2", 3);
        construct_disconnected_rank2_rank2("r2", "prop_vv", "r3");
        perm("r3", "(12|34)");
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q3a", t_start_Q3a, print_level);
    }

    // Q3b
    double t_start_Q3b = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("r1", "t2c+_v1", "r2", 2);
    mult("r2", "prop_vp", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update(result_name, 0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3b", t_start_Q3b, print_level);

    // Q3c
    double t_start_Q3c = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("t2c_v1", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c+_v12", "r3", "r4", 2);
    perm("r4", "(34)");
    update(result_name, 0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3c", t_start_Q3c, print_level);

    // Q3d
    double t_start_Q3d = abs_time();
    reorder("t2c_v12", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("t2c+_v12", "r2", "r3", 2);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3d", t_start_Q3d, print_level);

    // Q4a - disconnected term
    if (disconnected) {
        double t_start_Q4a = abs_time();
        reorder("s2c_v1", "r1", "1342");
        mult("r1", "t1c+", "r2", 2);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, +1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q4a", t_start_Q4a, print_level);
    }

    // Q4b - disconnected term
    if (disconnected) {
        double t_start_Q4b = abs_time();
        reorder("s2c+_v1", "r1", "1342");
        mult("r1", "t1c", "r2", 2);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, +1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q4b", t_start_Q4b, print_level);
    }

    // Q4c
    double t_start_Q4c = abs_time();
    reorder("t1c+", "r1", "21");
    reorder("s2c_v12", "r2", "3412");
    mult("prop_vp", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    closed("r5", "r6");
    perm("r6", "(12)");
    update(result_name, -1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4c", t_start_Q4c, print_level);

    // Q4d
    double t_start_Q4d = abs_time();
    reorder("prop_pv", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("s2c+_v12", "r2", "r3", 1);
    closed("r3", "r4");
    perm("r4", "(34)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4d", t_start_Q4d, print_level);

    // Q4e
    double t_start_Q4e = abs_time();
    reorder("s2c_v2", "r1", "1423");
    reorder("prop_pv", "r2", "21");
    mult("r2", "r1", "r3", 1);
    mult("r3", "t1c+_v", "r4", 1);
    reorder("r4", "r5", "2413");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4e", t_start_Q4e, print_level);

    // Q4f
    double t_start_Q4f = abs_time();
    reorder("s2c+_v2", "r1", "2341");
    reorder("t1c_v", "r2", "21");
    mult("prop_vp", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(12|34)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4f", t_start_Q4f, print_level);

    // Q4g
    double t_start_Q4g = abs_time();
    reorder("s2c_v12", "r1", "3412");
    reorder("prop_hh", "r2", "21");
    mult("t1c+_v", "r2", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4g", t_start_Q4g, print_level);

    // Q4h
    double t_start_Q4h = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("s2c+_v12", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4h", t_start_Q4h, print_level);

    // Q5a
    double t_start_Q5a = abs_time();
    reorder("t2c_v1", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("s2c+_v12", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5a", t_start_Q5a, print_level);

    // Q5b
    double t_start_Q5b = abs_time();
    reorder("s2c_v12", "r1", "3412");
    reorder("t2c+_v1", "r2", "1342");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5b", t_start_Q5b, print_level);

    // Q5c
    double t_start_Q5c = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    reorder("t2c_v2", "r2", "4213");
    mult("r1", "r2", "r3", 2);
    mult("r3", "prop_vh", "r4", 1);
    reorder("r4", "r5", "1423");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5c", t_start_Q5c, print_level);

    // Q5d
    double t_start_Q5d = abs_time();
    reorder("prop_hv", "r1", "21");
    reorder("t2c+_v2", "r2", "2314");
    reorder("s2c_v1", "r3", "1324");
    mult("r1", "r2", "r4", 1);
    mult("r3", "r4", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(12|34)");
    update(result_name, -1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q5d", t_start_Q5d, print_level);

    // Q6a
    double t_start_Q6a = abs_time();
    reorder("t1c_v", "r1", "21");
    mult("r1", "prop_ph", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6a", t_start_Q6a, print_level);

    // Q6b
    double t_start_Q6b = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    mult("r1", "prop_hp", "r2", 1);
    mult("r2", "t1c+_v", "r3", 1);
    reorder("r3", "r4", "3412");
    perm("r4", "(12)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6b", t_start_Q6b, print_level);

    // Q6c
    double t_start_Q6c = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6c", t_start_Q6c, print_level);

    // Q6d
    double t_start_Q6d = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    reorder("t1c", "r2", "21");
    mult("prop_vh", "r2", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q6d", t_start_Q6d, print_level);

    // Q7 - disconnected term
    if (disconnected) {
        double t_start_Q7 = abs_time();
        reorder("s1c+", "r1", "21");
        mult("s1c", "r1", "r2", 1);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7", t_start_Q7, print_level);
    }

    // Q8a
    double t_start_Q8a = abs_time();
    reorder("prop_ph", "r1", "21");
    reorder("s2c_v12", "r2", "3412");
    mult("s1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8a", t_start_Q8a, print_level);

    // Q8b
    double t_start_Q8b = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("s2c+_v12", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8b", t_start_Q8b, print_level);

    // Q8c
    double t_start_Q8c = abs_time();
    reorder("s1c+", "r1", "21");
    reorder("s2c_v2", "r2", "4123");
    mult("r1", "r2", "r3", 1);
    mult("r3", "prop_vh", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8c", t_start_Q8c, print_level);

    // Q8d
    double t_start_Q8d = abs_time();
    reorder("s2c+_v2", "r1", "2341");
    reorder("prop_hv", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(12|34)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q8d", t_start_Q8d, print_level);

    // Q9a - disconnected term
    if (disconnected) {
        double t_start_Q9a = abs_time();
        reorder("s2c+", "r1", "3412");
        mult("s2c", "r1", "r2", 3);
        closed("r2", "r3");
        construct_disconnected_rank2_rank2("r3", "prop_vv", "r5");
        perm("r5", "(12|34)");
        update(result_name, +0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9a", t_start_Q9a, print_level);
    }

    // Q9b
    double t_start_Q9b = abs_time();
    reorder("s2c+_v1", "r1", "1342");
    reorder("s2c", "s2c_21", "2143");
    reorder("s2c_21", "r2", "2134");
    reorder("prop_pv", "r3", "21");
    mult("r3", "r2", "r4", 1);
    mult("r1", "r4", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(12|34)");
    update(result_name, 1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9b", t_start_Q9b, print_level);

    // Q9c
    double t_start_Q9c = abs_time();
    reorder("s2c_v1", "r1", "1324");
    reorder("s2c+", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "4312");
    mult("prop_vp", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9c", t_start_Q9c, print_level);

    // Q9d
    double t_start_Q9d = abs_time();
    reorder("s2c_v1", "r1", "1342");
    reorder("s2c+_v1", "s2c+_21", "2143");
    reorder("s2c+_21", "r2", "2413");
    mult("r1", "prop_hh", "r3", 1);
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9d", t_start_Q9d, print_level);

    // Q9e
    double t_start_Q9e = abs_time();
    reorder("s2c+_v1", "s2c+_21", "2143");
    reorder("s2c+_21", "r1", "2431");
    reorder("s2c_v1", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12|34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q9e", t_start_Q9e, print_level);

    // Q10a
    double t_start_Q10a = abs_time();
    reorder("prop_pv", "r1", "21");
    reorder("s1c+", "r2", "21");
    mult("r1", "x2c", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(34)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10a", t_start_Q10a, print_level);

    // Q10b
    double t_start_Q10b = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("prop_vp", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    perm("r3", "(12)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10b", t_start_Q10b, print_level);

    // Q10c
    double t_start_Q10c = abs_time();
    reorder("s1c+", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10c", t_start_Q10c, print_level);

    // Q10d
    double t_start_Q10d = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    reorder("prop_pp", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q10d", t_start_Q10d, print_level);

    // Q11a
    double t_start_Q11a = abs_time();
    reorder("x2c+_v1", "r1", "3412");
    reorder("s2c", "r2", "1342");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11a", t_start_Q11a, print_level);

    // Q11b
    double t_start_Q11b = abs_time();
    reorder("s2c+", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("x2c_v1", "r2", "r3", 1);
    perm("r3", "(34)");
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11b", t_start_Q11b, print_level);

    // Q11c
    double t_start_Q11c = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("s2c", "r1", "r2", 2);
    reorder("r2", "r3", "3412");
    mult("r3", "prop_vh", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update(result_name, -0.5, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11c", t_start_Q11c, print_level);

    // Q11d
    double t_start_Q11d = abs_time();
    reorder("prop_hv", "r1", "21");
    mult("s2c+", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("x2c", "r3", "r4", 2);
    perm("r4", "(34)");
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q11d", t_start_Q11d, print_level);

    // Q12
    double t_start_Q12 = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("r1", "prop_pp", "r2", 1);
    mult("x2c", "r2", "r3", 2);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);
    print_time_mem_usage("Q12", t_start_Q12, print_level);

    if (triples_enabled()) {

        // Q13a
        double t_start_Q13a = abs_time();
        reorder("t3c_v12", "r1", "451263");
        mult("r1", "prop_ph", "r2", 2);
        mult("t2c+_v12", "r2", "r3", 2);
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13a", t_start_Q13a, print_level);

        // Q13b
        double t_start_Q13b = abs_time();
        reorder("t2c_v12", "r1", "3412");
        reorder("t3c+_v12", "r2", "124563");
        mult("r2", "prop_hp", "r3", 2);
        mult("r3", "r1", "r4", 2);
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13b", t_start_Q13b, print_level);

        // Q13c
        double t_start_Q13c = abs_time();
        reorder("t3c_v12", "r1", "451623");
        mult("t2c+_v1", "r1", "r2", 3);
        mult("prop_vh", "r2", "r3", 1);
        perm("r3", "(12)");
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13c", t_start_Q13c, print_level);

        // Q13d
        double t_start_Q13d = abs_time();
        reorder("t2c_v1", "r1", "3412");
        reorder("t3c+_v12", "r2", "124356");
        reorder("prop_hv", "r3", "21");
        mult("r1", "r2", "r4", 3);
        mult("r3", "r4", "r5", 1);
        reorder("r5", "r6", "3412");
        perm("r6", "(34)");
        update(result_name, 0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13d", t_start_Q13d, print_level);

        // Q14a
        double t_start_Q14a = abs_time();
        reorder("t3c_v2", "r1", "456123");
        reorder("prop_pv", "r2", "21");
        mult("t3c+_v12", "r1", "r3", 4);
        reorder("r3", "r4", "1243");
        mult("r4", "r2", "r5", 1);
        reorder("r5", "r6", "1243");
        perm("r6", "(34)");
        update(result_name, 1.0 / 6.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14a", t_start_Q14a, print_level);

        // Q14b
        double t_start_Q14b = abs_time();
        reorder("t3c_v12", "r1", "456123");
        mult("t3c+_v2", "r1", "r2", 4);
        reorder("r2", "r3", "2341");
        mult("prop_vp", "r3", "r4", 1);
        perm("r4", "(12)");
        update(result_name, 1.0 / 6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14b", t_start_Q14b, print_level);

        // Q14c
        double t_start_Q14c = abs_time();
        reorder("prop_pp", "r1", "21");
        mult("t3c_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("t3c+_v12", "r3", "r4", 4);
        update(result_name, 1.0 / 6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14c", t_start_Q14c, print_level);

        // Q14d
        double t_start_Q14d = abs_time();
        reorder("t3c_v12", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        mult("t3c+_v12", "r2", "r3", 4);
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14d", t_start_Q14d, print_level);

        // Q14e - disconnected term
        if (disconnected) {
            reorder("t3c", "r1", "456123");
            mult("t3c+", "r1", "r2", 5);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, -1.0 / 12.0, "r5");
            restore_stack_pos(pos);
        }

        // Q15a
        double t_start_Q15a = abs_time();
        reorder("s3c_v12", "r1", "145263");
        mult("r1", "t1c+", "r2", 2);
        mult("r2", "prop_vh", "r3", 1);
        reorder("r3", "r4", "1423");
        perm("r4", "(12)");
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15a", t_start_Q15a, print_level);

        // Q15b
        double t_start_Q15b = abs_time();
        reorder("s3c+_v12", "r1", "124563");
        reorder("prop_hv", "r2", "21");
        mult("r1", "t1c", "r3", 2);
        mult("r3", "r2", "r4", 1);
        perm("r4", "(34)");
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15b", t_start_Q15b, print_level);

        // Q15c
        double t_start_Q15c = abs_time();
        reorder("s3c_v12", "r1", "145263");
        mult("r1", "prop_ph", "r2", 2);
        mult("r2", "t1c+_v", "r3", 1);
        reorder("r3", "r4", "1423");
        perm("r4", "(12)");
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15c", t_start_Q15c, print_level);

        // Q15d
        double t_start_Q15d = abs_time();
        reorder("s3c+_v12", "r1", "124563");
        reorder("t1c_v", "r2", "21");
        mult("r1", "prop_hp", "r3", 2);
        mult("r3", "r2", "r4", 1);
        perm("r4", "(34)");
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15d", t_start_Q15d, print_level);

        // Q16a
        double t_start_Q16a = abs_time();
        reorder("s3c_v12", "r1", "145623");
        reorder("t2c+_v1", "r2", "3412");
        mult("r2", "prop_pp", "r3", 1);
        reorder("r3", "r4", "3412");
        mult("r1", "r4", "r5", 3);
        reorder("r5", "r6", "1423");
        perm("r6", "(12)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16a", t_start_Q16a, print_level);

        // Q16b
        double t_start_Q16b = abs_time();
        reorder("s3c+_v12", "r1", "124356");
        reorder("prop_pp", "r2", "21");
        mult("t2c_v1", "r2", "r3", 1);
        reorder("r3", "r4", "3412");
        mult("r1", "r4", "r5", 3);
        perm("r5", "(34)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16b", t_start_Q16b, print_level);

        // Q16c
        double t_start_Q16c = abs_time();
        reorder("prop_hh", "r1", "21");
        reorder("s3c_v12", "r2", "145623");
        mult("t2c+_v1", "r1", "r3", 1);
        mult("r2", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        perm("r5", "(12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16c", t_start_Q16c, print_level);

        // Q16d
        double t_start_Q16d = abs_time();
        reorder("s3c+_v12", "r1", "124356");
        reorder("t2c_v1", "r2", "3412");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 3);
        perm("r4", "(34)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16d", t_start_Q16d, print_level);

        // Q16e
        double t_start_Q16e = abs_time();
        reorder("s3c_v2", "r1", "154623");
        reorder("prop_pv", "r2", "21");
        mult("r1", "t2c+_v1", "r3", 3);
        reorder("r3", "r4", "1243");
        mult("r4", "r2", "r5", 1);
        reorder("r5", "r6", "1342");
        perm("r6", "(12|34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16e", t_start_Q16e, print_level);

        // Q16f
        double t_start_Q16f = abs_time();
        reorder("s3c+_v2", "r1", "241356");
        reorder("t2c_v1", "r2", "3412");
        mult("r1", "r2", "r3", 3);
        reorder("r3", "r4", "1243");
        mult("prop_vp", "r4", "r5", 1);
        perm("r5", "(12|34)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16f", t_start_Q16f, print_level);

        // Q16g - disconnected term
        if (disconnected) {
            reorder("s3c", "r1", "145623");
            mult("r1", "t2c+", "r2", 4);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, 0.25, "r5");
            restore_stack_pos(pos);
        }

        // Q16h - disconnected term
        if (disconnected) {
            reorder("s3c+", "r1", "145623");
            mult("r1", "t2c", "r2", 4);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, 0.25, "r5");
            restore_stack_pos(pos);
        }

        // Q17a
        double t_start_Q17a = abs_time();
        reorder("t3c+_v1", "r1", "124563");
        reorder("s3c_v12", "r2", "145623");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        perm("r5", "(12)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17a", t_start_Q17a, print_level);

        // Q17b
        double t_start_Q17b = abs_time();
        reorder("s3c+_v12", "r1", "124356");
        reorder("t3c_v1", "r2", "451263");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 3);
        perm("r4", "(34)");
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17b", t_start_Q17b, print_level);

        // Q17c
        double t_start_Q17c = abs_time();
        reorder("t3c+_v2", "r1", "245613");
        reorder("s3c_v2", "r2", "152346");
        reorder("prop_hv", "r3", "21");
        mult("r2", "r1", "r4", 4);
        mult("r4", "r3", "r5", 1);
        reorder("r5", "r6", "1342");
        perm("r6", "(12|34)");
        update(result_name, 0.25, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17c", t_start_Q17c, print_level);

        // Q17d
        double t_start_Q17d = abs_time();
        reorder("s3c+_v2", "r1", "245613");
        reorder("t3c_v2", "r2", "512346");
        mult("r1", "r2", "r3", 4);
        mult("prop_vh", "r3", "r4", 1);
        perm("r4", "(12|34)");
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17d", t_start_Q17d, print_level);

        // Q18a
        double t_start_Q18a = abs_time();
        reorder("s3c_v1", "r1", "142563");
        reorder("s2c+_v1", "r2", "1342");
        mult("r1", "prop_ph", "r3", 2);
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1324");
        perm("r5", "(12|34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18a", t_start_Q18a, print_level);

        // Q18b
        double t_start_Q18b = abs_time();
        reorder("s3c+_v1", "r1", "145263");
        reorder("s2c_v1", "r2", "1324");
        mult("r1", "prop_hp", "r3", 2);
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1324");
        perm("r5", "(12|34)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18b", t_start_Q18b, print_level);

        // Q18c
        double t_start_Q18c = abs_time();
        reorder("s3c_v2", "r1", "152346");
        reorder("s2c+", "r2", "3412");
        mult("r1", "r2", "r3", 3);
        reorder("r3", "r4", "4213");
        mult("r4", "prop_vh", "r5", 1);
        reorder("r5", "r6", "3412");
        perm("r6", "(12|34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18c", t_start_Q18c, print_level);

        // Q18d
        double t_start_Q18d = abs_time();
        reorder("s3c+_v2", "r1", "245613");
        reorder("prop_hv", "r2", "21");
        mult("r1", "s2c", "r3", 3);
        reorder("r3", "r4", "4123");
        mult("r4", "r2", "r5", 1);
        perm("r5", "(12|34)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18d", t_start_Q18d, print_level);

        // Q19a
        double t_start_Q19a = abs_time();
        reorder("prop_pv", "r1", "21");
        reorder("s3c+_v2", "r2", "245613");
        reorder("s3c", "r3", "152346");
        mult("r2", "r3", "r4", 4);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "3124");
        perm("r6", "(12|34)");
        update(result_name, -0.25, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q19a", t_start_Q19a, print_level);

        // Q19b
        double t_start_Q19b = abs_time();
        reorder("s3c+", "r1", "425613");
        reorder("s3c_v2", "r2", "152346");
        mult("r2", "r1", "r3", 4);
        mult("r3", "prop_vp", "r4", 1);
        reorder("r4", "r5", "1432");
        perm("r5", "(12|34)");
        update(result_name, -0.25, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q19b", t_start_Q19b, print_level);

        // Q19c
        double t_start_Q19c = abs_time();
        reorder("s3c_v1", "r1", "142356");
        reorder("s3c+_v1", "r2", "145623");
        reorder("prop_pp", "r3", "21");
        mult("r1", "r3", "r4", 1);
        mult("r4", "r2", "r5", 4);
        reorder("r5", "r6", "1324");
        perm("r6", "(12|34)");
        update(result_name, 0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q19c", t_start_Q19c, print_level);

        // Q19d
        double t_start_Q19d = abs_time();
        reorder("s3c_v1", "r1", "145623");
        reorder("prop_hh", "r2", "21");
        mult("s3c+_v1", "r2", "r3", 1);
        reorder("r3", "r4", "142356");
        mult("r1", "r4", "r5", 4);
        reorder("r5", "r6", "1324");
        perm("r6", "(12|34)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q19d", t_start_Q19d, print_level);

        // Q19e - disconnected term
        if (disconnected) {
            reorder("s3c+", "r1", "456123");
            mult("s3c", "r1", "r2", 5);
            closed("r2", "r3");
            copy("prop_pp", "r0");
            closed("r0", "r4");
            construct_disconnected_rank2_rank2("r3", "r4", "r5");
            perm("r5", "(12|34)");
            update(result_name, +1.0 / 12.0, "r5");
            restore_stack_pos(pos);
        }

        // Q20a
        double t_start_Q20a = abs_time();
        reorder("x3c_v1", "r1", "124356");
        reorder("s3c+", "r2", "451263");
        mult("r2", "prop_hp", "r3", 2);
        mult("r1", "r3", "r4", 3);
        perm("r4", "(34)");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q20a", t_start_Q20a, print_level);

        // Q20b
        double t_start_Q20b = abs_time();
        reorder("x3c+_v1", "r1", "451623");
        reorder("s3c", "r2", "124563");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 3);
        reorder("r4", "r5", "3412");
        perm("r5", "(12)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q20b", t_start_Q20b, print_level);

        // Q20c
        double t_start_Q20c = abs_time();
        reorder("s3c+", "r1", "456123");
        reorder("prop_hv", "r2", "21");
        mult("x3c", "r1", "r3", 4);
        mult("r3", "r2", "r4", 1);
        perm("r4", "(34)");
        update(result_name, -1.0 / 6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q20c", t_start_Q20c, print_level);

        // Q20d
        double t_start_Q20d = abs_time();
        reorder("x3c+", "r1", "456123");
        mult("r1", "s3c", "r2", 4);
        mult("r2", "prop_vh", "r3", 1);
        reorder("r3", "r4", "3412");
        perm("r4", "(12)");
        update(result_name, -1.0 / 6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q20d", t_start_Q20d, print_level);

        // Q21a
        double t_start_Q21a = abs_time();
        reorder("x3c_v1", "r1", "124563");
        reorder("prop_pv", "r2", "21");
        mult("r1", "t1c+", "r3", 2);
        mult("r3", "r2", "r4", 1);
        perm("r4", "(34)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21a", t_start_Q21a, print_level);

        // Q21b
        double t_start_Q21b = abs_time();
        reorder("x3c+_v1", "r1", "145263");
        mult("r1", "t1c", "r2", 2);
        mult("r2", "prop_vp", "r3", 1);
        reorder("r3", "r4", "1423");
        perm("r4", "(12)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21b", t_start_Q21b, print_level);

        // Q21c
        double t_start_Q21c = abs_time();
        reorder("prop_hh", "r1", "21");
        reorder("x3c_v12", "r2", "124536");
        mult("r1", "t1c+", "r3", 1);
        mult("r2", "r3", "r4", 2);
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21c", t_start_Q21c, print_level);

        // Q21d
        double t_start_Q21d = abs_time();
        reorder("x3c+_v12", "r1", "124536");
        reorder("t1c", "r2", "21");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 2);
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21d", t_start_Q21d, print_level);

        // Q21e
        double t_start_Q21e = abs_time();
        reorder("x3c_v12", "r1", "124536");
        reorder("t1c+", "r2", "21");
        mult("r2", "prop_pp", "r3", 1);
        mult("r1", "r3", "r4", 2);
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21e", t_start_Q21e, print_level);

        // Q21f
        double t_start_Q21f = abs_time();
        reorder("x3c+_v12", "r1", "124563");
        reorder("prop_pp", "r2", "21");
        mult("t1c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 2);
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21f", t_start_Q21f, print_level);

        // Q22a
        double t_start_Q22a = abs_time();
        reorder("x3c_v12", "r1", "124536");
        reorder("t2c+", "r2", "3142");
        mult("r2", "prop_hp", "r3", 2);
        mult("r1", "r3", "r4", 2);
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q22a", t_start_Q22a, print_level);

        // Q22b
        double t_start_Q22b = abs_time();
        reorder("x3c+_v12", "r1", "124563");
        reorder("t2c", "r2", "1342");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 2);
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q22b", t_start_Q22b, print_level);

        // Q23a
        double t_start_Q23a = abs_time();
        reorder("x3c_v1", "r1", "124563");
        reorder("s1c+", "r2", "21");
        mult("r1", "prop_ph", "r3", 2);
        mult("r3", "r2", "r4", 1);
        perm("r4", "(34)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q23a", t_start_Q23a, print_level);

        // Q23b
        double t_start_Q23b = abs_time();
        reorder("x3c+_v1", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "3412");
        perm("r4", "(12)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q23b", t_start_Q23b, print_level);

        // Q24a
        double t_start_Q24a = abs_time();
        reorder("x3c_v1", "r1", "124356");
        reorder("s2c+", "r2", "3412");
        mult("r2", "prop_pp", "r3", 1);
        mult("r1", "r3", "r4", 3);
        perm("r4", "(34)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24a", t_start_Q24a, print_level);

        // Q24b
        double t_start_Q24b = abs_time();
        reorder("x3c+_v1", "r1", "145623");
        reorder("prop_pp", "r2", "21");
        mult("s2c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        perm("r5", "(12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24b", t_start_Q24b, print_level);

        // Q24c
        double t_start_Q24c = abs_time();
        reorder("x3c_v1", "r1", "124356");
        reorder("prop_hh", "r2", "21");
        mult("s2c+", "r2", "r3", 1);
        reorder("r3", "r4", "3412");
        mult("r1", "r4", "r5", 3);
        perm("r5", "(34)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24c", t_start_Q24c, print_level);

        // Q24d
        double t_start_Q24d = abs_time();
        reorder("x3c+_v1", "r1", "145236");
        reorder("s2c", "r2", "1342");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 3);
        reorder("r4", "r5", "1423");
        perm("r5", "(12)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24d", t_start_Q24d, print_level);

        // Q24e
        double t_start_Q24e = abs_time();
        reorder("x3c", "r1", "124356");
        reorder("s2c+", "r2", "3412");
        reorder("prop_pv", "r3", "21");
        mult("r1", "r2", "r4", 3);
        reorder("r4", "r5", "4123");
        mult("r3", "r5", "r6", 1);
        reorder("r6", "r7", "3412");
        perm("r7", "(34)");
        update(result_name, 0.5, "r7");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24e", t_start_Q24e, print_level);

        // Q24f
        double t_start_Q24f = abs_time();
        reorder("x3c+", "r1", "145623");
        mult("r1", "s2c", "r2", 3);
        reorder("r2", "r3", "4231");
        mult("prop_vp", "r3", "r4", 1);
        perm("r4", "(12)");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24f", t_start_Q24f, print_level);

        // Q25a
        double t_start_Q25a = abs_time();
        reorder("x2c+", "r1", "3412");
        reorder("x3c", "r2", "124563");
        mult("r2", "prop_ph", "r3", 2);
        mult("r3", "r1", "r4", 2);
        closed("r4", "r5");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q25a", t_start_Q25a, print_level);

        // Q25b
        double t_start_Q25b = abs_time();
        reorder("x3c+", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("x2c", "r2", "r3", 2);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q25b", t_start_Q25b, print_level);

        // Q26a
        double t_start_Q26a = abs_time();
        reorder("x3c+", "r1", "456123");
        mult("r1", "prop_pp", "r2", 1);
        mult("x3c", "r2", "r3", 4);
        closed("r3", "r4");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q26a", t_start_Q26a, print_level);

        // Q26b
        double t_start_Q26b = abs_time();
        reorder("prop_hh", "r1", "21");
        mult("x3c+", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("x3c", "r3", "r4", 4);
        closed("r4", "r5");
        update(result_name, -1.0 / 6.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q26b", t_start_Q26b, print_level);
    }

    /*
     * for the "connected" scheme only: connected diagrams ((T^+T)_conn O)_conn.
     * In fact, "overlap" (or metric) operator (T^+ T)_conn is contracted with
     * the property operator. The sign standing before this term is "minus".
     */
    if (disconnected == 0) {
        pos = get_stack_pos();
        sector_0h2p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
        reorder("TT_conn", "r1", "3412");
        mult("r1", "prop_vv", "r2", 1);
        reorder("r2", "r3", "3412");
        perm("r3", "(12)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}

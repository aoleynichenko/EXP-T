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

#include "finite_order_prop_0h3p.h"

#include <stdlib.h>

#include "cc_properies.h"
#include "engine.h"
#include "spinors.h"

#include "finite_order_prop.h"
#include "finite_order_overlap.h"


void finite_order_property_0h3p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "pppppp", "111111", "123456", NOT_PERM_UNIQUE, operator_symmetry);

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
        restrict_valence("t3c", "t3c_v123", "000111", 0);
        restrict_valence("t3c+", "t3c+_v123", "111000", 0);

        restrict_valence("s3c", "s3c_v1", "100100", 0);
        restrict_valence("s3c+", "s3c+_v1", "100100", 0);
        restrict_valence("s3c", "s3c_v2", "100010", 0);
        restrict_valence("s3c+", "s3c+_v2", "010100", 0);
        restrict_valence("s3c", "s3c_v12", "100110", 0);
        restrict_valence("s3c+", "s3c+_v12", "110100", 0);
        restrict_valence("s3c", "s3c_v23", "100011", 0);
        restrict_valence("s3c+", "s3c+_v23", "011100", 0);
        restrict_valence("s3c", "s3c_v123", "100111", 0);
        restrict_valence("s3c+", "s3c+_v123", "111100", 0);

        restrict_valence("x3c", "x3c_v1", "110100", 0);
        restrict_valence("x3c+", "x3c+_v1", "100110", 0);
        restrict_valence("x3c", "x3c_v3", "110001", 0);
        restrict_valence("x3c+", "x3c+_v3", "001110", 0);
        restrict_valence("x3c", "x3c_v12", "110110", 0);
        restrict_valence("x3c+", "x3c+_v12", "110110", 0);
        restrict_valence("x3c", "x3c_v13", "110101", 0);
        restrict_valence("x3c+", "x3c+_v13", "101110", 0);
        restrict_valence("x3c", "x3c_v123", "110111", 0);
        restrict_valence("x3c+", "x3c+_v123", "111110", 0);

        restrict_valence("z3c", "z3c_v1", "111100", 0);
        restrict_valence("z3c+", "z3c+_v1", "100111", 0);
        restrict_valence("z3c", "z3c_v12", "111110", 0);
        restrict_valence("z3c+", "z3c+_v12", "110111", 0);
    }

    dg_stack_pos_t pos = get_stack_pos();
    int print_level = cc_opts->print_level;

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    if (print_level >= CC_PRINT_HIGH) {
        diagram_stack_print();
        printf(" diagram                 time (sec)  max mem (mb)  max mem (gb)\n");
    }

    // linear terms
    // (no linear terms for the CCSD approximation)

    if (triples_enabled()) {

        // L1a
        double t_start_L1a = abs_time();
        reorder("prop_pv", "r1", "21");
        mult("z3c_v12", "r1", "r2", 1);
        perm("r2", "(6/45)");
        update(result_name, 1.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("L1a", t_start_L1a, print_level);

        // L1b
        double t_start_L1b = abs_time();
        reorder("z3c+_v12", "r1", "456123");
        mult("r1", "prop_vp", "r2", 1);
        diagram_stack_erase("r1");
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        perm("r3", "(3/12)");
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("L1b", t_start_L1b, print_level);

        // L2a
        double t_start_L2a = abs_time();
        reorder("x3c_v123", "r1", "456123");
        mult("r1", "prop_vh", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        perm("r3", "(3/12)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("L2a", t_start_L2a, print_level);

        // L2b
        double t_start_L2b = abs_time();
        reorder("prop_hv", "r1", "21");
        mult("x3c+_v123", "r1", "r2", 1);
        perm("r2", "(6/45)");
        update(result_name, -1.0, "r2");
        restore_stack_pos(pos);
        print_time_mem_usage("L2b", t_start_L2b, print_level);
    }

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // quadratic terms

    // Q1a
    double t_start_Q1a = abs_time();
    reorder("s2c_v12", "s2c_v12_r", "2143");
    reorder("s2c_v12_r", "r1", "2341");
    reorder("prop_ph", "r2", "21");
    mult("x2c_v1", "r2", "r3", 1);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "124356");
    perm("r5", "(3/12|4/56)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1a", t_start_Q1a, print_level);

    // Q1b
    double t_start_Q1b = abs_time();
    reorder("s2c+_v12", "r1", "2134");
    reorder("x2c+_v1", "r2", "1342");
    mult("r2", "prop_hp", "r3", 1);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "145236");
    perm("r5", "(1/23|6/45)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q1b", t_start_Q1b, print_level);

    // Q2a
    double t_start_Q2a = abs_time();
    reorder("s2c+_v12", "r1", "2134");
    reorder("s2c_v12", "r2", "1342");
    mult("r2", "prop_hh", "r3", 1);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "145236");
    perm("r5", "(1/23|6/45)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2a", t_start_Q2a, print_level);

    // Q2b
    double t_start_Q2b = abs_time();
    reorder("t2c_v12", "r1", "3421");
    reorder("s2c+_v12", "r2", "2134");
    mult("prop_vh", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "145236");
    perm("r5", "(1/23|6/45)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2b", t_start_Q2b, print_level);

    // Q2c
    double t_start_Q2c = abs_time();
    reorder("s2c_v12", "r1", "1342");
    reorder("prop_hv", "r2", "21");
    mult("r2", "t2c+_v12", "r3", 1);
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "156234");
    perm("r5", "(1/23|6/45)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2c", t_start_Q2c, print_level);

    // Q2d
    double t_start_Q2d = abs_time();
    reorder("s2c_v2", "r1", "1423");
    reorder("s2c+_v12", "r2", "2134");
    reorder("prop_pv", "r3", "21");
    mult("r1", "r3", "r4", 1);
    reorder("r4", "r5", "1243");
    mult("r5", "r2", "r6", 1);
    reorder("r6", "r7", "145326");
    perm("r7", "(1/23|456)");
    update(result_name, -1.0, "r7");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2d", t_start_Q2d, print_level);

    // Q2e
    double t_start_Q2e = abs_time();
    reorder("s2c+_v2", "s2c+_v2_r", "2143");
    reorder("s2c+_v2_r", "r1", "1432");
    reorder("s2c_v12", "r2", "1342");
    mult("prop_vp", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "154236");
    perm("r5", "(123|6/45)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q2e", t_start_Q2e, print_level);

    // Q3a
    double t_start_Q3a = abs_time();
    reorder("x2c+_v1", "x2c+_v2", "2143");
    reorder("x2c+_v2", "r2", "2341");
    reorder("s2c_v1", "r1", "1342");
    mult("r1", "prop_vh", "r3", 1);
    reorder("r3", "r4", "1243");
    mult("r4", "r2", "r5", 1);
    reorder("r5", "r6", "134256");
    diagram_stack_erase("r5");
    perm("r6", "(123|4/56)");
    update(result_name, -1.0, "r6");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3a", t_start_Q3a, print_level);

    // Q3b
    double t_start_Q3b = abs_time();
    reorder("s2c+_v1", "s2c+_v1_r", "2143");
    reorder("s2c+_v1_r", "r1", "2413");
    reorder("prop_hv", "r2", "21");
    mult("r2", "r1", "r3", 1);
    mult("x2c_v1", "r3", "r4", 1);
    reorder("r4", "r5", "125346");
    diagram_stack_erase("r4");
    perm("r5", "(3/12|456)");
    update(result_name, -1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q3b", t_start_Q3b, print_level);

    // Q4a
    double t_start_Q4a = abs_time();
    reorder("x2c+_v1", "x2c+_v2", "2143");
    reorder("x2c+_v2", "r1", "2341");
    reorder("x2c", "r2", "1243");
    reorder("prop_pv", "r3", "21");
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1243");
    mult("r5", "r1", "r6", 1);
    reorder("r6", "r7", "124356");
    diagram_stack_erase("r6");
    perm("r7", "(3/12|4/56)");
    update(result_name, 1.0, "r7");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4a", t_start_Q4a, print_level);

    // Q4b
    double t_start_Q4b = abs_time();
    reorder("x2c+", "r1", "3412");
    mult("prop_vp", "r1", "r2", 1);
    mult("x2c_v1", "r2", "r3", 1);
    reorder("r3", "r4", "124356");
    diagram_stack_erase("r3");
    perm("r4", "(3/12|4/56)");
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4b", t_start_Q4b, print_level);

    // Q4c
    double t_start_Q4c = abs_time();
    reorder("x2c+_v1", "x2c+_v2", "2143");
    reorder("x2c+_v2", "r1", "2341");
    reorder("prop_pp", "r2", "21");
    mult("x2c_v1", "r2", "r3", 1);
    mult("r3", "r1", "r4", 1);
    reorder("r4", "r5", "124356");
    diagram_stack_erase("r4");
    perm("r5", "(3/12|4/56)");
    update(result_name, 1.0, "r5");
    restore_stack_pos(pos);
    print_time_mem_usage("Q4c", t_start_Q4c, print_level);

    if (triples_enabled()) {

        // Q5a
        double t_start_Q5a = abs_time();
        reorder("s3c_v123", "r1", "456123");
        mult("r1", "prop_vh", "r2", 1);
        reorder("r2", "r3", "123465");
        diagram_stack_erase("r2");
        mult("r3", "t1c+_v", "r4", 1);
        reorder("r4", "r5", "465123");
        diagram_stack_erase("r4");
        perm("r5", "(123)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q5a", t_start_Q5a, print_level);

        // Q5b
        double t_start_Q5b = abs_time();
        reorder("prop_hv", "r1", "21");
        reorder("t1c_v", "r2", "21");
        mult("s3c+_v123", "r1", "r3", 1);
        reorder("r3", "r4", "123465");
        diagram_stack_erase("r3");
        mult("r4", "r2", "r5", 1);
        reorder("r5", "r6", "123465");
        diagram_stack_erase("r5");
        perm("r6", "(456)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q5b", t_start_Q5b, print_level);

        // Q6a
        double t_start_Q6a = abs_time();
        reorder("prop_pv", "r1", "21");
        mult("x3c_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("r3", "t1c+_v", "r4", 1);
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12|6/45)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q6a", t_start_Q6a, print_level);

        // Q6b
        double t_start_Q6b = abs_time();
        reorder("t1c_v", "r1", "21");
        mult("x3c+_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("r3", "prop_vp", "r4", 1);
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12|6/45)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q6b", t_start_Q6b, print_level);

        // Q6c
        double t_start_Q6c = abs_time();
        reorder("prop_hh", "r1", "21");
        reorder("x3c_v123", "r2", "456123");
        mult("t1c+_v", "r1", "r3", 1);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q6c", t_start_Q6c, print_level);

        // Q6d
        double t_start_Q6d = abs_time();
        reorder("t1c_v", "r1", "21");
        mult("r1", "prop_hh", "r2", 1);
        mult("x3c+_v123", "r2", "r3", 1);
        perm("r3", "(6/45)");
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q6d", t_start_Q6d, print_level);

        // Q7a
        double t_start_Q7a = abs_time();
        reorder("x3c_v123", "r1", "456123");
        reorder("prop_ph", "r2", "21");
        mult("s1c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 1);
        diagram_stack_erase("r1");
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7a", t_start_Q7a, print_level);

        // Q7b
        double t_start_Q7b = abs_time();
        reorder("s1c+", "r1", "21");
        mult("r1", "prop_hp", "r2", 1);
        mult("x3c+_v123", "r2", "r3", 1);
        perm("r3", "(6/45)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7b", t_start_Q7b, print_level);

        // Q7c
        double t_start_Q7c = abs_time();
        reorder("s1c+", "r1", "21");
        mult("x3c_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("r3", "prop_vh", "r4", 1);
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12|6/45)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7c", t_start_Q7c, print_level);

        // Q7d
        double t_start_Q7d = abs_time();
        reorder("prop_hv", "r1", "21");
        mult("x3c+_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("r3", "s1c", "r4", 1);
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12|6/45)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q7d", t_start_Q7d, print_level);

        // Q8a
        double t_start_Q8a = abs_time();
        reorder("t1c_v", "r1", "21");
        mult("r1", "prop_ph", "r2", 1);
        mult("z3c_v12", "r2", "r3", 1);
        perm("r3", "(6/45)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q8a", t_start_Q8a, print_level);

        // Q8b
        double t_start_Q8b = abs_time();
        reorder("prop_hp", "r1", "21");
        reorder("z3c+_v12", "r2", "456123");
        mult("t1c+_v", "r1", "r3", 1);
        mult("r2", "r3", "r4", 1);
        diagram_stack_erase("r2");
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q8b", t_start_Q8b, print_level);

        // Q9a
        double t_start_Q9a = abs_time();
        reorder("prop_pv", "r1", "21");
        reorder("s1c+", "r2", "21");
        mult("z3c_v1", "r1", "r3", 1);
        reorder("r3", "r4", "123465");
        diagram_stack_erase("r3");
        mult("r4", "r2", "r5", 1);
        reorder("r5", "r6", "123465");
        diagram_stack_erase("r5");
        perm("r6", "(456)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9a", t_start_Q9a, print_level);

        // Q9b
        double t_start_Q9b = abs_time();
        reorder("z3c+_v1", "r1", "456123");
        mult("r1", "prop_vp", "r2", 1);
        diagram_stack_erase("r1");
        reorder("r2", "r3", "123465");
        diagram_stack_erase("r2");
        mult("r3", "s1c", "r4", 1);
        reorder("r4", "r5", "465123");
        diagram_stack_erase("r4");
        perm("r5", "(123)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9b", t_start_Q9b, print_level);

        // Q9c
        double t_start_Q9c = abs_time();
        reorder("s1c+", "r1", "21");
        mult("r1", "prop_pp", "r2", 1);
        mult("z3c_v12", "r2", "r3", 1);
        perm("r3", "(6/45)");
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9c", t_start_Q9c, print_level);

        // Q9d
        double t_start_Q9d = abs_time();
        reorder("z3c+_v12", "r1", "456123");
        reorder("prop_pp", "r2", "21");
        mult("s1c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 1);
        diagram_stack_erase("r1");
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q9d", t_start_Q9d, print_level);

        // Q10a
        double t_start_Q10a = abs_time();
        reorder("t3c_v123", "r1", "456123");
        mult("r1", "prop_vh", "r2", 1);
        diagram_stack_erase("r1");
        reorder("r2", "r3", "123645");
        diagram_stack_erase("r2");
        mult("r3", "t2c+_v12", "r4", 2);
        reorder("r4", "r5", "564123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10a", t_start_Q10a, print_level);

        // Q10b
        double t_start_Q10b = abs_time();
        reorder("prop_hv", "r1", "21");
        reorder("t2c_v12", "r2", "3412");
        mult("t3c+_v123", "r1", "r3", 1);
        diagram_stack_erase("r1");
        reorder("r3", "r4", "123645");
        diagram_stack_erase("r3");
        mult("r4", "r2", "r5", 2);
        diagram_stack_erase("r4");
        reorder("r5", "r6", "123564");
        diagram_stack_erase("r5");
        perm("r6", "(6/45)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10b", t_start_Q10b, print_level);

        // Q11a
        double t_start_Q11a = abs_time();
        reorder("prop_pv", "r1", "21");
        mult("s3c_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("r3", "t2c+_v12", "r4", 2);
        diagram_stack_erase("r3");
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(1/23|6/45)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11a", t_start_Q11a, print_level);

        // Q11b
        double t_start_Q11b = abs_time();
        reorder("t2c_v12", "r1", "3412");
        mult("s3c+_v12", "r1", "r2", 2);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("r3", "prop_vp", "r4", 1);
        diagram_stack_erase("r3");
        reorder("r4", "r5", "456123");
        diagram_stack_erase("r4");
        perm("r5", "(3/12|4/56)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11b", t_start_Q11b, print_level);

        // Q11c
        double t_start_Q11c = abs_time();
        reorder("s3c_v123", "r1", "456123");
        mult("r1", "t2c+_v1", "r2", 2);
        diagram_stack_erase("r1");
        mult("r2", "prop_vp", "r3", 1);
        reorder("r3", "r4", "456123");
        diagram_stack_erase("r3");
        perm("r4", "(123)");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11c", t_start_Q11c, print_level);

        // Q11d
        double t_start_Q11d = abs_time();
        reorder("t2c_v1", "r1", "3412");
        reorder("prop_pv", "r2", "21");
        mult("s3c+_v123", "r1", "r3", 2);
        mult("r3", "r2", "r4", 1);
        diagram_stack_erase("r3");
        perm("r4", "(456)");
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11d", t_start_Q11d, print_level);

        // Q11e
        double t_start_Q11e = abs_time();
        reorder("s3c_v123", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        diagram_stack_erase("r1");
        mult("r2", "t2c+_v12", "r3", 2);
        reorder("r3", "r4", "456123");
        diagram_stack_erase("r3");
        perm("r4", "(1/23)");
        update(result_name, -1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11e", t_start_Q11e, print_level);

        // Q11f
        double t_start_Q11f = abs_time();
        reorder("t2c_v12", "r1", "3412");
        mult("r1", "prop_hh", "r2", 1);
        mult("s3c+_v123", "r2", "r3", 2);
        perm("r3", "(4/56)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11f", t_start_Q11f, print_level);

        // Q12a
        double t_start_Q12a = abs_time();
        reorder("s3c_v12", "r1", "145263");
        reorder("s2c+_v2", "s2c+_v2_r", "2143");
        reorder("s2c+_v2_r", "r2", "1432");
        mult("r1", "prop_vh", "r3", 1);
        reorder("r3", "r4", "123645");
        mult("r4", "r2", "r5", 2);
        reorder("r5", "r6", "145236");
        perm("r6", "(123|6/45)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12a", t_start_Q12a, print_level);

        // Q12b
        double t_start_Q12b = abs_time();
        reorder("prop_hv", "r1", "21");
        reorder("s2c_v2", "s2c_v2_r", "2143");
        reorder("s2c_v2_r", "r2", "2314");
        mult("s3c+_v12", "r1", "r3", 1);
        reorder("r3", "r4", "124653");
        mult("r4", "r2", "r5", 2);
        reorder("r5", "r6", "125364");
        perm("r6", "(3/12|456)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12b", t_start_Q12b, print_level);

        // Q13a
        double t_start_Q13a = abs_time();
        reorder("t2c+_v1", "r1", "1342");
        reorder("x3c_v123", "r2", "456123");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "456123");
        perm("r5", "(3/12)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13a", t_start_Q13a, print_level);

        // Q13b
        double t_start_Q13b = abs_time();
        reorder("t2c_v1", "r1", "3142");
        mult("r1", "prop_ph", "r2", 2);
        mult("x3c+_v123", "r2", "r3", 1);
        perm("r3", "(6/45)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13b", t_start_Q13b, print_level);

        // Q13c
        double t_start_Q13c = abs_time();
        reorder("prop_hv", "r1", "21");
        reorder("x3c_v12", "r2", "124536");
        mult("t2c+_v2", "r1", "r3", 1);
        reorder("r3", "r4", "2431");
        mult("r2", "r4", "r5", 2);
        reorder("r5", "r6", "125346");
        perm("r6", "(3/12|6/45)");
        update(result_name, -1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13c", t_start_Q13c, print_level);

        // Q13d
        double t_start_Q13d = abs_time();
        reorder("x3c+_v12", "r1", "124563");
        reorder("t2c_v2", "r2", "4132");
        mult("prop_vh", "r2", "r3", 1);
        mult("r1", "r3", "r4", 2);
        reorder("r4", "r5", "125346");
        perm("r5", "(3/12|6/45)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13d", t_start_Q13d, print_level);

        // Q14a
        double t_start_Q14a = abs_time();
        reorder("prop_pv", "r1", "21");
        reorder("s2c+_v2", "r2", "2341");
        mult("x3c_v1", "r1", "r3", 1);
        reorder("r3", "r4", "124635");
        mult("r4", "r2", "r5", 2);
        reorder("r5", "r6", "125364");
        perm("r6", "(3/12|456)");
        update(result_name, -1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14a", t_start_Q14a, print_level);

        // Q14b
        double t_start_Q14b = abs_time();
        reorder("x3c+_v1", "r1", "145362");
        reorder("s2c_v2", "r2", "1423");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "123564");
        mult("r4", "prop_vp", "r5", 1);
        reorder("r5", "r6", "146235");
        perm("r6", "(123|6/45)");
        update(result_name, -1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14b", t_start_Q14b, print_level);

        // Q14c
        double t_start_Q14c = abs_time();
        reorder("s2c+", "s2c+_r", "2143");
        reorder("s2c+_r", "r1", "4312");
        reorder("x3c_v12", "r2", "124536");
        mult("prop_vp", "r1", "r3", 1);
        mult("r2", "r3", "r4", 2);
        reorder("r4", "r5", "125346");
        perm("r5", "(3/12|6/45)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14c", t_start_Q14c, print_level);

        // Q14d
        double t_start_Q14d = abs_time();
        reorder("x3c+_v12", "r1", "124563");
        reorder("s2c", "r2", "1243");
        reorder("prop_pv", "r3", "21");
        mult("r3", "r2", "r4", 1);
        mult("r1", "r4", "r5", 2);
        reorder("r5", "r6", "126345");
        perm("r6", "(3/12|6/45)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14d", t_start_Q14d, print_level);

        // Q14e
        double t_start_Q14e = abs_time();
        reorder("x3c_v13", "r1", "124635");
        reorder("prop_hh", "r2", "21");
        mult("s2c+_v2", "r2", "r3", 1);
        reorder("r3", "r4", "2341");
        mult("r1", "r4", "r5", 2);
        reorder("r5", "r6", "125364");
        perm("r6", "(3/12|5/46)");
        update(result_name, 1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14e", t_start_Q14e, print_level);

        // Q14f
        double t_start_Q14f = abs_time();
        reorder("x3c+_v13", "r1", "134526");
        reorder("s2c_v2", "r2", "1432");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 2);
        reorder("r4", "r5", "152346");
        perm("r5", "(2/13|6/45)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14f", t_start_Q14f, print_level);

        // Q14g
        double t_start_Q14g = abs_time();
        reorder("s2c+_v2", "r1", "2341");
        reorder("x3c_v13", "r2", "124635");
        mult("r1", "prop_pp", "r3", 1);
        mult("r2", "r3", "r4", 2);
        reorder("r4", "r5", "125364");
        perm("r5", "(3/12|5/46)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14g", t_start_Q14g, print_level);

        // Q14h
        double t_start_Q14h = abs_time();
        reorder("s2c_v2", "r1", "1423");
        reorder("prop_pp", "r2", "21");
        reorder("x3c+_v13", "r3", "134562");
        mult("r1", "r2", "r4", 1);
        mult("r3", "r4", "r5", 2);
        reorder("r5", "r6", "152346");
        perm("r6", "(2/13|6/45)");
        update(result_name, -1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14h", t_start_Q14h, print_level);

        // Q15a
        double t_start_Q15a = abs_time();
        reorder("x3c_v3", "r1", "126345");
        reorder("x2c+", "r2", "3412");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "123564");
        mult("r4", "prop_vh", "r5", 1);
        reorder("r5", "r6", "126453");
        perm("r6", "(3/12|6/45)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15a", t_start_Q15a, print_level);

        // Q15b
        double t_start_Q15b = abs_time();
        reorder("prop_hv", "r1", "21");
        mult("x3c+_v3", "r1", "r2", 1);
        reorder("r2", "r3", "345612");
        mult("r3", "x2c", "r4", 2);
        reorder("r4", "r5", "561234");
        perm("r5", "(3/12|6/45)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15b", t_start_Q15b, print_level);

        // Q16a
        double t_start_Q16a = abs_time();
        reorder("prop_hv", "r0", "21");
        mult("s2c+", "r0", "r1", 1);
        reorder("r1", "r2", "3412");
        mult("z3c_v1", "r2", "r3", 2);
        perm("r3", "(456)");
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16a", t_start_Q16a, print_level);

        // Q16b
        double t_start_Q16b = abs_time();
        reorder("z3c+_v1", "r0", "456123");
        reorder("s2c", "r1", "3412");
        mult("r1", "prop_vh", "r2", 1);
        reorder("r2", "r3", "3412");
        mult("r0", "r3", "r4", 2);
        reorder("r4", "r5", "456123");
        perm("r5", "(123)");
        update(result_name, -0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16b", t_start_Q16b, print_level);

        // Q16c
        double t_start_Q16c = abs_time();
        reorder("s2c+", "r1", "3142");
        mult("r1", "prop_hp", "r2", 2);
        mult("z3c_v12", "r2", "r3", 1);
        perm("r3", "(6/45)");
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16c", t_start_Q16c, print_level);

        // Q16d
        double t_start_Q16d = abs_time();
        reorder("s2c", "r1", "1342");
        reorder("z3c+_v12", "r2", "456123");
        mult("r1", "prop_ph", "r3", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "456123");
        perm("r5", "(3/12)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16d", t_start_Q16d, print_level);

        // Q17a
        double t_start_Q17a = abs_time();
        reorder("x2c+", "r1", "3412");
        reorder("prop_pv", "r2", "21");
        mult("z3c", "r2", "r3", 1);
        reorder("r3", "r4", "123645");
        mult("r4", "r1", "r5", 2);
        reorder("r5", "r6", "123564");
        perm("r6", "(6/45)");
        update(result_name, 0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17a", t_start_Q17a, print_level);

        // Q17b
        double t_start_Q17b = abs_time();
        reorder("z3c+", "r1", "456123");
        mult("r1", "prop_vp", "r2", 1);
        reorder("r2", "r3", "123645");
        mult("r3", "x2c", "r4", 2);
        reorder("r4", "r5", "564123");
        perm("r5", "(3/12)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17b", t_start_Q17b, print_level);

        // Q17c
        double t_start_Q17c = abs_time();
        reorder("prop_pp", "r1", "21");
        reorder("x2c+", "r2", "3412");
        mult("z3c_v1", "r1", "r3", 1);
        mult("r3", "r2", "r4", 2);
        perm("r4", "(4/56)");
        update(result_name, 1.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17c", t_start_Q17c, print_level);

        // Q17d
        double t_start_Q17d = abs_time();
        reorder("prop_pp", "r1", "21");
        reorder("z3c+_v1", "r2", "456123");
        mult("x2c", "r1", "r3", 1);
        mult("r2", "r3", "r4", 2);
        reorder("r4", "r5", "456123");
        perm("r5", "(1/23)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q17d", t_start_Q17d, print_level);

        // Q18a
        double t_start_Q18a = abs_time();
        reorder("prop_pv", "r1", "21");
        mult("t3c_v12", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("t3c+_v123", "r3", "r4", 3);
        perm("r4", "(6/45)");
        update(result_name, -1.0/6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18a", t_start_Q18a, print_level);

        // Q18b
        double t_start_Q18b = abs_time();
        reorder("t3c_v123", "r1", "456123");
        mult("r1", "t3c+_v12", "r2", 3);
        mult("r2", "prop_vp", "r3", 1);
        reorder("r3", "r4", "456123");
        perm("r4", "(3/12)");
        update(result_name, -1.0/6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18b", t_start_Q18b, print_level);

        // Q18c
        double t_start_Q18c = abs_time();
        reorder("t3c_v123", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        mult("t3c+_v123", "r2", "r3", 3);
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q18c", t_start_Q18c, print_level);

        // Q19a
        double t_start_Q19a = abs_time();
        reorder("t3c+_v12", "r1", "124563");
        reorder("s3c_v123", "r2", "456123");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 2);
        reorder("r4", "r5", "456123");
        perm("r5", "(1/23)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q19a", t_start_Q19a, print_level);

        // Q19b
        double t_start_Q19b = abs_time();
        reorder("t3c_v12", "r1", "451263");
        mult("r1", "prop_ph", "r2", 2);
        mult("s3c+_v123", "r2", "r3", 2);
        perm("r3", "(4/56)");
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q19b", t_start_Q19b, print_level);

        // Q20a
        double t_start_Q20a = abs_time();
        reorder("s3c+_v23", "r1", "234561");
        reorder("s3c_v23", "r2", "561234");
        reorder("prop_pp", "r3", "21");
        mult("r2", "r3", "r4", 1);
        mult("r4", "r1", "r5", 3);
        reorder("r5", "r6", "345612");
        perm("r6", "(1/23|4/56)");
        update(result_name, 0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q20a", t_start_Q20a, print_level);

        // Q20b
        double t_start_Q20b = abs_time();
        reorder("s3c+_v23", "r1", "234156");
        reorder("s3c_v23", "r2", "156423");
        mult("r2", "prop_hh", "r3", 1);
        mult("r3", "r1", "r4", 3);
        reorder("r4", "r5", "145623");
        perm("r5", "(1/23|4/56)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q20b", t_start_Q20b, print_level);

        // Q21a
        double t_start_Q21a = abs_time();
        reorder("x3c+_v1", "x3c+_v1_321", "321654");
        reorder("x3c+_v1_321", "r1", "356421");
        reorder("s3c_v1", "r2", "143265");
        mult("r1", "r2", "r3", 3);
        mult("r3", "prop_vh", "r4", 1);
        reorder("r4", "r5", "461523");
        perm("r5", "(123|4/56)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21a", t_start_Q21a, print_level);

        // Q21b
        double t_start_Q21b = abs_time();
        reorder("x3c_v1", "x3c_v1_321", "321654");
        reorder("x3c_v1_321", "r1", "236145");
        reorder("prop_hv", "r2", "21");
        mult("s3c+_v1", "r2", "r3", 1);
        reorder("r3", "r4", "146523");
        mult("r4", "r1", "r5", 3);
        reorder("r5", "r6", "145236");
        perm("r6", "(1/23|456)");
        update(result_name, -1.0, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21b", t_start_Q21b, print_level);

        // Q21c
        double t_start_Q21c = abs_time();
        reorder("x3c+_v12", "x3c+_v12_321", "321654");
        reorder("x3c+_v12_321", "r1", "235641");
        reorder("s3c_v1", "r2", "143652");
        mult("r2", "prop_ph", "r3", 2);
        mult("r1", "r3", "r4", 2);
        reorder("r4", "r5", "512634");
        perm("r5", "(1/23|4/56)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21c", t_start_Q21c, print_level);

        // Q21d
        double t_start_Q21d = abs_time();
        reorder("s3c+_v1", "r1", "146352");
        reorder("x3c_v12", "x3c_v12_321", "321654");
        reorder("x3c_v12_321", "r2", "235614");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 2);
        reorder("r4", "r5", "512634");
        perm("r5", "(1/23|4/56)");
        update(result_name, 1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q21d", t_start_Q21d, print_level);

        // Q22a
        double t_start_Q22a = abs_time();
        reorder("prop_pv", "r1", "21");
        reorder("x3c+_v3", "r2", "345612");
        mult("x3c", "r1", "r3", 1);
        reorder("r3", "r4", "612345");
        mult("r4", "r2", "r5", 3);
        reorder("r5", "r6", "234561");
        perm("r6", "(3/12|6/45)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q22a", t_start_Q22a, print_level);

        // Q22b
        double t_start_Q22b = abs_time();
        reorder("x3c+", "r1", "456123");
        reorder("x3c_v3", "r2", "612345");
        mult("r1", "prop_vp", "r3", 1);
        reorder("r3", "r4", "612345");
        mult("r2", "r4", "r5", 3);
        reorder("r5", "r6", "234561");
        perm("r6", "(3/12|6/45)");
        update(result_name, -0.5, "r6");
        restore_stack_pos(pos);
        print_time_mem_usage("Q22b", t_start_Q22b, print_level);

        // Q22c
        double t_start_Q22c = abs_time();
        reorder("x3c+_v3", "r1", "345621");
        reorder("x3c_v3", "r2", "126354");
        mult("r1", "prop_pp", "r3", 1);
        mult("r2", "r3", "r4", 3);
        reorder("r4", "r5", "124563");
        perm("r5", "(3/12|6/45)");
        update(result_name, -1.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q22c", t_start_Q22c, print_level);

        // Q22d
        double t_start_Q22d = abs_time();
        reorder("x3c+_v3", "r1", "345126");
        reorder("x3c_v3", "r2", "126453");
        mult("r2", "prop_hh", "r3", 1);
        mult("r3", "r1", "r4", 3);
        reorder("r4", "r5", "124563");
        perm("r5", "(3/12|6/45)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q22d", t_start_Q22d, print_level);

        // Q23a
        double t_start_Q23a = abs_time();
        reorder("x3c", "r2", "456123");
        mult("r2", "prop_vh", "r3", 1);
        diagram_stack_erase("r2");
        reorder("r3", "r4", "456123");
        diagram_stack_erase("r3");
        reorder("z3c+", "r1", "456123");
        mult("r4", "r1", "r5", 3);
        diagram_stack_erase("r1");
        diagram_stack_erase("r4");
        perm("r5", "(3/12)");
        update(result_name, -1.0/6.0, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q23a", t_start_Q23a, print_level);

        // Q23b
        double t_start_Q23b = abs_time();
        reorder("prop_hv", "r1", "21");
        mult("x3c+", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        diagram_stack_erase("r2");
        mult("z3c", "r3", "r4", 3);
        diagram_stack_erase("r3");
        perm("r4", "(6/45)");
        update(result_name, -1.0/6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q23b", t_start_Q23b, print_level);

        // Q23c
        double t_start_Q23c = abs_time();
        reorder("x3c", "r1", "124563");
        reorder("z3c+_v1", "r2", "145623");
        mult("r1", "prop_ph", "r3", 2);
        mult("r2", "r3", "r4", 2);
        reorder("r4", "r5", "156234");
        perm("r5", "(1/23)");
        update(result_name, 0.5, "r5");
        restore_stack_pos(pos);
        print_time_mem_usage("Q23c", t_start_Q23c, print_level);

        // Q23d
        double t_start_Q23d = abs_time();
        reorder("x3c+", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("z3c_v1", "r2", "r3", 2);
        perm("r3", "(4/56)");
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q23d", t_start_Q23d, print_level);

        // Q24
        double t_start_Q24 = abs_time();
        reorder("z3c+", "r1", "456123");
        mult("r1", "prop_pp", "r2", 1);
        mult("z3c", "r2", "r3", 3);
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q24", t_start_Q24, print_level);
    }

    /*
     * for the "connected" scheme only: connected diagrams ((T^+T)_conn O)_conn.
     * In fact, "overlap" (or metric) operator (T^+ T)_conn is contracted with
     * the property operator. The sign standing before this term is "minus".
     */
    if (disconnected == 0) {
        pos = get_stack_pos();

        // renormalization term
        double t_start_Q_renorm = abs_time();
        sector_0h3p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
        reorder("TT_conn", "r1", "456123");
        mult("r1", "prop_vv", "r2", 1);
        reorder("r2", "r3", "456123");
        perm("r3", "(3/12)");
        update(result_name, -1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q_renorm", t_start_Q_renorm, print_level);
    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}

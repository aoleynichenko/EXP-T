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

#include "finite_order_prop_1h0p.h"

#include "cc_properies.h"
#include "engine.h"
#include <stdlib.h>

#include "finite_order_overlap.h"

void print_time_mem_usage(char *diagram_name, double t_start, int print_level);

/*
 * TODO:
 * <> connected scheme
 * <> triples contributions
 */

void finite_order_property_1h0p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "hh", "11", "12", NOT_PERM_UNIQUE, operator_symmetry);

    /*
     * prepare diagrams
     */
    restrict_valence("t1c", "t1c_g", "10", 0);
    restrict_valence("t1c+", "t1c+_g", "01", 0);

    restrict_valence("t2c", "t2c_g1", "1000", 0);
    restrict_valence("t2c+", "t2c+_g1", "0010", 0);

    restrict_valence("h2c", "h2c_g1", "1010", 0);
    restrict_valence("h2c+", "h2c+_g1", "1010", 0);

    restrict_valence("prop_hp", "prop_gp", "10", 0);
    restrict_valence("prop_ph", "prop_pg", "01", 0);
    restrict_valence("prop_hh", "prop_gh", "10", 0);
    restrict_valence("prop_hh", "prop_hg", "01", 0);
    restrict_valence("prop_hh", "prop_gg", "11", 0);

    if (triples_enabled()) {
        restrict_valence("t3c", "t3c_g1", "100000", 0);
        restrict_valence("t3c+", "t3c+_g1", "000100", 0);
        restrict_valence("h3c", "h3c_g1", "100100", 0);
        restrict_valence("h3c+", "h3c+_g1", "100100", 0);
    }


    dg_stack_pos_t pos = get_stack_pos();
    int print_level = cc_opts->print_level;

    // O1
    copy("prop_gg", "r1");
    closed("r1", "r2");
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    // L1a
    reorder("prop_pg", "r1", "21");
    mult("t1c_g", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L1b
    reorder("t1c+_g", "r1", "21");
    mult("prop_gp", "r1", "r2", 1);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L2a
    reorder("h1c", "r1", "21");
    mult("prop_gh", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L2b
    reorder("prop_hg", "r1", "21");
    mult("h1c+", "r1", "r2", 1);
    update(result_name, -1.0, "r2");
    restore_stack_pos(pos);

    // L3a
    reorder("h2c_g1", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    // L3b
    reorder("h2c+_g1", "r1", "1342");
    mult("r1", "prop_hp", "r2", 2);
    update(result_name, 1.0, "r2");
    restore_stack_pos(pos);

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return;
    }

    // Q1a
    reorder("t1c+_g", "r1", "21");
    mult("r1", "prop_pp", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q1b
    reorder("t1c+_g", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("prop_gh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q1c
    reorder("prop_hg", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q1d - disconnected term
    if (disconnected) {
        double t1_t1 = scalar_product("C", "N", "t1c", "t1c");
        copy("prop_hh", "r0");
        closed("r0", "r1");
        update(result_name, t1_t1, "r1");
        restore_stack_pos(pos);
    }

    // Q2a
    reorder("t2c_g1", "r1", "1342");
    reorder("t1c+_g", "r2", "21");
    mult("r1", "prop_ph", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q2b
    reorder("t2c+_g1", "r1", "3142");
    mult("r1", "prop_hp", "r2", 2);
    mult("t1c_g", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q2c
    reorder("t2c_g1", "r1", "1342");
    reorder("prop_pg", "r2", "21");
    mult("r1", "t1c+", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q2d
    reorder("t2c+_g1", "r1", "3142");
    mult("r1", "t1c", "r2", 2);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q3a - disconnected term - to be removed
    if (disconnected) {
        double t2_t2 = 0.25 * scalar_product("C", "N", "t2c", "t2c");
        copy("prop_hh", "r0");
        closed("r0", "r1");
        update(result_name, t2_t2, "r1");
        restore_stack_pos(pos);
    }

    // Q3a
    reorder("prop_hg", "r1", "21");
    reorder("t2c+", "r2", "4123");
    mult("r1", "r2", "r3", 1);
    mult("t2c_g1", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q3b
    reorder("t2c+_g1", "r1", "3412");
    reorder("t2c", "r2", "2341");
    mult("prop_gh", "r2", "r3", 1);
    mult("r3", "r1", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q3c
    reorder("prop_hh", "r1", "21");
    mult("t2c+_g1", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("t2c_g1", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q3d
    reorder("t2c+_g1", "r1", "3412");
    mult("r1", "prop_pp", "r2", 1);
    mult("t2c_g1", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q4a
    reorder("prop_ph", "r1", "21");
    reorder("h1c", "r2", "21");
    mult("t1c_g", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q4b
    reorder("t1c+_g", "r1", "21");
    mult("r1", "prop_hp", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4c
    reorder("prop_pg", "r1", "21");
    mult("r1", "t1c", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q4d
    reorder("h1c", "r1", "21");
    mult("r1", "t1c+", "r2", 1);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q5a
    reorder("h2c_g1", "r1", "1342");
    reorder("prop_hh", "r2", "21");
    mult("t1c+", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5b
    reorder("h2c+_g1", "r1", "1342");
    reorder("t1c", "r2", "21");
    mult("prop_hh", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q5c
    reorder("t1c+", "r1", "21");
    reorder("h2c_g1", "r2", "1324");
    mult("r1", "prop_pp", "r3", 1);
    mult("r2", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q5d
    reorder("h2c+_g1", "r1", "1324");
    reorder("prop_pp", "r2", "21");
    mult("r2", "t1c", "r3", 1);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6a
    reorder("h2c", "r1", "3412");
    mult("r1", "t2c+", "r2", 3);
    mult("prop_gp", "r2", "r3", 1);
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);

    // Q6b
    reorder("t2c", "r1", "3412");
    reorder("prop_pg", "r2", "21");
    mult("h2c+", "r1", "r3", 3);
    mult("r3", "r2", "r4", 1);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q6c
    reorder("h2c_g1", "r1", "1324");
    reorder("t2c+", "r2", "3142");
    mult("r2", "prop_hp", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q6d
    reorder("h2c+_g1", "r1", "1342");
    reorder("t2c", "r2", "1342");
    mult("r2", "prop_ph", "r3", 2);
    mult("r1", "r3", "r4", 2);
    update(result_name, 1.0, "r4");
    restore_stack_pos(pos);

    // Q7
    reorder("h1c", "r1", "21");
    mult("r1", "prop_hh", "r2", 1);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    // Q8a
    reorder("h2c", "r1", "3142");
    mult("r1", "prop_ph", "r2", 2);
    mult("h1c+", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q8b
    reorder("h2c+", "r1", "1342");
    reorder("h1c", "r2", "21");
    mult("r1", "prop_hp", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q8c
    reorder("h2c", "r1", "3142");
    mult("r1", "t1c+", "r2", 2);
    mult("prop_gh", "r2", "r3", 1);
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // Q8d
    reorder("h2c+", "r1", "1342");
    reorder("prop_hg", "r2", "21");
    mult("r1", "t1c", "r3", 2);
    mult("r3", "r2", "r4", 1);
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // Q9a
    reorder("prop_pp", "r1", "21");
    mult("h2c", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    mult("h2c+", "r3", "r4", 3);
    update(result_name, -0.5, "r4");
    restore_stack_pos(pos);

    // Q9b
    reorder("h2c", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    mult("h2c+", "r2", "r3", 3);
    update(result_name, 1.0, "r3");
    restore_stack_pos(pos);

    if (triples_enabled()) {

        // Q10a
        double t_start_Q10a = abs_time();
        reorder("t3c_g1", "r1", "145623");
        reorder("prop_pg", "r2", "21");
        mult("r1", "t2c+", "r3", 4);
        mult("r3", "r2", "r4", 1);
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10a", t_start_Q10a, print_level);

        // Q10b
        double t_start_Q10b = abs_time();
        reorder("t3c+_g1", "r1", "415623");
        mult("r1", "t2c", "r2", 4);
        mult("prop_gp", "r2", "r3", 1);
        update(result_name, 0.25, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10b", t_start_Q10b, print_level);

        // Q10c
        double t_start_Q10c = abs_time();
        reorder("t3c_g1", "r1", "124563");
        reorder("t2c+_g1", "r2", "3412");
        mult("r1", "prop_ph", "r3", 2);
        mult("r3", "r2", "r4", 3);
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10c", t_start_Q10c, print_level);

        // Q10d
        double t_start_Q10d = abs_time();
        reorder("t3c+_g1", "r1", "451263");
        mult("r1", "prop_hp", "r2", 2);
        mult("t2c_g1", "r2", "r3", 3);
        update(result_name, 0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q10d", t_start_Q10d, print_level);

        // Q11a
        double t_start_Q11a = abs_time();
        reorder("t3c+_g1", "r1", "456123");
        mult("r1", "t3c", "r2", 5);
        mult("prop_gh", "r2", "r3", 1);
        update(result_name, -1.0/12.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11a", t_start_Q11a, print_level);

        // Q11b
        double t_start_Q11b = abs_time();
        reorder("t3c+", "r1", "456123");
        reorder("prop_hg", "r2", "21");
        mult("t3c_g1", "r1", "r3", 5);
        mult("r2", "r3", "r4", 1);
        update(result_name, -1.0/12.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11b", t_start_Q11b, print_level);

        // Q11c
        double t_start_Q11c = abs_time();
        reorder("t3c+_g1", "r1", "456123");
        reorder("prop_pp", "r2", "21");
        mult("t3c_g1", "r2", "r3", 1);
        mult("r3", "r1", "r4", 5);
        update(result_name, 1.0/4.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11c", t_start_Q11c, print_level);

        // Q11d
        double t_start_Q11d = abs_time();
        reorder("prop_hh", "r1", "21");
        mult("t3c+_g1", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("t3c_g1", "r3", "r4", 5);
        update(result_name, -1.0/6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q11d", t_start_Q11d, print_level);

        // Q12a
        double t_start_Q12a = abs_time();
        reorder("h3c_g1", "r1", "145263");
        mult("r1", "prop_ph", "r2", 2);
        mult("r2", "t1c+", "r3", 2);
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12a", t_start_Q12a, print_level);

        // Q12b
        double t_start_Q12b = abs_time();
        reorder("h3c+_g1", "r1", "145263");
        mult("r1", "prop_hp", "r2", 2);
        mult("r2", "t1c", "r3", 2);
        update(result_name, 1.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q12b", t_start_Q12b, print_level);

        // Q13a
        double t_start_Q13a = abs_time();
        reorder("h3c", "r1", "415623");
        mult("r1", "t2c+", "r2", 4);
        mult("prop_gh", "r2", "r3", 1);
        update(result_name, -0.25, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13a", t_start_Q13a, print_level);

        // Q13b
        double t_start_Q13b = abs_time();
        reorder("h3c+", "r1", "145623");
        reorder("prop_hg", "r2", "21");
        mult("r1", "t2c", "r3", 4);
        mult("r3", "r2", "r4", 1);
        update(result_name, -0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13b", t_start_Q13b, print_level);

        // Q13c
        double t_start_Q13c = abs_time();
        reorder("h3c_g1", "r1", "142356");
        reorder("t2c+", "r2", "3412");
        mult("r2", "prop_pp", "r3", 1);
        mult("r1", "r3", "r4", 4);
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13c", t_start_Q13c, print_level);

        // Q13d
        double t_start_Q13d = abs_time();
        reorder("h3c+_g1", "r1", "145623");
        reorder("prop_pp", "r2", "21");
        mult("t2c", "r2", "r3", 1);
        mult("r1", "r3", "r4", 4);
        update(result_name, 0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13d", t_start_Q13d, print_level);

        // Q13e
        double t_start_Q13e = abs_time();
        reorder("prop_hh", "r1", "21");
        reorder("h3c_g1", "r2", "145623");
        mult("t2c+", "r1", "r3", 1);
        mult("r2", "r3", "r4", 4);
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13e", t_start_Q13e, print_level);

        // Q13f
        double t_start_Q13f = abs_time();
        reorder("h3c+_g1", "r1", "142356");
        reorder("t2c", "r2", "3412");
        mult("r2", "prop_hh", "r3", 1);
        mult("r1", "r3", "r4", 4);
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q13f", t_start_Q13f, print_level);

        // Q14a
        double t_start_Q14a = abs_time();
        reorder("t3c+", "r1", "451263");
        reorder("h3c_g1", "r2", "142356");
        mult("r1", "prop_hp", "r3", 2);
        mult("r2", "r3", "r4", 4);
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14a", t_start_Q14a, print_level);

        // Q14b
        double t_start_Q14b = abs_time();
        reorder("t3c", "r1", "124563");
        reorder("h3c+_g1", "r2", "145623");
        mult("r1", "prop_ph", "r3", 2);
        mult("r2", "r3", "r4", 4);
        update(result_name, 0.25, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14b", t_start_Q14b, print_level);

        // Q14c
        double t_start_Q14c = abs_time();
        reorder("t3c", "r1", "456123");
        reorder("prop_pg", "r2", "21");
        mult("h3c+", "r1", "r3", 5);
        mult("r3", "r2", "r4", 1);
        update(result_name, -1.0/12.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14c", t_start_Q14c, print_level);

        // Q14d
        double t_start_Q14d = abs_time();
        reorder("h3c", "r1", "456123");
        mult("r1", "t3c+", "r2", 5);
        mult("prop_gp", "r2", "r3", 1);
        update(result_name, -1.0/12.0, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q14d", t_start_Q14d, print_level);

        // Q15a
        double t_start_Q15a = abs_time();
        reorder("h3c", "r1", "451263");
        mult("r1", "prop_ph", "r2", 2);
        mult("h2c+", "r2", "r3", 3);
        update(result_name, -0.5, "r3");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15a", t_start_Q15a, print_level);

        // Q15b
        double t_start_Q15b = abs_time();
        reorder("h2c", "r1", "3412");
        reorder("h3c+", "r2", "124563");
        mult("r2", "prop_hp", "r3", 2);
        mult("r3", "r1", "r4", 3);
        update(result_name, -0.5, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q15b", t_start_Q15b, print_level);

        // Q16a
        double t_start_Q16a = abs_time();
        reorder("prop_pp", "r1", "21");
        mult("h3c", "r1", "r2", 1);
        reorder("r2", "r3", "456123");
        mult("h3c+", "r3", "r4", 5);
        update(result_name, -1.0/6.0, "r4");
        restore_stack_pos(pos);
        print_time_mem_usage("Q16a", t_start_Q16a, print_level);

        // Q16b
        double t_start_Q16b = abs_time();
        reorder("h3c", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        mult("h3c+", "r2", "r3", 5);
        update(result_name, 0.25, "r3");
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
        sector_1h0p_overlap("TT_conn", CC_PROPERTIES_APPROX_QUADRATIC);
        reorder("prop_gg", "r1", "21");
        mult("TT_conn", "r1", "r2", 1);
        update(result_name, 1.0, "r2");
        restore_stack_pos(pos);
    }

    if (approximation <= CC_PROPERTIES_APPROX_QUADRATIC) {
        return;
    }
}


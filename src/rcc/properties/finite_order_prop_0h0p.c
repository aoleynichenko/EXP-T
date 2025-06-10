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

#include "finite_order_prop_0h0p.h"

#include "cc_properies.h"
#include "engine.h"
#include "spinors.h"
#include "symmetry.h"


double complex finite_order_property_0h0p(int approximation)
{
    dg_stack_pos_t pos = get_stack_pos();

    double complex d = 0.0 + 0.0 * I;

    if (approximation <= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return d;
    }

    // L1a + L2a
    d += 2.0 * scalar_product("C", "N", "prop_hp", "t1c");

    if (approximation <= CC_PROPERTIES_APPROX_LINEAR) {
        return d;
    }

    // Q1a
    reorder("t1c", "r1", "21");
    mult("prop_hh", "r1", "r2", 1);
    d -= scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    // Q1b
    reorder("prop_pp", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    d += scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    // Q2a + Q2b
    reorder("t2c", "r1", "1342");
    mult("r1", "prop_ph", "r2", 2);
    d += 2.0 * scalar_product("C", "N", "t1c", "r2");
    restore_stack_pos(pos);

    // Q3a
    reorder("t2c", "r1", "3412");
    mult("r1", "prop_hh", "r2", 1);
    reorder("r2", "r3", "3412");
    d += -0.5 * scalar_product("C", "N", "t2c", "r3");
    restore_stack_pos(pos);

    // Q3b
    reorder("prop_pp", "r1", "21");
    mult("t2c", "r1", "r2", 1);
    d += 0.5 * scalar_product("C", "N", "t2c", "r2");
    restore_stack_pos(pos);

    if (triples_enabled()) {

        // Q4a + Q4b
        reorder("t3c", "r1", "145623");
        mult("r1", "t2c+", "r2", 4);
        d += (0.25 + 0.25) * scalar_product("C", "N", "prop_hp", "r2");
        restore_stack_pos(pos);

        // Q5a
        reorder("t3c", "r1", "456123");
        mult("r1", "prop_hh", "r2", 1);
        reorder("r2", "r3", "456123");
        d -= (1.0 / 12.0) * scalar_product("C", "N", "t3c", "r3");
        restore_stack_pos(pos);

        // Q5b
        reorder("prop_pp", "r1", "21");
        mult("t3c", "r1", "r2", 1);
        d += (1.0 / 12.0) * scalar_product("C", "N", "t3c", "r2");
        restore_stack_pos(pos);
    }

    return d;
}

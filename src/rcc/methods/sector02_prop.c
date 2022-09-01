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


#include "cc_properies.h"
#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ccutils.h"
#include "engine.h"
#include "datamodel.h"
#include "diis.h"
#include "heff.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"


void property_0h2p_linear_D1(char *result_name);
void property_0h2p_linear_D2(char *result_name);


void sector_0h2p_prop(char *result_name, int approximation)
{
    tmplt(result_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE);

    if (approximation >= CC_PROPERTIES_APPROX_LINEAR) {
        property_0h2p_linear_D1(result_name);
        property_0h2p_linear_D2(result_name);
    }
}


void property_0h2p_linear_D1(char *result_name)
{
    dg_stack_pos_t pos = get_stack_pos();

    // D1a
    reorder("s2c", "r1", "1342");
    mult("r1", "op_ph", "r2", 1);
    reorder("r2", "r3", "1423");
    closed("r3", "r4");
    perm("r4", "(12)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // D1b
    reorder("op_hp", "r1", "21");
    mult("s2c+", "r1", "r2", 1);
    closed("r2", "r3");
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);
}


void property_0h2p_linear_D2(char *result_name)
{
    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("op_pp", "r1", "21");
    mult("x2c", "r1", "r2", 1);
    closed("r2", "r3");
    perm("r3", "(34)");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // D2b
    reorder("x2c+", "r1", "3412");
    mult("r1", "op_pp", "r2", 1);
    reorder("r2", "r3", "3412");
    closed("r3", "r4");
    perm("r4", "(12)");
    update(result_name, +1.0, "r4");
    restore_stack_pos(pos);
}


void sector_0h2p_norm(char *result_name, int approximation)
{
    tmplt(result_name, "pppp", "1111", "1234", NOT_PERM_UNIQUE);

    if (approximation == CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // D1a
    reorder("s2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 1);
    reorder("r2", "r3", "1423");
    closed("r3", "r4");
    perm("r4", "(12)");
    update(result_name, -1.0, "r4");
    restore_stack_pos(pos);

    // D1b
    reorder("t1c", "r1", "21");
    mult("s2c+", "r1", "r2", 1);
    closed("r2", "r3");
    perm("r3", "(34)");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // D2a
    reorder("s1c+", "r1", "21");
    mult("x2c", "r1", "r2", 1);
    closed("r2", "r3");
    perm("r3", "(34)");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // D2b
    reorder("x2c+", "r1", "3412");
    mult("r1", "s1c", "r2", 1);
    reorder("r2", "r3", "3412");
    closed("r3", "r4");
    perm("r4", "(12)");
    update(result_name, +1.0, "r4");
    restore_stack_pos(pos);

    // D3
    reorder("t2c", "r1", "3412");
    mult("t2c+", "r1", "r2", 2);
    closed("r2", "r3");
    update(result_name, +0.5, "r3");
    restore_stack_pos(pos);

    // D4
    reorder("s2c+", "s2c+_21", "2143"); // interchange electrons 1 <-> 2
    reorder("s2c+_21", "r1", "2431");
    reorder("s2c", "r2", "1324");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1324");
    closed("r4", "r5");
    perm("r5", "(12|34)");
    update(result_name, +1.0, "r5");
    restore_stack_pos(pos);

    // D5
    reorder("x2c+", "r1", "3412");
    mult("x2c", "r1", "r2", 2);
    closed("r2", "r3");
    update(result_name, +0.5, "r3");
    restore_stack_pos(pos);
}


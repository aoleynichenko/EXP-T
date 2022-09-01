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


int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

void guess_operator_symmetry(int nspinors, double complex *prop_matrix);

double complex calculate_vacuum_expectation_value(int nspinors, double complex *prop_matrix);

void diagram_conjugate(char *source_name, char *target_name);


void property_0h1p_zero_D1(char *result);
void property_0h1p_linear_D2(char *result);
void property_0h1p_linear_D3(char *result);
void property_0h1p_linear_D4(char *result);


void sector_0h1p_prop(char *result_name, int approximation)
{
    tmplt(result_name, "pp", "11", "12", NOT_PERM_UNIQUE);

    if (approximation >= CC_PROPERTIES_APPROX_MODEL_SPACE) {
        property_0h1p_zero_D1(result_name);
    }
    if (approximation >= CC_PROPERTIES_APPROX_LINEAR) {
        property_0h1p_linear_D2(result_name);
        property_0h1p_linear_D3(result_name);
        property_0h1p_linear_D4(result_name);
    }
}


void property_0h1p_zero_D1(char *result)
{
    dg_stack_pos_t pos = get_stack_pos();

    copy("op_pp", "r1");
    closed("r1", "r2");
    update(result, 1.0, "r2");

    restore_stack_pos(pos);
}


void property_0h1p_linear_D2(char *result)
{
    dg_stack_pos_t pos = get_stack_pos();

    // D2a
    reorder("t1c", "r1", "21");
    mult("op_ph", "r1", "r2", 1);
    closed("r2", "r3");
    update(result, -1.0, "r3");
    restore_stack_pos(pos);

    // D2b
    reorder("op_hp", "r1", "21");
    mult("t1c+", "r1", "r2", 1);
    closed("r2", "r3");
    update(result, -1.0, "r3");
    restore_stack_pos(pos);
}


void property_0h1p_linear_D3(char *result)
{
    dg_stack_pos_t pos = get_stack_pos();

    // D3a
    reorder("op_pp", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    closed("r2", "r3");
    update(result, +1.0, "r3");
    restore_stack_pos(pos);

    // D3b
    reorder("s1c+", "r1", "21");
    mult("op_pp", "r1", "r2", 1);
    closed("r2", "r3");
    update(result, +1.0, "r3");
    restore_stack_pos(pos);
}


void property_0h1p_linear_D4(char *result)
{
    dg_stack_pos_t pos = get_stack_pos();

    // D4a
    reorder("s2c", "r1", "1342");
    mult("r1", "op_ph", "r2", 2);
    closed("r2", "r3");
    update(result, +1.0, "r3");
    restore_stack_pos(pos);

    // D4a
    reorder("s2c+", "r1", "1342");
    mult("r1", "op_hp", "r2", 2);
    closed("r2", "r3");
    update(result, +1.0, "r3");
    restore_stack_pos(pos);
}


void sector_0h1p_norm(char *result_name, int approximation)
{
    tmplt(result_name, "pp", "11", "12", NOT_PERM_UNIQUE);

    if (approximation == CC_PROPERTIES_APPROX_MODEL_SPACE) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // D1
    reorder("t1c", "r1", "21");
    mult("t1c+", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, -1.0, "r3");
    restore_stack_pos(pos);

    // D2a
    reorder("s2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    closed("r2", "r3");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // D2b
    reorder("s2c+", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    closed("r2", "r3");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // D3
    reorder("s2c+", "r1", "3412");
    mult("s2c", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, +0.5, "r3");
    restore_stack_pos(pos);

    // D4
    reorder("s1c+", "r1", "21");
    mult("s1c", "r1", "r2", 1);
    closed("r2", "r3");
    update(result_name, +1.0, "r3");
    restore_stack_pos(pos);

    // D5
    reorder("t2c", "r1", "3412");
    mult("t2c+", "r1", "r2", 3);
    closed("r2", "r3");
    update(result_name, -0.5, "r3");
    restore_stack_pos(pos);
}

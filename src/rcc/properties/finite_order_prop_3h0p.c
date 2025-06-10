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

#include "finite_order_prop_3h0p.h"

#include <stdlib.h>

#include "cc_properies.h"
#include "engine.h"
#include "spinors.h"

#include "finite_order_prop.h"
#include "finite_order_overlap.h"


void finite_order_property_3h0p(char *result_name, int operator_symmetry, int approximation, int disconnected)
{
    tmplt_sym(result_name, "hhhhhh", "111111", "123456", NOT_PERM_UNIQUE, operator_symmetry);
    return;
}



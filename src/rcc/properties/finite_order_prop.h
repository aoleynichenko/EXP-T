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


#ifndef CC_FINITE_ORDER_PROP_INCLUDED
#define CC_FINITE_ORDER_PROP_INCLUDED

void calculate_properties_finite_order_method(int sect_h, int sect_p);

void print_time_mem_usage(char *diagram_name, double t_start, int print_level);

#include "finite_order_overlap.h"
#include "finite_order_prop_0h0p.h"
#include "finite_order_prop_1h0p.h"
#include "finite_order_prop_0h1p.h"
#include "finite_order_prop_2h0p.h"
#include "finite_order_prop_0h2p.h"
#include "finite_order_prop_1h1p.h"

#endif // CC_FINITE_ORDER_PROP_INCLUDED


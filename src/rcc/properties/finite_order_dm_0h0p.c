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

/*
 * finite-order construction of the density matrix for the vacuum state (the 0h0p sector).
 * only terms quadratic in cluster operators are retained.
 */

#include "finite_order_dm_0h0p.h"

#include "engine.h"
#include "symmetry.h"


void finite_order_density_matrix_block_hp_0h0p(int dm_sym)
{
    if (dm_sym != get_totally_symmetric_irrep()) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // DM_00_L1a
    update("dm_hp", 1.0, "t1c");
    restore_stack_pos(pos);

    // DM_00_Q2a
    reorder("t2c", "r1", "1342");
    mult("r1", "t1c+", "r2", 2);
    update("dm_hp", 1.0, "r2");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_ph_0h0p(int dm_sym)
{
    if (dm_sym != get_totally_symmetric_irrep()) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // DM_00_L1b
    update("dm_ph", 1.0, "t1c+");
    restore_stack_pos(pos);

    // DM_00_Q2b
    reorder("t2c+", "r1", "1342");
    mult("r1", "t1c", "r2", 2);
    update("dm_ph", 1.0, "r2");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_hh_0h0p(int dm_sym)
{
    if (dm_sym != get_totally_symmetric_irrep()) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // DM_00_Q1a
    reorder("t1c+", "r1", "21");
    mult("t1c", "r1", "r2", 1);
    update("dm_hh", -1.0, "r2");
    restore_stack_pos(pos);

    // DM_00_Q3a
    reorder("t2c+", "r1", "3412");
    mult("t2c", "r1", "r2", 3);
    update("dm_hh", -0.5, "r2");
    restore_stack_pos(pos);
}


void finite_order_density_matrix_block_pp_0h0p(int dm_sym)
{
    if (dm_sym != get_totally_symmetric_irrep()) {
        return;
    }

    dg_stack_pos_t pos = get_stack_pos();

    // DM_00_Q1b
    reorder("t1c", "r1", "21");
    mult("t1c+", "r1", "r2", 1);
    update("dm_pp", 1.0, "r2");
    restore_stack_pos(pos);

    // DM_00_Q3b
    reorder("t2c", "r1", "3412");
    mult("t2c+", "r1", "r2", 3);
    update("dm_pp", 0.5, "r2");
    restore_stack_pos(pos);
}

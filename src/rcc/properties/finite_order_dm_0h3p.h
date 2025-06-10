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

#ifndef CC_FINITE_ORDER_DM_0H3P_H_INCLUDED
#define CC_FINITE_ORDER_DM_0H3P_H_INCLUDED

void construct_n_part_active_density_matrix_diagram_0h3p(
    char *diagram_name,
    char *bra_irrep_name, int bra_state_number,
    char *ket_irrep_name, int ket_state_number,
    int n_particles
);

void finite_order_density_matrix_block_hp_0h3p(int dm_sym);

void finite_order_density_matrix_block_ph_0h3p(int dm_sym);

void finite_order_density_matrix_block_hh_0h3p(int dm_sym);

void finite_order_density_matrix_block_pp_0h3p(int dm_sym);

#endif // CC_FINITE_ORDER_DM_0H3P_H_INCLUDED

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

#ifndef FINITE_ORDER_OVERLAP_H_INCLUDED
#define FINITE_ORDER_OVERLAP_H_INCLUDED

#include <complex.h>
#include "../heff/slater_det.h"

void calculate_wavefunction_norms_and_overlaps(int sect_h, int sect_p);

double calculate_wavefunction_norm(int sect_h, int sect_p, char *irrep_name, int state_number);

double sector_0h0p_overlap(int approximation);

void sector_1h0p_overlap(char *result_name, int approximation);

void sector_0h1p_overlap(char *result_name, int approximation);

void sector_1h1p_overlap(char *result_name, int approximation);

void sector_2h0p_overlap(char *result_name, int approximation);

void sector_0h2p_overlap(char *result_name, int approximation);

void sector_3h0p_overlap(char *result_name, int approximation);

void sector_0h3p_overlap(char *result_name, int approximation);

#endif // FINITE_ORDER_OVERLAP_H_INCLUDED

/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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
 * Tools for construction and diagonalization of the effective Hamiltonian
 * matrix and the subsequent analysis of its eigenvalues and eigenvectors.
 */

#ifndef CC_HEFF_H_INCLUDED
#define CC_HEFF_H_INCLUDED

#include <complex.h>
#include <stdarg.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

#include "engine.h"
#include "options.h"
#include "../heff/slater_rules.h"
#include "../heff/ms_prop.h"
#include "../heff/direct_prop.h"
#include "../heff/ms_nat_spinors.h"
#include "../heff/cc_wfn_overlap.h"


void heff_analysis(int sect_h, int sect_p, ...);

void construct_model_space_natural_spinors(int sect_h, int sect_p, int rep1, int state1, int rep2, int state2);

void dipole_length_tdms(int sect_h, int sect_p);

void model_space_property(cc_ms_prop_query_t *prop_query, int approximation);

void model_space_properties_and_natural_orbitals();

void write_formatted_heff_0h0p(double total_energy);

#endif /* CC_HEFF_H_INCLUDED */

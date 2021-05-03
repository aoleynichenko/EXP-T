/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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

/*******************************************************************************
 * renorm_omega.h
 * ==============
 *
 * Restoration of the intermediate normalization of the wave operator
 * Omega = { exp(T) }.
 * This code is used primarily for the FF calculation of off-diagonal matrix
 * elements between the 0h0p ground and 1h1p excited states.

 * For details, see:
 * A. Zaitsevskii, A. V. Oleynichenko, E. Eliav,
 * Finite-Field Calculations of Transition Properties by the Fock Space
 * Relativistic Coupled Cluster Method: Transitions between Different
 * Fock Space Sectors.
 * Symmetry 2020, 12(11), 1845; https://doi.org/10.3390/sym12111845
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_RENORM_OMEGA_H_INCLUDED
#define CC_RENORM_OMEGA_H_INCLUDED

#include "comdef.h"
#include "slater_det.h"
#include "symmetry.h"

void restore_intermediate_normalization(size_t *block_dims, slater_det_t **det_list, double complex **heff,
                                        double complex **eigvalues);


#endif /* CC_RENORM_OMEGA_H_INCLUDED */

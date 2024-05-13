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

#ifndef EXPT_TRANSFORM_H
#define EXPT_TRANSFORM_H

#include <complex.h>

void transform(int n_nat_spinors, int n_act_spinors, int basis_dim,
               int *map_act_to_all, int *kramers_conj,
               double complex *nat_coeffs_spinor_basis,
               double complex *cmo_alpha, double complex *cmo_beta,
               double complex *cno_alpha, double complex *cno_beta
);

#endif //EXPT_TRANSFORM_H

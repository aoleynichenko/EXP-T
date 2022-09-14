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

/*
 * Slater rules: evaluation of matrix elements in the basis of Slater
 * determinants.
 *
 * 2019-2021 Alexander Oleynichenko
 */

#ifndef CC_SLATER_RULES_H_INCLUDED
#define CC_SLATER_RULES_H_INCLUDED

#include "comdef.h"
#include "engine.h"
#include "slater_det.h"

typedef double complex (*matrix_getter_fun)(void *source, void *indices);

extern double complex (*slater_rule)(slater_det_t *d1, slater_det_t *d2);

void setup_slater(void *source, matrix_getter_fun getter,
                  int bra_sect_h, int bra_sect_p, int ket_sect_h, int ket_sect_p, int npart);

#endif /* CC_SLATER_RULES_H_INCLUDED */

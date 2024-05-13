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
 * Tools for operating with Slater determinants.
 */

#ifndef CC_SLATER_DET_H_INCLUDED
#define CC_SLATER_DET_H_INCLUDED

#include <stdio.h>

#include "comdef.h"

typedef struct {
    moindex_t indices[MAX_SECTOR_RANK];
    int16_t sym;
} slater_det_t;

int is_vacuum_det(slater_det_t *det);

int set_vacuum_det(slater_det_t *det);

double model_det_energy(int sect_h, int sect_p, slater_det_t *det);

void str_print_slater_det(char *str, int sect_h, int sect_p, slater_det_t *det);

void print_slater_det(FILE *f, int sect_h, int sect_p, slater_det_t *det);

double kronecker_delta(int sect_h, int sect_p, slater_det_t *det_1, slater_det_t *det_2);

#endif /* CC_SLATER_DET_H_INCLUDED */

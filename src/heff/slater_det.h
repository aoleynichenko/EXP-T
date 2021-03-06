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
 * slater_det.h
 * ============
 *
 * Tools for operating with Slater determinants.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_SLATER_DET_H_INCLUDED
#define CC_SLATER_DET_H_INCLUDED

#include <stdio.h>

#include "comdef.h"
#include "engine.h"

typedef struct {
    moindex_t indices[MAX_SECTOR_RANK];
    int8_t sym;
} slater_det_t;

int is_vacuum_det(slater_det_t *det);

int set_vacuum_det(slater_det_t *det);

int detcmp(const void *_d1, const void *_d2);

void print_slater_det(FILE *f, int sect_h, int sect_p, slater_det_t *det);

#endif /* CC_SLATER_DET_H_INCLUDED */

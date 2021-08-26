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
 * intham1.h
 * =========
 *
 * Simple version of the intermediate Hamiltonian-like technique
 *
 * 2021 Alexander Oleynichenko
 ******************************************************************************/

#ifndef CC_INTHAM1_H_INCLUDED
#define CC_INTHAM1_H_INCLUDED

#include "comdef.h"

enum {
    IH1_SHIFT_TO_GIVEN_BOUNDS = 0,
    IH1_SHIFT_TO_MAIN_BOUNDS,
    IH1_SHIFT_TO_BARYCENTER
};

enum {
    IH1_SUM_FORMULA_BOX = 0,
    IH1_SUM_FORMULA_LINE
};

typedef struct {
    int sectors[MAX_SECTOR_RANK][MAX_SECTOR_RANK];
    int main_holes;
    int main_particles;
    int shift_to_formula;
    double holes_shift_to;
    double particles_shift_to;
    int npower;
    int shift_type;
    int sum_formula;
} ih1_options_t;

void intham1_calculate_shifts(int sect_h, int sect_p);
double intham1_get_shift_value(int sect_h, int sect_p, int rank, int *indices, int *valence_labels);

#endif /* CC_INTHAM1_H_INCLUDED */

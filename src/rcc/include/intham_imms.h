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
 * Simple Intermediate Hamiltonian for incomplete main model spaces.
 *
 * 2021-2022 Alexander Oleynichenko
 */

#ifndef CC_INTHAM_IMMS_H_INCLUDED
#define CC_INTHAM_IMMS_H_INCLUDED

#include "comdef.h"

#define IH_IMMS_MAX_SPINOR_SUBSPACES 50
#define IH_IMMS_MAX_MAIN_SUBSPACES   50

typedef struct {
    int sectors[MAX_SECTOR_RANK][MAX_SECTOR_RANK];
    double det_shift_to[MAX_SECTOR_RANK][MAX_SECTOR_RANK];
    int det_shift_auto[MAX_SECTOR_RANK][MAX_SECTOR_RANK];
    int npower;
    int shift_type;
    double scale_shift;
    int n_spinor_subspaces;
    double subspace_energy_ranges[IH_IMMS_MAX_SPINOR_SUBSPACES][2];
    int n_main_subspaces;
    int main_occ[IH_IMMS_MAX_MAIN_SUBSPACES][IH_IMMS_MAX_SPINOR_SUBSPACES];
} ih_imms_options_t;

int intham_imms_in_sector(int sect_h, int sect_p);

void intham_imms_setup(int sect_h, int sect_p);

double intham_imms_get_shift_value(int sect_h, int sect_p, int rank, int *indices, int *valence_labels);

#endif /* CC_INTHAM_IMMS_H_INCLUDED */
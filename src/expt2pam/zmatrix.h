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
 * Tools for conversion of coordinates: zmatrix -> cartesian
 *
 * 2023 Alexander Oleynichenko
 */

#ifndef EXPT_ZMATRIX_H
#define EXPT_ZMATRIX_H

#include "molecule.h"

void convert_zmatrix_to_xyz(int n_atoms, int *charges, int *rconnect, double *rlist,
                            int *aconnect, double *alist, int *dconnect, double *dlist, molecule_t *mol);

#endif //EXPT_ZMATRIX_H

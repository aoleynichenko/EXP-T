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
 * Sorting of integrals (creating basic diagrams from the raw integral arrays).
 *
 * 2018-2021 Alexander Oleynichenko
 */

#ifndef CC_SORT_H_INCLUDED
#define CC_SORT_H_INCLUDED

#include <complex.h>
#include <stdint.h>

#include "datamodel.h"

#define CC_SORTING_IO_BUF_SIZE 16384

// leave a request for sorting of diagram named 'name' with given 'qparts',
// 'valence', 'order' characteristic
void request_sorting(char *name, char *qparts, char *valence, char *order);

// performs sorting for all the requests leaved
void perform_sorting();

void sort_prop(int nspinors, double complex *oper_ints);

#endif /* CC_SORT_H_INCLUDED */

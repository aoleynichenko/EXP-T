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

/**
* Electronic structure models.
*
* 2018-2021 Alexander Oleynichenko
*/

#ifndef CC_METHODS_H_INCLUDED
#define CC_METHODS_H_INCLUDED

#include "options.h"

// orders of perturbation theory
enum {
    PT_2 = 2,
    PT_3 = 3,
    PT_4 = 4,
    PT_5 = 5,
    PT_INF = 100
};

int sector00(cc_options_t *opts);

int sector01(cc_options_t *opts);

int sector10(cc_options_t *opts);

int sector11(cc_options_t *opts);

int sector02(cc_options_t *opts);

int sector20(cc_options_t *opts);

int sector03(cc_options_t *opts);

void mixed_00_11(cc_options_t *opts);

#endif /* CC_METHODS_H_INCLUDED */

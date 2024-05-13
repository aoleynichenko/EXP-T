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
 * Electronic structure models.
 */

#ifndef CC_METHODS_H_INCLUDED
#define CC_METHODS_H_INCLUDED

#include "options.h"

// orders of perturbation theory
enum {
    PT_1 = 1,
    PT_2 = 2,
    PT_3 = 3,
    PT_4 = 4,
    PT_5 = 5,
    PT_INF = 100
};

int sector_0h0p(cc_options_t *opts);

int sector_0h1p(cc_options_t *opts);

int sector_1h0p(cc_options_t *opts);

int sector_1h1p(cc_options_t *opts);

int sector_0h2p(cc_options_t *opts);

int sector_2h0p(cc_options_t *opts);

int sector_0h3p(cc_options_t *opts);

int sector_3h0p(cc_options_t *opts);

int sector_1h2p(cc_options_t *opts);

int sector_2h1p(cc_options_t *opts);

#endif /* CC_METHODS_H_INCLUDED */

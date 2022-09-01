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

#ifndef CC_PROPERIES_H_INCLUDED
#define CC_PROPERIES_H_INCLUDED

#include <complex.h>

#include "../heff/slater_det.h"

enum {
    CC_PROPERTIES_APPROX_MODEL_SPACE = 0,
    CC_PROPERTIES_APPROX_LINEAR,
    CC_PROPERTIES_APPROX_QUADRATIC,
    CC_PROPERTIES_APPROX_CUBIC,
    CC_PROPERTIES_APPROX_QUARTIC
};

void sector_0h0p_calculate_properties();

double complex sector_0h0p_norm(int approximation);

void sector_0h1p_norm(char *result_name, int approximation);

void sector_0h2p_norm(char *result_name, int approximation);

double complex sector_0h0p_prop(int approximation);

void sector_0h1p_prop(char *result_name, int approximation);

void sector_0h2p_prop(char *result_name, int approximation);

#endif /* CC_PROPERIES_H_INCLUDED */

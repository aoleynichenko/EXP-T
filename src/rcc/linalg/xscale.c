/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
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

#include "linalg.h"

#include <assert.h>
#include <math.h>


/*
 * scale vector element-wise by a number
 */
void xscale(data_type_t type, size_t n, void *x, void *multiplier)
{
    assert(x != NULL);
    assert(multiplier != NULL);
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);

    if (type == CC_DOUBLE) {
        double *dbuf = x;
        double a = *((double *) multiplier);

        for (size_t i = 0; i < n; i++) {
            dbuf[i] *= a;
        }
    }
    else { // CC_DOUBLE_COMPLEX
        double complex *zbuf = x;
        double complex a = *((double complex *) multiplier);

        for (size_t i = 0; i < n; i++) {
            zbuf[i] *= a;
        }
    }
}

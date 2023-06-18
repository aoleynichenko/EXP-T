/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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

#include "datamodel.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "error.h"
#include "memory.h"
#include "spinors.h"
#include "templates.h"


/*
 * if data to be compressed are of type 'double'
 */
#define TYPENAME_UNCOMPRESSED double

#define TYPENAME_COMPRESSED double
#include "compress_triples_template.c"
#undef TYPENAME_COMPRESSED
#define TYPENAME_COMPRESSED float
#include "compress_triples_template.c"
#undef TYPENAME_COMPRESSED

#undef TYPENAME_UNCOMPRESSED

/*
 * if data to be compressed are of type 'double complex'
 */
#define TYPENAME_UNCOMPRESSED double_complex_t

#define TYPENAME_COMPRESSED double_complex_t
#include "compress_triples_template.c"
#undef TYPENAME_COMPRESSED
#define TYPENAME_COMPRESSED float_complex_t
#include "compress_triples_template.c"
#undef TYPENAME_COMPRESSED

#undef TYPENAME_UNCOMPRESSED


void compress_triples_rank6(block_t *block)
{
    double zero_thresh = cc_opts->compress_triples_thresh;

    if (arith == CC_ARITH_COMPLEX) {
        if (cc_opts->compress_triples_data_type == CC_DOUBLE) {
            TEMPLATE2(compress_amplitudes, double_complex_t, double_complex_t)(block, zero_thresh);
        }
        else {
            TEMPLATE2(compress_amplitudes, double_complex_t, float_complex_t)(block, zero_thresh);
        }
    }
    else {
        if (cc_opts->compress_triples_data_type == CC_DOUBLE) {
            TEMPLATE2(compress_amplitudes, double, double)(block, zero_thresh);
        }
        else {
            TEMPLATE2(compress_amplitudes, double, float)(block, zero_thresh);
        }
    }
}


void decompress_triples_rank6(block_t *block)
{
    if (arith == CC_ARITH_COMPLEX) {
        if (cc_opts->compress_triples_data_type == CC_DOUBLE) {
            TEMPLATE2(decompress_amplitudes, double_complex_t, double_complex_t)(block);
        }
        else {
            TEMPLATE2(decompress_amplitudes, double_complex_t, float_complex_t)(block);
        }
    }
    else {
        if (cc_opts->compress_triples_data_type == CC_DOUBLE) {
            TEMPLATE2(decompress_amplitudes, double, double)(block);
        }
        else {
            TEMPLATE2(decompress_amplitudes, double, float)(block);
        }
    }
}



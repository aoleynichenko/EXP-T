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

#ifdef TYPENAME_UNCOMPRESSED
#ifdef TYPENAME_COMPRESSED

#include "templates.h"


/**
 * Functions used to calculate absolute value are different for real and complex numbers,
 * 'fabs' and 'cabs', respectively. 'cabs' is much slower than 'fabs'
 */
#if TYPENAME_UNCOMPRESSED == double
    #undef ABS_FUN
    #define ABS_FUN(x) fabs(x)
#else
    #undef ABS_FUN
    #define ABS_FUN(x) cabs(x)
#endif

/**
 * Counts non-zero amplitudes.
 */
size_t TEMPLATE2(count_nonzero_amplitudes, TYPENAME_UNCOMPRESSED, TYPENAME_COMPRESSED)(block_t *block, double thresh)
{
    size_t n_nonzero = 0;
    TYPENAME_UNCOMPRESSED *buf = (TYPENAME_UNCOMPRESSED *) block->buf;

    for (size_t i = 0; i < block->size; i++) {
        if (ABS_FUN(buf[i]) >= thresh) {
            n_nonzero++;
        }
    }

    return n_nonzero;
}


/**
 * compression.
 * TYPENAME_UNCOMPRESSED = data type of uncompressed amplitudes ('double' or 'double complex')
 * TYPENAME_COMPRESSED = data type of compressed amplitudes ('double', 'float', 'double complex' or 'float complex')
 */
void TEMPLATE2(compress_amplitudes, TYPENAME_UNCOMPRESSED, TYPENAME_COMPRESSED)(block_t *block, double zero_thresh)
{
    if (block->is_compressed == 1) {
        return;
    }

    size_t n_nonzero = TEMPLATE2(count_nonzero_amplitudes, TYPENAME_UNCOMPRESSED, TYPENAME_COMPRESSED)(block, zero_thresh);
    if (n_nonzero == 0) {
        cc_free(block->buf);
        // 'block->buf' contains only number of compressed elements (zero)
        block->buf = (double complex *) cc_malloc(sizeof(size_t) * 1);
        *((size_t *) block->buf) = 0;
        block->is_compressed = 1;
        return;
    }

    size_t dst_size = sizeof(size_t) + (sizeof(TYPENAME_COMPRESSED) + sizeof(size_t)) * n_nonzero;
    char *dst = (char *) cc_malloc(dst_size);
    *((size_t *) dst) = n_nonzero;

    char *dst_ptr = dst + sizeof(size_t); // skip 8 bytes for the length of the compressed data

    for (size_t i = 0; i < block->size; i++) {
        TYPENAME_UNCOMPRESSED amplt_value = ((TYPENAME_UNCOMPRESSED *) block->buf)[i];
        if (ABS_FUN(amplt_value) < zero_thresh) {
            continue;
        }

        // write pair: linear index (size_t) + amplitude (double)
        *((size_t *) dst_ptr) = i;
        *((TYPENAME_COMPRESSED *) (dst_ptr + sizeof(size_t))) = (TYPENAME_COMPRESSED) amplt_value;

        dst_ptr += sizeof(size_t) + sizeof(TYPENAME_COMPRESSED);
    }

    size_t n_compressed = (dst_ptr - dst - sizeof(size_t)) / (sizeof(size_t) + sizeof(TYPENAME_COMPRESSED));
    if (n_compressed != n_nonzero) {
        printf("n_compressed = %ld\n", n_compressed);
        printf("n_nonzero = %ld\n", n_nonzero);
        errquit("compress_triples: n_compressed != n_nonzero %d", *((char *) NULL));
    }

    cc_free(block->buf);
    block->buf = (double complex *) dst;
    block->is_compressed = 1;
}


/**
 * decompression.
 */
void TEMPLATE2(decompress_amplitudes, TYPENAME_UNCOMPRESSED, TYPENAME_COMPRESSED)(block_t *block)
{
    if (block->is_compressed == 0) {
        return;
    }

    char *compressed_data = (char *) block->buf;
    size_t src_size = *((size_t *) compressed_data);
    compressed_data += sizeof(size_t); // skip bytes containing 'src_size'

    double complex *dst = (double complex *) cc_calloc(block->size, SIZEOF_WORKING_TYPE);

    for (size_t i = 0; i < src_size; i++) {

        char *ptr_index = compressed_data;
        char *ptr_value = compressed_data + sizeof(size_t);

        size_t index = *((size_t *) ptr_index);
        TYPENAME_UNCOMPRESSED amplt_value = *((TYPENAME_COMPRESSED *) ptr_value);

        ((TYPENAME_UNCOMPRESSED *) dst)[index] = amplt_value;

        compressed_data += sizeof(size_t) + sizeof(TYPENAME_COMPRESSED);
    }

    cc_free(block->buf);
    block->buf = dst;
    block->is_compressed = 0;
}


#endif /* TYPENAME_COMPRESSED defined */
#endif /* TYPENAME_UNCOMPRESSED defined */
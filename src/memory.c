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
 * memory.c
 * ========
 *
 * Memory allocator.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t n_allocated = 0;
static size_t max_allocated = 0;
static size_t max_available = 100 * 1024u * 1024u; // start with 100 Mb

void cc_memory_usage();


void cc_init_allocator(size_t max_mem)
{
    max_available = max_mem;
}


/**
 * Wrapper for malloc() from libc.
 * Allocates memory:
 * [sizeof(size_t) block-size-nbytes] [useful space]
 *                                    ^
 *                                 returns
 */
void *cc_malloc(size_t nbytes)
{
    void *mem;

    if (n_allocated + nbytes > max_available) {
        printf("bytes allocated        = %ld\n", n_allocated);
        printf("max memory usage limit = %ld\n", max_available);
        printf("bytes available (free) = %ld\n", max_available - n_allocated);
        printf("bytes required for the current allocation = %ld\n", nbytes);
        printf("cc_malloc(): cannot allocate memory (not enough memory)\n");
        abort();
        return NULL;
    }

    mem = malloc(sizeof(size_t) + nbytes);
    if (mem == NULL) {
        printf("cc_malloc(): cannot allocate memory (%ld bytes): %s\n", nbytes,
               strerror(errno));
        return NULL;
    }

    // first entry = number of "useful" bytes
    *((size_t *) mem) = nbytes;

    // update counters
    n_allocated += nbytes;
    if (n_allocated > max_allocated) {
        max_allocated = n_allocated;
    }

    return mem + sizeof(size_t);
}


/**
 * Allocates a block of memory for an array of num elements,
 * each of them size bytes long, and initializes all its bits to zero.
 * @param num
 * @param size
 * @return pointer to memory if success, NULL if error
 */
void *cc_calloc(size_t num, size_t size)
{
    void *mem = cc_malloc(num * size);
    if (mem != NULL) {
        memset(mem, 0, num * size);
    }
    return mem;
}


/**
 * Wrapper for free() from libc
 */
void cc_free(void *p)
{
    size_t nbytes;
    void *mem;

    if (p == NULL) { return; }

    mem = p - sizeof(size_t);
    nbytes = *((size_t *) mem);

    // update counter
    n_allocated -= nbytes;

    // free memory
    free(mem);
}


void cc_finalize_allocator()
{
    cc_memory_usage();
}


void cc_memory_usage()
{
    const double b2mb = 1.0 / (1024.0 * 1024.0);
    const double b2gb = 1.0 / (1024.0 * 1024.0 * 1024.0);

    printf("\n");
    printf(" memory in usage  = %ld bytes = %.1f Mb = %.2f Gb\n",
           n_allocated, b2mb * n_allocated, b2gb * n_allocated);
    printf(" max memory usage = %ld bytes = %.1f Mb = %.2f Gb\n",
           max_allocated, b2mb * max_allocated, b2gb * max_allocated);
}


size_t cc_get_current_memory_usage()
{
    return n_allocated;
}


size_t cc_get_peak_memory_usage()
{
    return max_allocated;
}

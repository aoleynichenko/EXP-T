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

/*******************************************************************************
* memory.h
* ========
*
* Memory allocator.
*
* 2018-2021 Alexander Oleynichenko
******************************************************************************/

#ifndef CC_MEMORY_H_INCLUDED
#define CC_MEMORY_H_INCLUDED

// initialization of memory allocation subsystem
void cc_init_allocator(size_t max_mem);

// wrapper for malloc() from libc
void *cc_malloc(size_t nbytes);

// allocate memory and set all bytes to zero
void *cc_calloc(size_t num, size_t size);

// wrapper for free() from libc
void cc_free(void *p);

// finalization of memory allocation subsystem
void cc_finalize_allocator();

// prints memory usage statistics to stdout
void cc_memory_usage();

// returns information about memory usage
size_t cc_get_current_memory_usage();

size_t cc_get_peak_memory_usage();

char *cc_strdup(const char *src);

#endif /* CC_MEMORY_H_INCLUDED */

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

/*
 * Common definitions in types used in the EXP-T project.
 */

#ifndef CC_COMDEF_H_INCLUDED
#define CC_COMDEF_H_INCLUDED

#include <complex.h>
#include <stdint.h>
#include <stddef.h>

/* sometimes CUDA compiler knows nothing about the "complex" literal */
#define complex _Complex

// types
typedef int16_t moindex_t;

typedef int bool_t;

typedef double complex double_complex_t;
typedef float complex float_complex_t;

// limits (max size of static arrays)
#define MAX_SECTOR_RANK 4
#define CC_MAX_TITLE 256
#define CC_MAX_FILE_NAME_LENGTH 256
#define CC_DIAGRAM_MAX_NAME 64
#define CC_DIAGRAM_MAX_RANK 16
#define CC_MAX_SPINORS 2048
#define CC_MAX_PATH_LENGTH 1024
#define CC_MAX_DPD_STACK_SIZE 1024
#define CC_MAX_NUM_IRREPS 256

#endif /* CC_COMDEF_H_INCLUDED */

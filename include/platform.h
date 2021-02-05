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
* platform.h
* ==========
*
* Functions and definitions which make the code more compatible with a variety
* of compilers and platforms (platform- & compiler-specific code).
*
* 2018-2021 Alexander Oleynichenko
******************************************************************************/

#ifndef CC_PLATFORM_H_INCLUDED
#define CC_PLATFORM_H_INCLUDED

#include <stddef.h>
#include <stdint.h>

#if defined __ICC
#define COMPILER_INTEL
#elif defined __GNUC__
#define COMPILER_GNU
#else
#define COMPILER_OTHER
#endif

// detection of BLAS/LAPACK implementation
#if defined BLAS_Intel10_64ilp \
 || defined BLAS_Intel10_64lp  \
 || defined BLAS_Intel10_64ilp_seq  \
 || defined BLAS_Intel10_64lp_seq  \
 || defined BLAS_Intel10_32  \
 || defined BLAS_Intel
#define BLAS_MKL

#include "mkl.h"

#elif defined BLAS_OpenBLAS
#define BLAS_OPENBLAS
#include "cblas.h"
#else // other BLAS implementations
// macro __has_include is implemented in modern GNU and Intel compilers
#ifdef __has_include
#if __has_include("cblas.h")
#include "cblas.h"
#else
#define CBLAS_USE_OWN
#include "cblas-decl.h"
#endif
#define BLAS_OTHER
// __has_include is not defined
// try to include something...
#else
#define CBLAS_USE_OWN
#define LAPACKE_USE_OWN

#include "cblas-decl.h"

#define BLAS_OTHER
#endif
#endif

void print_sizeof_types();

void get_host_name(char *name, size_t len);

void print_blas_info();

int run_program(char *cmd, ...);

#ifdef COMPILER_INTEL
char *strdup(const char *s);
#endif /* COMPILER_INTEL */

#endif /* CC_PLATFORM_H_INCLUDED */

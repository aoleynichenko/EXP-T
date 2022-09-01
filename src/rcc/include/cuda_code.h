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

/*
 * Interface to the CUDA code.
 * Will be compiled only in case CUDA device is located and drivers
 * and libraries are configured properly.
 *
 * 2020-2021 Alexander Oleynichenko
 */

#ifndef CC_CUDA_H_INCLUDED
#define CC_CUDA_H_INCLUDED

#include <complex.h>

#ifdef CUDA_FOUND

int mulblocks_cuda(int carith,
                    double complex *A, int ma, int na,
                    double complex *B, int mb, int nb,
                    double beta, double complex *C);

int reorder_cuda(int carith, int rank, int *perm, int *dims1, int *dims2, double complex *v1, double complex *v2);

int cuda_device_query();

#endif /* CUDA_FOUND */

#endif /* CC_CUDA_H_INCLUDED */

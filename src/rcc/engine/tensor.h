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

#ifndef CC_TENSOR_H_INCLUDED
#define CC_TENSOR_H_INCLUDED

#include <complex.h>
#include <stddef.h>

size_t tensor_num_elements(int rank, const int *shape);

size_t tensor_index_to_linear(int rank, const int *shape, const int *idx);

void tensor_index_to_compound(int rank, const int *shape, size_t linear_idx, int *idx);

int tensor_index_in_range(int rank, const int *shape, const int *idx);

double real_tensor_get_element(int rank, double *tensor, int *shape, int *idx);

double complex complex_tensor_get_element(int rank, double complex *tensor, int *shape, int *idx);

void real_tensor_set_element(int rank, double *tensor, int *shape, int *idx, double value);

void complex_tensor_set_element(int rank, double complex *tensor, int *shape, int *idx, double complex value);

#endif // CC_TENSOR_H_INCLUDED

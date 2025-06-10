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

#include "tensor.h"

#include <complex.h>
#include <stddef.h>

/*
 * TODO: check of index range
 */


/**
 * returns total number of elements in a tensor
 */
size_t tensor_num_elements(int rank, const int *shape)
{
    size_t num = 1;

    for (int i = 0; i < rank; i++) {
        num *= shape[i];
    }

    return num;
}


/**
 * converts compound index into linear one
 * all indices must be >= 0
 * rank - number of dimensions in the tensor
 *
 * Example: 4-dimensional tensor with shape [0:L-1,0:M-1,0:N-1,0:K-1].
 * For element (l,i,j,k) linear index is equal to
 * (((M*l+i))*N+j*K+k ==
 *  == MNK*l + NK*i + K*j + k
 */
size_t tensor_index_to_linear(int rank, const int *shape, const int *idx)
{
    size_t linear_index = 0;

    for (int i = 0; i < rank; i++) {
        size_t coef = 1;
        for (int j = i + 1; j < rank; j++) {
            coef *= shape[j];
        }
        linear_index += coef * idx[i];
    }

    return linear_index;
}


/**
 * converts linear index into compound one (is inverse of the 'tensor_index_to_linear')
 * all indices must be >= 0.
 *
 * rank - number of dimensions in the tensor
 * lin_idx -- linear index (input)
 * idx[rank] -- compound index (output)
 */
void tensor_index_to_compound(int rank, const int *shape, size_t linear_idx, int *idx)
{
    int coefs[rank];

    // calculate transformation coefficients
    for (int i = 0; i < rank; i++) {
        coefs[i] = 1;
        for (int j = i + 1; j < rank; j++) {
            coefs[i] *= shape[j];
        }
    }

    for (int i = 0; i < rank; i++) {
        idx[i] = (int) (linear_idx / coefs[i]);
        linear_idx = linear_idx % coefs[i];
    }
}


/**
 * checks if compound index is in ranges or not
 */
int tensor_index_in_range(int rank, const int *shape, const int *idx)
{
    for (int i = 0; i < rank; i++) {
        if (idx[i] >= shape[i]) {
            return 0;
        }
    }

    return 1;
}


/**
 * get the value of the matrix element with compound index 'idx'.
 * (version for real tensors)
 */
double real_tensor_get_element(int rank, double *tensor, int *shape, int *idx)
{
    size_t i = tensor_index_to_linear(rank, shape, idx);
    return tensor[i];
}


/**
 * get the value of the matrix element with compound index 'idx'.
 * (version for complex tensors)
 */
double complex complex_tensor_get_element(int rank, double complex *tensor, int *shape, int *idx)
{
    size_t i = tensor_index_to_linear(rank, shape, idx);
    return tensor[i];
}


/**
 * set the value of the matrix element with compound index 'idx'.
 * (version for complex tensors)
 */
void real_tensor_set_element(int rank, double *tensor, int *shape, int *idx, double value)
{
    size_t i = tensor_index_to_linear(rank, shape, idx);
    tensor[i] = value;
}


/**
 * set the value of the matrix element with compound index 'idx'.
 * (version for real tensors)
 */
void complex_tensor_set_element(int rank, double complex *tensor, int *shape, int *idx, double complex value)
{
    size_t i = tensor_index_to_linear(rank, shape, idx);
    tensor[i] = value;
}

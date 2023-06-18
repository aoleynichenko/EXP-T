/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

/**
 * Multidimensional tensor transposition on GPU.
 * Current implementation can work only with 4-dim tensors
 * (extension is rather straightforward)
 */

#include "cuda_code.h"

#include <stdio.h>

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "engine.h"

extern "C" {

// must be descreased for Fermi (comp cap = 2.x)
#define BLOCK_SIZE 1024

// constant memory
__constant__ int64_t
        d_nsize;
__constant__ int32_t
        d_perm[CC_DIAGRAM_MAX_RANK];
__constant__ int64_t
        d_coef1[CC_DIAGRAM_MAX_RANK];
__constant__ int64_t
        d_coef2[CC_DIAGRAM_MAX_RANK];

/*******************************************************************************
 * reorder_kernel_rank4
 *
 * reorder -- kernel function for 4-dim tensors.
 * instructions are grouped in order to have as less register variables
 * as possible (in order to avoid calls to the local memory storage).
 * this can be of extreme importance for old architectures (Tesla, Fermi).
 * Note that the cached constant memory is intensively used.
 *
 * Do not try to understand anything, better see reorder.c (the same algorithm,
 * but with normal number of variables, much more readable).
 ******************************************************************************/
__global__ void reorder_kernel_rank4_complex(cuDoubleComplex *v1, cuDoubleComplex *v2)
{
    __shared__
    int16_t idx4[4 * BLOCK_SIZE];
    int64_t coef;

    int64_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= d_nsize) { return; }

    int64_t offs = 4 * threadIdx.x;

    // dim 0
    coef = d_coef1[0];
    idx4[offs + 0] = i / coef;
    i = i % coef;

    // dim 1
    coef = d_coef1[1];
    idx4[offs + 1] = i / coef;
    i = i % coef;

    // dim 2
    coef = d_coef1[2];
    idx4[offs + 2] = i / coef;
    i = i % coef;

    // dim 3
    idx4[offs + 3] = i;

    // use 'coef' variable instead of 'index2' variable
    coef = d_coef2[0] * idx4[offs + d_perm[0]] + d_coef2[1] * idx4[offs + d_perm[1]] +
           d_coef2[2] * idx4[offs + d_perm[2]] + idx4[offs + d_perm[3]];

    // copy matrix element from global memory to global
    i = threadIdx.x + blockIdx.x * blockDim.x;
    v2[coef /* index2 */] = v1[i];
}

__global__ void reorder_kernel_rank4_real(double *v1, double *v2)
{
    __shared__
    int16_t idx4[4 * BLOCK_SIZE];
    int64_t coef;

    int64_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= d_nsize) { return; }

    int64_t offs = 4 * threadIdx.x;

    // dim 0
    coef = d_coef1[0];
    idx4[offs + 0] = i / coef;
    i = i % coef;

    // dim 1
    coef = d_coef1[1];
    idx4[offs + 1] = i / coef;
    i = i % coef;

    // dim 2
    coef = d_coef1[2];
    idx4[offs + 2] = i / coef;
    i = i % coef;

    // dim 3
    idx4[offs + 3] = i;

    // use 'coef' variable instead of 'index2' variable
    coef = d_coef2[0] * idx4[offs + d_perm[0]] + d_coef2[1] * idx4[offs + d_perm[1]] +
           d_coef2[2] * idx4[offs + d_perm[2]] + idx4[offs + d_perm[3]];

    // copy matrix element from global memory to global
    i = threadIdx.x + blockIdx.x * blockDim.x;
    v2[coef /* index2 */] = v1[i];
}

int reorder_cuda(int carith, int rank, int *perm, int *dims1, int *dims2, double _Complex *v1, double _Complex *v2)
{
    int64_t coef1[CC_DIAGRAM_MAX_RANK];
    int64_t coef2[CC_DIAGRAM_MAX_RANK];
    int64_t nsize;
    //int block_size = 32;//engine_config.n_cuda_threads;
    int i, j, idim;
    cudaError_t cuerr;
    cuDoubleComplex *d_v1, *d_v2;
    int SZ_WORKING_TYPE = carith ? sizeof(double _Complex) : sizeof(double);

    if (rank != 4) {
        fprintf(stderr, "reorder on GPU is now available only for 4-dimensional tensors\n");
        fprintf(stderr, "please, contact alexvoleynichenko@mail.com the new extended code is needed\n");
        fprintf(stderr, "see %s for details\n", __FILE__);
        return 1;
    }

    // calculate total number of elements to be reordered
    nsize = 1;
    for (idim = 0; idim < rank; idim++) {
        nsize *= dims1[idim];
    }

    // prepare coefficients for recalculation: compound index <-> linear index
    for (i = 0; i < rank; i++) {
        coef1[i] = 1;
        coef2[i] = 1;
    }
    for (i = 0; i < rank - 1; i++) {
        coef1[i] = 1;
        coef2[i] = 1;
        for (j = i + 1; j < rank; j++) {
            coef1[i] *= dims1[j];
            coef2[i] *= dims2[j];
        }
    }

    // copy all constant data to the constant memory
    cuerr = cudaMemcpyToSymbol(d_nsize, &nsize, sizeof(int64_t) * 1, 0, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy 'nsize' from host to device: %s\n", cudaGetErrorString(cuerr));
        return 1;
    }
    cuerr = cudaMemcpyToSymbol(d_perm, perm, sizeof(int32_t) * rank, 0, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy 'perm' from host to device: %s\n", cudaGetErrorString(cuerr));
        return 1;
    }
    cuerr = cudaMemcpyToSymbol(d_coef1, coef1, sizeof(int64_t) * rank, 0, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy 'coef1' from host to device: %s\n", cudaGetErrorString(cuerr));
        return 1;
    }
    cuerr = cudaMemcpyToSymbol(d_coef2, coef2, sizeof(int64_t) * rank, 0, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy 'coef2' from host to device: %s\n", cudaGetErrorString(cuerr));
        return 1;
    }

    // alloc arrays on device
    cuerr = cudaMalloc((void **) &d_v1, nsize * SZ_WORKING_TYPE);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot allocate GPU memory for d_v1: %s\n",
                cudaGetErrorString(cuerr));
        fflush(stderr);
        return 1;
    }
    cuerr = cudaMalloc((void **) &d_v2, nsize * SZ_WORKING_TYPE);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot allocate GPU memory for d_v2: %s\n",
                cudaGetErrorString(cuerr));
        fflush(stderr);
        return 1;
    }

    // copy v1 to GPU
    cuerr = cudaMemcpy(d_v1, v1, nsize * SZ_WORKING_TYPE, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy data (v1) from host to device: %s\n",
                cudaGetErrorString(cuerr));
        fflush(stderr);
        return 1;
    }

    // launch kernel
    dim3 threads = dim3(BLOCK_SIZE, 1, 1);
    dim3 blocks = dim3((nsize + BLOCK_SIZE - 1) / BLOCK_SIZE, 1, 1);

    if (carith) {
        reorder_kernel_rank4_complex <<< blocks, threads >>>(d_v1, d_v2);
    }
    else {
        reorder_kernel_rank4_real <<< blocks, threads >>>((double *) d_v1, (double *) d_v2);
    }
    cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot launch kernel: %s\n",
                cudaGetErrorString(cuerr));
        fflush(stderr);
        return 1;
    }

    cuerr = cudaDeviceSynchronize();
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot synchronize device: %s\n",
                cudaGetErrorString(cuerr));
        fflush(stderr);
        return 1;
    }

    // get result (v2) from GPU
    cuerr = cudaMemcpy(v2, d_v2, nsize * SZ_WORKING_TYPE, cudaMemcpyDeviceToHost);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy result (v2) from device to host: %s\n",
                cudaGetErrorString(cuerr));
        fflush(stderr);
        return 1;
    }

    cudaFree(d_v1);
    cudaFree(d_v2);

    return 0;
}

} /* extern "C" */

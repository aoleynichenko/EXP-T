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
 * mulblocks_cuda.c
 *
 * Matrix multiplication on GPU (for tensor contractions).
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "cuda_code.h"

#include <complex.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cublas_v2.h>


/*******************************************************************************
 * mulblocks_cuda
 *
 * Performs matrix multiplication C = A * B^T on GPU using CUBLAS z/d-gemm.
 ******************************************************************************/
int mulblocks_cuda(int carith,
                   double complex *A, int ma, int na,
                   double complex *B, int mb, int nb,
                   double beta, double complex *C)
{
    int m = ma;
    int n = mb;
    int k = na;
    int size_A = m * k;
    int size_B = n * k;
    int size_C = m * n;
    int SZ_WORKING_TYPE = (carith ? sizeof(double complex) : sizeof(double));
    int nbytes_A = SZ_WORKING_TYPE * size_A;
    int nbytes_B = SZ_WORKING_TYPE * size_B;
    int nbytes_C = SZ_WORKING_TYPE * size_C;
    double complex zalpha = 1.0 + 0.0 * I;
    double complex zbeta = beta + 0.0 * I;
    cudaError_t cuerr;

    // allocate memory on device
    double complex *d_A, *d_B, *d_C;
    cuerr = cudaMalloc((void **) &d_A, nbytes_A);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot allocate GPU memory for the first operand (block): %s\n",
                cudaGetErrorString(cuerr));
        return 1;
    }
    cuerr = cudaMalloc((void **) &d_B, nbytes_B);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot allocate GPU memory for the second operand (block): %s\n",
                cudaGetErrorString(cuerr));
        return 1;
    }
    cuerr = cudaMalloc((void **) &d_C, nbytes_C);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot allocate GPU memory for the resulting block: %s\n",
                cudaGetErrorString(cuerr));
        return 1;
    }

    // copy A, B, C matrices to device
    cuerr = cudaMemcpy(d_A, A, nbytes_A, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy data (block A) from host to device: %s\n",
                cudaGetErrorString(cuerr));
    }
    cuerr = cudaMemcpy(d_B, B, nbytes_B, cudaMemcpyHostToDevice);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy data (block B) from host to device: %s\n",
                cudaGetErrorString(cuerr));
    }
    if (fabs(beta) > 1e-13) {    // copy 'C' only in case 'C' is really needed
        cuerr = cudaMemcpy(d_C, C, nbytes_C, cudaMemcpyHostToDevice);
        if (cuerr != cudaSuccess) {
            fprintf(stderr, "Cannot copy data (block C) from host to device: %s\n",
                    cudaGetErrorString(cuerr));
        }
    }

    // init CUBLAS handle
    cublasStatus_t cberr;
    cublasHandle_t handle;
    cberr = cublasCreate(&handle);
    if (cberr != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "Cannot create CUBLAS handle (err code = %d)\n", cberr);
        return 1;
    }

    // multiply C = A*B^T

    // setup timer
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, NULL);

    if (carith) {
        cberr = cublasZgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            (cuDoubleComplex * ) & zalpha, (cuDoubleComplex *) d_B, k, (cuDoubleComplex *) d_A, k,
                            (cuDoubleComplex * ) & zbeta, (cuDoubleComplex *) d_C, n);
    }
    else {
        cberr = cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, k,
                            (double *) &zalpha, (double *) d_B, k, (double *) d_A, k,
                            (double *) &zbeta, (double *) d_C, n);
    }
    if (cberr != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "Error launching cublasZgemm (err code = %s)\n", cberr);
        return 1;
    }

    cudaEventRecord(stop, NULL);

    // Wait for the stop event to complete
    cuerr = cudaEventSynchronize(stop);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot synchronize device: %s\n", cudaGetErrorString(cuerr));
        return 1;
    }
    float msecTotal = 0.0f;
    cudaEventElapsedTime(&msecTotal, start, stop);

    // copy result from device to host
    cuerr = cudaMemcpy(C, d_C, nbytes_C, cudaMemcpyDeviceToHost);
    if (cuerr != cudaSuccess) {
        fprintf(stderr, "Cannot copy data (block C) from device to host: %s\n",
                cudaGetErrorString(cuerr));
    }

    // finalize CUBLAS
    cberr = cublasDestroy(handle);
    if (cberr != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "Cannot destroy CUBLAS handle (err code = %d)\n", cberr);
        return 1;
    }

    // finalize
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}

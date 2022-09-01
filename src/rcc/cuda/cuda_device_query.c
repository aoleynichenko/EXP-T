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
 * device_query.c
 * ==============
 *
 * Extracts info about NVIDIA GPU(s) and prints them to stdout.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "cuda_code.h"

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "engine.h"


/******************************************************************************
 * device_query
 *
 * This code is based on example deviceQuery.cpp from CUDA samples.
 *****************************************************************************/
int cuda_device_query()
{
    int device_count;
    int device;
    int driver_version = 0;
    int runtime_version = 0;
    int block_size;
    struct cudaDeviceProp dev_prop;
    cudaError_t err_id;

    // get number of CUDA devices
    err_id = cudaGetDeviceCount(&device_count);
    if (err_id != cudaSuccess) {
        printf("cudaGetDeviceCount returned %d\n-> %s\n", err_id, cudaGetErrorString(err_id));
        printf("Result = FAIL\n");
        return 1;
    }
    if (device_count == 0) {
        printf("There are no available device(s) that support CUDA\n");
    }
    else {
        printf("Detected %d CUDA Capable device(s)\n", device_count);
    }

    // for each device: get device properties and print them to stdout

    for (device = 0; device < device_count; device++) {
        cudaSetDevice(device);
        cudaGetDeviceProperties(&dev_prop, device);

        printf("\nDevice %d: \"%s\"\n", device, dev_prop.name);

        cudaDriverGetVersion(&driver_version);
        cudaRuntimeGetVersion(&runtime_version);
        printf("  CUDA Driver Version / Runtime Version           : %d.%d / %d.%d\n",
               driver_version / 1000, (driver_version % 100) / 10,
               runtime_version / 1000, (runtime_version % 100) / 10);

        printf("  CUDA Capability Major/Minor version number      : %d.%d\n",
               dev_prop.major, dev_prop.minor);
        printf("  Total amount of global memory                   : %.0f MBytes \n",
               dev_prop.totalGlobalMem / 1048576.0f);
        printf("  Number of multiprocessors                       : %d\n",
               dev_prop.multiProcessorCount);
        printf("  GPU Max Clock rate (frequency)                  : %.0f MHz\n",
               dev_prop.clockRate * 1e-3);
        printf("  Maximum Texture Dimension Size (x,y,z)          : "
               "1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
               dev_prop.maxTexture1D, dev_prop.maxTexture2D[0],
               dev_prop.maxTexture2D[1], dev_prop.maxTexture3D[0],
               dev_prop.maxTexture3D[1], dev_prop.maxTexture3D[2]);
        printf("  Maximum Layered 1D Texture Size, (num) layers   : 1D=(%d), %d layers\n",
               dev_prop.maxTexture1DLayered[0], dev_prop.maxTexture1DLayered[1]);
        printf("  Maximum Layered 2D Texture Size, (num) layers   : 2D=(%d,%d), %d layers\n",
               dev_prop.maxTexture2DLayered[0], dev_prop.maxTexture2DLayered[1],
               dev_prop.maxTexture2DLayered[2]);
        printf("  Total amount of constant memory                 : %lu bytes\n",
               dev_prop.totalConstMem);
        printf("  Total amount of shared memory per block         : %lu bytes\n",
               dev_prop.sharedMemPerBlock);
        printf("  Total number of registers available per block   : %d\n",
               dev_prop.regsPerBlock);
        printf("  Warp size                                       : %d\n",
               dev_prop.warpSize);
        printf("  Maximum number of threads per multiprocessor    : %d\n",
               dev_prop.maxThreadsPerMultiProcessor);
        printf("  Maximum number of threads per block             : %d\n",
               dev_prop.maxThreadsPerBlock);
        printf("  Max dimension size of a thread block (x,y,z)    : (%d,%d,%d)\n",
               dev_prop.maxThreadsDim[0], dev_prop.maxThreadsDim[1],
               dev_prop.maxThreadsDim[2]);
        printf("  Max dimension size of a grid size    (x,y,z)    : (%d,%d,%d)\n",
               dev_prop.maxGridSize[0], dev_prop.maxGridSize[1],
               dev_prop.maxGridSize[2]);
        printf("  Maximum memory pitch                            : %lu bytes\n",
               dev_prop.memPitch);
        printf("  Texture alignment                               : %lu bytes\n",
               dev_prop.textureAlignment);
        printf("  Concurrent copy and kernel execution            : %s with %d copy "
               "engine(s)\n", (dev_prop.deviceOverlap ? "Yes" : "No"),
               dev_prop.asyncEngineCount);
        printf("  Run time limit on kernels                       : %s\n",
               dev_prop.kernelExecTimeoutEnabled ? "Yes" : "No");
        printf("  Integrated GPU sharing Host Memory              : %s\n",
               dev_prop.integrated ? "Yes" : "No");
        printf("  Support host page-locked memory mapping         : %s\n",
               dev_prop.canMapHostMemory ? "Yes" : "No");
        printf("  Alignment requirement for Surfaces              : %s\n",
               dev_prop.surfaceAlignment ? "Yes" : "No");
        printf("  Device has ECC support                          : %s\n",
               dev_prop.ECCEnabled ? "Enabled" : "Disabled");
        printf("  TCC (Tesla Compute Cluster Driver)              : %s\n",
               dev_prop.tccDriver ? "Yes" : "No");
        printf("  Device supports Unified Addressing (UVA)        : %s\n",
               dev_prop.unifiedAddressing ? "Yes" : "No");
        printf("  Device supports Compute Preemption              : %s\n",
               dev_prop.computePreemptionSupported ? "Yes" : "No");
        printf("  Supports Cooperative Kernel Launch              : %s\n",
               dev_prop.cooperativeLaunch ? "Yes" : "No");
        printf("  Supports MultiDevice Co-op Kernel Launch        : %s\n",
               dev_prop.cooperativeMultiDeviceLaunch ? "Yes" : "No");
        printf("  Device PCI Domain ID / Bus ID / location ID     : %d / %d / %d\n",
               dev_prop.pciDomainID, dev_prop.pciBusID, dev_prop.pciDeviceID);

        const char *sComputeMode[] = {
                "Default (multiple host threads can use ::cudaSetDevice() with device "
                "simultaneously)",
                "Exclusive (only one host thread in one process is able to use "
                "::cudaSetDevice() with this device)",
                "Prohibited (no host thread can use ::cudaSetDevice() with this "
                "device)",
                "Exclusive Process (many threads in one process is able to use "
                "::cudaSetDevice() with this device)",
                "Unknown",
                NULL};
        printf("  Compute Mode:\n");
        printf("     < %s >\n", sComputeMode[dev_prop.computeMode]);

        block_size = dev_prop.maxThreadsPerBlock;
    }

    printf("\n");

    // configure CC engine
    //engine_config.n_cuda_threads = block_size;

    return 0;
}

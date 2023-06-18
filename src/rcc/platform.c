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
 * Functions and definitions which make the code more compatible with a variety
 * of compilers and platforms (platform- & compiler-specific code).
 */

#include "platform.h"

#include <complex.h>
#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>

#include "error.h"
#include "timer.h"


/**
 * Returns string with hame of the host (server etc)
 * NOTE: array 'name' must be pre-allocated
 */
#if defined(__linux__) || defined(__APPLE__)

#include <unistd.h>

void get_host_name(char *name, size_t len)
{
    int status = gethostname(name, len);
    if (status != 0) {
        strcpy(name, "\n");
    }
}


#endif /* linux or macosx */


/**
 * Detects BLAS/LAPACK implementation and its version.
 */
void print_blas_info()
{
    printf("BLAS/LAPACK implementation: ");

#ifdef BLAS_MKL
    //printf("Intel MKL v %d.%d.%d", __INTEL_MKL__, __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__);
    printf("Intel MKL");
#elif defined BLAS_OPENBLAS
    printf("OpenBLAS");
#else /* BLAS_OTHER */
    printf("undetected");
#endif

#ifdef CBLAS_USE_OWN
    printf(" (own cblas.h file is used!)");
#endif

    printf("\n");
}


/**
 * Executes external command (rm, cp, mv, etc).
 *
 * NOTE: whis function doesn't check care about buffer cmd_buf overflow
 * TODO: fix it (not very necessary actually for our purposes)
 */
int execute_external_program(char *cmd, ...)
{
    // construct command line from arguments
    const size_t cmd_len = 1024;
    char cmd_buf[cmd_len];
    va_list vl;
    char *arg = NULL;

    strcpy(cmd_buf, cmd);
    va_start(vl, cmd);

    while ((arg = va_arg(vl, char *)) != NULL) {
        strcat(cmd_buf, " ");
        strcat(cmd_buf, arg);
    }

    va_end(vl);

    // flush all buffered output to stdout
    // this is done to guarantee right order of outputs of driver and subprograms
    fflush(stdout);

    // execute command
    int s = system(cmd_buf);
    if (s != 0) {
        errquit("execution of '%s' finished with non-zero status", cmd_buf);
        return -1;
    }
    printf("\n");
    return 0;
}




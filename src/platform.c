/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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
 * platform.c
 *
 * Functions and definitions which make the code more compatible with a variety
 * of compilers and platforms (platform- & compiler-specific code).
 *
 * 2018-2020 Alexander Oleynichenko
 ******************************************************************************/

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


/*******************************************************************************
 * print_sizeof_types
 *
 * Prints sizes (in bytes) of the most commonly used data types.
 ******************************************************************************/
void print_sizeof_types()
{
    printf("\n");
    printf("------------------------- size of data types (bytes) -------------------------\n");
    printf("                               integer types                                  \n");
    printf(" char           %2d  unsigned char  %2d   signed char    %2d   short          %2d\n",
           sizeof(char), sizeof(unsigned char), sizeof(signed char), sizeof(short));
    printf(" unsigned short %2d  int            %2d   unsigned int   %2d   long           %2d\n",
           sizeof(unsigned short), sizeof(int), sizeof(unsigned int), sizeof(long));
    printf(" unsigned long  %2d  long long      %2d   unsigned long long                 %2d\n",
           sizeof(unsigned long), sizeof(long long), sizeof(unsigned long long));
    printf(" int8_t         %2d  int16_t        %2d   int32_t        %2d   int64_t        %2d\n",
           sizeof(int8_t), sizeof(int16_t), sizeof(int32_t), sizeof(int64_t));
    printf(" uint8_t        %2d  uint16_t       %2d   uint32_t       %2d   uint64_t       %2d\n",
           sizeof(uint8_t), sizeof(uint16_t), sizeof(uint32_t), sizeof(uint64_t));
    printf("                      size, pointer, difference types                         \n");
    printf(" size_t         %2d  ptrdiff_t      %2d   void *         %2d\n",
           sizeof(size_t), sizeof(ptrdiff_t), sizeof(void *));

    printf("                                 real types                                   \n");
    printf(" float          %2d  double         %2d   long double    %2d\n",
           sizeof(float), sizeof(double), sizeof(long double));
    printf("                               complex types                                  \n");
    printf(" float complex  %2d  double complex %2d   long double complex                %2d\n",
           sizeof(float _Complex), sizeof(double _Complex), sizeof(long double _Complex));
    printf("------------------------------------------------------------------------------\n");
    printf("\n");
}


/*******************************************************************************
 * get_host_name
 *
 * Returns string with hame of the host (server etc)
 * NOTE: array 'name' must be pre-allocated
 ******************************************************************************/
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


/*******************************************************************************
 * blas_info
 *
 * Detects BLAS/LAPACK versions.
 ******************************************************************************/
void print_blas_info()
{
    printf("BLAS/LAPACK implementation: ");

#ifdef BLAS_MKL
    printf("Intel MKL v %d.%d.%d", __INTEL_MKL__, __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__);
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


/*******************************************************************************
 * run_program
 ******************************************************************************/

// NOTE: whis function doesn't check care about buffer cmd_buf overflow
// TODO: fix it (not very necessary actually for our purposes)
int run_program(char *cmd, ...)
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




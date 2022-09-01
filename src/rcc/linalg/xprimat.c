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

#include "linalg.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

void dprimat(const double *A, int n, int m, const char *comment);

void zprimat(const double complex *A, int n, int m, const char *comment);


/**
 * Prints matrix. Only one dot is printed instead of zeroes.
 *
 * @param type CC_DOUBLE or CC_COMPLEX
 * @param A    n by m matrix stored as 1D array
 * @param n    number of rows
 * @param m    number of columns
 * @param comment
 */
void xprimat(data_type_t type, const void *A, int n, int m, const char *comment)
{
    assert(type == CC_DOUBLE || type == CC_DOUBLE_COMPLEX);
    assert(A != NULL);
    assert(comment != NULL);
    assert(n > 0);
    assert(m > 0);

    if (type == CC_DOUBLE) {
        dprimat((double *) A, n, m, comment);
    }
    else { // type == CC_DOUBLE_COMPLEX
        zprimat((double complex *) A, n, m, comment);
    }
}


void dprimat(const double *A, int n, int m, const char *comment)
{
    int i, j;
    const int maxcol = 10;
    double thresh = 1e-14;
    int ncol;
    int n_remain;
    int n_done;
    double a_ij;

    printf(" [ %d x %d ] %s\n", n, m, comment);

    n_remain = m;
    n_done = 0;
    while (n_remain > 0) {
        ncol = (n_remain < maxcol) ? n_remain : maxcol;
        // print block's header
        printf("     ");
        for (i = 0; i < ncol; i++) {
            printf("              [%3d]           ", n_done + i);
        }
        printf("\n");
        // print block
        for (j = 0; j < n; j++) {
            printf(" [%3d] ", j);
            for (i = 0; i < ncol; i++) {
                a_ij = A[j * n + n_done + i];
                fabs(a_ij) < thresh ? printf("       .      ") : printf("%14.5e", a_ij);
                printf("  ");
            }
            printf("\n");
        }
        n_remain -= ncol;
        n_done += ncol;
    }
}


void zprimat(const double complex *A, int n, int m, const char *comment)
{
    int i, j;
    const int maxcol = 5;
    double thresh = 1e-14;
    int ncol;
    int n_remain;
    int n_done;
    double re, im;

    printf(" [ %d x %d ] %s\n", n, m, comment);

    n_remain = m;
    n_done = 0;
    while (n_remain > 0) {
        ncol = (n_remain < maxcol) ? n_remain : maxcol;
        // print block's header
        printf("     ");
        for (i = 0; i < ncol; i++) {
            printf("              [%3d]           ", n_done + i);
        }
        printf("\n");
        // print block
        for (j = 0; j < n; j++) {
            printf(" [%3d] ", j);
            for (i = 0; i < ncol; i++) {
                re = creal(A[j * n + n_done + i]);
                im = cimag(A[j * n + n_done + i]);
                fabs(re) < thresh ? printf("       .      ") : printf("%14.5e", re);
                fabs(im) < thresh ? printf("       .      ") : printf("%14.5e", im);
                printf("  ");
            }
            printf("\n");
        }
        n_remain -= ncol;
        n_done += ncol;
    }
}

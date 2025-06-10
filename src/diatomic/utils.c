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

#include "utils.h"

#include <stdlib.h>
#include <stdio.h>


double *new_double_array(int n, double *values)
{
    double *arr = (double *) calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        arr[i] = values[i];
    }

    return arr;
}


double **new_2d_array(int n, int m)
{
    double **A = (double **) calloc(n, sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *) calloc(m, sizeof(double));
    }

    return A;
}


void delete_2d_array(double **A, int n, int m)
{
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
}


double ***new_3d_array(int n, int m, int k)
{
    double ***A = (double ***) calloc(n, sizeof(double **));
    for (int i = 0; i < n; i++) {
        A[i] = (double **) calloc(m, sizeof(double *));
        for (int j = 0; j < m; j++) {
            A[i][j] = (double *) calloc(k, sizeof(double));
        }
    }

    return A;
}


void delete_3d_array(double ***A, int n, int m, int k)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            free(A[i][j]);
        }
        free(A[i]);
    }
    free(A);
}


/*
 * array := factor * array
 */
void rescale_array(int len, double *array, double factor)
{
    for (int i = 0; i < len; i++) {
        array[i] *= factor;
    }
}


/**
 * allocates square matrix filled with zeros
 */
double *new_zero_matrix(int n)
{
    return (double *) calloc(n * n, sizeof(double));
}


/**
 * allocates square tridiagonal n x n matrix.
 *
 * elements on the main diagonal are equal to 'diag_values',
 * and elements on the 1st and -1st diagonals are equal to 'offdiag_values'.
 * other elements are zero.
 */
double *new_tridiagonal_matrix(int n, double diag_values, double offdiag_values)
{
    double *A = (double *) calloc(n * n, sizeof(double));

    for (int i = 0; i < n; i++) {
        A[i * n + i] = diag_values;
    }
    for (int i = 1; i < n; i++) {
        A[(i - 1) * n + i] = offdiag_values;
        A[i * n + (i - 1)] = offdiag_values;
    }

    return A;
}


/**
 * prints n x n matrix to stdout
 */
void print_matrix(int n, double *A, char *comment)
{
    printf(" Matrix %dx%d: %s\n", n, n, comment);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.6f", A[i * n + j]);
        }
        printf("\n");
    }
}


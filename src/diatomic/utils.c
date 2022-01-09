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

#include "utils.h"

#include <stdlib.h>


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

/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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


#include "mdprop.h"

#include <complex.h>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "libunf.h"
#include "mrconee.h"

void analyze_complex_matrix(int dim, double _Complex *matrix,
                            int *re_zero, int *im_zero, int *re_symmetric, int *im_symmetric);

void analyze_nonzero_blocks(int dim, double _Complex *matrix, mrconee_data_t *mrconee_data);


void read_mdprop(char *path, mrconee_data_t *mrconee_data)
{
    printf("\n");
    printf("                **********************************\n");
    printf("                ** DIRAC integrals: MDPROP file **\n");
    printf("                **********************************\n");
    printf("\n");

    unf_file_t *file = unf_open(path, "r", UNF_ACCESS_SEQUENTIAL);
    if (file == NULL) {
        printf(" cannot open MDPPROP unformatted file\n");
        printf(" exp-t will be continue without any property calculations\n");
        return;
    }

    for (int prop_count = 1; unf_next_rec_size(file) != 0; prop_count++) {
        /*
         * name of the property
         */
        char oper_name[32];
        unf_read(file, "c32", oper_name);
        memmove(oper_name, oper_name + 24, 8);
        oper_name[8] = '\0';

        if (strcmp(oper_name, "EOFLABEL") == 0) {
            break;
        }

        printf("\n[%d]  %s\n\n", prop_count, oper_name);

        /*
         * get property matrix elements
         */
        int record_size = unf_next_rec_size(file);
        int nspinors = round(sqrt(record_size / (sizeof(double _Complex))));
        double _Complex *oper_matrix = (double _Complex *) calloc(nspinors * nspinors, sizeof(double _Complex));
        int n_matrix_elements = nspinors * nspinors;
        unf_read(file, "z8[i4]", oper_matrix, &n_matrix_elements);

        /*
         * analysis of a property matrix
         */
        int re_zero = 1;
        int im_zero = 1;
        int re_symm = 1;
        int im_symm = 1;

        analyze_complex_matrix(nspinors, oper_matrix, &re_zero, &im_zero, &re_symm, &im_symm);
        printf(" real part: ");
        if (re_zero) {
            printf("zero\n");
        }
        else {
            printf("non-zero %s\n", re_symm ? "symmetric" : "antisymmetric");
        }
        printf(" imag part: ");
        if (im_zero) {
            printf("zero\n");
        }
        else {
            printf("non-zero %s\n", im_symm ? "symmetric" : "antisymmetric");
        }

        if (mrconee_data) {
            analyze_nonzero_blocks(nspinors, oper_matrix, mrconee_data);
        }

        /*
         * cleanup
         */
        free(oper_matrix);
    }

    unf_close(file);

    printf("\n");
}


void analyze_complex_matrix(int dim, double _Complex *matrix,
                            int *re_zero, int *im_zero, int *re_symmetric, int *im_symmetric)
{
    const double zero_thresh = 1e-14;

    *re_zero = 1;
    *im_zero = 1;
    *re_symmetric = 1;
    *im_symmetric = 1;

    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            double _Complex a_ij = matrix[i * dim + j];
            double _Complex a_ji = matrix[j * dim + i];

            if (fabs(creal(a_ij)) > zero_thresh || fabs(creal(a_ji)) > zero_thresh) {
                *re_zero = 0;
            }

            if (fabs(cimag(a_ij)) > zero_thresh || fabs(cimag(a_ji)) > zero_thresh) {
                *im_zero = 0;
            }

            if (fabs(creal(a_ij) - creal(a_ji)) > zero_thresh) {
                *re_symmetric = 0;
            }
            if (fabs(cimag(a_ij) - cimag(a_ji)) > zero_thresh) {
                *im_symmetric = 0;
            }
        }
    }
}


void analyze_nonzero_blocks(int dim, double _Complex *matrix, mrconee_data_t *mrconee_data)
{
    const double zero_thresh = 1e-14;

    printf(" non-zero blocks:\n");

    for (int irep = 0; irep < mrconee_data->num_irreps; irep++) {
        for (int jrep = irep; jrep < mrconee_data->num_irreps; jrep++) {
            int non_zero = 0;

            for (int i = 0; i < mrconee_data->num_spinors; i++) {
                for (int j = 0; j < mrconee_data->num_spinors; j++) {
                    if (mrconee_data->spinor_irreps[i] != irep) {
                        continue;
                    }
                    if (mrconee_data->spinor_irreps[j] != jrep) {
                        continue;
                    }

                    if (cabs(matrix[i * dim + j]) > zero_thresh) {
                        non_zero = 1;
                    }
                }
            }

            if (non_zero) {
                printf(" %s - %s\n", mrconee_data->irrep_names[irep], mrconee_data->irrep_names[jrep]);
            }
        }
    }
}

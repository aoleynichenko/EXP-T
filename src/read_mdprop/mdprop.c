//
// Created by Alexander Oleynichenko on 29.12.2023.
//

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

#include "unformatted.h"
#include "mrconee.h"

void analyze_complex_matrix(int dim, double _Complex *matrix,
                            int *re_zero, int *im_zero, int *re_symmetric, int *im_symmetric);

void analyze_nonzero_blocks(int dim, double _Complex *matrix, mrconee_data_t *mrconee_data);


void read_mdprop(char *path, mrconee_data_t *mrconee_data)
{
    unformatted_file_t *file = unformatted_open(path, "r");
    if (file == NULL) {
        printf(" MDPROP file not found\n");
        return;
    }

    for (int prop_count = 1; unformatted_next_record_size(file) != 0; prop_count++) {
        /*
         * name of the property
         */
        char oper_name[32];
        unformatted_read(file, "c32", oper_name);
        memmove(oper_name, oper_name + 24, 8);
        oper_name[8] = '\0';

        if (strcmp(oper_name, "EOFLABEL") == 0) {
            break;
        }

        printf("\n[%d]  %s\n\n", prop_count, oper_name);

        /*
         * get property matrix elements
         */
        int record_size = unformatted_next_record_size(file);
        int nspinors = round(sqrt(record_size / (sizeof(double _Complex))));
        double _Complex *oper_matrix = (double _Complex *) calloc(nspinors * nspinors, sizeof(double _Complex));
        int n_matrix_elements = nspinors * nspinors;
        unformatted_read(file, "z8[i4]", oper_matrix, &n_matrix_elements);

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

    unformatted_close(file);

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

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

#include "read_dirac_output.h"


#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "errquit.h"

#define MAX_BUF_SIZE 1024


void fgets_skip_line(FILE *inp_file);

void parse_num_basis_functions(char *buf, int *n_basis_fun);

void parse_electronic_eigenvalue(char *buf, double *eigenvalue);

void parse_quaternion_coeffs_block(FILE *inp_file, double complex *coeffs_alpha, double complex *coeffs_beta,
                                   int n_basis_fun, char **basis_fun_labels);

void extract_integers(char *buf, int *integer_list, int *n_numbers);

int sum_integers(int *integer_list, int n_numbers);

char *my_strdup(char *s);


/*
 * extracts molecular spinors expanded in the basis of atomic orbitals
 * from the DIRAC output file.
 *
 * the following arrays are allocated inside the function:
 *  - eigenvalues
 *  - coeffs_alpha
 *  - coeffs_beta
 *  - basis_fun_labels
 */
void read_dirac_output(
        char *path, int *n_spinors, int *n_basis_fun, double **eigenvalues,
        double complex **coeffs_alpha, double complex **coeffs_beta,
        char ***basis_fun_labels
)
{
    char buf[MAX_BUF_SIZE];
    int integer_list[16];
    int n_numbers;

    const char *pattern_basis_dim = "  Number of large orbitals in each symmetry:";
    const char *pattern_n_spinors_axial = "No. of positive energy orbitals (NESH):"; // for axial symmetries
    const char *pattern_n_spinors_inversion = "No. of electronic orbitals (NESH):"; // for point groups with inversions
    const char *pattern_electronic = "* Electronic eigenvalue no.";

    printf(" > reading DIRAC output file:\n");
    printf("   %s\n", path);

    FILE *dirac_out = fopen(path, "r");
    if (dirac_out == NULL) {
        errquit("cannot open DIRAC output '%s'", path);
    }

    *n_basis_fun = 0;
    *n_spinors = 0;
    *eigenvalues = NULL;
    *coeffs_alpha = NULL;
    *coeffs_beta = NULL;

    int count_spinors = 0;

    while (fgets(buf, MAX_BUF_SIZE, dirac_out) != NULL) {

        // get total number of basis functions
        if (strstr(buf, pattern_basis_dim)) {
            extract_integers(buf, integer_list, &n_numbers);
            *n_basis_fun = sum_integers(integer_list, n_numbers);
            printf(" > total number of basis functions: %d\n", *n_basis_fun);
        }

            // get total number of spinors:
            // point groups without inversion
        else if (strstr(buf, pattern_n_spinors_axial)) {
            extract_integers(buf, integer_list, &n_numbers);
            *n_spinors = integer_list[0];
            printf(" > total number of spinors: %d\n", *n_spinors);
        }

            // get total number of spinors:
            // point groups with inversion
        else if (strstr(buf, pattern_n_spinors_inversion)) {
            extract_integers(buf, integer_list, &n_numbers);
            *n_spinors = integer_list[0] + integer_list[1];
            printf(" > number of gerade spinors: %d\n", integer_list[0]);
            printf(" > number of ungerade spinors: %d\n", integer_list[1]);
            printf(" > total number of spinors: %d\n", *n_spinors);
        }

            // read spinor (block of quaternion coefficients)
        else if (strstr(buf, pattern_electronic)) {

            // allocate memory for coefficients if needed
            if (*eigenvalues == NULL) {
                *eigenvalues = (double *) calloc(*n_spinors, sizeof(double));
                *coeffs_alpha = (double complex *) calloc((*n_basis_fun) * (*n_spinors), sizeof(double complex));
                *coeffs_beta = (double complex *) calloc((*n_basis_fun) * (*n_spinors), sizeof(double complex));
                *basis_fun_labels = (char **) calloc((*n_basis_fun), sizeof(char *));
            }

            // get spinor one-electron energy
            parse_electronic_eigenvalue(buf, *eigenvalues + count_spinors);
            fgets_skip_line(dirac_out);

            // get coefficients
            parse_quaternion_coeffs_block(
                    dirac_out,
                    *coeffs_alpha + (*n_basis_fun) * count_spinors,
                    *coeffs_beta + (*n_basis_fun) * count_spinors,
                    *n_basis_fun, *basis_fun_labels
            );

            count_spinors++;
        }
    }

    fclose(dirac_out);
}


/*
 * extract list of trailing integers from the string 'buf'
 */
void extract_integers(char *buf, int *integer_list, int *n_numbers)
{
    // go to the list of numbers
    buf = strchr(buf, ':') + 1;

    *n_numbers = 0;

    // parse list of integers
    char *pch = strtok(buf, " ");
    while (pch != NULL) {
        integer_list[*n_numbers] = atoi(pch);
        *n_numbers = *n_numbers + 1;

        pch = strtok(NULL, " ");
    }
}


/*
 * Example:
 * "* Electronic eigenvalue no.  1: -41.187345433878"
 */
void parse_electronic_eigenvalue(char *buf, double *eigenvalue)
{
    buf = strchr(buf, ':') + 1;
    *eigenvalue = atof(buf);
}


void parse_quaternion_coeffs_block(FILE *inp_file, double complex *coeffs_alpha, double complex *coeffs_beta,
                                   int n_basis_fun, char **basis_fun_labels)
{
    char buf[MAX_BUF_SIZE];

    while (fgets(buf, MAX_BUF_SIZE, inp_file) != NULL) {

        if (strcmp(buf, "\n") == 0) {
            // end of quaternion expansion
            break;
        }

        // index of the basis function
        char *pch = strtok(buf, " ");
        int index = atoi(pch) - 1;
        if (index >= n_basis_fun) {
            errquit("index of basis function exceeds the total size of the large basis\n"
                    "use the *PRIVEC/.AOLAB flag in the DIRAC input file while generating spinors");
        }

        // buffer to which the label of the basis function will be stored
        char label[MAX_BUF_SIZE];
        label[0] = '\0';

        // parse quaternion number
        double q[4];
        int iq = 0;
        pch = strtok(NULL, " ");

        while (pch != NULL) {

            if (strchr(pch, '.')) {
                // floating-point number
                q[iq++] = atof(pch);
            } else {
                // append this substring to the label of the basis function
                strcat(label, pch);
                strcat(label, " ");
            }

            pch = strtok(NULL, " ");
        }

        // save quaternion coefficient
        coeffs_alpha[index] = q[0] + I * q[1];
        coeffs_beta[index] = -q[2] + I * q[3];

        // save the label of the basis function (if needed)
        if (basis_fun_labels[index] == NULL) {
            basis_fun_labels[index] = my_strdup(label);
        }
    }
}


int sum_integers(int *integer_list, int n_numbers)
{
    int sum = 0;

    for (int i = 0; i < n_numbers; i++) {
        sum += integer_list[i];
    }

    return sum;
}


void fgets_skip_line(FILE *inp_file)
{
    char buf[MAX_BUF_SIZE];

    fgets(buf, MAX_BUF_SIZE, inp_file);
}


char *my_strdup(char *s)
{
    char *p = (char *) malloc(strlen(s) + 1); /* +1 for the '\0' */

    if (p != NULL) {
        strcpy(p, s);
    }

    return p;
}

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

#include "read_input.h"

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "errquit.h"
#include "units.h"
#include "input_data.h"


#define MAX_LINE_LEN 1024
#define MAX_NUM_POINTS 1024

void read_potential_curve(FILE *inp, int *n_points, double **r, double **pot1, double **pot2, double **prop);

void parse_pec_units(char *line, double *distance_conv_factor, double *energy_conv_factor);

int parse_doubles(char *line, double *values);

void rescale_array(int len, double *array, double factor);

int starts_with(char *str, char *pattern);

void remove_trailing_whitespaces(char *str);

char *my_strdup(char *s);

double *new_double_array(int n, double *values);


input_data_t *read_input(char *input_file)
{
    char line[MAX_LINE_LEN];
    char buf[MAX_LINE_LEN];
    int nread;
    double mass1, mass2;
    double distance_conv_factor = 1.0;
    double energy_conv_factor = 1.0;

    input_data_t *input_data = (input_data_t *) calloc(1, sizeof(input_data_t));

    /*
     * defaults
     */
    input_data->grid_size = 300;

    /*
     * open file and read it line-by-line
     */
    FILE *inp = fopen(input_file, "r");
    if (inp == NULL) {
        errquit("cannot open input file '%s'", input_file);
    }

    printf(" > echo of the input file\n\n");
    while (fgets(line, MAX_LINE_LEN, inp) != NULL) {
        printf("   %s", line);
    }
    printf("\n");

    rewind(inp);

    while (fgets(line, MAX_LINE_LEN, inp) != NULL) {

        /*
         * skip comments
         */
        if (starts_with(line, "#")) {
            continue;
        }

        /*
         * atomic masses
         */
        else if (starts_with(line, "masses")) {
            nread = sscanf(line, "%s%lf%lf", buf, &mass1, &mass2);
            if (nread != 3) {
                errquit("syntax error in 'masses' keyword");
            }
            input_data->reduced_mass = mass1 * mass2 / (mass1 + mass2) * AMU_TO_ELECTRON_MASS;
        }

        /*
         * lower and upper vibrational quantum number (v)
         */
        else if (starts_with(line, "viblevels")) {
            int lower, upper;

            nread = sscanf(line, "%s%d%d", buf, &lower, &upper);
            if (nread != 3) {
                errquit("syntax error in the 'viblevels' keyword");
            }
            if (lower < 0 || upper < 0) {
                errquit("vibrational numbers should be non-negative");
            }
            if (lower > upper) {
                errquit("wrong order of vibrational numbers");
            }

            input_data->v_min = lower;
            input_data->v_max = upper;
        }

        /*
         * lower and upper rotational quantum number (J)
         */
        else if (starts_with(line, "rotlevels")) {
            int lower, upper;

            nread = sscanf(line, "%s%d%d", buf, &lower, &upper);
            if (nread != 3) {
                errquit("syntax error in the 'rotlevels' keyword");
            }
            if (lower < 0 || upper < 0) {
                errquit("rotational numbers should be non-negative");
            }
            if (lower > upper) {
                errquit("wrong order of rotational numbers");
            }

            input_data->J_min = lower;
            input_data->J_max = upper;
        }

        /*
         * number of points in the integration grid
         */
        else if (starts_with(line, "grid_size")) {
            int grid_size;

            nread = sscanf(line, "%s%d", buf, &grid_size);
            if (nread != 2) {
                errquit("syntax error in the 'grid_size' keyword");
            }
            if (grid_size <= 0) {
                errquit("grid size must be positive");
            }

            input_data->grid_size = grid_size;
        }

        /*
         * potential energy curve U(r)
         */
        else if (starts_with(line, "potential")) {

            int n_points = 0;
            double *r = NULL;
            double *pot1 = NULL;
            double *pot2 = NULL;
            double *prop = NULL;

            // check if units are specified
            parse_pec_units(line, &distance_conv_factor, &energy_conv_factor);

            // read points
            read_potential_curve(inp, &n_points, &r, &pot1, &pot2, &prop);

            // transform curve to atomic units (bohrs, Hartrees)
            rescale_array(n_points, r, distance_conv_factor);
            rescale_array(n_points, pot1, energy_conv_factor);
            if (pot2) {
                rescale_array(n_points, pot2, energy_conv_factor);
            }

            // potentials and property curves are stored as splines
            input_data->n_points = n_points;
            input_data->r = r;
            input_data->pot1 = construct_cubic_spline(n_points, r, pot1);
            if (pot2) {
                input_data->pot2 = construct_cubic_spline(n_points, r, pot2);
            }
            if (prop) {
                input_data->prop = construct_cubic_spline(n_points, r, prop);
            }

            free(pot1);
            free(pot2);
            free(prop);
        }
    }

    fclose(inp);

    return input_data;
}


void read_potential_curve(FILE *inp, int *n_points, double **r, double **pot1, double **pot2, double **prop)
{
    char line[MAX_LINE_LEN];
    double buf_r[MAX_NUM_POINTS];
    double buf_pot1[MAX_NUM_POINTS];
    double buf_pot2[MAX_NUM_POINTS];
    double buf_prop[MAX_NUM_POINTS];
    double double_buf[10];

    /*
     * by default, we have only 2 columns: radial grid and potential
     * optional: 3rd column - propetry matrix element (expectation values will be calculated)
     * optional: 4th column - excited state potential (transition moments will be calculated)
     */
    int ncols = 2;

    /*
     * read PEC points until the 'end' keyword
     */
    *n_points = 0;

    while (fgets(line, MAX_LINE_LEN, inp) != NULL) {

        if (starts_with(line, "end")) {
            break;
        }

        int count = parse_doubles(line, double_buf);

        // the first line determines the number of columns
        if (*n_points == 0) {
            ncols = count;
            if (ncols < 2 || ncols > 4) {
                errquit("while reading potential energy curve");
            }
        }

        if (count != ncols) {
            errquit("while reading potential energy curve");
        }

        buf_r[*n_points] = double_buf[0];
        buf_pot1[*n_points] = double_buf[1];
        buf_prop[*n_points] = double_buf[2];
        buf_pot2[*n_points] = double_buf[3];

        *n_points = *n_points + 1;
    }

    /*
     * store potential curves
     */
    *r = new_double_array(*n_points, buf_r);
    *pot1 = new_double_array(*n_points, buf_pot1);
    if (ncols >= 3) {
        *prop = new_double_array(*n_points, buf_prop);
    }
    if (ncols == 4) {
        *pot2 = new_double_array(*n_points, buf_pot2);
    }
}


void parse_pec_units(char *line, double *distance_conv_factor, double *energy_conv_factor)
{
    *distance_conv_factor = ATOMIC_TO_ATOMIC;
    *energy_conv_factor = ATOMIC_TO_ATOMIC;

    int count = 0;
    char *pch = strtok(line, " \t");
    while (pch != NULL) {

        remove_trailing_whitespaces(pch);

        /*
         * units of distance
         */
        if (count == 1) {
            if (strcmp(pch, "angstrom") == 0) {
                *distance_conv_factor = ANGSTROM_TO_ATOMIC;
            }
            else if (strcmp(pch, "atomic") == 0 ||
                     strcmp(pch, "bohr") == 0 ||
                     strcmp(pch, "au") == 0) {
                // nothing to do
            }
            else {
                errquit("wrong units of distance: '%s'\n"
                        "allowed units are 'angstrom', 'atomic', 'bohr', 'au'",
                        pch);
            }
        }

        /*
         * units of energy
         */
        else if (count == 2) {
            if (strcmp(pch, "cm-1") == 0 ||
                strcmp(pch, "cm") == 0) {
                *energy_conv_factor = CM_TO_ATOMIC;
            }
            else if (strcmp(pch, "atomic") == 0) {
                // nothing to do
            }
            else {
                errquit("wrong units of energy: '%s'\n"
                        "allowed units are 'cm-1', 'atomic'",
                        pch);
            }
        }

        pch = strtok(NULL, " ");
        count++;
    }
}


int parse_doubles(char *line, double *values)
{
    int count = 0;
    char *buf = my_strdup(line);

    char *pch = strtok(buf, " \t");
    while (pch != NULL) {

        int nread = sscanf(pch, "%lf", values + count);
        if (nread != 1) {
            return -1;
        }

        pch = strtok(NULL, " ");
        count++;
    }

    free(buf);
    return count;
}


void rescale_array(int len, double *array, double factor)
{
    for (int i = 0; i < len; i++) {
        array[i] *= factor;
    }
}


int starts_with(char *str, char *pattern)
{
    char *p = strstr(str, pattern);
    if (p != NULL && p == str) {
        return 1;
    }

    return 0;
}


void remove_trailing_whitespaces(char *str)
{
    for (int i = strlen(str) - 1; i >= 0; i--) {
        if (isspace(str[i])) {
            str[i] = '\0';
        }
    }
}


char *my_strdup(char *s)
{
    char *p = (char *) malloc(strlen(s) + 1); /* +1 for the '\0' */

    if (p != NULL) {
        strcpy(p, s);
    }

    return p;
}


double *new_double_array(int n, double *values)
{
    double *arr = (double *) calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        arr[i] = values[i];
    }

    return arr;
}
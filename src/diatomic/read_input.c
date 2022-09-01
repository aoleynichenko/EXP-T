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
#include "utils.h"


#define MAX_LINE_LEN 1024
#define MAX_NUM_POINTS 1024

void parse_option_masses(input_data_t *input_data, char *line);

void parse_option_charge(input_data_t *input_data, char *line);

void parse_option_viblevels(input_data_t *input_data, char *line);

void parse_option_rotlevels(input_data_t *input_data, char *line);

void parse_option_grid_size(input_data_t *input_data, char *line);

void parse_option_solver(input_data_t *input_data, char *line);

void parse_option_mapping(input_data_t *input_data, char *line);

void parse_option_potential(input_data_t *input_data, FILE *file, char *line);

void read_potential_curve(FILE *inp, int *n_points, double **r, double **pot1, double **pot2, double **prop);

void parse_pec_units(char *line, double *distance_conv_factor, double *energy_conv_factor);

int parse_doubles(char *line, double *values);

int starts_with(char *str, char *pattern);

void remove_trailing_whitespaces(char *str);

char *my_strdup(char *s);

double charge_modified_reduced_mass(double mass1, double mass2, int charge);


/**
 * reads user-specified options (masses, potentials, etc) from the formatted text file.
 * returns structure containing all the necessary information.
 */
input_data_t *read_input(char *input_file)
{
    char line[MAX_LINE_LEN];

    input_data_t *input_data = (input_data_t *) calloc(1, sizeof(input_data_t));

    /*
     * defaults
     */
    input_data->write_psi = 0;
    input_data->grid_size = 300;
    input_data->solver = SOLVER_NUMEROV;
    input_data->mapping = new_mapping(MAPPING_IDENTITY, NULL);

    /*
     * open file and read it line-by-line
     */
    FILE *inp_file = fopen(input_file, "r");
    if (inp_file == NULL) {
        errquit("cannot open input file '%s'", input_file);
    }

    printf(" > echo of the input file\n\n");
    while (fgets(line, MAX_LINE_LEN, inp_file) != NULL) {
        printf("   %s", line);
    }
    printf("\n");

    rewind(inp_file);

    while (fgets(line, MAX_LINE_LEN, inp_file) != NULL) {

        /*
         * skip comments
         */
        if (starts_with(line, "#")) {
            continue;
        }

        /*
         * atomic masses
         */
        if (starts_with(line, "masses")) {
            parse_option_masses(input_data, line);
        }

        /*
         * net charge of the molecule
         */
        if (starts_with(line, "charge")) {
            parse_option_charge(input_data, line);
        }

        /*
         * lower and upper vibrational quantum number (v)
         */
        if (starts_with(line, "viblevels")) {
            parse_option_viblevels(input_data, line);
        }

        /*
         * lower and upper rotational quantum number (J)
         */
        if (starts_with(line, "rotlevels")) {
            parse_option_rotlevels(input_data, line);
        }

        /*
         * number of points in the integration grid
         */
        if (starts_with(line, "grid_size")) {
            parse_option_grid_size(input_data, line);
        }

        /*
         * solver: Numerov or finite-difference
         */
        if (starts_with(line, "solver")) {
            parse_option_solver(input_data, line);
        }

        /*
         * mapping of the integration variable
         */
        if (starts_with(line, "mapping")) {
            parse_option_mapping(input_data, line);
        }

        /*
         * write vibrational wavefunctions to the formatted files or not
         */
        if (starts_with(line, "write_psi")) {
            input_data->write_psi = 1;
        }

        /*
         * potential energy curve U(r)
         */
        if (starts_with(line, "potential")) {
            parse_option_potential(input_data, inp_file, line);
        }
    }

    /*
     * final validation of the options
     */
    if (input_data->solver == SOLVER_NUMEROV && input_data->mapping->type != MAPPING_IDENTITY) {
        errquit("mapping is currently not available for the Numerov integrator, use the FD2 integrator instead");
    }

    /*
     * post-processing of input parameters
     */
    input_data->reduced_mass = charge_modified_reduced_mass(input_data->mass1, input_data->mass2, input_data->charge);

    fclose(inp_file);

    return input_data;
}


/**
 * atomic masses
 *
 * syntax:
 * masses <real mass1> <real mass2>
 */
void parse_option_masses(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    double mass1, mass2;

    int nread = sscanf(line, "%s%lf%lf", buf, &mass1, &mass2);
    if (nread != 3) {
        errquit("syntax error in 'masses' keyword");
    }

    input_data->mass1 = mass1;
    input_data->mass2 = mass2;
}


/**
 * net charge of the molecule
 *
 * syntax:
 * charge <integer net-charge>
 */
void parse_option_charge(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    int charge;

    int nread = sscanf(line, "%s%d", buf, &charge);
    if (nread != 2) {
        errquit("syntax error in the 'charge' keyword");
    }

    input_data->charge = charge;
}


/**
 * lower and upper vibrational quantum number (v)
 *
 * syntax:
 * viblevels <integer vmin> <integer vmax>
 */
void parse_option_viblevels(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    int lower, upper;

    int nread = sscanf(line, "%s%d%d", buf, &lower, &upper);
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


/**
 * lower and upper rotational quantum number (J)
 *
 * syntax:
 * rotlevels <integer Jmin> <integer Jmax>
 */
void parse_option_rotlevels(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    int lower, upper;

    int nread = sscanf(line, "%s%d%d", buf, &lower, &upper);
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


/**
 * number of points in the integration grid
 *
 * syntax:
 * grid_size <integer number-of-points>
 */
void parse_option_grid_size(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    int grid_size;

    int nread = sscanf(line, "%s%d", buf, &grid_size);
    if (nread != 2) {
        errquit("syntax error in the 'grid_size' keyword");
    }
    if (grid_size <= 0) {
        errquit("grid size must be positive");
    }

    input_data->grid_size = grid_size;
}


/**
 * solver: Numerov or finite-difference
 *
 * syntax:
 * solver [ numerov | fd2 ]
 */
void parse_option_solver(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    char solver_name[MAX_LINE_LEN];

    int nread = sscanf(line, "%s%s", buf, solver_name);
    if (nread != 2) {
        errquit("syntax error in the 'solver' keyword");
    }
    if (strcmp(solver_name, "numerov") == 0) {
        input_data->solver = SOLVER_NUMEROV;
    } else if (strcmp(solver_name, "fd2") == 0) {
        input_data->solver = SOLVER_FD2;
    } else {
        errquit("unknown solver: %s", solver_name);
    }
}


/**
 * mapping of the integration variable
 *
 * syntax:
 * mapping meshkov08 <double ra (angstroms)> <double beta>
 */
void parse_option_mapping(input_data_t *input_data, char *line)
{
    char buf[MAX_LINE_LEN];
    char mapping_type[MAX_LINE_LEN];

    int nread = sscanf(line, "%s%s", buf, mapping_type);
    if (nread != 2) {
        errquit("syntax error in the 'mapping' keyword");
    }
    if (strcmp(mapping_type, "meshkov08") == 0) {
        double params[2];

        nread = sscanf(line, "%s%s%lf%lf", buf, mapping_type, &params[0], &params[1]);
        if (nread != 4) {
            errquit("syntax error in the 'mapping/meshkov08' option");
        }

        params[0] *= ANGSTROM_TO_ATOMIC;

        if (input_data->mapping != NULL) {
            free(input_data->mapping);
        }
        input_data->mapping = new_mapping(MAPPING_MESHKOV_08, params);
    } else {
        errquit("unknown mapping type: %s", mapping_type);
    }
}


/**
 * potential energy curve U(r)
 *
 * syntax:
 * potential [<units of r>] [<units of energy>]
 * <r>  <U1> [<property(R)>] [<U2>]
 * ...
 * end
 */
void parse_option_potential(input_data_t *input_data, FILE *file, char *line)
{
    double distance_conv_factor = 1.0;
    double energy_conv_factor = 1.0;

    int n_points = 0;
    double *r = NULL;
    double *pot1 = NULL;
    double *pot2 = NULL;
    double *prop = NULL;

    // check if units are specified
    parse_pec_units(line, &distance_conv_factor, &energy_conv_factor);

    // read points
    read_potential_curve(file, &n_points, &r, &pot1, &pot2, &prop);

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


/**
 * extracts units of distance and energy from the 'line' string.
 *
 * allowed units are:
 *   distance: angstrom, atomic = bohr = au
 *   energy:   cm-1 = cm, atomic = hartree
 */
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
            } else if (strcmp(pch, "atomic") == 0 ||
                       strcmp(pch, "bohr") == 0 ||
                       strcmp(pch, "au") == 0) {
                // nothing to do
            } else {
                errquit("wrong units of distance: '%s'\n"
                        "allowed units are 'angstrom', 'atomic', 'bohr', 'au'",
                        pch);
            }
        }

        /*
         * units of energy
         */
        if (count == 2) {
            if (strcmp(pch, "cm-1") == 0 ||
                strcmp(pch, "cm") == 0) {
                *energy_conv_factor = CM_TO_ATOMIC;
            } else if (strcmp(pch, "atomic") == 0 ||
                       strcmp(pch, "hartree") == 0) {
                // nothing to do
            } else {
                errquit("wrong units of energy: '%s'\n"
                        "allowed units are 'cm-1', 'cm', 'atomic', 'hartree'",
                        pch);
            }
        }

        pch = strtok(NULL, " ");
        count++;
    }
}


/**
 * calculates the Watson's "charge-modified reduced mass".
 * masses are assumed to be in amu.
 * returns the value in the electron mass units.
 *
 * see also:
 * [1] J. K. G. Watson, The isotope dependence of diatomic Dunham coefficients.
 *     J. Mol. Spectrosc. 80, 411 (1980)
 *     doi: 10.1016/0022-2852(80)90152-6
 * [2] R. J. Le Roy, LEVEL: A computer program for solving the radial
 *     Schrödinger equation for bound and quasibound levels.
 *     JQSRT, 186, 167 (2017)
 *     doi: 10.1016/j.jqsrt.2016.05.028
 */
double charge_modified_reduced_mass(double mass1, double mass2, int charge)
{
    mass1 = mass1 * AMU_TO_ELECTRON_MASS;
    mass2 = mass2 * AMU_TO_ELECTRON_MASS;

    return mass1 * mass2 / (mass1 + mass2 - charge);
}


int parse_doubles(char *line, double *values)
{
    int count = 0;
    char *buf = my_strdup(line);
    remove_trailing_whitespaces(buf);

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
        } else {
            return;
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



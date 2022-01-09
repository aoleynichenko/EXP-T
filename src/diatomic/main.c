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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matrix_element.h"
#include "read_input.h"
#include "cubic_spline.h"
#include "harmonic.h"
#include "morse.h"
#include "numerov.h"
#include "errquit.h"
#include "units.h"
#include "string.h"
#include "utils.h"


#define MAX(a, b) (((a)>(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))

void print_header();

void print_footer();

void parse_arguments(int argc, char **argv, char **path_input_file);

void print_energy_levels(input_data_t *input_data, cubic_spline_t *pot, int *nroots, double **energies,
                         double ***wavefunctions);

void calc_transition_moments(input_data_t *input_data,
                             int *nroots1, double **energies1, double ***wavefunctions1,
                             int *nroots2, double **energies2, double ***wavefunctions2);

int main(int argc, char **argv)
{
    /*
     * parse command line arguments
     */

    char *path_input_file = NULL;
    parse_arguments(argc, argv, &path_input_file);

    print_header();

    /*
     * read input file
     */
    input_data_t *input_data = read_input(path_input_file);
    print_input_data(input_data);

    /*
     * solve rovibrational problem for the ground state potential
     */
    int *nroots1 = NULL;
    double emin1, re1;
    double **energies1 = NULL;
    double ***wavefunctions1 = NULL;
    cubic_spline_t *pot1 = input_data->pot1;

    printf(" > begin ground state potential\n\n");

    printf(" > equilibrium parameters\n\n");
    find_spline_minimum(pot1, &re1, &emin1, 0);
    input_data->min_energy = emin1;
    printf("   e(min) = %.12f a.u.\n", emin1);
    printf("   re     = %.4f a.u. = %.4f angstrom\n", re1, re1 * ATOMIC_TO_ANGSTROM);
    printf("\n");

    harmonic_analysis(input_data, pot1, re1);
    morse_analysis(input_data, pot1, re1);

    solve(input_data, pot1, &nroots1, &energies1, &wavefunctions1);
    print_energy_levels(input_data, pot1, nroots1, energies1, wavefunctions1);

    /*
     * solve rovibrational problem for the excited state potential
     * (if required)
     */
    if (input_data->pot2 != NULL) {
        int *nroots2 = NULL;
        double emin2, re2;
        double **energies2 = NULL;
        double ***wavefunctions2 = NULL;
        cubic_spline_t *pot2 = input_data->pot2;

        printf(" > begin excited state potential\n\n");

        printf(" > equilibrium parameters\n\n");
        find_spline_minimum(pot2, &re2, &emin2, 0);
        printf("   e(min) = %.12f a.u.\n", emin2);
        printf("   re     = %.4f a.u. = %.4f angstrom\n", re2, re2 * ATOMIC_TO_ANGSTROM);
        printf("   Te     = %.4f cm^-1\n", (emin2 - emin1) * ATOMIC_TO_CM);
        printf("\n");

        harmonic_analysis(input_data, pot2, re2);
        morse_analysis(input_data, pot2, re2);

        solve(input_data, pot2, &nroots2, &energies2, &wavefunctions2);
        print_energy_levels(input_data, pot2, nroots2, energies2, wavefunctions2);

        /*
         * transition moments between ground and excited states
         */
        calc_transition_moments(input_data,
                                nroots1, energies1, wavefunctions1,
                                nroots2, energies2, wavefunctions2);
    }

    print_footer();
}


/*
 * prints beautiful table with energies of the J-v states
 * and expectation values of properties
 */
void print_energy_levels(input_data_t *input_data, cubic_spline_t *pot, int *nroots, double **energies,
                         double ***wavefunctions)
{
    int v_min = input_data->v_min;
    int v_max = input_data->v_max;
    int J_min = input_data->J_min;
    int J_max = input_data->J_max;

    double r_min = pot->x[0];
    double r_max = pot->x[pot->n - 1];
    int n_rot = J_max - J_min + 1;
    int n_vib = v_max - v_min + 1;
    int ngrid = input_data->grid_size;
    double emin = input_data->min_energy;

    cubic_spline_t *r_spline = construct_cubic_spline(input_data->n_points, input_data->r, input_data->r);

    printf(" > rovibrational energy levels\n\n");
    printf("     J   v       energy, cm^-1        < r >, A        < prop >\n\n");

    for (int J = J_min; J <= J_max; J++) {
        int index_J = J - J_min;

        for (int v = v_min; v <= MIN(v_max, nroots[index_J]); v++) {
            int index_v = v - v_min;

            double E_Jv = energies[index_J][index_v];
            double *psi = wavefunctions[index_J][index_v];

            double e_cm = (E_Jv - emin) * ATOMIC_TO_CM;
            double r_v = matrix_element_spline(ngrid, r_min, r_max, psi, psi, r_spline);
            double prp = 0.0;
            if (input_data->prop != NULL) {
                prp = matrix_element_spline(ngrid, r_min, r_max, psi, psi, input_data->prop);
            }

            printf("@ %4d%4d%20.4f%16.6f%16.6f\n", J, v, e_cm, r_v * ATOMIC_TO_ANGSTROM, prp);
        }
        printf("\n");
    }

    delete_spline(r_spline);
}


void print_header()
{
    char buf[100];

    sprintf(buf, "%s %s", __DATE__, __TIME__);

    printf("\n");
    printf("\t**********************************************************************************\n"
           "\t**                                                                              **\n"
           "\t**           program for rotational-vibrational levels and properties           **\n"
           "\t**                              of diatomic molecules                           **\n"
           "\t**                                                                              **\n"
           "\t**                              2022 A. Oleynichenko                            **\n"
           "\t**                       compiled %-30s                **\n"
           "\t**                                                                              **\n"
           "\t**********************************************************************************\n\n",
           buf);

    time_t lt = time(NULL);
    struct tm *ptr = localtime(&lt);
    printf(" > started at %s\n", asctime(ptr));
}


void print_footer()
{
    time_t lt = time(NULL);
    struct tm *ptr = localtime(&lt);
    printf(" > finished at %s\n", asctime(ptr));
}


void parse_arguments(int argc, char **argv, char **path_input_file)
{
    if (argc != 2) {
        printf("\n");
        printf(" usage: %s <input-file>\n", argv[0]);
        printf("\n");
        exit(1);
    }

    *path_input_file = argv[1];
}



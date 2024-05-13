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

#include "harmonic.h"

#include <math.h>
#include <stdio.h>

#include "cubic_spline.h"
#include "input_data.h"
#include "units.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double rotational_constant(double r_e, double reduced_mass);

double distortion_constant(double r_e, double reduced_mass, double we);


/**
 * based on the given potential, calculates molecular constants
 * in the harmonic + rigid rotor approximation
 */
void harmonic_analysis(input_data_t *data, cubic_spline_t *pot, double rmin)
{
    printf(" > analysis in the harmonic oscillator + rigid rotor approximation\n\n");

    double we;
    calculate_harmonic_frequency(pot, rmin, data->reduced_mass, &we);
    printf("   we = %.4f cm^-1\n", we * ATOMIC_TO_CM);

    double Be = rotational_constant(rmin, data->reduced_mass);
    printf("   Be = %.6f cm^-1\n", Be * ATOMIC_TO_CM);

    double De = distortion_constant(rmin, data->reduced_mass, we);
    printf("   De = %.6f cm^-1\n", De * ATOMIC_TO_CM);

    printf("\n");
}


/**
 * returns harmonic frequency w_e (in atomic units of energy)
 */
void calculate_harmonic_frequency(cubic_spline_t *pot, double rmin, double reduced_mass, double *we)
{
    double k = spline_second_derivative(pot, rmin);
    *we = sqrt(k / reduced_mass);
}


/**
 * returns rigid-rotor value of the rotational constant (in atomic units of energy)
 * B = 1/(2*mu*r^2)
 */
double rotational_constant(double r_e, double reduced_mass)
{
    double I = reduced_mass * r_e * r_e;
    return 1.0 / (2 * I);
}


/**
 * returns the centrifugal distortion constant (in atomic units of energy)
 * D = 4 B^3 / w^2
 */
double distortion_constant(double r_e, double reduced_mass, double we)
{
    double B = rotational_constant(r_e, reduced_mass);
    return 4.0 * B * B * B / (we * we);
}

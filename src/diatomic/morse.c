/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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

#include "morse.h"

#include <math.h>
#include <stdio.h>

#include "units.h"


/**
 * Morse potential energy function.
 *
 * @param r radial variable
 * @param re equilibrium internuclear separation
 * @param a alpha parameter
 * @param De dissociation energy
 */
double morse_potential(double r, double re, double a, double De)
{
    return De * (exp(-2.0*a*(r-re)) - 2.0 * exp(-a*(r-re)));
}


/**
 * estimates we, wexe, De, vmax parameters of the potential using
 * the analytic formulas for the Morse oscillator
 */
void morse_analysis(input_data_t *data, cubic_spline_t *pot, double rmin)
{
    printf(" > analysis in the morse approximation\n\n");

    double y2 = spline_second_derivative(pot, rmin);
    double y3 = spline_third_derivative(pot, rmin);

    double alpha = -0.5 * y3 / y2;
    double De = y2 / (2.0 * alpha * alpha);
    double we = alpha * sqrt(2.0 * De / data->reduced_mass);
    double wexe = alpha * alpha / (2.0 * data->reduced_mass);
    int vmax = floor(we / (2.0 * wexe) - 0.5);

    printf("   (parameters are estimated at the minimum point only)\n");
    printf("   alpha = %.4f bohr^-1\n", alpha);
    printf("   De    = %.6f a.u. = %.0f cm^-1\n", De, De * ATOMIC_TO_CM);
    printf("   we    = %.2f cm^-1\n", we * ATOMIC_TO_CM);
    printf("   wexe  = %.2f cm^-1\n", wexe * ATOMIC_TO_CM);
    printf("   vmax  = %d\n", vmax);

    printf("\n");
}

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

#include "input_data.h"

#include <stdlib.h>
#include <stdio.h>

#include "units.h"

void print_input_data(input_data_t *data)
{
    printf(" > input data:\n");
    printf("\n");

    printf("   atomic mass 1       : %.8f (amu)\n", data->mass1);
    printf("   atomic mass 2       : %.8f (amu)\n", data->mass2);
    printf("   net charge          : %+d\n", data->charge);
    printf("   reduced mass        : %.8f (amu) = %.8f (electron masses)\n",
           data->reduced_mass * ELECTRON_MASS_TO_AMU, data->reduced_mass);
    printf("   vibrational levels  : from %d to %d\n", data->v_min, data->v_max);
    printf("   rotational levels   : from %d to %d\n", data->J_min, data->J_max);
    printf("   elec transition     : %s\n", (data->pot2 != NULL) ? "yes" : "no");
    printf("   write wavefun files : %s\n", data->write_psi ? "yes" : "no");
    printf("   grid size           : %d points\n", data->grid_size);
    printf("   solver (integrator) : %s\n", (data->solver == SOLVER_NUMEROV) ? "numerov" : "2nd order finite difference (FD2)");
    printf("   mapping of r grid   : ");
    if (data->mapping->type == MAPPING_IDENTITY) {
        printf("identity = off\n");
    }
    else if (data->mapping->type == MAPPING_MESHKOV_08) {
        double re = data->mapping->params[0];
        double beta = data->mapping->params[1];
        printf("meshkov08 (re=%g angstroms, beta=%g)\n", re, beta);
    }
    else {
        printf("other\n");
    }
    printf("\n");

    // potential 1
    if (data->pot1 != NULL) {
        printf("           r [a.u.]          U1(r) [a.u.]  \n");
        for (int i = 0; i < data->n_points; i++) {
            printf("   %20.12f%20.12f\n", data->r[i], data->pot1->y[i]);
        }
        printf("\n");
    }

    // potential 2
    if (data->pot2 != NULL) {
        printf("           r [a.u.]          U2(r) [a.u.]  \n");
        for (int i = 0; i < data->n_points; i++) {
            printf("   %20.12f%20.12f\n", data->r[i], data->pot2->y[i]);
        }
        printf("\n");
    }

    // property matrix element
    if (data->prop != NULL) {
        printf("           r [a.u.]           property(r)  \n");
        for (int i = 0; i < data->n_points; i++) {
            printf("   %20.12f%20.12f\n", data->r[i], data->prop->y[i]);
        }
        printf("\n");
    }
}


/**
 * deallocates memory used by the 'data' object
 */
void delete_input_data(input_data_t *data)
{
    free(data->r);
    if (data->pot1 != NULL) {
        delete_spline(data->pot1);
    }
    if (data->pot2 != NULL) {
        delete_spline(data->pot2);
    }
    if (data->prop != NULL) {
        delete_spline(data->prop);
    }
    free(data);
}

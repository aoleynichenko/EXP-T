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

#include "write_psi.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "input_data.h"
#include "units.h"


#define MAX(a, b) (((a)>(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))

void write_wavefunction_file(int v, int J, int n_points, double *radial_grid, mapping_t *mapping, double *psi);


/**
 * flushes wavefunctions into the formatted files.
 */
void write_wavefunctions(input_data_t *input_data, int *nroots, double **energies, double ***wavefunctions)
{
    int v_min = input_data->v_min;
    int v_max = input_data->v_max;
    int J_min = input_data->J_min;
    int J_max = input_data->J_max;

    cubic_spline_t *pot = input_data->pot1;
    double r_min = pot->x[0];
    double r_max = pot->x[pot->n - 1];
    int n_rot = J_max - J_min + 1;
    int n_vib = v_max - v_min + 1;
    int ngrid = input_data->grid_size;
    double emin = input_data->min_energy;
    double *radial_grid = construct_radial_grid(ngrid, r_min, r_max, input_data->mapping);

    for (int J = J_min; J <= J_max; J++) {
        int index_J = J - J_min;

        for (int v = v_min; v <= MIN(v_max, nroots[index_J]-1); v++) {
            int index_v = v - v_min;

            double *psi = wavefunctions[index_J][index_v];
            write_wavefunction_file(v, J, ngrid, radial_grid, input_data->mapping, psi);
        }
    }

    free(radial_grid);
}


void write_wavefunction_file(int v, int J, int n_points, double *radial_grid, mapping_t *mapping, double *phi)
{
    // create the "wavefunctions" directory if needed
    struct stat st = {0};
    if (stat("wavefunctions", &st) == -1) {
        mkdir("wavefunctions", 0700);
    }

    // construct filename
    char file_name[256];
    sprintf(file_name, "wavefunctions/psi_v%d_J%d.dat", v, J);

    // if mapping: restore wavefunction psi(r) = sqrt(g(y)) phi(y)
    double *psi = (double *) calloc(n_points, sizeof(double));
    restore_wavefunction(n_points, radial_grid, mapping, phi, psi);

    // write file
    FILE *f = fopen(file_name, "w");
    fprintf(f, "%24s%24s%24s%24s\n", "y", "phi(y)", "r", "psi(r)");

    for (int i = 0; i < n_points; i++) {
        double ri = radial_grid[i];
        double yi = mapping->r2y(mapping, ri);

        fprintf(f, "%24.16f%24.16f%24.16f%24.16f\n", yi, phi[i], ri, psi[i]);
    }

    free(psi);
    fclose(f);
}

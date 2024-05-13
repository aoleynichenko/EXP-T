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

#include "read_natural_spinors.h"

#include <complex.h>
#include <stdio.h>

#include "errquit.h"

#define MAX_BUF_SIZE 1024


int skip_line(FILE *fp);


/*
 * reads formatted file with natural spinors expanded in the basis of molecular spinors.
 * such formatted files are produced by the EXP-T program.
 */
void read_natural_spinors(
        char *file_name,
        int *n_mol_spinors, double **spinor_energies,
        int *n_nat_spinors, double **occ_numbers, double complex **nto_coeffs
        )
{
    char buf[MAX_BUF_SIZE];
    double eps;
    double occ_number;
    double coef_re, coef_im;
    int n_scan;

    printf(" > reading file with natural spinors\n");

    FILE *nto_file = fopen(file_name, "r");
    if (nto_file == NULL) {
        errquit("cannot open file with (quasi)natural orbitals");
    }

    // read number of active molecular spinors
    n_scan = fscanf(nto_file, "%s%d", buf, n_mol_spinors);
    if (n_scan != 2) {
        errquit("cannot read number of active molecular spinors");
    }
    skip_line(nto_file);
    skip_line(nto_file);

    // read energies of active molecular spinors
    *spinor_energies = (double *) calloc(*n_mol_spinors, sizeof(double));
    for (int ispinor = 0; ispinor < *n_mol_spinors; ispinor++) {
        n_scan = fscanf(nto_file, "%s%lf%s", buf, &eps, buf);
        if (n_scan != 3) {
            errquit("cannot read energies of active molecular spinors");
        }
        (*spinor_energies)[ispinor] = eps;
    }

    // read number of natural spinors
    n_scan = fscanf(nto_file, "%s%d", buf, n_nat_spinors);
    if (n_scan != 2) {
        errquit("cannot read number of natural spinors");
    }
    skip_line(nto_file);

    // read natural spinors
    *occ_numbers = (double *) calloc(*n_nat_spinors, sizeof(double));
    *nto_coeffs = (double complex *) calloc((*n_mol_spinors) * (*n_nat_spinors), sizeof(double complex));

    for (int i_nts = 0; i_nts < *n_nat_spinors; i_nts++) {

        // occupation number
        n_scan = fscanf(nto_file, "%s%lf", buf, &occ_number);
        if (n_scan != 2) {
            errquit("cannot read occ number of natural spinor %d", i_nts + 1);
        }
        (*occ_numbers)[i_nts] = occ_number;

        // coefficients of expansion in the basis of active molecular spinors
        for (int ispinor = 0; ispinor < *n_mol_spinors; ispinor++) {
            n_scan = fscanf(nto_file, "%s%lf%lf", buf, &coef_re, &coef_im);
            if (n_scan != 3) {
                errquit("cannot read coefficient of natural spinor %d", i_nts + 1);
            }
            (*nto_coeffs)[(*n_mol_spinors) * i_nts + ispinor] = coef_re + coef_im * I;
        }
    }

    fclose(nto_file);
}


int skip_line(FILE *fp)
{
    int c;

    while (c = fgetc(fp), c != '\n' && c != EOF);

    return c;
}

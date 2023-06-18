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

/*
 * Operations with formatted files containing effective Hamiltonians.
 * The format of formatted HEFFF files was primarily introduced by A. Zaitsevskii.
 */

#include "formatted_heff.h"

#include <complex.h>
#include <stdio.h>
#include <string.h>

#include "options.h"
#include "symmetry.h"


/**
 * Writes blocks of the effective Hamiltonian to the formatted file 'HEFF'
 *
 * @param sect_h
 * @param sect_p
 * @param block_dims
 * @param heff
 */
void write_formatted_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff)
{
    int first_irrep = first_nonzero_irrep(block_dims);

    // open file with formatted Heff
    FILE *hefff = open_formatted_heff_file(sect_h, sect_p, NULL);

    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_dims[irep] == 0) {
            continue;
        }
        double complex *heff_block = heff[irep];
        size_t ms_size = block_dims[irep];
        hefff_write_block(hefff, arith == CC_ARITH_COMPLEX, irep - first_irrep + 1, ms_size, heff_block);
    }

    close_formatted_heff_file(hefff);
}


/**
 * Writes total energy to the file with formatted Heff.
 * file will be truncated
 */
void write_formatted_heff_0h0p(double total_energy)
{
    FILE *hefff = fopen("HEFF", "w");
    if (arith == CC_ARITH_COMPLEX) {
        fprintf(hefff, "complex      # arithmetic\n");
    }
    else {
        fprintf(hefff, "real         # arithmetic\n");
    }
    fprintf(hefff, "0h0p         # sector\n");
    fprintf(hefff, "   1     1   # rep No & heff size\n");
    if (arith == CC_ARITH_COMPLEX) {
        fprintf(hefff, "%21.12E%21.12E\n", total_energy, 0.0);
    }
    else {
        fprintf(hefff, "%21.12E\n", total_energy);
    }
    fclose(hefff);
}


/**
 * Opens formatted HEFF file.
 * Writes header lines.
 * Returns pointer to the FILE structure.
 */
FILE *open_formatted_heff_file(int sect_h, int sect_p, char *label)
{
    FILE *hefff;

    hefff = fopen("HEFF", "a");
    if (hefff == NULL) {
        return NULL;
    }

    if (label == NULL) {
        fprintf(hefff, "%dh%dp         # sector\n", sect_h, sect_p);
    }
    else {
        fprintf(hefff, "%s         # sector\n", label);
    }

    return hefff;
}


/**
 * Closes formatted HEFF file.
 */
void close_formatted_heff_file(FILE *hefff)
{
    fclose(hefff);
}


/**
 * Writes a block of Heff to the formatted file.
 */
void hefff_write_block(FILE *hefff, int carith, int rep_no, size_t dim, double complex *heff)
{
    fprintf(hefff, "%4d%6d   # rep No & heff size\n", rep_no, dim);
    for (int i = 0; i < dim * dim; i++) {
        if (carith) {
            fprintf(hefff, "%21.12E%21.12E", creal(heff[i]), cimag(heff[i]));
            if (i > 0 && i % 2 != 0) { fprintf(hefff, "\n"); }
        }
        else {
            fprintf(hefff, "%21.12E", creal(heff[i]));
            if (i > 0 && (i + 1) % 4 == 0) { fprintf(hefff, "\n"); }
        }
    }
    if (dim * dim % 2 != 0) { fprintf(hefff, "\n"); }
}



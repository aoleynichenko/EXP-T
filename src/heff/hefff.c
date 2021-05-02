/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
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

/**
 * Operations with formatted files containing effective Hamiltonians.
 * The format of formatted HEFFF files was primarily introduced by
 * A. Zaitsevskii.
 * 2019-2021 Alexander Oleynichenko
 */

#include <complex.h>
#include <stdio.h>
#include <string.h>


/**
 * Opens formatted HEFF file.
 * Writes header lines.
 * Returns pointer to the FILE structure.
 */
FILE *hefff_open(int sect_h, int sect_p, char *label)
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
void hefff_close(FILE *hefff)
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

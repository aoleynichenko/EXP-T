/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2021 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

/*******************************************************************************
 * tm2c_interface.c
 * ================
 *
 * Implementation of the interface to the tm2c integral code by R. Berger,
 * K. Gaul, S. Giesen and other our colleagues from Marburg.
 * Reads info about spinors and transformed MO integrals from the formatted
 * files.
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "interfaces.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "engine.h"
#include "platform.h"
#include "error.h"
#include "io.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"

#define MAX_TITLE 512


// reads info about spinors and transformed MO integrals from formatted file
// NOTE: no error checking at all!
void tm2c_interface(char *moints_file_1, char *moints_file_2, char *moints_file_prop, cc_options_t *opts)
{
    char title[MAX_TITLE];
    int occ_numbers[CC_MAX_SPINORS];
    double orb_energies[CC_MAX_SPINORS];
    FILE *inp_file;
    int tile_size;

    timer_new_entry("dirac", "tm2c interface");
    timer_start("dirac");

    inp_file = fopen("../expt_moldata", "r");
    if (inp_file == NULL) {
        errquit("tm2c interface: file 'expt_moldata' not found");
    }

    fgets(title, MAX_TITLE, inp_file);
    title[strlen(title) - 1] = '\0'; // trim '\n'

    fscanf(inp_file, "%d", &NSPINORS);
    fscanf(inp_file, "%lf", &opts->enuc);

    for (int i = 0; i < NSPINORS; i++) {
        fscanf(inp_file, "%d", &occ_numbers[i]);
    }
    for (int i = 0; i < NSPINORS; i++) {
        fscanf(inp_file, "%lf", &orb_energies[i]);
    }

    fclose(inp_file);

    //======================== core hamiltonian matrix =========================

    FILE *hcore_file = fopen("../expt_hcore", "r");
    double complex *hcore;
    double h_re, h_im;
    int i1, i2;
    if (hcore_file == NULL) {
        errquit("tm2c interface: file 'expt_hcore' not found");
    }
    size_t hcore_nbytes = sizeof(double complex) * NSPINORS * NSPINORS;
    hcore = (double complex *) cc_malloc(hcore_nbytes);
    memset(hcore, 0, hcore_nbytes);

    while (fscanf(hcore_file, "%d%d%lf%lf", &i1, &i2, &h_re, &h_im) == 4) {
        hcore[i1 * NSPINORS + i2] = h_re + h_im * I;
        hcore[i2 * NSPINORS + i1] = h_re - h_im * I;
    }

    int fd = io_open("HINT", "w");
    io_write_compressed(fd, hcore, hcore_nbytes);
    io_close(fd);

    cc_free(hcore);
    fclose(hcore_file);

    //==========================================================================

    printf("Comment: %s\n", title);
    printf("NSPINORS = %d\n", NSPINORS);

    opts->escf = 0.0;
    printf("Nuclear repulsion energy = %20.12f\n", opts->enuc);
    printf("Total SCF energy = %20.12f\n", opts->escf);

    // only C1 group with complex arithmetics
    arith = CC_ARITH_COMPLEX;
    carith = 1;
    printf("Arithmetic: complex\n");
    printf("Spin-orbit: always ON\n");

    // Number of double group irreps
    nsym = 1;

    // Representation names
    rep_names = (char **) cc_malloc(nsym * sizeof(char *));
    rep_names[0] = (char *) cc_malloc(sizeof(char) * 8);
    strcpy(rep_names[0], "A");
    irrep_a1 = 0;
    printf("Point symmetry group: C1\n");
    printf("Totally symmetric irrep: %s\n", rep_names[irrep_a1]);

    // "Multiplication table" ( = dir prod table for Abelian groups)
    // in fact, 3d array

    // allocate table
    dir_prod_table = (int ***) cc_malloc(nsym * sizeof(int **));
    for (int i = 0; i < nsym; i++) {
        dir_prod_table[i] = (int **) cc_malloc(nsym * sizeof(int *));
        for (int j = 0; j < nsym; j++) {
            dir_prod_table[i][j] = (int *) cc_malloc((nsym + 1) * sizeof(int));
            for (int k = 0; k < nsym + 1; k++) {
                dir_prod_table[i][j][k] = -1;
            }
            // get values from the dirac_data structure (from MRCONEE files)
        }
    }
    is_abelian = 1;
    dir_prod_table[0][0][0] = 1;
    dir_prod_table[0][0][1] = 0; // A x A = A

    tile_size = opts->tile_size;

    // fill spinor info array
    //  *** Spinor, irrep, occupation, energy ***
    spinor_info = (spinor_attr_t *) cc_malloc(sizeof(spinor_attr_t) * NSPINORS);
    for (int i = 0; i < NSPINORS; i++) {
        spinor_info[i].repno = 0;
        spinor_info[i].eps = orb_energies[i];
        spinor_info[i].space_flags = 0;
        if (occ_numbers[i] == 1) {
            spinor_info[i].space_flags |= CC_FLAG_OCCUPIED;
        }
    }

    create_spinor_blocks(tile_size);
    setup_occupation_numbers(opts, NSPINORS, spinor_info);
    setup_active_space(opts);
    setup_fast_access_spinor_lists();

    print_symmetry_info();
    print_spinor_info_table();

    //========================= two-electron integrals =========================

    // todo: tiles -> VINT-*-*-*-* files

    printf("\nstart 2-el integrals...\n");

    FILE *twoel_file = fopen("../expt_twoel", "r");
    const int max_buf = 4096;
    int32_t nint = 0;
    int i3, i4;
    double v_re, v_im;

    static int spinor_to_block[CC_MAX_SPINORS];
    for (int isb = 0; isb < n_spinor_blocks; isb++) {
        spinor_block_t *spb = &spinor_blocks[isb];
        for (int i = 0; i < spb->size; i++) {
            spinor_to_block[spb->indices[i]] = isb;
        }
    }

    typedef struct {
        double complex vint[max_buf];
        int16_t idx4[4 * max_buf];
        int32_t nint;
    } output_buffer_t;
    output_buffer_t *output_buffers;

    int n_out_buf = n_spinor_blocks * n_spinor_blocks * n_spinor_blocks * n_spinor_blocks;
    output_buffers = (output_buffer_t *) cc_malloc(sizeof(output_buffer_t) * n_out_buf);
    for (int i = 0; i < n_out_buf; i++) {
        output_buffers[i].nint = 0;
    }

    while (fscanf(twoel_file, "%d%d%d%d%lf%lf", &i1, &i2, &i3, &i4, &v_re, &v_im) == 6) {

        int isb1 = spinor_to_block[i1];
        int isb2 = spinor_to_block[i2];
        int isb3 = spinor_to_block[i3];
        int isb4 = spinor_to_block[i4];
        int nsb = n_spinor_blocks;

        output_buffer_t *out_buf = &output_buffers[isb1 * nsb * nsb * nsb + isb2 * nsb * nsb + isb3 * nsb + isb4];

        nint = out_buf->nint;
        out_buf->vint[nint] = v_re + v_im * I;
        out_buf->idx4[4 * nint] = i1 + 1;
        out_buf->idx4[4 * nint + 1] = i2 + 1;
        out_buf->idx4[4 * nint + 2] = i3 + 1;
        out_buf->idx4[4 * nint + 3] = i4 + 1;
        out_buf->nint++;

        if (out_buf->nint == max_buf) {
            char vint_file[256];
            sprintf(vint_file, "VINT-%d-%d-%d-%d", isb1 + 1, isb2 + 1, isb3 + 1, isb4 + 1);
            fd = io_open(vint_file, "a");
            io_write_compressed(fd, &out_buf->nint, sizeof(int32_t));
            io_write_compressed(fd, out_buf->idx4, 4 * sizeof(int16_t) * out_buf->nint);
            io_write_compressed(fd, out_buf->vint, sizeof(double complex) * out_buf->nint);
            io_close(fd);
            out_buf->nint = 0;
        }
    }

    // flush remaining integrals
    for (int isb1 = 0; isb1 < n_spinor_blocks; isb1++) {
        for (int isb2 = 0; isb2 < n_spinor_blocks; isb2++) {
            for (int isb3 = 0; isb3 < n_spinor_blocks; isb3++) {
                for (int isb4 = 0; isb4 < n_spinor_blocks; isb4++) {
                    int nsb = n_spinor_blocks;
                    output_buffer_t *out_buf = &output_buffers[isb1 * nsb * nsb * nsb + isb2 * nsb * nsb + isb3 * nsb +
                                                               isb4];

                    char vint_file[256];
                    sprintf(vint_file, "VINT-%d-%d-%d-%d", isb1 + 1, isb2 + 1, isb3 + 1, isb4 + 1);
                    fd = io_open(vint_file, "a");

                    if (out_buf->nint > 0) {
                        io_write_compressed(fd, &out_buf->nint, sizeof(int32_t));
                        io_write_compressed(fd, out_buf->idx4, 4 * sizeof(int16_t) * out_buf->nint);
                        io_write_compressed(fd, out_buf->vint, sizeof(double complex) * out_buf->nint);
                        out_buf->nint = 0;
                    }
                    io_write_compressed(fd, &out_buf->nint, sizeof(int32_t));
                    io_close(fd);
                }
            }
        }
    }

    fclose(twoel_file);
    cc_free(output_buffers);

    printf("end of 2-el integrals\n");

    timer_stop("dirac");
}

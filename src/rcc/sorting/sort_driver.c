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

/**
 * Sorting of integrals (creating basic diagrams from the raw integral arrays).
 *
 * Algorithm in brief:
 * (1) Split MDCINT into smaller files
 *     (one 4-dim integral block <IJ|KL> = one "VINT-I-J-K-L" file)
 *     (file dirac_binary.f90)
 * (2) Put diagrams to be sorted into the queue
 * (3) Read "VINT-I-J-K-L" files one-by-one; for each file find corresponding
 *     blocks in 2-particle diagrams in the queue; fill these blocks with
 *     integrals from "VINT-I-J-K-L" file.
 * (4) Read files again, but add integrals with interchanged K,L indices with
 *     the opposite sign (antisymmetrzation): <ij||kl> = <ij|kl> - <ij|lk>
 * (5) Construct Fock matrix (using the "hhhh" diagram and the "HINT" file)
 * (6) Construct 1-particle diagrams (requests are again in the queue)
 */

#include "sort.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "sorting_request.h"

#include "datamodel.h"
#include "engine.h"
#include "error.h"
#include "io.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "timer.h"
#include "utils.h"
#include "symmetry.h"


// TODO: separate file with sorting stats ?

void sorting_print_configuration();

void sort_onel();

void sort_twoel();


/**
 * prints sorting configuration:
 * size characteristics of the system & memory requirements
 */
void sorting_print_configuration()
{
    const size_t io_integrals_buf_size = sizeof(double complex) * CC_SORTING_IO_BUF_SIZE;
    const size_t io_indices_buf_size = 4 * sizeof(int16_t) * CC_SORTING_IO_BUF_SIZE;
    const size_t io_total_buf_size = io_integrals_buf_size + io_indices_buf_size;
    const double io_total_buf_size_mb = io_total_buf_size / (1024.0 * 1024.0);
    const double twoelec_buf_size_mb = pow(get_max_spinor_block_size(), 4) / (1024.0 * 1024.0);

    printf(" number of spinors                              %d\n", get_num_spinors());
    printf(" number of spinor blocks                        %d\n", n_spinor_blocks);
    printf(" tile size                                      %d\n", cc_opts->tile_size);
    printf(" max spinor block size                          %d\n", get_max_spinor_block_size());
    printf(" size of i/o buffer for integrals and indices   %.3f MB\n", io_total_buf_size_mb);
    printf(" size of the buffer for two-elec integrals      %.3f MB\n", twoelec_buf_size_mb);
}


/**
 * Leaves requests for sorting of 1- and 2-particle diagrams with given
 * "holes/particles" and "valence" characteristics from the raw integrals.
 * The new (empty) diagram will be added to the diagram stack.
 * Diagram will have the 'standard' order (12) or (1234)
 *
 * Arguments:
 *   name     symbolic name of a new diagram
 *   qparts   holes/particles list (hh, hp, ph or pp)
 *   valence  which indices are valence (1/0) ("00", "10", "01" or "11")
 *   order    which order must be assigned to the diagram
 *            (NO ACTUAL REORDERING IS PERFORMED)
 *            allowed orders are only 12, 21, 1234, 3412
 */
void request_sorting_sym(char *name, char *qparts, char *valence, char *order, int operator_symmetry)
{
    char file_name[CC_MAX_FILE_NAME_LENGTH];
    int rank;

    if (n_requests + 1 == CC_MAX_SORTING_REQUESTS) {
        errquit("number of sorting requests exceeds CC_MAX_SORTING_REQUESTS (=%d)",
                CC_MAX_SORTING_REQUESTS);
    }

    if (strlen(qparts) == 2 && strlen(valence) == 2 && strlen(order) == 2) {
        rank = 2;
    }
        // 2. two-particle diagram
    else if (strlen(qparts) == 4 && strlen(valence) == 4 && strlen(order) == 4) {
        rank = 4;
    }
    else {
        errquit("in request sorting(): cannot determine rank of the diagram "
                "'%s' to be sorted: [%s/%s/%s]. Rank must be equal to 2 or 4\n",
                name, qparts, valence, order);
    }

    // ??? здесь считывать уже готовые ?
    // try to read from disk
    sprintf(file_name, "%s.dg", name);
    if (rank == 2 && cc_opts->reuse_integrals_1) {
        if (diagram_read(file_name) != NULL) {
            printf(" Reuse 1-electron integrals file '%s'\n", name);
            return;
        }
        else {
            errquit(" 1-electron integrals file '%s' not found\n", name);
        }
    }
    if (rank == 4 && cc_opts->reuse_integrals_2) {
        if (diagram_read(file_name) != NULL) {
            printf(" Reuse 2-electron integrals file '%s'\n", name);
            return;
        }
        else {
            errquit(" 2-electron integrals file '%s' not found\n", name);
        }
    }

    // create the diagram template:
    // 1. one-particle diagram
    if (rank == 2) {
        tmplt_sym(name, qparts, valence, "12", NOT_PERM_UNIQUE, operator_symmetry);
    }
        // 2. two-particle diagram
    else if (rank == 4) {
        tmplt_sym(name, qparts, valence, "1234", IS_PERM_UNIQUE, operator_symmetry); // здесь нужно писать 'order'
    }
    else {
        errquit("in request sorting(): only diagrams with rank = 2 or 4 are allowed: %s[%s/%s/%s]\n",
                name, qparts, valence, order);
    }

    append_sorting_request(sorting_requests, &n_requests, name, qparts, valence, order);
}


void request_sorting(char *name, char *qparts, char *valence, char *order)
{
    request_sorting_sym(name, qparts, valence, order, get_totally_symmetric_irrep());
}


/**
 * main function of the sorting module
 */
void perform_sorting()
{
    int ireq;
    sorting_request_t *req;
    int sect_h, sect_p;

    sect_h = cc_opts->curr_sector_h;
    sect_p = cc_opts->curr_sector_p;

    timer_new_entry("sort", "Sorting of integrals");
    timer_start("sort");

    printf("\n");
    printf(" Integral sorting for the %dh%dp sector\n", sect_h, sect_p);
    printf(" started at");
    print_asctime();

    // try to use already sorted diagrams (read them from disk)
    // если получится без проблем ВСЕ требуемые диаграммы, значит,
    // можно спокойно проскочить сортировки и выполнять лишь расчеты ESCF и т.д.
    // для успешно считанных не надо делать reorder (точнее, нельзя!)
    // повторно писать их на диск тоже не нужно

    sorting_print_configuration();

    sort_twoel();
    sort_onel();

    // только для реально отсортированных в этом запуске запросов!
    // reorder some diagrams (if required)
    for (ireq = 0; ireq < n_requests; ireq++) {
        req = &sorting_requests[ireq];
        if (strcmp(req->order, "12") == 0 || strcmp(req->order, "1234") == 0) {
            // nothing to reorder, already in the required ("normal") order
            continue;
        }
        // for PPPP diagram only: complex conjugate instead of reordering (already done, SKIP)
        if (strcmp(req->hp, "pppp") == 0 &&
            strcmp(req->valence, "0000") == 0 &&
            strcmp(req->order, "3412") == 0) {
            // set order 3412 instead of 1234
            req->dg->order[0] = 3;
            req->dg->order[1] = 4;
            req->dg->order[2] = 1;
            req->dg->order[3] = 2;
            continue;
        }

        reorder(req->dg_name, req->dg_name, req->order);
    }

    // write diagrams to disk
    // только реально обработанные запросы!
    for (ireq = 0; ireq < n_requests; ireq++) {
        char dg_file_name[CC_MAX_PATH_LENGTH];
        req = &sorting_requests[ireq];
        sprintf(dg_file_name, "%s.dg", req->dg_name);
        diagram_write(diagram_stack_find(req->dg_name), dg_file_name);
        //printf("%s ", req->dg_name);
    }

    /*
     * finalize
     */
    n_requests = 0;

    timer_stop("sort");
    //printf("   number of blocks read from disk: %d\n", n_blocks_read);
    //printf("   total number of integrals read from disk: %ld\n", n_integrals_read);
    //printf("   total number of bytes read from disk: %ld (%.2f GB)\n", n_bytes_read,
    //       (double) n_bytes_read / (1024.0 * 1024.0 * 1024.0));
    //printf("   sorting performance, Mb/sec: %.2f\n", (double) n_bytes_read / (1024.0 * 1024.0) / timer_get("sort"));
    printf(" time for 2-e integrals sorting, sec: %.2f\n", timer_get("sort"));
    printf(" time for DIRAC interface (integral extraction & write), sec: %.2f\n", timer_get("dirac"));
    printf(" total time for sorting operations, sec: %.2f\n", timer_get("dirac") + timer_get("sort"));

    printf(" finished at at");
    print_asctime();
}


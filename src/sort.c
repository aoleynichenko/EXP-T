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

/*******************************************************************************
 * sort.c
 * ======
 *
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
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "datamodel.h"
#include "engine.h"
#include "error.h"
#include "io.h"
#include "memory.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"
#include "utils.h"

// arrays with "raw" integrals -- sparse and with zeros (but fast)
double complex *h_ints;
double complex *f_ints;

void init_sorting();

void sort_onel();

void finalize_sorting();

int read_prop_two_files(int nspinors, char *file_re, char *file_im, double complex *prop_mat);
int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat);

#define CC_MAX_SORTING_REQUESTS 64
#define CC_SORTING_IO_BUF_SIZE 16384

enum {
    INTS_COULOMB = 1,
    INTS_GAUNT = 2,
    INTS_BREIT = 3
};

typedef struct {
    diagram_t *dg;
    char dg_name[256];
    char hp[CC_DIAGRAM_MAX_RANK];
    char valence[CC_DIAGRAM_MAX_RANK];
    char order[CC_DIAGRAM_MAX_RANK];
    int done;
} sorting_request_t;
static sorting_request_t sorting_requests[CC_MAX_SORTING_REQUESTS];
static int n_requests = 0;

static int max_spinor_block_size;
static int32_t *symblock_indices = NULL;
static double complex ****vint = NULL;


/*******************************************************************************
 * init_sorting
 *
 * Allocates memory and initializes data structures used by sorting routines.
 ******************************************************************************/
void init_sorting()
{
    int tile_size;
    const size_t integrals_buf_size = sizeof(double complex) * CC_SORTING_IO_BUF_SIZE;
    const size_t indices_buf_size = 4 * sizeof(int16_t) * CC_SORTING_IO_BUF_SIZE;
    size_t symblock_indices_size;
    size_t vint_size;
    size_t sz;
    int i, j, k, l;

    tile_size = cc_opts->tile_size;

    max_spinor_block_size = 0;
    for (i = 0; i < n_spinor_blocks; i++) {
        if (spinor_blocks[i].size > max_spinor_block_size) {
            max_spinor_block_size = spinor_blocks[i].size;
        }
    }
    sz = max_spinor_block_size;

    // allocate memory for the temporary arrays
    // alloc arrays for integrals
    h_ints = (double complex *) cc_malloc(nspinors * nspinors * sizeof(double complex));
    f_ints = (double complex *) cc_malloc(nspinors * nspinors * sizeof(double complex));

    symblock_indices_size = sizeof(int32_t) * 4 * sz * sz * sz * sz;
    symblock_indices = (int32_t *) cc_malloc(symblock_indices_size);
    vint = (double complex ****) cc_malloc(sizeof(double complex ***) * sz);
    for (i = 0; i < sz; i++) {
        vint[i] = (double complex ***) cc_malloc(sizeof(double complex **) * sz);
        for (j = 0; j < sz; j++) {
            vint[i][j] = (double complex **) cc_malloc(sizeof(double complex *) * sz);
            for (k = 0; k < sz; k++) {
                vint[i][j][k] = (double complex *) cc_malloc(sizeof(double complex) * sz);
            }
        }
    }
    vint_size = sizeof(double complex ****) + sz * sizeof(double complex ***) +
                sz * sz * sizeof(double complex **) + sz * sz * sz * sizeof(double complex *) +
                sz * sz * sz * sz * sizeof(double complex);

    if (cc_opts->print_level >= CC_PRINT_HIGH) {
        printf("   diagrams to be sorted: ");
        for (i = 0; i < n_requests; i++) {
            printf("[%s/%s] ", sorting_requests[i].hp, sorting_requests[i].valence);
        }
        printf("\n");
    }
    printf("   number of spinor blocks: %d\n", n_spinor_blocks);
    printf("   tile size: %d\n", tile_size);
    printf("   max spinor block length: %d\n", max_spinor_block_size);
    printf("   i/o buffers size, bytes: %d (indices) + %d (integrals) [%.2f KB total]\n",
           integrals_buf_size, indices_buf_size,
           (integrals_buf_size + indices_buf_size) / 1024.0);
    printf("   estimated max memory required for the straightforward algorithm (not used!): %.2f GB\n",
           sizeof(double complex) * nspinors * nspinors * nspinors * nspinors / (1024.0 * 1024.0 * 1024.0));
    printf("   size of working arrays, bytes: %d (indices) + %d (integrals) [%.2f MB total]\n",
           symblock_indices_size, vint_size,
           (symblock_indices_size + vint_size) / (1024.0 * 1024.0));
}


/*******************************************************************************
 * request_sorting
 *
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
 ******************************************************************************/
void request_sorting(char *name, char *qparts, char *valence, char *order)
{
    diagram_t *dg;
    sorting_request_t *req;
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
        tmplt(name, qparts, valence, "12", NOT_PERM_UNIQUE);
    }
        // 2. two-particle diagram
    else if (rank == 4) {
        tmplt(name, qparts, valence, "1234", IS_PERM_UNIQUE); // здесь нужно писать 'order'
    }
    else {
        errquit("in request sorting(): only diagrams with rank = 2 or 4 are allowed: %s[%s/%s/%s]\n",
                name, qparts, valence, order);
    }
    dg = diagram_stack_find(name);

    req = &sorting_requests[n_requests];
    req->dg = dg;
    strcpy(req->dg_name, name);
    strcpy(req->hp, qparts);
    strcpy(req->valence, valence);
    strcpy(req->order, order);
    req->done = 0;
    n_requests++;
}


/**
 * Returns number of integrals read
 * vint_array = factor1 * vint_array + factor2 * new_integrals
 */
int read_twoel_block_unformatted(char *vint_file_name, double complex ****vint_array, double complex factor1,
                                 double complex factor2, size_t *n_bytes_read)
{
    int32_t nint;
    int32_t iint;
    int i, j, k, l;
    size_t n_integrals_read = 0;
    size_t nbr = 0;
    int fd;
    int nloop = 0;
    int empty_file = 0;
    double dfactor1 = creal(factor1);
    double dfactor2 = creal(factor2);

    static double complex buf_integrals[CC_SORTING_IO_BUF_SIZE];
    static int16_t buf_indices[4 * CC_SORTING_IO_BUF_SIZE];

    fd = io_open(vint_file_name, "r");

    while (1) {
        int s = io_read_compressed(fd, &nint, sizeof(int32_t));
        nloop++;

        if (nint == 0) {
            if (nloop == 1) { empty_file = 1; }
            break;
        }

        s = io_read_compressed(fd, buf_indices, nint * 4 * sizeof(int16_t));
        if (carith) {
            s = io_read_compressed(fd, buf_integrals, nint * sizeof(double complex));
        }
        else {
            s = io_read_compressed(fd, buf_integrals, nint * sizeof(double));
        }
        for (iint = 0; iint < nint; iint++) {
            // absolute spinor indices
            i = buf_indices[iint * 4] - 1; // -1 since numeration from 0 in C and from 1 in F90
            j = buf_indices[iint * 4 + 1] - 1;
            k = buf_indices[iint * 4 + 2] - 1;
            l = buf_indices[iint * 4 + 3] - 1;
            i = spinor_index_global2local[i];
            j = spinor_index_global2local[j];
            k = spinor_index_global2local[k];
            l = spinor_index_global2local[l];

            // "local" spinor indices (inside the spinor block under consideration)
            if (carith) {
                vint_array[i][j][k][l] = factor1 * vint_array[i][j][k][l] + factor2 * buf_integrals[iint];
            }
            else {
                double d = ((double *) buf_integrals)[iint];
                vint_array[i][j][k][l] = dfactor1 * vint_array[i][j][k][l] + dfactor2 * d + 0.0 * I;
            }
        }
        n_integrals_read += nint;
        if (arith == CC_ARITH_REAL) {
            nbr += nint * sizeof(double) + nint * 4 * sizeof(int16_t) + sizeof(int32_t);
        }
        else {  // complex arith
            nbr += nint * sizeof(double complex) + nint * 4 * sizeof(int16_t) + sizeof(int32_t);
        }
    }
    io_close(fd);
    *n_bytes_read = nbr;

    if (empty_file) {
        return 0;
    }

    return n_integrals_read;
}


/**
 * main function of the sorting module
 */
void perform_sorting()
{
    int spinor_block_1, spinor_block_2, spinor_block_3, spinor_block_4;
    char vint_file_name[256];
    int fd;
    int sign;
    int antisym_inplace = 0;
    int ireq, isb;
    diagram_t *dg;
    block_t *sb;
    sorting_request_t *req;
    int *p_spinor_blocks;
    int i1, i2, i3, i4;
    int n_blocks_read = 0;
    size_t n_integrals_read = 0;
    size_t n_bytes_read = 0;
    int sect_h, sect_p;

    sect_h = cc_opts->curr_sector_h;
    sect_p = cc_opts->curr_sector_p;

    timer_new_entry("sort", "Sorting of integrals");
    timer_start("sort");

    // count numbers of 1-el and 2-el diagrams to be sorted
    int n_requests_1el = 0;
    int n_requests_2el = 0;
    for (int ireq = 0; ireq < n_requests; ireq++) {
        req = &sorting_requests[ireq];
        if (req->dg->rank == 2) {
            n_requests_1el++;
        }
        else {
            n_requests_2el++;
        }
    }

    printf("\n");
    printf("  Integral sorting routine: integrals for the %dh%dp sector\n", sect_h, sect_p);
    printf(" ------------------------------------------------------------------------------------------\n");
    printf("  ");
    print_asctime();

    // try to use already sorted diagrams (read them from disk)
    // если получится без проблем ВСЕ требуемые диаграммы, значит,
    // можно спокойно проскочить сортировки и выполнять лишь расчеты ESCF и т.д.
    // для успешно считанных не надо делать reorder (точнее, нельзя!)
    // повторно писать их на диск тоже не нужно

    init_sorting();

    if (n_requests_2el > 0) {
        printf("   sorting two-electron integrals ...\n");

        for (sign = -1; sign <= 1; sign += 2) {
            for (spinor_block_1 = 0; spinor_block_1 < n_spinor_blocks; spinor_block_1++) {
                for (spinor_block_2 = 0; spinor_block_2 < n_spinor_blocks; spinor_block_2++) {
                    for (spinor_block_3 = 0; spinor_block_3 < n_spinor_blocks; spinor_block_3++) {
                        for (spinor_block_4 = 0; spinor_block_4 < n_spinor_blocks; spinor_block_4++) {

                            if (spinor_block_3 == spinor_block_4) {
                                antisym_inplace = 1;
                            }
                            else {
                                antisym_inplace = 0;
                            }

                            if (antisym_inplace && sign == -1) {
                                // skip reading the second block which is equal to the first
                                // (it is redundant)
                                continue;
                            }

                            // sizes of spinor blocks
                            int spb1_sz = spinor_blocks[spinor_block_1].size;
                            int spb2_sz = spinor_blocks[spinor_block_2].size;
                            int spb3_sz = spinor_blocks[spinor_block_3].size;
                            int spb4_sz = spinor_blocks[spinor_block_4].size;
                            // init tmp array with integrals
                            for (int i = 0; i < spb1_sz; i++) {
                                for (int j = 0; j < spb2_sz; j++) {
                                    for (int k = 0; k < (sign == 1 ? spb3_sz : spb4_sz); k++) {
                                        for (int l = 0; l < (sign == 1 ? spb4_sz : spb3_sz); l++) {
                                            vint[i][j][k][l] = 0.0 + 0.0 * I;
                                        }
                                    }
                                }
                            }

                            if (sign == 1) {
                                sprintf(vint_file_name, "VINT-%d-%d-%d-%d",
                                        spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_3 + 1, spinor_block_4 + 1);
                            }
                            else {
                                sprintf(vint_file_name, "VINT-%d-%d-%d-%d",
                                        spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_4 + 1, spinor_block_3 + 1);
                            }
                            if (io_file_exists(vint_file_name) == 0) {
                                continue;
                            }
                            n_blocks_read++;

                            // load Coulomb integrals
                            size_t nbr = 0;
                            n_integrals_read = read_twoel_block_unformatted(vint_file_name, vint, 0.0 + 0.0 * I,
                                                                            1.0 + 0.0 * I, &nbr);
                            if (n_integrals_read == 0) {
                                continue;
                            }
                            n_bytes_read += nbr;

                            // load Gaunt integrals (if needed)
                            if (cc_opts->gaunt_defined) {
                                // read GINT* files
                                if (sign == 1) {
                                    sprintf(vint_file_name, "GINT-%d-%d-%d-%d",
                                            spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_3 + 1,
                                            spinor_block_4 + 1);
                                }
                                else {
                                    sprintf(vint_file_name, "GINT-%d-%d-%d-%d",
                                            spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_4 + 1,
                                            spinor_block_3 + 1);
                                }
                                n_integrals_read = read_twoel_block_unformatted(vint_file_name, vint, 1.0 + 0.0 * I,
                                                                                1.0 + 0.0 * I, &nbr);
                                if (n_integrals_read == 0) {
                                    continue;
                                }
                                n_bytes_read += nbr;
                            }

                            // load two-electron property integrals (if needed)
                            for (int iprop = 0; iprop < cc_opts->n_twoprop; iprop++) {
                                // read TWOPROP* files
                                if (sign == 1) {
                                    sprintf(vint_file_name, "TWOPROP%d-%d-%d-%d-%d", iprop + 1,
                                            spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_3 + 1,
                                            spinor_block_4 + 1);
                                }
                                else {
                                    sprintf(vint_file_name, "TWOPROP%d-%d-%d-%d-%d", iprop + 1,
                                            spinor_block_1 + 1, spinor_block_2 + 1, spinor_block_4 + 1,
                                            spinor_block_3 + 1);
                                }
                                double complex lambda = cc_opts->twoprop_lambda[iprop];
                                n_integrals_read = read_twoel_block_unformatted(vint_file_name, vint, 1.0 + 0.0 * I,
                                                                                lambda, &nbr);
                                /*if (n_integrals_read == 0) {
                                    continue;
                                }*/
                                n_bytes_read += nbr;
                            }
                            if (cc_opts->n_twoprop > 0 && n_integrals_read == 0) {
                                continue;
                            }

                            // process requests -- fill all diagrams
                            for (ireq = 0; ireq < n_requests; ireq++) {
                                req = &sorting_requests[ireq];
                                dg = req->dg;
                                if (dg->rank != 4) { // one-electron diagrams will be sorted later
                                    continue;
                                }
                                for (isb = 0; isb < dg->n_blocks; isb++) {
                                    sb = dg->blocks[isb];
                                    if (sb->is_unique == 0) {  // skip permutationally unique block
                                        continue;
                                    }
                                    p_spinor_blocks = sb->spinor_blocks;
                                    if (p_spinor_blocks[0] != spinor_block_1 ||
                                        p_spinor_blocks[1] != spinor_block_2 ||
                                        p_spinor_blocks[2] != spinor_block_3 ||
                                        p_spinor_blocks[3] != spinor_block_4) {
                                        continue; // go to the next block
                                    }
                                    // block was found. fill it!
                                    symblock_load(sb);
                                    double *dbuf = (double *) sb->buf;
                                    //if (sign==-1) memset(sb->buf, 0, sizeof(double complex)*sb->size);
                                    //if (sign==1 && antisym_inplace) memset(sb->buf, 0, sizeof(double complex)*sb->size);
                                    symblock_gen_indices(sb, symblock_indices);
                                    if (antisym_inplace) {
                                        for (size_t i = 0; i < sb->size; i++) {
                                            i1 = symblock_indices[sb->rank * i];  // sb->rank == 4
                                            i2 = symblock_indices[sb->rank * i + 1];
                                            i3 = symblock_indices[sb->rank * i + 2];
                                            i4 = symblock_indices[sb->rank * i + 3];
                                            // global -> local
                                            i1 = spinor_index_global2local[i1];
                                            i2 = spinor_index_global2local[i2];
                                            i3 = spinor_index_global2local[i3];
                                            i4 = spinor_index_global2local[i4];
                                            if (carith) {
                                                sb->buf[i] = vint[i1][i2][i3][i4] - vint[i1][i2][i4][i3];
                                            }
                                            else {
                                                dbuf[i] = creal(vint[i1][i2][i3][i4] - vint[i1][i2][i4][i3]);
                                            }
                                        }
                                    }
                                    else {
                                        for (size_t i = 0; i < sb->size; i++) {
                                            i1 = symblock_indices[sb->rank * i];  // sb->rank == 4
                                            i2 = symblock_indices[sb->rank * i + 1];
                                            i3 = symblock_indices[sb->rank * i + 2];
                                            i4 = symblock_indices[sb->rank * i + 3];
                                            // global -> local
                                            i1 = spinor_index_global2local[i1];
                                            i2 = spinor_index_global2local[i2];
                                            i3 = spinor_index_global2local[i3];
                                            i4 = spinor_index_global2local[i4];
                                            if (carith) {
                                                sb->buf[i] += (sign == 1) ? vint[i1][i2][i3][i4]
                                                                          : -vint[i1][i2][i4][i3];
                                            }
                                            else {
                                                dbuf[i] += (sign == 1) ? creal(vint[i1][i2][i3][i4])
                                                                       : creal(-vint[i1][i2][i4][i3]);
                                            }
                                        }
                                    }
                                    // for PPPP diagram only: complex conjugate instead of reordering
                                    if (strcmp(req->hp, "pppp") == 0 &&
                                        strcmp(req->valence, "0000") == 0 &&
                                        strcmp(req->order, "3412") == 0 && sign == 1 && carith) {
                                        for (size_t i = 0; i < sb->size; i++) {
                                            sb->buf[i] = conj(sb->buf[i]);
                                        }
                                    }
                                    symblock_store(sb);
                                    break;
                                } // end of loop over blocks
                                //restore_diagram(dg);
                            } // end of loop over requests
                        }
                    }
                }
            }
        } // end of 4 nested 'for' loops over spinor blocks
    } // end if (n_requests_2el > 0)

    printf("   sorting one-electron integrals ...\n");

    // тут нужно обязательно войти внутрь, сколько бы запросов ни было
    sort_onel();

    // только для реально отсортированных в этом запуске запросов!
    // reorder some diagrams (if required)
    printf("   reorder diagrams (if required): ");
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
        printf("%s[%s] ", req->dg_name, req->order);
        reorder(req->dg_name, req->dg_name, req->order);
    }
    printf("done\n");

    // write diagrams to disk
    // только реально обработанные запросы!
    printf("   save sorted diagrams to disk: ");
    for (ireq = 0; ireq < n_requests; ireq++) {
        char dg_file_name[CC_MAX_PATH_LENGTH];
        req = &sorting_requests[ireq];
        sprintf(dg_file_name, "%s.dg", req->dg_name);
        diagram_write(diagram_stack_find(req->dg_name), dg_file_name);
        printf("%s ", req->dg_name);
    }
    printf("done\n");

    timer_stop("sort");
    printf("   number of blocks read from disk: %d\n", n_blocks_read);
    printf("   total number of integrals read from disk: %ld\n", n_integrals_read);
    printf("   total number of bytes read from disk: %ld (%.2f GB)\n", n_bytes_read,
           (double) n_bytes_read / (1024.0 * 1024.0 * 1024.0));
    printf("   sorting performance, Mb/sec: %.2f\n", (double) n_bytes_read / (1024.0 * 1024.0) / timer_get("sort"));
    printf("   time for 2-e integrals sorting, sec: %.2f\n", timer_get("sort"));
    printf("   time for DIRAC interface (integral extraction & write), sec: %.2f\n", timer_get("dirac"));
    printf("   total time for sorting operations, sec: %.2f\n", timer_get("dirac") + timer_get("sort"));

    finalize_sorting();

    printf(" ------------------------------------------------------------------------------------------\n");
}


void reconstruct_fock()
{
    int idx4[4];

    // get pointer to the diagrams required for the construction of the
    // exchange-correlation part of the Fock matrix
    diagram_t *dg_hhhh = diagram_stack_find("hhhh");
    if (dg_hhhh == NULL) {
        errquit("in sort_onel(): diagram 'hhhh' is required to construct Fock matrix, but it was not found");
    }
    diagram_t *dg_hphh = diagram_stack_find("hphh");
    if (dg_hphh == NULL) {
        errquit("in sort_onel(): diagram 'hphh' is required to construct Fock matrix, but it was not found");
    }
    diagram_t *dg_hhhp = diagram_stack_find("hhhp");
    if (dg_hhhp == NULL) {
        errquit("in sort_onel(): diagram 'hhhp' is required to construct Fock matrix, but it was not found");
    }
    diagram_t *dg_hphp = diagram_stack_find("hphp");
    if (dg_hphp == NULL) {
        errquit("in sort_onel(): diagram 'hphp' is required to construct Fock matrix, but it was not found");
    }

    printf("     Fock matrix reconstruction ...\n");
    for (int i = 0; i < nspinors; i++) {
        for (int j = 0; j < nspinors; j++) {
            f_ints[i * nspinors + j] = h_ints[i * nspinors + j];
            int nocc = spinor_types.holes[0];   // size of list of occ spinors
            for (int p = 0; p < nocc; p++) {
                int idx_p = spinor_types.holes[p + 1];
                idx4[0] = idx_p;
                idx4[1] = i;
                idx4[2] = idx_p;
                idx4[3] = j;
                if (is_hole(i)) {
                    if (is_hole(j)) {
                        f_ints[i * nspinors + j] += diagram_get(dg_hhhh, idx4);
                    }
                    else { // 'j is particle'
                        f_ints[i * nspinors + j] += diagram_get(dg_hhhp, idx4);
                    }
                }
                else { // 'i' is particle
                    if (is_hole(j)) {
                        f_ints[i * nspinors + j] += diagram_get(dg_hphh, idx4);
                    }
                    else { // 'j is particle'
                        f_ints[i * nspinors + j] += diagram_get(dg_hphp, idx4);
                    }
                }
            }
        }
    }
}


double recalculate_scf_energy()
{
    int idx4[4];
    double new_escf = cc_opts->enuc;

    diagram_t *dg_hhhh = diagram_stack_find("hhhh");
    if (dg_hhhh == NULL) {
        errquit("in sort_onel(): diagram 'hhhh' is required to construct Fock matrix, but it was not found");
    }

    for (int i = 0; i < nspinors; i++) {
        if (!is_hole(i)) { continue; }
        new_escf += h_ints[i * nspinors + i];
        for (int j = 0; j < nspinors; j++) {
            if (!is_hole(j)) { continue; }
            idx4[0] = i;
            idx4[1] = j;
            idx4[2] = i;
            idx4[3] = j;
            new_escf += 0.5 * diagram_get(dg_hhhh, idx4);
        }
    }

    return new_escf;
}


void sort_onel()
{
    int fd, i, j, ireq;
    diagram_t *dg;
    size_t isb;
    int idx1, idx2;
    block_t *sb;
    int need_update_eps;

    // read one-electron integrals -- core Fock operator
    // ? а что если HINT нет?
    fd = io_open("HINT", "r");
    io_read_compressed(fd, h_ints, sizeof(double complex) * nspinors * nspinors);
    io_close(fd);

    // read one-electron operators from the OneProp code by Leonid V. Skripnikov
    if (cc_opts->oneprop_on) {
        double complex *oper = (double complex *) cc_malloc(nspinors * nspinors * sizeof(double complex));

        for (int ioper = 0; ioper < cc_opts->n_oneprop; ioper++) {
            double complex lambda = cc_opts->oneprop_lambda[ioper];
            printf("     Reading one-electron property from Oneprop "
                   "(lambda = %e %e, file_re = %s, file_im = %s)... \n",
                   creal(lambda), cimag(lambda),
                   cc_opts->oneprop_file_re[ioper], cc_opts->oneprop_file_im[ioper]);

            // read files with (a) real part of the prop matrix (b) imaginary part
            int status = read_prop_two_files(nspinors, cc_opts->oneprop_file_re[ioper], cc_opts->oneprop_file_im[ioper], oper);
            if (status == EXIT_FAILURE) {
                errquit("Cannot open OneProp files '%s' (re) '%s' (im)",
                        cc_opts->oneprop_file_re[ioper], cc_opts->oneprop_file_im[ioper]);
            }

            // add operator matrix to the Fock matrix:
            // F = F + lambda*Oper (lambda = perturbation parameter)
            for (i = 0; i < nspinors; i++) {
                for (j = 0; j < nspinors; j++) {
                    h_ints[i * nspinors + j] += lambda * oper[i * nspinors + j];
                }
            }
        }

        // all done
        cc_free(oper);
    }

    if (cc_opts->n_mdprop > 0) {
        double complex *oper = (double complex *) cc_malloc(nspinors * nspinors * sizeof(double complex));

        for (int ioper = 0; ioper < cc_opts->n_mdprop; ioper++) {
            double complex lambda = cc_opts->mdprop_lambda[ioper];
            printf("     Reading one-electron property %s from MDPROP "
                   "(lambda = %e %e): ",
                   cc_opts->mdprop_file[ioper],
                   creal(lambda), cimag(lambda));

            // read file containing the prop matrix
            int status = read_prop_single_file(nspinors, cc_opts->mdprop_file[ioper], oper);
            if (status == EXIT_FAILURE) {
                errquit("Cannot open property matrix file '%s'", cc_opts->mdprop_file[ioper]);
            }

            // add operator matrix to the Fock matrix:
            // F = F + lambda*Oper (lambda = perturbation parameter)
            for (i = 0; i < nspinors; i++) {
                for (j = 0; j < nspinors; j++) {
                    h_ints[i * nspinors + j] += lambda * oper[i * nspinors + j];
                }
            }
        }

        // all done
        cc_free(oper);
    }

    // construct Fock matrix
    if (!cc_opts->x2cmmf) {
        reconstruct_fock();
    }
    else {  // x2cmmf hamiltonian, no recompute fock matrix
        printf("     Fock matrix reconstruction will be skipped\n");
        memset(f_ints, 0, sizeof(double complex) * nspinors * nspinors);
        for (int i = 0; i < nspinors; i++) {
            f_ints[i * nspinors + i] = spinor_info[i].eps;
        }
    }

    // recalculate SCF energy (and update it if needed)
    if (!cc_opts->x2cmmf) {
        double new_escf = recalculate_scf_energy();
        if (fabs(cc_opts->escf - new_escf) > 1e-10) {
            printf("     SCF energy (energy of reference determinant) was updated:\n");
            printf("       old energy = %20.12f a.u.\n", cc_opts->escf);
            printf("       new energy = %20.12f a.u.\n", new_escf);
            cc_opts->escf = new_escf;
        }
    }

    // oneprop was here

    // compare diagonal elements of the reconstructed Fock matrix with orbital
    // energies read from the MRCONEE file. if they don't coincide then
    // substitute eps[] with new Fock matrix diagonal elements for the subsequent
    // use in perturbation expressions
    need_update_eps = 0;
    double eps_diff = 0.0;
    double max_eps_diff = 0.0;
    for (i = 0; i < nspinors; i++) {
        double complex f_ii = f_ints[i * nspinors + i];
        if (cimag(f_ii) > 1e-13) {
            printf("     (!) in sort_onel(): imaginary value of the [%d,%d] diagonal "
                   "element of the Fock matrix = %.8f %.3E (threshold = 1e-13)\n",
                   i, i, creal(f_ii), cimag(f_ii));
        }
        eps_diff = cabs(f_ii - spinor_info[i].eps);
        if (eps_diff > max_eps_diff) {
            max_eps_diff = eps_diff;
        }
        if (eps_diff > 1e-10) {
            need_update_eps = 1;
        }
    }
    // update diagonal matrix elements
    if (need_update_eps) {
        if (max_eps_diff > 1e-6) {
            printf("     The diagonal elements of the reconstructed Fock matrix don't "
                   "coincide with the orbital energies!\n");
            printf("     NOTE: The diagonal elements of the recomputed Fock matrix "
                   "(right column) are used in perturbation expressions.\n");

            printf("\n");
            printf("     no    rep         occ    active     one-el energy        recalc energy            delta    \n");
            printf("    ---------------------------------------------------------------------------------------------\n");
            for (i = 0; i < nspinors; i++) {
                double eps = spinor_info[i].eps;
                double f_ii = creal(f_ints[i * nspinors + i]);
                printf("    %4d%4d \"%-6s\"%4d       %1s     %16.10f     %16.10f     %16.6e\n",
                       spinor_info[i].seqno, spinor_info[i].repno,
                       rep_names[spinor_info[i].repno], spinor_info[i].occ,
                       (spinor_info[i].active == 1) ? "a" : "i", eps, f_ii, f_ii - eps);
            }
            printf("    ---------------------------------------------------------------------------------------------\n");
        }
        printf("     max deviation of the diagonal elements of the reconstructed"
               " Fock matrix and orbital energies = %.6e\n", max_eps_diff);
        // recalculate energies
        for (i = 0; i < nspinors; i++) {
            spinor_info[i].eps = creal(f_ints[i * nspinors + i]);
        }
    }

    // только запросы ?
    printf("     fill 1-electron diagrams ... ");
    for (ireq = 0; ireq < n_requests; ireq++) {
        dg = sorting_requests[ireq].dg;
        if (dg->rank != 2) { // only one-electron diagrams
            continue;
        }
        printf("%s ", dg->name);
        // assign values of matrix elements
        for (isb = 0; isb < dg->n_blocks; isb++) {
            sb = dg->blocks[isb];
            symblock_load(sb);
            symblock_gen_indices(sb, symblock_indices);
            double *dbuf = (double *) sb->buf;

            for (i = 0; i < sb->size; i++) {
                idx1 = symblock_indices[sb->rank * i];  // sb->rank == 2
                idx2 = symblock_indices[sb->rank * i + 1];
                if (idx1 == idx2) {
                    if (carith) {
                        sb->buf[i] = 0.0 + 0.0 * I;
                    }
                    else {
                        dbuf[i] = 0.0;
                    }
                }
                else {
                    if (carith) {
                        sb->buf[i] = f_ints[idx1 * nspinors + idx2];
                    }
                    else {
                        dbuf[i] = creal(f_ints[idx1 * nspinors + idx2]);
                    }
                }
            }
            symblock_unload(sb);
        }
    }
    printf("done\n");
}


/**
 * Reads pair of formatted files containing real square matrices of properties.
 * NOTE: transposition is included to be able to use files from DIRAC.
 *
 * @param nspinors number of spinors
 * @param file_re real part
 * @param file_im imag part
 * @param prop_mat (is returned) nspinors x nspinors complex matrix
 */
int read_prop_two_files(int nspinors, char *file_re, char *file_im, double complex *prop_mat)
{
    FILE *f_prop;
    int idx_i, idx_j;
    double value;

    memset(prop_mat, 0, nspinors * nspinors * sizeof(double complex));

    // read files with (a) real part of the prop matrix
    f_prop = fopen(file_re, "r");
    if (f_prop == NULL) {
        return EXIT_FAILURE;
    }
    while (fscanf(f_prop, "%d%d%lf", &idx_i, &idx_j, &value) == 3) {
        prop_mat[(idx_j - 1) * nspinors + (idx_i - 1)] = value + 0.0 * I;
    }
    fclose(f_prop);

    // (b) imaginary part
    f_prop = fopen(file_im, "r");
    if (f_prop == NULL) {
        return EXIT_FAILURE;
    }
    while (fscanf(f_prop, "%d%d%lf", &idx_i, &idx_j, &value) == 3) {
        prop_mat[(idx_j - 1) * nspinors + (idx_i - 1)] += 0.0 + value * I;
    }
    fclose(f_prop);

    return EXIT_SUCCESS;
}


/**
 * Reads property matrix from file 'prop_name' (extracted from MDPROP)
 *
 * @param nspinors number of spinors
 * @param prop_name name of the property under consideration in the MDPROP file
 * @param prop_mat (is returned) nspinors x nspinors complex matrix
 */
int read_prop_single_file(int nspinors, char *prop_name, double complex *prop_mat)
{
    FILE *f_prop;
    int idx_i, idx_j;
    double re_value, im_value;

    memset(prop_mat, 0, nspinors * nspinors * sizeof(double complex));

    // read file with the prop matrix
    f_prop = fopen(prop_name, "r");
    if (f_prop == NULL) {
        return EXIT_FAILURE;
    }
    while (fscanf(f_prop, "%d%d%lf%lf", &idx_i, &idx_j, &re_value, &im_value) == 4) {
        prop_mat[(idx_i - 1) * nspinors + (idx_j - 1)] = re_value + im_value * I;
    }
    fclose(f_prop);

    return EXIT_SUCCESS;
}


void finalize_sorting()
{
    size_t i, j, k;

    n_requests = 0;

    cc_free(h_ints);
    cc_free(f_ints);
    cc_free(symblock_indices);

    for (i = 0; i < max_spinor_block_size; i++) {
        for (j = 0; j < max_spinor_block_size; j++) {
            for (k = 0; k < max_spinor_block_size; k++) {
                cc_free(vint[i][j][k]);
            }
            cc_free(vint[i][j]);
        }
        cc_free(vint[i]);
    }
    cc_free(vint);
}

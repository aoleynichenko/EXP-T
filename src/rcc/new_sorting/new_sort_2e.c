/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
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

#include "new_sort_2e.h"

#include <spinors.h>

#include "error.h"
#include "options.h"
#include "libunf.h"
#include "../sorting/sorting_request.h"
#include "../include/timer.h"
#include "mrconee.h"
#include "../include/engine.h"

#include <stdio.h>
#include <stdlib.h>

enum {
    INT_CLASS_NO_BARS = 0,
    INT_CLASS_ONE_BAR,
    INT_CLASS_TWO_BARS
};

static int sign(int x);

void expand_ints(diagram_t *dg, int ikr, int jkr, int nonzr, void *ind_buf, double *val_buf, int32_t *kr, int nkr);

int int_class(int ikr, int jkr, int kkr, int lkr);

void put_integral(diagram_t *dg, int i, int j, int k, int l, double complex val);

void perm_symm(diagram_t *dg, int i, int j, int k, int l, double complex val, int nkr, int32_t *kr);

int kr2abs(int ikr, int32_t *kr, int nkr);


void new_sort_2e()
{
    printf("new sort 2e\n");
    mrconee_data_t *mrconee_data = cc_opts->mrconee_data;

    if (num_twoelec_requests(sorting_requests, n_requests) == 0) {
        // nothing to sort
        return;
    }

    unf_file_t *mdcint = unf_open(cc_opts->integral_file_2, "r", UNF_ACCESS_SEQUENTIAL);
    if (mdcint == NULL) {
        errquit("cannot open MDCINT unformatted file\n");
        return;
    }

    int use_int4 = mrconee_data->dirac_int_size == 4;
    int use_int8 = mrconee_data->dirac_int_size == 8;

    printf(" two-electron integrals file\n");

    /*
     * read date and time, total number of Kramers pairs
     */
    char date_time[100];
    int32_t nkr = 0;
    int64_t nkr_8 = 0;

    int nread;
    if (use_int4) {
        nread = unf_read(mdcint, "c18,i4", date_time, &nkr);
    }
    else {
        nread = unf_read(mdcint, "c18,i8", date_time, &nkr_8);
        nkr = (int32_t) nkr_8;
    }
    if (nread != 2 || unf_error(mdcint)) {
        perror(" error while reading MDCINT file");
        return;
    }

    /*
     * also read indices of Kramers pairs
     */
    unf_backspace(mdcint);

    int32_t num_spinors = 2 * nkr;
    int32_t *kr = (int32_t *) calloc(num_spinors, sizeof(int32_t *));
    int64_t *kr_8 = (int64_t *) calloc(num_spinors, sizeof(int64_t *));

    if (use_int4) {
        nread = unf_read(mdcint, "c18,i4,i4[i4]", date_time, &nkr, kr, &num_spinors);
    }
    else {
        nread = unf_read(mdcint, "c18,i8,i8[i4]", date_time, &nkr_8, kr_8, &num_spinors);
    }
    if (nread != 3 || unf_error(mdcint)) {
        perror(" error while reading MDCINT file");
        return;
    }

    if (use_int8) {
        for (int i = 0; i < num_spinors; i++) {
            kr[i] = (int32_t) kr_8[i];
        }
    }
    free(kr_8);

    date_time[18] = '\0';
    printf(" date and time              %s\n", date_time);
    printf(" number of Kramers pairs    %d\n", nkr);
    for (int i = 0; i < nkr; i++) {
        printf("%4d%4d    %4d%4d\n", kr[2*i], kr[2*i + 1], kr2abs(i + 1, kr, nkr), kr2abs(-(i + 1), kr, nkr));
    }

    /*
     * read chunks of non-zero two-electron integrals
     *
     * read (luint, end = 1301, err = 1302) ikr, jkr, nonzr, &
     * (indk(inz), indl(inz), inz = 1, nonzr), &
     * (cbuf(1, inz), inz = 1, nonzr)
     */
    int64_t count_non_zero = 0;

    char *ind_buf = (char *) calloc(num_spinors * num_spinors, 2 * sizeof(int32_t));
    char *ind_buf_8 = (char *) calloc(num_spinors * num_spinors, 2 * sizeof(int64_t));
    double *val_buf_real = (double *) calloc(num_spinors * num_spinors, sizeof(double));
    double *val_buf_complex = (double *) calloc(num_spinors * num_spinors, sizeof(double));

    int32_t ikr = 0;
    int32_t jkr = 0;
    int32_t nonzr = 0;
    int64_t ikr8 = 0;
    int64_t jkr8 = 0;
    int64_t nonzr8 = 0;
    double time_start = abs_time();

    while (true) {
        if (mrconee_data->group_arith == 1 || mrconee_data->is_spinfree == 1) {
            if (use_int4) {
                nread = unf_read(mdcint, "3i4,c8[i4],r8[i4]", &ikr, &jkr, &nonzr, ind_buf, &nonzr, val_buf_real,
                                 &nonzr);

                for (int ireq = 0; ireq < n_requests; ireq++) {
                    sorting_request_t *req = &sorting_requests[ireq];
                    diagram_t *dg = req->dg;
                    if (dg->rank != 4) { // one-electron diagrams will be sorted later
                        continue;
                    }

                    expand_ints(dg, ikr, jkr, nonzr, ind_buf, val_buf_real, kr, nkr);
                }

            }
            else {
                nread = unf_read(mdcint, "3i8,c16[i8],r8[i8]", &ikr8, &jkr8, &nonzr8, ind_buf_8, &nonzr8, val_buf_real,
                                 &nonzr8);
                ikr = (int32_t) ikr8;
                jkr = (int32_t) jkr8;
                nonzr = (int32_t) nonzr8;
            }
        }
        else {
            if (use_int4) {
                nread = unf_read(mdcint, "3i4,c8[i4],z8[i4]", &ikr, &jkr, &nonzr, ind_buf, &nonzr, val_buf_complex,
                                 &nonzr);
            }
            else {
                nread = unf_read(mdcint, "3i8,c16[i8],z8[i8]", &ikr8, &jkr8, &nonzr8, ind_buf_8, &nonzr8,
                                 val_buf_complex, &nonzr8);
                ikr = (int32_t) ikr8;
                jkr = (int32_t) jkr8;
                nonzr = (int32_t) nonzr8;
            }
        }

        if (nread != 5 || unf_error(mdcint)) {

            perror(" error while reading MDCINT file");
            break;
        }

        if (ikr == 0 && jkr == 0) {
            break;
        }

        count_non_zero += nonzr;
    }

    double time_finish = abs_time();
    printf(" number of non-zero ints    %lld\n", count_non_zero);
    printf(" time for reading 2e ints   %.2f sec\n\n", time_finish - time_start);


    /*
     * cleanup
     */

    unf_close(mdcint);

    free(ind_buf);
    free(ind_buf_8);
    free(val_buf_real);
    free(val_buf_complex);
    free(kr);

    print_amplitude_distribution_analysis("hhpp");
    //print_amplitude_distribution_analysis("pphh");

    //exit(0);
}


// ikr = [1 ... nkr]
int kr2abs(int ikr, int32_t *kr, int nkr)
{
    if (ikr < 0) {
        return kr[2 * (-ikr-1) + 1];
    }
    else {
        return kr[2 * (ikr-1)];
    }
}


void expand_ints(diagram_t *dg, int ikr, int jkr, int nonzr, void *ind_buf, double *val_buf, int32_t *kr, int nkr)
{
    int32_t *indk = (int32_t *) calloc(nonzr, sizeof(int32_t));
    int32_t *indl = (int32_t *) calloc(nonzr, sizeof(int32_t));

    for (int i = 0; i < nonzr; i++) {
        indk[i] = ((int32_t *) ind_buf)[2 * i];
        indl[i] = ((int32_t *) ind_buf)[2 * i + 1];
    }

    int kclass = int_class(ikr, indk[0], jkr, indl[0]);

    for (int idx = 0; idx < nonzr; idx++) {
        int kkr = indk[idx];
        int lkr = indl[idx];
        double complex val = val_buf[idx] + 0.0 * I;

        // pass to Dirac notation: (ij|kl) -> <ik|jl>
        int ikrd = ikr;
        int jkrd = kkr;
        int kkrd = jkr;
        int lkrd = lkr;

        //printf("%3d%3d%3d%3d\n", ikrd, jkrd, kkrd, lkrd);

        if (kclass == INT_CLASS_NO_BARS || kclass == INT_CLASS_TWO_BARS) {
            perm_symm(dg, ikrd, jkrd, kkrd, lkrd, val, nkr, kr);
            if (cc_opts->mrconee_data->is_spinfree == 1) {
                perm_symm(dg, ikrd, -lkrd, kkrd, -jkrd, val, nkr, kr);  // MUST be in NR case
                perm_symm(dg, -kkrd, jkrd, -ikrd, lkrd, val, nkr, kr);  // MUST be in NR case
            }
            perm_symm(dg, -kkrd, -lkrd, -ikrd, -jkrd, val, nkr, kr);
        }
        else {
            perm_symm(dg, ikrd, jkrd, kkrd, lkrd, val, nkr, kr);
            perm_symm(dg, -kkrd, -lkrd, -ikrd, -jkrd, -val, nkr, kr);
        }


        /*
            ! construct integral
            cint = cbuf(1, inz) + cbuf(2, inz) * (0.0d0, 1.0d0)

            if (kclass == INTCLASS_NO_BARS .or. kclass == INTCLASS_TWO_BARS) then
                call perm_symm(cint, kr(ikrd), kr(jkrd), kr(kkrd), kr(lkrd))
                if (is_spinfree == 1) then
                    call perm_symm(cint, kr(ikrd), kr(-lkrd), kr(kkrd), kr(-jkrd))   ! MUST be in NR case
                    call perm_symm(cint, kr(-kkrd), kr(jkrd), kr(-ikrd), kr(lkrd))   ! MUST be in NR case
                end if
                call perm_symm(cint, kr(-kkrd), kr(-lkrd), kr(-ikrd), kr(-jkrd))

            else if (kclass == INTCLASS_ONE_BAR) then
                call perm_symm(cint, kr(ikrd), kr(jkrd), kr(kkrd), kr(lkrd))
                !call perm_symm(cint,kr(ikrd),kr(-lkrd),kr(kkrd),kr(-jkrd))
                !call perm_symm(-cint,kr(-kkrd),kr(jkrd),kr(-ikrd),kr(lkrd))
                call perm_symm(-cint, kr(-kkrd), kr(-lkrd), kr(-ikrd), kr(-jkrd))
            end if
         */
    }

    free(indk);
    free(indl);
}


void perm_symm(diagram_t *dg, int ikr, int jkr, int kkr, int lkr, double complex val, int nkr, int32_t *kr)
{
    int i = kr2abs(ikr, kr, nkr) - 1;
    int j = kr2abs(jkr, kr, nkr) - 1;
    int k = kr2abs(kkr, kr, nkr) - 1;
    int l = kr2abs(lkr, kr, nkr) - 1;

    put_integral(dg, i, j, k, l, val);
    put_integral(dg, j, i, l, k, val);
    put_integral(dg, k, l, i, j, conj(val));
    put_integral(dg, l, k, j, i, conj(val));

    put_integral(dg, i, j, l, k, -val);
    put_integral(dg, j, i, k, l, -val);
    put_integral(dg, k, l, j, i, -conj(val));
    put_integral(dg, l, k, i, j, -conj(val));
}


void put_integral(diagram_t *dg, int i, int j, int k, int l, double complex val)
{
    int idx[4];
    idx[0] = i;
    idx[1] = j;
    idx[2] = k;
    idx[3] = l;
    /*idx[0] = spinor_index_global2local[i];
    idx[1] = spinor_index_global2local[j];
    idx[2] = spinor_index_global2local[k];
    idx[3] = spinor_index_global2local[l];*/

    //printf("set %d/%d/%d/%d   %.6g\n", i, j, k, l, creal(val));


    diagram_set(dg, val, idx);
}


/*
    subroutine perm_symm(cint, i, j, k, l)
        complex(8), intent(in) :: cint
        integer(4), intent(in) :: i, j, k, l

        call put_integral(int(i, 2), int(j, 2), int(k, 2), int(l, 2), cint)
        call put_integral(int(j, 2), int(i, 2), int(l, 2), int(k, 2), cint)
        call put_integral(int(k, 2), int(l, 2), int(i, 2), int(j, 2), conjg(cint))
        call put_integral(int(l, 2), int(k, 2), int(j, 2), int(i, 2), conjg(cint))

    end subroutine perm_symm
 */



int int_class(int ikr, int jkr, int kkr, int lkr)
{
    int sign_product = sign(ikr) * sign(jkr) * sign(kkr) * sign(lkr);
    if (sign_product > 0) {
        // INT_CLASS_TWO_BARS is the same for Coulomb integrals
        return INT_CLASS_NO_BARS;
    }
    else {
        return INT_CLASS_ONE_BAR;
    }
}


static int sign(int x)
{
    if (x < 0) {
        return -1;
    }
    else if (x > 0) {
        return 1;
    }
    else {
        return 0;
    }
}
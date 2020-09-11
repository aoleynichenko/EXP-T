/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2020 The EXP-T developers.
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
* Utility functions for FSCC models.
*
* 2020 Alexander Oleynichenko
*/

#include "ccutils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "engine.h"
#include "datamodel.h"
#include "options.h"
#include "symmetry.h"

/**
 * flushes to file sizes of all blocks in the diagram
 * (in order to build distribution)
 */
void flush_block_sizes(char *diag_name, char *file_name)
{
    diagram_t *diag = diagram_stack_find(diag_name);
    if (diag == NULL) {
        errquit("flush_block_sizes(): diagram not found");
    }
    FILE *f = fopen(file_name, "w");
    for (size_t i = 0; i < diag->n_blocks; i++) {
        fprintf(f, "%d\n", diag->blocks[i]->size);
    }
    fclose(f);
}


/**
 * Damping of cluster amplitudes
 */
void damping(int h, int p, char *old_ampl, char *new_ampl, int iter)
{
    assert(h >= 0 && h < MAX_SECTOR_RANK);
    assert(p >= 0 && p < MAX_SECTOR_RANK);

    if (cc_opts->damping[h][p].enabled &&
        iter <= cc_opts->damping[h][p].stop) {
        double alpha = cc_opts->damping[h][p].factor;
        if (cc_opts->print_level >= CC_PRINT_DEBUG) {
            printf("\n damping: %.3f %s[iter=%d] + %.3f %s[iter=%d]\n",
                   1.0 - alpha, new_ampl, iter, alpha, old_ampl, iter - 1);
        }
        add(1.0 - alpha, new_ampl, alpha, old_ampl, new_ampl);
    }
}


/**
 * Prints norm of the operator represented by the diagram 'dg_name'.
 * Norm = sqrt ( sum |a_i|^2 for a_i in diagram 'a' )
 */
void print_ampl_norm(char *label, char *dg_name)
{
    double complex norm2z = scalar_product("C", "N", dg_name, dg_name);
    double norm2 = creal(norm2z);
    double norm = sqrt(norm2);
    printf(" %s = ", label);
    if (norm > 1e-4) {
        printf("%.8f\n", norm);
    }
    else {
        printf("%.4E\n", norm);
    }
}


/**
 * Prints max (by abs value) element of diagram 'dg_name'
 */
void print_max_ampl(char *label, char *dg_name)
{
    double max_t;
    int max_idx[CC_DIAGRAM_MAX_RANK];
    int i;
    int rk;
    int occ, active, repno;
    double eps;
    char rep_name[256];

    rk = rank(dg_name);
    findmax(dg_name, &max_t, max_idx);
    printf(" %s [ ", label);
    for (i = 0; i < rk / 2; i++){
        printf("%d ", max_idx[i]+1);
    }
    printf("-> ");
    for (i = rk / 2; i < rk; i++){
        printf("%d ", max_idx[i]+1);
    }
    printf("] = ");
    if (max_t > 1e-4) {
        printf("%.8f\n", max_t);
    }
    else {
        printf("%.4E\n", max_t);
    }

    // print more information about these spinors in order to simplify analysis
    for (i = 0; i < rk; i++){
        get_spinor_info(max_idx[i], &repno, &occ, &active, &eps);
        printf("   [%4d] eps=%16.8f rep=%s %s %s\n",
                max_idx[i]+1, eps, get_irrep_name(repno),
               (active == 1) ? "active" : "inactive",
               (occ == 1) ? "occ" : "virt");
    }
}


double t1_diagnostic(char *dg_name, int nelec)
{
    double t1diagn = 0.0;

    diagram_t *dg = diagram_stack_find(dg_name);
    if (dg == NULL) {
        errquit("t1_diagnostic(): diagram '%s' not found", dg_name);
    }

    for (size_t isb = 0; isb < dg->n_blocks; isb++){
        block_t *sb = dg->blocks[isb];
        symblock_load(sb);
        for (size_t i = 0; i < sb->size; i++){
            double abs_val = carith ? cabs(sb->buf[i]) : fabs(((double *)sb->buf)[i]);
            t1diagn += abs_val * abs_val;
        }
        symblock_unload(sb);
    }

    t1diagn = sqrt(t1diagn / nelec);
    return t1diagn;
}


/**
 * Checks if one of cluster amplitudes |t| > 1
 * If |t| > 1 returns 1, if all is OK returns 0
 */
int check_divergence(char *dg_name)
{
    double max_t;
    int max_idx[CC_DIAGRAM_MAX_RANK];

    findmax(dg_name, &max_t, max_idx);

    if (max_t > 1.0) {
        return 1;
    }
    else{
        return 0;
    }
}


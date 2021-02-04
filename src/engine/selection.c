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

/*******************************************************************************
 * selection.c
 *
 * Tools for selection of cluster amplitudes which satisfy some condition.
 * All other amplitudes will be set to zero.
 * These subroutines not affect performance; they are designed just for fast
 * construction and testing of new approximate schemes.
 *
 * 2020 Alexander Oleynichenko
 ******************************************************************************/

#include "selection.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "datamodel.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"


void set_zero(double complex *t)
{
    if (carith) {
        *t = 0.0 + 0.0*I;
    }
    else {
        *((double *)t) = 0.0;
    }
}


void selection_all(ampl_selection_t *rule, int *idx, double complex *t)
{
    if (rule->task == CC_SELECTION_SET_ZERO_EXCEPT) {
        return;
    }

    set_zero(t);
}


void selection_exc_window(ampl_selection_t *rule, int *idx, double complex *data)
{
    double exc_energy = 0;  // sum eps{part} - sum eps{hole}, normally positive
    for (int i = 0; i < rule->rank/2; i++) {
        exc_energy -= spinor_info[idx[i]].eps;
    }
    for (int i = rule->rank/2; i < rule->rank; i++) {
        exc_energy += spinor_info[idx[i]].eps;
    }

    if (rule->e1 <= exc_energy && exc_energy <= rule->e2) {
        if (rule->task == CC_SELECTION_SET_ZERO) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
    else {
        if (rule->task == CC_SELECTION_SET_ZERO_EXCEPT) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
}


void selection_eps_window(ampl_selection_t *rule, int *idx, double complex *data)
{
    int in_window = 1;
    for (int i = 0; i < rule->rank; i++) {
        double eps = spinor_info[idx[i]].eps;
        if (!(rule->e1 <= eps && eps <= rule->e2)) {
            in_window = 0;
        }
    }


    if (in_window) {
        if (rule->task == CC_SELECTION_SET_ZERO) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
    else {
        if (rule->task == CC_SELECTION_SET_ZERO_EXCEPT) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
}


// what about epv? (two or more coinciding indices)
void selection_spectator(ampl_selection_t *rule, int *idx, double complex *data)
{
    int n_spect_lines = 0;

    for (int i = 0; i < rule->rank/2; i++) {
        if (!is_active(idx[i])) {
            continue;
        }
        for (int j = rule->rank/2; j < rule->rank; j++) {
            if (idx[i] == idx[j])  {
                n_spect_lines++;
                break;
            }
        }
    }

    if (n_spect_lines == (rule->sect_h + rule->sect_p)) {
        if (rule->task == CC_SELECTION_SET_ZERO) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
    else {
        if (rule->task == CC_SELECTION_SET_ZERO_EXCEPT) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
}


// at least one "bra" index should be valence
void selection_act_to_act(ampl_selection_t *rule, int *idx, double complex *data)
{
    int n_val_indices = 0;
    for (int i = rule->rank/2; i < rule->rank; i++) {
        if (is_active(idx[i])) {
            n_val_indices++;
        }
    }

    if (n_val_indices > 0) {
        if (rule->task == CC_SELECTION_SET_ZERO) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
    else {
        if (rule->task == CC_SELECTION_SET_ZERO_EXCEPT) {
            set_zero(data);
        }
        else { /* do nothing*/ }
    }
}


void apply_selections(int sect_h, int sect_p, char *diag_name)
{
    diagram_t *dg;

    dg = diagram_stack_find(diag_name);
    if (dg == NULL) {
        errquit("apply_selection(): diagram '%s' not found", diag_name);
    }

    int rank = dg->rank;

    // for each selection rule
    for (int irule = 0; irule < cc_opts->n_select; irule++) {
        if (cc_opts->selects[irule].sect_h != sect_h ||
            cc_opts->selects[irule].sect_p != sect_p ||
            cc_opts->selects[irule].rank != rank) {
            continue;
        }

        ampl_selection_t *rule = cc_opts->selects + irule;

        for (size_t ib = 0; ib < dg->n_blocks; ib++) {
            if (dg->blocks[ib]->is_unique == 0) {
                restore_block(dg, dg->blocks[ib]);
            }

            block_t *block = dg->blocks[ib];

            size_t i, j;
            int *indices = (int *) cc_malloc(block->size * block->rank * sizeof(int));

            symblock_load(block);
            symblock_gen_indices(block, indices);

            for (i = 0; i < block->size; i++) {
                int *idx = indices + block->rank * i;
                double complex *data = carith ? (block->buf + i) : (double complex *)((double*)block->buf + i);

                switch (rule->rule) {
                    case CC_SELECTION_ALL:
                        selection_all(rule, idx, data);
                        break;
                    case CC_SELECTION_SPECTATOR:
                        selection_spectator(rule, idx, data);
                        break;
                    case CC_SELECTION_ACT_TO_ACT:
                        selection_act_to_act(rule, idx, data);
                        break;
                    case CC_SELECTION_EXC_WINDOW:
                        selection_exc_window(rule, idx, data);
                        break;
                    case CC_SELECTION_EPS_WINDOW:
                        selection_eps_window(rule, idx, data);
                        break;
                    default:
                        break;
                }
            }

            symblock_unload(block);

            cc_free(indices);

            if (dg->blocks[ib]->is_unique == 0) {
                destroy_block(dg->blocks[ib]);
            }
        }
    }
}

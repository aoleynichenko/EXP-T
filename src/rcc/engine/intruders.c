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

/*
 * Analysis of possible intruder states.
 */

#include <math.h>
#include <stdio.h>

#include "engine.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"
#include "utils.h"

typedef struct {
    int rank;
    int idx[CC_DIAGRAM_MAX_RANK];
    double denom;
} indexed_denom_t;


static void save_denominator(indexed_denom_t *denoms, int len_denoms, int dim, int *idx, double new_denom);

static void print_denominator(indexed_denom_t *denom);

static bool contain_repeating_numbers(int n, const int *arr);


/**
 * Denominator-based prediction of possible intruder states.
 *
 * @param diagram_name symbolic name of the diagram containing
 *      some cluster operator
 * @param nmax print information about only 'nmax' intruders
 *      with largest energy denominators
 */
void predict_intruders_for_diagram(char *diagram_name, int nmax)
{
    int idx[CC_DIAGRAM_MAX_RANK];

    timer_new_entry("intruders", "Prediction of intruder states");
    timer_start("intruders");

    assert_diagram_exists(diagram_name);
    diagram_t *dg = diagram_stack_find(diagram_name);
    int rank = dg->rank;

    indexed_denom_t *max_denoms = (indexed_denom_t *) cc_calloc(nmax, sizeof(indexed_denom_t));

    /*
     * extract one-electron energies
     */
    int nspinors = get_num_spinors();
    double *eps = (double *) cc_malloc(nspinors * sizeof(double));
    get_spinor_energies(nspinors, eps);

    /*
     * collect information about energy denominators
     */
    for (size_t iblock = 0; iblock < dg->n_blocks; iblock++) {
        block_t *block = dg->blocks[iblock];
        if (block->is_unique == 0) {
            continue;
        }

        // loop over matrix elements
        for (int i = 0; i < block->size; i++) {
            tensor_index_to_compound(rank, block->shape, i /* linear index */, idx);

            // relative indices -> spinor indices
            for (int j = 0; j < rank; j++) {
                idx[j] = block->indices[j][idx[j]];
            }

            // amplitudes with all active indices must be excluded
            int n_active = 0;
            for (int j = 0; j < rank; j++) {
                if (is_active(idx[j])) {
                    n_active++;
                }
            }
            if (n_active == rank) {
                continue;
            }

            // exclude epv amplitudes (which are always zero), example: t_ii^ab
            if (contain_repeating_numbers(rank / 2, idx) ||
                contain_repeating_numbers(rank / 2, idx + rank / 2)) {
                continue;
            }

            // construct energy denominator
            double denom = 0.0;
            for (int j = 0; j < rank / 2; j++) {
                denom += eps[idx[j]];
            }
            for (int j = rank / 2; j < rank; j++) {
                denom -= eps[idx[j]];
            }

            save_denominator(max_denoms, nmax, rank, idx, denom);
        }
    }

    /*
     * print max energy denominators
     */
    int n_positive_denoms = 0;
    for (int iden = 0; iden < nmax; iden++) {
        if (max_denoms[iden].denom > 0) {
            n_positive_denoms++;
        }
    }

    if (n_positive_denoms > 0) {
        printf("\n Warning: intruder states are possible!\n");

        for (int iden = 0; iden < nmax && n_positive_denoms > 0; iden++) {
            if (max_denoms[iden].denom <= 0) {
                continue;
            }
            print_denominator(max_denoms + iden);
        }

        printf("\n");
    }

    cc_free(max_denoms);
    cc_free(eps);

    timer_stop("intruders");
}


/*
 * puts information about energy denominator into the 'denoms' array
 * if it exceeds (by absolute value) one of the denominators already stored.
 */
static void save_denominator(indexed_denom_t *denoms, int len_denoms, int dim, int *idx, double new_denom)
{
    for (int i = 0; i < len_denoms; i++) {
        if (new_denom > denoms[i].denom) {
            denoms[i].rank = dim;
            denoms[i].denom = new_denom;
            intcpy(denoms[i].idx, idx, dim);
            return;
        }
        else if (fabs(new_denom - denoms[i].denom) < 1e-6) { // omit repeating max values
            return;
        }
    }
}


static void print_denominator(indexed_denom_t *denom)
{
    int rk = denom->rank;

    printf(" denominator [ ");
    for (int i = 0; i < rk / 2; i++) {
        printf("%d ", denom->idx[i] + 1);
    }
    printf("-> ");
    for (int i = rk / 2; i < rk; i++) {
        printf("%d ", denom->idx[i] + 1);
    }
    printf("] = ");
    printf("%.4E\n", denom->denom);

    // print more information about these spinors in order to simplify analysis
    for (int i = 0; i < rk; i++) {
        int ispinor = denom->idx[i];
        double spinor_energy;
        int repno;
        int occ;
        int active;

        get_spinor_info(ispinor, &repno, &occ, &active, &spinor_energy);

        printf("     [%4d] eps=%16.8f rep=%s %s %s",
               ispinor + 1, spinor_energy, get_irrep_name(repno),
               (active == 1) ? "active" : "inactive",
               (occ == 1) ? "occ" : "virt");

        if (cc_opts->spinor_labels[ispinor] != NULL) {
            printf(" (%s)", cc_opts->spinor_labels[ispinor]);
        }

        printf("\n");
    }
}


static bool contain_repeating_numbers(int n, const int *arr)
{
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (arr[i] == arr[j]) {
                return true;
            }
        }
    }

    return false;
}

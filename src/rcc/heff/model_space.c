/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
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
 * Model space in FS-CC consists of Slater determinants.
 * Current implementation deals with only the quasi-complete model spaces.
 *
 * 2019-2021 Alexander Oleynichenko
 */

#include "model_space.h"

#include <string.h>

#include "memory.h"
#include "spinors.h"
#include "symmetry.h"


/**
 * Returns total number of model-space determinants for the given FS sector
 */
size_t get_model_space_size(int sect_h, int sect_p, int nspinors, spinor_attr_t *spinor_info)
{
    int nacth, nactp;
    size_t ndeth, ndetp;

    get_active_space_size(&nacth, &nactp);

    // count number of independent determinants
    // 1 added electron:  Ndet = nactp
    // 2 added electrons: Ndet = nactp * (nactp-1)/2
    // ... (independent off-diagonal elements of n-dim matrix)

    ndeth = 1;
    ndetp = 1;

    // contribution from holes
    if (sect_h == 1) {
        ndeth = nacth;
    }
    else if (sect_h == 2) {
        ndeth = nacth * (nacth - 1) / 2;
    }
    else if (sect_h == 3) { // ?? / n! ?
        ndeth = nacth * (nacth - 1) * (nacth - 2) / 6;
    }
    else if (sect_h == 4) {
        ndeth = nacth * (nacth - 1) * (nacth - 2) * (nacth - 3) / 24;
    }

    // contibution from particles
    if (sect_p == 1) {
        ndetp = nactp;
    }
    else if (sect_p == 2) {
        ndetp = nactp * (nactp - 1) / 2;
    }
    else if (sect_p == 3) { // ?? / n! ?
        ndetp = nactp * (nactp - 1) * (nactp - 2) / 6;
    }
    else if (sect_p == 4) {
        ndetp = nactp * (nactp - 1) * (nactp - 2) * (nactp - 3) / 24;
    }

    return ndeth * ndetp;
}


/**
 * Returns list of model-space Slater determinants and number of determinants
 * belonging to each symmetry (via the 'ms_rep_sizes' array).
 * NOTE: the detlist[] array MUST be pre-allocated
 */
slater_det_t **construct_model_space(int sect_h, int sect_p, size_t *ms_rep_sizes)
{
    moindex_t *active_holes_indices;
    moindex_t *active_parts_indices;
    int nacth, nactp;
    size_t count;
    size_t total_ms_size = get_model_space_size(sect_h, sect_p, get_num_spinors(), spinor_info);
    slater_det_t *det_buf = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * total_ms_size);

    // set counters to zero
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        ms_rep_sizes[irep] = 0;
    }

    int nspinors = get_num_spinors();
    active_holes_indices = (moindex_t *) cc_malloc(nspinors * sizeof(moindex_t));
    active_parts_indices = (moindex_t *) cc_malloc(nspinors * sizeof(moindex_t));

    get_active_holes_particles(&nacth, &nactp, active_holes_indices, active_parts_indices);

    // in order to avoid dynamically-nested loops we shall consider
    // different Fock-space sectors simply as different cases

    if (sect_h == 0 && sect_p == 0) {
        set_vacuum_det(&det_buf[0]);
        int vac_irrep = get_vacuum_irrep();
        ms_rep_sizes[vac_irrep]++;
        count = 1;
    }
    else if (sect_h == 1 && sect_p == 0) {
        count = 0;
        for (int k = 0; k < nacth; k++) {
            moindex_t idx_k = active_holes_indices[k];
            int rep_k = spinor_info[idx_k].repno;
            det_buf[count].indices[0] = idx_k;
            det_buf[count].sym = rep_k;
            ms_rep_sizes[rep_k]++;
            count++;
        }
    }
    else if (sect_h == 0 && sect_p == 1) {
        count = 0;
        for (int a = 0; a < nactp; a++) {
            moindex_t idx_a = active_parts_indices[a];
            int rep_a = spinor_info[idx_a].repno;
            det_buf[count].indices[0] = idx_a;
            det_buf[count].sym = rep_a;
            ms_rep_sizes[rep_a]++;
            count++;
        }
    }
    else if (sect_h == 0 && sect_p == 2) {
        count = 0;
        for (int a = 0; a < nactp; a++) {
            for (int b = a + 1; b < nactp; b++) {
                moindex_t idx_a = active_parts_indices[a];
                moindex_t idx_b = active_parts_indices[b];
                det_buf[count].indices[0] = idx_a;
                det_buf[count].indices[1] = idx_b;

                int rep_a = spinor_info[idx_a].repno;
                int rep_b = spinor_info[idx_b].repno;
                int rep_ab = mulrep2_abelian(rep_a, rep_b);
                det_buf[count].sym = rep_ab;
                ms_rep_sizes[rep_ab]++;
                count++;
            }
        }
    }
    else if (sect_h == 2 && sect_p == 0) {
        count = 0;
        for (int k = 0; k < nacth; k++) {
            for (int m = k + 1; m < nacth; m++) {
                moindex_t idx_k = active_holes_indices[k];
                moindex_t idx_m = active_holes_indices[m];
                det_buf[count].indices[0] = idx_k;
                det_buf[count].indices[1] = idx_m;

                int rep_k = spinor_info[idx_k].repno;
                int rep_m = spinor_info[idx_m].repno;
                int rep_km = mulrep2_abelian(rep_k, rep_m);
                det_buf[count].sym = rep_km;
                ms_rep_sizes[rep_km]++;
                count++;
            }
        }
    }
    else if (sect_h == 3 && sect_p == 0) {
        count = 0;
        for (int k = 0; k < nacth; k++) {
            for (int m = k + 1; m < nacth; m++) {
                for (int l = m + 1; l < nacth; l++) {
                    moindex_t idx_k = active_holes_indices[k];
                    moindex_t idx_m = active_holes_indices[m];
                    moindex_t idx_l = active_holes_indices[l];
                    det_buf[count].indices[0] = idx_k;
                    det_buf[count].indices[1] = idx_m;
                    det_buf[count].indices[2] = idx_l;

                    int rep_k = spinor_info[idx_k].repno;
                    int rep_m = spinor_info[idx_m].repno;
                    int rep_l = spinor_info[idx_l].repno;
                    int rep_km = mulrep2_abelian(rep_k, rep_m);
                    int rep_kml = mulrep2_abelian(rep_km, rep_l);
                    det_buf[count].sym = rep_kml;

                    ms_rep_sizes[rep_kml]++;
                    count++;
                }
            }
        }
    }
    else if (sect_h == 1 && sect_p == 1) {
        count = 0;
        // 1h1p determinants
        for (int i = 0; i < nacth; i++) {
            for (int a = 0; a < nactp; a++) {
                moindex_t idx_i = active_holes_indices[i];
                moindex_t idx_a = active_parts_indices[a];
                det_buf[count].indices[0] = idx_i;
                det_buf[count].indices[1] = idx_a;

                int rep_i = inverse_irrep_abelian(spinor_info[idx_i].repno);
                int rep_a = spinor_info[idx_a].repno;
                int rep_ia = mulrep2_abelian(rep_i, rep_a);
                det_buf[count].sym = rep_ia;
                ms_rep_sizes[rep_ia]++;
                count++;
            }
        }
    }
    else if (sect_h == 1 && sect_p == 2) {
        count = 0;
        for (int i = 0; i < nacth; i++) {
            for (int a = 0; a < nactp; a++) {
                for (int b = a + 1; b < nactp; b++) {
                    moindex_t idx_i = active_holes_indices[i];
                    moindex_t idx_a = active_parts_indices[a];
                    moindex_t idx_b = active_parts_indices[b];
                    det_buf[count].indices[0] = idx_i;
                    det_buf[count].indices[1] = idx_a;
                    det_buf[count].indices[2] = idx_b;

                    int rep_i = inverse_irrep_abelian(spinor_info[idx_i].repno);
                    int rep_a = spinor_info[idx_a].repno;
                    int rep_b = spinor_info[idx_b].repno;
                    int rep_ia = mulrep2_abelian(rep_i, rep_a);
                    int rep_iab = mulrep2_abelian(rep_ia, rep_b);
                    det_buf[count].sym = rep_iab;
                    ms_rep_sizes[rep_iab]++;
                    count++;
                }
            }
        }
    }
    else if (sect_h == 0 && sect_p == 3) {
        count = 0;
        for (int a = 0; a < nactp; a++) {
            for (int b = a + 1; b < nactp; b++) {
                for (int c = b + 1; c < nactp; c++) {
                    moindex_t idx_a = active_parts_indices[a];
                    moindex_t idx_b = active_parts_indices[b];
                    moindex_t idx_c = active_parts_indices[c];
                    det_buf[count].indices[0] = idx_a;
                    det_buf[count].indices[1] = idx_b;
                    det_buf[count].indices[2] = idx_c;

                    int rep_a = spinor_info[idx_a].repno;
                    int rep_b = spinor_info[idx_b].repno;
                    int rep_c = spinor_info[idx_c].repno;
                    int rep_ab = mulrep2_abelian(rep_a, rep_b);
                    int rep_abc = mulrep2_abelian(rep_ab, rep_c);
                    det_buf[count].sym = rep_abc;

                    ms_rep_sizes[rep_abc]++;
                    count++;
                }
            }
        }
    }

    // sort determinants by irrep (in ascending order)
    //qsort(det_buf, count, sizeof(slater_det_t), detcmp);

    // allocate memory
    slater_det_t **det_basis = (slater_det_t **) cc_malloc(sizeof(slater_det_t *) * get_num_irreps());
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (ms_rep_sizes[irrep] == 0) {
            det_basis[irrep] = NULL;
        }
        else {
            det_basis[irrep] = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * ms_rep_sizes[irrep]);
        }
    }

    // copy determinants from the buffer to the target arrays
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t count = 0;
        for (size_t idet = 0; idet < total_ms_size; idet++) {
            if (det_buf[idet].sym == irrep) {
                det_basis[irrep][count] = det_buf[idet];
                count++;
            }
        }
    }

    cc_free(det_buf);
    cc_free(active_holes_indices);
    cc_free(active_parts_indices);

    return det_basis;
}


/**
 * Prints information about model space determinants.
 */
void print_model_space_info(int sect_h, int sect_p, size_t *block_dims, slater_det_t **det_basis)
{
    if (cc_opts->print_model_space) {
        printf(" List of model space determinants:\n");
        for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
            for (int i = 0; i < block_dims[irrep]; i++) {
                print_slater_det(stdout, sect_h, sect_p, det_basis[irrep] + i);
            }
        }
        printf("\n");
    }

    printf(" Model space dimensions:\n");
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] != 0) {
            printf("  [%3d] dim = %-6d %-4s\n", irrep, block_dims[irrep], get_irrep_name(irrep));
        }
    }
    printf("\n");
}


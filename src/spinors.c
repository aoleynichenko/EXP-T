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
 * spinors.c
 * =========
 *
 * All information about one-particle functions (spin-orbitals or spinors)
 * (and their symmetry).
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "engine.h"
#include "spinors.h"
#include "symmetry.h"
#include "options.h"
#include "utils.h"

// number of one-particle functions (e.g. spinors)
int nspinors;

// array of spinor attributes
spinor_attr_t *spinor_info;

// list of spinor indices -- of each symmetry
size_t n_spinor_blocks;
spinor_block_t *spinor_blocks;

// quasiparticle type
spinor_types_t spinor_types;

spinor_block_t *spinor_blocks_h0;
spinor_block_t *spinor_blocks_h1;
spinor_block_t *spinor_blocks_p0;
spinor_block_t *spinor_blocks_p1;
spinor_block_t *spb_hp01[2][2];

// mapping "global spinor index -> local spinor index (in spinor block)"
int spinor_index_global2local[CC_MAX_SPINORS];


// info functions ("getters")
inline int is_active(int idx)
{
    return (spinor_info[idx].active == 1) ? 1 : 0;
    //return spinor_info[idx].active;
}


inline int is_hole(int idx)
{
    return spinor_info[idx].occ /*== 1*/;
}


inline int is_act_hole(int idx)
{
    return is_hole(idx) && is_active(idx);
}


inline int is_inact_hole(int idx)
{
    return is_hole(idx) && !is_active(idx);
}


inline int is_part(int idx)
{
    return spinor_info[idx].occ == 0;
}


inline int is_act_part(int idx)
{
    return is_part(idx) && is_active(idx);
}


inline int is_inact_part(int idx)
{
    return is_part(idx) && !is_active(idx);
}


// returns irrep to which the vacuum (0h0p) electronic state belongs
int get_vacuum_irrep()
{
    return get_totally_symmetric_irrep();
}


void get_active_space_size(int *nacth, int *nactp)
{
    *nacth = 0;
    *nactp = 0;

    for (int i = 0; i < nspinors; i++) {
        if (is_act_hole(i)) {
            (*nacth)++;
        }
        else if (is_act_part(i)) {
            (*nactp)++;
        }
    }
}


void get_active_space(int *nacth, int *nactp, moindex_t *active_holes_indices, moindex_t *active_parts_indices)
{
    *nacth = 0;
    *nactp = 0;

    for (int i = 0; i < nspinors; i++) {
        if (is_act_hole(i)) {
            active_holes_indices[*nacth] = i;
            (*nacth)++;
        }
        else if (is_act_part(i)) {
            active_parts_indices[*nactp] = i;
            (*nactp)++;
        }
    }
}


double get_eps(int idx)
{
    return spinor_info[idx].eps;
}


int get_num_electrons()
{
    int nelec = 0;

    for (size_t i = 0; i < nspinors; i++) {
        nelec += spinor_info[i].occ;
    }

    return nelec;
}


void get_spinor_info(int idx, int *repno, int *occ, int *active, double *eps)
{
    *repno = spinor_info[idx].repno;
    *occ = spinor_info[idx].occ;
    *active = spinor_info[idx].active;
    *eps = spinor_info[idx].eps;
}


// at the moment: only symmetry blocks are supported (no subblocks)
void create_spinor_blocks(int tilesize)
{
    int *sizes;
    int i, j, isp;
    int prt_lvl;

    prt_lvl = cc_opts->print_level;

    // count number of spinors in each irrep
    sizes = (int *) cc_malloc(nsym * sizeof(int));
    memset(sizes, 0, nsym * sizeof(int));
    for (i = 0; i < nspinors; i++) {
        sizes[spinor_info[i].repno]++;
    }

    // count number of nonzero blocks
    n_spinor_blocks = 0;
    for (i = 0; i < nsym; i++) {
        if (sizes[i] > 0) {
            for (j = 0; j < sizes[i]; j += tilesize) {
                n_spinor_blocks++;
            }
        }
    }
    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        printf("nsym = %d\n", nsym);
        printf("sizes = ");
        for (i = 0; i < nsym; i++) { printf("%d ", sizes[i]); }
        printf("\n");
        printf("number of nonzero blocks = %d\n", n_spinor_blocks);
    }

    // allocate memory
    spinor_blocks = (spinor_block_t *) cc_malloc(n_spinor_blocks * sizeof(spinor_block_t));
    int k;
    j = 0;
    for (i = 0; i < nsym; i++) {
        if (sizes[i] > 0) {
            for (k = 0; k < sizes[i]; k += tilesize) {
                int size = (k + tilesize > sizes[i]) ? (sizes[i] - k) : tilesize;
                if (cc_opts->print_level >= CC_PRINT_DEBUG) {
                    printf("i = %d j = %d repno = %d size = %d\n", i, j, i, size);
                }
                spinor_blocks[j].size = 0;
                spinor_blocks[j].repno = i;
                spinor_blocks[j].indices = (int *) cc_malloc(size * sizeof(int));
                j++;
            }
        }
    }

    // assign indices
    for (i = 0; i < nspinors; i++) {
        for (j = 0; j < n_spinor_blocks; j++) {
            if (spinor_info[i].repno == spinor_blocks[j].repno) {
                if (spinor_blocks[j].size == tilesize) {
                    continue;
                }
                spinor_blocks[j].indices[spinor_blocks[j].size++] = i;
                spinor_info[i].blockno = j;
                break;
            }
        }
    }

    int field_width = 7;
    for (int i = 0; i < n_spinor_blocks; i++) {
        int len = strlen(get_irrep_name(spinor_info[i].repno));
        if (len > field_width) {
            field_width = len;
        }
    }

    /* beautiful table with information about spinor blocks */
    printf(" Blocks of molecular spinors:\n");
    print_hyphens(1, 80);
    printf("   #  size%*s%s\n", field_width, "irrep", "  spinor indices");
    print_hyphens(1, 80);
    for (int i = 0; i < n_spinor_blocks; i++) {
        size_t sz = spinor_blocks[i].size;
        char *rep_name = get_irrep_name(spinor_blocks[i].repno);
        int *indices = spinor_blocks[i].indices;
        printf(" [%2d] %4d%*s  ", i, sz, field_width, rep_name);
        for (int j = 0; j < sz; j++) {
            if (j != 0) {
                printf(",");
            }
            printf("%d", indices[j] + 1);

            int nskip = 0;
            for (int k = j + 1; k < sz; k++) {
                if (indices[k] == indices[k - 1] + 1) {
                    nskip++;
                }
                else {
                    break;
                }
            }
            if (nskip > 0) {
                printf("-%d", indices[j + nskip] + 1);
                j += nskip;
            }
        }
        printf("\n");
    }
    print_hyphens(1, 80);
    printf("\n");

    // construct mapping: "global spinor index -> local spinor index"
    for (isp = 0; isp < nspinors; isp++) {
        for (i = 0; i < n_spinor_blocks; i++) {
            for (j = 0; j < spinor_blocks[i].size; j++) {
                if (spinor_blocks[i].indices[j] == isp) {
                    spinor_index_global2local[isp] = j;
                    goto contin;
                }
            }
        }
        contin:
        continue;
    }

    if (prt_lvl >= CC_PRINT_DEBUG) {
        printf("spinor indices, global to local mapping:\n");
        for (i = 0; i < nspinors; i++) {
            printf("%3d->%3d ", i, spinor_index_global2local[i]);
            if ((i + 1) % 10 == 0) {
                printf("\n");
            }
        }
        printf("\n");
    }

    // cleanup
    cc_free(sizes);

    // spinor blocks for difference combinations of flags:
    // (hole,general)
    // (particle,general)
    // (hole,active)
    // (particle,active)
}


int is_symblock_zero_rank2(int *spinor_blocks, int *qparts, int *valence)
{
    int hp, val;
    spinor_block_t *arr;

    hp = (qparts[0] == 'h') ? 0 : 1;
    val = valence[0] == 1;
    arr = spb_hp01[hp][val];
    int sz1 = arr[spinor_blocks[0]].size;
    if (sz1 == 0) { return 1; }

    hp = (qparts[1] == 'h') ? 0 : 1;
    val = valence[1] == 1;
    arr = spb_hp01[hp][val];
    int sz2 = arr[spinor_blocks[1]].size;
    if (sz2 == 0) { return 1; }

    return 0;
}


int is_symblock_zero_rank4(int *spinor_blocks, int *qparts, int *valence)
{
    int hp, val;
    spinor_block_t *arr;
    int sz;

    hp = (qparts[0] == 'h') ? 0 : 1;
    val = valence[0] == 1;
    arr = spb_hp01[hp][val];
    sz = arr[spinor_blocks[0]].size;
    if (sz == 0) { return 1; }

    hp = (qparts[1] == 'h') ? 0 : 1;
    val = valence[1] == 1;
    arr = spb_hp01[hp][val];
    sz = arr[spinor_blocks[1]].size;
    if (sz == 0) { return 1; }

    hp = (qparts[2] == 'h') ? 0 : 1;
    val = valence[2] == 1;
    arr = spb_hp01[hp][val];
    sz = arr[spinor_blocks[2]].size;
    if (sz == 0) { return 1; }

    hp = (qparts[3] == 'h') ? 0 : 1;
    val = valence[3] == 1;
    arr = spb_hp01[hp][val];
    sz = arr[spinor_blocks[3]].size;
    if (sz == 0) { return 1; }

    return 0;
}


int (*is_symblock_zerox[CC_DIAGRAM_MAX_RANK])(int *spinor_blocks, int *qparts, int *valence) = {
        is_symblock_zero_rank2, is_symblock_zero_rank4, NULL, NULL, NULL
};


int is_symblock_zero(int rank, int *spinor_blocks, int *qparts, int *valence)
{
    for (int i = 0; i < rank; i++) {
        int hp = (qparts[i] == 'h') ? 0 : 1;
        int val = (valence[i] == 1) ? 1 : 0;
        spinor_block_t *arr = spb_hp01[hp][val];
        int bsize = arr[spinor_blocks[i]].size;
        if (bsize == 0) {
            return 1;
        }
    }
    return 0;
}


int spinor_attr_cmp_by_eps(const void *i, const void *j)
{
    double eps1 = ((const spinor_attr_t *) i)->eps;
    double eps2 = ((const spinor_attr_t *) j)->eps;

    if (fabs(eps1 - eps2) < 1e-10) {
        return 0;
    }
    if (eps1 > eps2) {
        return 1;
    }
    else {
        return -1;
    }
}


/* classification of spinors: find all holes, active holes, active particles,
   particles.
   these data will be stored in the structure spinor_types (data type
   spinor_types_t).
   [from spinor.h]
   typedef struct {
       int *holes;      // including active holes
       int *inact_holes;
       int *act_holes;
       int *act_parts;
       int *inact_parts;
       int *parts;      // including active particles
   } spinor_types_t;
*/
void classify_spinors(double actsp_min, double actsp_max,
                      int nacth, int nactp)
{
    int *h, *ih, *ah, *p, *ip, *ap;
    int nh, nih, nah, np, nip, nap;
    size_t i;
    int occ;
    double eps;
    int act;
    int nocc;
    spinor_attr_t *sorted_spinor_info;

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        printf("\nspinor classification (define active space)...\n");
        printf("active space MIN energy = %g a.u.\n", actsp_min);
        printf("active space MAX energy = %g a.u.\n", actsp_max);
        printf("I/A = inactive/active\n");
        printf("H/P = hole/particle\n");
    }

    // set all counters to zero
    nh = 0;
    nih = 0;
    nah = 0;
    np = 0;
    nip = 0;
    nap = 0;

    // alloc tmp arrays (max possible size == nspinors)
    h = (int *) cc_malloc(nspinors * sizeof(int));
    ih = (int *) cc_malloc(nspinors * sizeof(int));
    ah = (int *) cc_malloc(nspinors * sizeof(int));
    p = (int *) cc_malloc(nspinors * sizeof(int));
    ip = (int *) cc_malloc(nspinors * sizeof(int));
    ap = (int *) cc_malloc(nspinors * sizeof(int));
    sorted_spinor_info = (spinor_attr_t *) cc_malloc(nspinors * sizeof(spinor_attr_t));

    if (!cc_opts->occ_defined && cc_opts->nelec_defined) {
        int nelec[CC_MAX_NUM_IRREPS];
        for (i = 0; i < CC_MAX_NUM_IRREPS; i++) {
            nelec[i] = cc_opts->nelec[i];
        }
        for (i = 0; i < nspinors; i++) {
            int repno = spinor_info[i].repno;
            if (nelec[repno] > 0) {
                spinor_info[i].occ = 1;
                nelec[repno]--;
            }
            else {
                spinor_info[i].occ = 0;
            }
        }
    }
    if (cc_opts->occ_defined) {
        for (i = 0; i < nspinors; i++) {
            spinor_info[i].occ = cc_opts->occ[i];
        }
    }

    // IF 'nactp/nacth first' DIRECTIVE IS SPECIFIED:
    // sort spinors by energy (for subsequent extraction of the active space)
    // and find range of "active" energies
    // (REDEFINE active_eps_min and active_eps_max !)
    if (nactp != 0 || nacth != 0) {
        memcpy(sorted_spinor_info, spinor_info, sizeof(spinor_attr_t) * nspinors);
        qsort(sorted_spinor_info, nspinors, sizeof(spinor_attr_t), spinor_attr_cmp_by_eps);
        nocc = 0;
        for (i = 0; i < nspinors; i++) {
            if (sorted_spinor_info[i].occ == 1) {
                nocc++;
            }
        }
        actsp_min = sorted_spinor_info[nocc - nacth].eps - 1e-7;
        actsp_max = sorted_spinor_info[nocc + nactp - 1].eps + 1e-7;
    }

    if (cc_opts->actsp_defined == 1 || cc_opts->actsp_defined == 2) { // by energy or by total number: nactp/nacth
        for (i = 0; i < nspinors; i++) {
            occ = spinor_info[i].occ;
            eps = spinor_info[i].eps;
            if (occ == 1) { // hole
                act = (eps >= actsp_min) ? 1 : 0;
            }
            else { // particle
                act = (eps <= actsp_max) ? 1 : 0;
            }
            spinor_info[i].active = act;
        }
    }
    else if (cc_opts->actsp_defined == 3) { // by irrep
        // here we assume that all spinors inside the irrep are sorted by energy
        // this is true for DIRAC
        for (size_t i = 0; i < CC_MAX_NREP; i++) {
            cc_active_spec_t *sp = cc_opts->active_specs + i;
            if (*(sp->rep_name) == '\0') {
                continue;
            }
            int irep = get_rep_number(sp->rep_name);
            if (irep == -1) {
                errquit("Error in the input file: wrong irrep name in active space specification: '%s'",
                        sp->rep_name);
            }
            // count number of spinors and occupied spinors in this irrep
            int n_spi = 0;
            int n_occ = 0;
            for (size_t i = 0; i < nspinors; i++) {
                if (spinor_info[i].repno == irep) {
                    n_spi++;
                    if (spinor_info[i].occ == 1) {
                        n_occ++;
                    }
                }
            }
            if (n_spi == 0) {
                printf(" Warning: there are no spinors belonging to the irrep '%s'\n", sp->rep_name);
                continue;
            }
            int n_inact_h = n_occ - sp->nacth;
            int count = 0;
            for (size_t ispinor = 0; ispinor < nspinors; ispinor++) {
                spinor_attr_t *info = &spinor_info[ispinor];
                if (info->repno != irep) {
                    continue;
                }
                if (count < n_inact_h) {
                    count++;
                    continue;
                }
                if (count < n_occ) {
                    info->active = 1;
                    count++;
                    continue;
                }
                if (count < n_occ + sp->nactp) {
                    info->active = 1;
                    count++;
                    continue;
                }
            }
        }
    }
    else if (cc_opts->actsp_defined == 4) {
        for (i = 0; i < nspinors; i++) {
            spinor_info[i].active = cc_opts->active_each[i];
        }
    }

    // count spinors of each type
    for (i = 0; i < nspinors; i++) {
        occ = spinor_info[i].occ;
        act = spinor_info[i].active;

        // hole
        if (occ == 1) {
            h[nh++] = i;
            if (act == 1) {
                ah[nah++] = i;
            }
            else {
                ih[nih++] = i;
            }
        }
            // particle
        else {
            p[np++] = i;
            if (act == 1) {
                ap[nap++] = i;
            }
            else {
                ip[nip++] = i;
            }
        }
        spinor_info[i].active = act;
    }

    // total count
    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        printf("NH = NIH + NAH = %d + %d = %d\n", nih, nah, nh);
        printf("NP = NIP + NAP = %d + %d = %d\n", nip, nap, np);
    }

    // fill struct spinor_types (global variable, see spinor.h)
    // alloc memory Q.S.
    spinor_types.holes = (int *) cc_malloc((nh + 1) * sizeof(int));
    spinor_types.inact_holes = (int *) cc_malloc((nih + 1) * sizeof(int));
    spinor_types.act_holes = (int *) cc_malloc((nah + 1) * sizeof(int));
    spinor_types.parts = (int *) cc_malloc((np + 1) * sizeof(int));
    spinor_types.inact_parts = (int *) cc_malloc((nip + 1) * sizeof(int));
    spinor_types.act_parts = (int *) cc_malloc((nap + 1) * sizeof(int));

    // first element [0] = size
    spinor_types.holes[0] = nh;
    spinor_types.inact_holes[0] = nih;
    spinor_types.act_holes[0] = nah;
    spinor_types.parts[0] = np;
    spinor_types.inact_parts[0] = nip;
    spinor_types.act_parts[0] = nap;

    // copy spinor indices from tmp arrays
    memcpy(spinor_types.holes + 1, h, nh * sizeof(int));
    memcpy(spinor_types.inact_holes + 1, ih, nih * sizeof(int));
    memcpy(spinor_types.act_holes + 1, ah, nah * sizeof(int));
    memcpy(spinor_types.parts + 1, p, np * sizeof(int));
    memcpy(spinor_types.inact_parts + 1, ip, nip * sizeof(int));
    memcpy(spinor_types.act_parts + 1, ap, nap * sizeof(int));

    // cleanup tmp arrays
    cc_free(h);
    cc_free(ih);
    cc_free(ah);
    cc_free(p);
    cc_free(ip);
    cc_free(ap);
    cc_free(sorted_spinor_info);

    spinor_blocks_h0 = (spinor_block_t *) cc_malloc(sizeof(spinor_block_t) * n_spinor_blocks);
    spinor_blocks_h1 = (spinor_block_t *) cc_malloc(sizeof(spinor_block_t) * n_spinor_blocks);
    spinor_blocks_p0 = (spinor_block_t *) cc_malloc(sizeof(spinor_block_t) * n_spinor_blocks);
    spinor_blocks_p1 = (spinor_block_t *) cc_malloc(sizeof(spinor_block_t) * n_spinor_blocks);

    for (size_t i = 0; i < n_spinor_blocks; i++) {
        spinor_blocks_h0[i].repno = spinor_blocks[i].repno;
        spinor_blocks_h0[i].size = 0;
        spinor_blocks_h0[i].indices = (int *) cc_malloc(spinor_blocks[i].size * sizeof(int));
        for (int j = 0; j < spinor_blocks[i].size; j++) {
            if (is_hole(spinor_blocks[i].indices[j])) {
                spinor_blocks_h0[i].indices[spinor_blocks_h0[i].size++] = spinor_blocks[i].indices[j];
            }
        }
    }
    for (size_t i = 0; i < n_spinor_blocks; i++) {
        spinor_blocks_p0[i].repno = spinor_blocks[i].repno;
        spinor_blocks_p0[i].size = 0;
        spinor_blocks_p0[i].indices = (int *) cc_malloc(spinor_blocks[i].size * sizeof(int));
        for (int j = 0; j < spinor_blocks[i].size; j++) {
            if (is_part(spinor_blocks[i].indices[j])) {
                spinor_blocks_p0[i].indices[spinor_blocks_p0[i].size++] = spinor_blocks[i].indices[j];
            }
        }
    }
    for (size_t i = 0; i < n_spinor_blocks; i++) {
        spinor_blocks_h1[i].repno = spinor_blocks[i].repno;
        spinor_blocks_h1[i].size = 0;
        spinor_blocks_h1[i].indices = (int *) cc_malloc(spinor_blocks[i].size * sizeof(int));
        for (int j = 0; j < spinor_blocks[i].size; j++) {
            if (is_act_hole(spinor_blocks[i].indices[j])) {
                spinor_blocks_h1[i].indices[spinor_blocks_h1[i].size++] = spinor_blocks[i].indices[j];
            }
        }
    }
    for (size_t i = 0; i < n_spinor_blocks; i++) {
        spinor_blocks_p1[i].repno = spinor_blocks[i].repno;
        spinor_blocks_p1[i].size = 0;
        spinor_blocks_p1[i].indices = (int *) cc_malloc(spinor_blocks[i].size * sizeof(int));
        for (int j = 0; j < spinor_blocks[i].size; j++) {
            if (is_act_part(spinor_blocks[i].indices[j])) {
                spinor_blocks_p1[i].indices[spinor_blocks_p1[i].size++] = spinor_blocks[i].indices[j];
            }
        }
    }

    spb_hp01[0][0] = spinor_blocks_h0;
    spb_hp01[0][1] = spinor_blocks_h1;
    spb_hp01[1][0] = spinor_blocks_p0;
    spb_hp01[1][1] = spinor_blocks_p1;
}


void print_spinor_info()
{
    int i;
    int occ_i[256], occ_a[256], virt_i[256], virt_a[256];
    int occ, act, rep;
    int ccs_only_space[CC_MAX_SPINORS];

    // set counters to zero
    for (i = 0; i < nsym; i++) {
        occ_i[i] = 0;
        occ_a[i] = 0;
        virt_a[i] = 0;
        virt_i[i] = 0;
    }

    // calculate CCSD only subspace of spinors
    setup_singles_only_space(ccs_only_space);

    printf("\n");
    printf("\t\t\t\tSpinors info\n");
    printf("\t\t\t\t------------\n");
    printf("\n");
    printf("     no    rep         occ    active        one-el energy\n");
    printf("    -------------------------------------------------------\n");
    for (i = 0; i < nspinors; i++) {
        printf("    %4d%4d \"%-6s\"%4d       %1s     %20.12f  %s\n",
               spinor_info[i].seqno + 1, spinor_info[i].repno,
               rep_names[spinor_info[i].repno], spinor_info[i].occ,
               (spinor_info[i].active == 1) ? "a" : "-", spinor_info[i].eps,
               ccs_only_space[i] ? "ccs only" : "");
        // count some statistics
        rep = spinor_info[i].repno;
        occ = spinor_info[i].occ;
        act = spinor_info[i].active;
        if (occ == 1 && act == 0) {
            occ_i[rep]++;
        }
        else if (occ == 1 && act == 1) {
            occ_a[rep]++;
        }
        else if (occ == 0 && act == 0) {
            virt_i[rep]++;
        }
        else if (occ == 0 && act == 1) {
            virt_a[rep]++;
        }
    }
    printf("    -------------------------------------------------------\n\n");

    // how many reps have non-zero number of spinors
    int to_be_printed[64];
    for (int irep = 0; irep < nsym; irep++) {
        if (occ_i[irep] == 0 && occ_a[irep] == 0 && virt_a[irep] == 0 && virt_i[irep] == 0) {
            to_be_printed[irep] = 0;
        }
        else {
            to_be_printed[irep] = 1;
        }
    }

    // print little statistics
    printf("    irreps            ");
    for (i = 0; i < nsym; i++) {
        if (to_be_printed[i]) { printf("%6s", rep_names[i]); }
    }
    printf("\n");
    printf("    occupied inactive ");
    for (i = 0; i < nsym; i++) {
        if (to_be_printed[i]) { printf("%6d", occ_i[i]); }
    }
    printf("\n");
    printf("    occupied active   ");
    for (i = 0; i < nsym; i++) {
        if (to_be_printed[i]) { printf("%6d", occ_a[i]); }
    }
    printf("\n");
    printf("    virtual active    ");
    for (i = 0; i < nsym; i++) {
        if (to_be_printed[i]) { printf("%6d", virt_a[i]); }
    }
    printf("\n");
    printf("    virtual inactive  ");
    for (i = 0; i < nsym; i++) {
        if (to_be_printed[i]) { printf("%6d", virt_i[i]); }
    }
    printf("\n\n");
}


int get_spinor_block_size(int spinor_block_number)
{
    assert(spinor_block_number >= 0 && spinor_block_number < n_spinor_blocks);

    return spinor_blocks[spinor_block_number].size;
}


int get_max_spinor_block_size()
{
    int max_sz = 0;

    for (size_t i = 0; i < n_spinor_blocks; i++) {
        if (spinor_blocks[i].size > max_sz) {
            max_sz = spinor_blocks[i].size;
        }
    }

    return max_sz;
}


void spinors_cleanup()
{
    size_t i;

    cc_free(spinor_info);

    for (i = 0; i < n_spinor_blocks; i++) {
        cc_free(spinor_blocks[i].indices);
    }
    cc_free(spinor_blocks);

    cc_free(spinor_types.holes);
    cc_free(spinor_types.inact_holes);
    cc_free(spinor_types.act_holes);
    cc_free(spinor_types.parts);
    cc_free(spinor_types.inact_parts);
    cc_free(spinor_types.act_parts);
}


/**
 * @param ccs_only_space pre-allocated array (dim=nspinors);
 * on exit, contains 0 or 1 for each spinor:
 *   0 -- to be used for relaxation only
 *   1 -- to be used for both relaxation and correlation
 * @return ccs_only_space
 * @note here we assume that all spinors inside the irrep are sorted by energy;
 * this is true for DIRAC
 */
void setup_singles_only_space(int *ccs_only_space)
{
    double core_emin = 0.0, core_emax = 0.0;
    double virt_emin = 0.0, virt_emax = 0.0;

    memset(ccs_only_space, 0, sizeof(int) * nspinors);

    if (!cc_opts->do_relax) {
        return;
    }

    cc_space_t *space_core = &cc_opts->relax_core;
    cc_space_t *space_virt = &cc_opts->relax_virtual;

    /* setup core space */

    if (space_core->deftype == CC_SPACE_BY_TOTAL) {
        spinor_attr_t *sorted_spinor_info = (spinor_attr_t *) cc_malloc(nspinors * sizeof(spinor_attr_t));
        memcpy(sorted_spinor_info, spinor_info, sizeof(spinor_attr_t) * nspinors);
        qsort(sorted_spinor_info, nspinors, sizeof(spinor_attr_t), spinor_attr_cmp_by_eps);
        if (space_core->total > nspinors) {
            space_core->total = nspinors;
        }
        core_emin = sorted_spinor_info[0].eps - 1e-7;
        core_emax = sorted_spinor_info[space_core->total - 1].eps + 1e-7;
    }
    if (space_core->deftype == CC_SPACE_BY_ENERGY) {
        core_emin = space_core->emin;
        core_emax = space_core->emax;
    }
    if (space_core->deftype == CC_SPACE_BY_ENERGY ||
        space_core->deftype == CC_SPACE_BY_TOTAL) {
        for (int i = 0; i < nspinors; i++) {
            if (spinor_info[i].eps >= core_emin && spinor_info[i].eps <= core_emax) {
                ccs_only_space[i] = 1;
            }
        }
    }

    if (space_core->deftype == CC_SPACE_BY_IRREPS) {
        for (size_t i = 0; i < CC_MAX_NREP; i++) {
            if (*(space_core->rep_names[i]) == '\0') {
                continue;
            }
            int irep = get_rep_number(space_core->rep_names[i]);
            if (irep == -1) { // no such irrep
                continue;
            }
            // mark first N spinors in this irrep
            int n_spi = 0;
            for (int j = 0; j < nspinors; j++) {
                if (spinor_info[j].repno == irep) {
                    n_spi++;
                    if (n_spi <= space_core->dim[i]) {
                        ccs_only_space[j] = 1;
                    }
                }
            }
        }
    }

    /* setup core virtual space */

    if (space_virt->deftype == CC_SPACE_BY_TOTAL) {
        spinor_attr_t *sorted_spinor_info = (spinor_attr_t *) cc_malloc(nspinors * sizeof(spinor_attr_t));
        memcpy(sorted_spinor_info, spinor_info, sizeof(spinor_attr_t) * nspinors);
        qsort(sorted_spinor_info, nspinors, sizeof(spinor_attr_t), spinor_attr_cmp_by_eps);
        if (space_virt->total > nspinors) {
            space_virt->total = nspinors;
        }
        virt_emin = sorted_spinor_info[nspinors - space_virt->total].eps - 1e-7;
        virt_emax = sorted_spinor_info[nspinors - 1].eps + 1e-7;
    }
    if (space_virt->deftype == CC_SPACE_BY_ENERGY) {
        virt_emin = space_virt->emin;
        virt_emax = space_virt->emax;
    }
    if (space_virt->deftype == CC_SPACE_BY_ENERGY ||
        space_virt->deftype == CC_SPACE_BY_TOTAL) {
        for (int i = 0; i < nspinors; i++) {
            if (spinor_info[i].eps >= virt_emin && spinor_info[i].eps <= virt_emax) {
                ccs_only_space[i] = 1;
            }
        }
    }

    if (space_virt->deftype == CC_SPACE_BY_IRREPS) {
        for (size_t i = 0; i < CC_MAX_NREP; i++) {
            if (*(space_virt->rep_names[i]) == '\0') {
                continue;
            }
            int irep = get_rep_number(space_virt->rep_names[i]);
            if (irep == -1) { // no such irrep
                continue;
            }
            // mark last N spinors in this irrep
            int n_spi = 0;
            for (int j = nspinors - 1; j >= 0; j--) {
                if (spinor_info[j].repno == irep) {
                    n_spi++;
                    if (n_spi <= space_virt->dim[i]) {
                        ccs_only_space[j] = 1;
                    }
                }
            }
        }
    }
}
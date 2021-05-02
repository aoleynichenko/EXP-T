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
 * heff.c
 * ======
 *
 * (1) construction and diagonalization of the effective Hamiltonian matrix;
 * (2) analysis of its eigenvalues and eigenvectors (model vectors);
 * (3) write matrix and vectors to disk.
 * NOTE: ALL THIS CODE IS DESIGNED FOR ABELIAN GROUPS ONLY!
 *
 * 2019-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "heff.h"

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codata.h"
#include "io.h"
#include "engine.h"
#include "datamodel.h"
#include "linalg.h"
#include "mvcoef.h"
#include "natorb.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"

void zero_order_heff(int sect_h, int sect_p, size_t dim, slater_det_t *det_list, double complex *heff);
void omega_0h0p_0h1p(size_t dim, slater_det_t *det_list_1h1p, double complex *omega,
                     char *exc_oper_name, char *deexc_oper_name);
FILE *hefff_open(int sect_h, int sect_p, char *label);
void hefff_close(FILE *hefff);
void hefff_write_block(FILE *hefff, int carith, int rep_no, size_t dim, double complex *heff);
void
print_model_vector(FILE *f_out, int sect_h, int sect_p, int rep_no, char *rep_name, int state_no, double complex energy,
                   size_t len, double complex *coeffs, slater_det_t *det_list, double coef_thresh);
int get_nroots_for_irrep(char *irrep_name);
int get_intham_main_subspace_for_irrep(char *irrep_name);
void renormalize_wave_operator_0h0p_0h1p(size_t dim, slater_det_t *det_list,
                                         double complex *heff, double complex *heff_prime, double complex *omega);

void construct_quasi_natural_orbitals(int sect_h, int sect_p, int ms_size, slater_det_t *det_list,
                                      double complex *coef_left, double complex *coef_right,
                                      double complex *natorb_left, double complex *natorb_right,
                                      double *nat_occ_numbers, double *config);


/*******************************************************************************
 * get_model_space_size
 *
 * Returns total number of model-space determinants for the given FS sector
 ******************************************************************************/
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


/*******************************************************************************
 * detcmp
 *
 * Comparator for Slater determinants; the lower the number of irrep =>
 * the older the determinant
 ******************************************************************************/
int detcmp(const void *_d1, const void *_d2)
{
    const slater_det_t *d1 = (const slater_det_t *) _d1;
    const slater_det_t *d2 = (const slater_det_t *) _d2;
    if (d1->sym != d2->sym) {
        return abs(d1->sym) - abs(d2->sym);
    }
    else {
        return d2->sym - d1->sym;
    }
}


/*******************************************************************************
 * create_model_dets
 *
 * Returns list of model-space Slater determinants and number of determinants
 * belonging to each symmetry (via the 'ms_rep_sizes' array).
 * NOTE: detlist[] array MUST be pre-allocated
 *
 * global variables used: NSPINORS, nsym
 ******************************************************************************/
void create_model_dets(int sect_h, int sect_p, size_t *ms_rep_sizes, slater_det_t *detlist)
{
    moindex_t *active_holes_indices;
    moindex_t *active_parts_indices;
    int nacth, nactp;
    size_t ndet;

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
        return;
    }
    else if (sect_h == 1 && sect_p == 0) {
        ndet = 0;
        for (int k = 0; k < nacth; k++) {
            moindex_t idx_k = active_holes_indices[k];
            int rep_k = spinor_info[idx_k].repno;
            detlist[ndet].indices[0] = idx_k;
            detlist[ndet].sym = rep_k;
            ms_rep_sizes[rep_k]++;
            ndet++;
        }
    }
    else if (sect_h == 0 && sect_p == 1) {
        ndet = 0;
        for (int a = 0; a < nactp; a++) {
            moindex_t idx_a = active_parts_indices[a];
            int rep_a = spinor_info[idx_a].repno;
            detlist[ndet].indices[0] = idx_a;
            detlist[ndet].sym = rep_a;
            ms_rep_sizes[rep_a]++;
            ndet++;
        }
    }
    else if (sect_h == 0 && sect_p == 2) {
        ndet = 0;
        for (int a = 0; a < nactp; a++) {
            for (int b = a + 1; b < nactp; b++) {
                moindex_t idx_a = active_parts_indices[a];
                moindex_t idx_b = active_parts_indices[b];
                detlist[ndet].indices[0] = idx_a;
                detlist[ndet].indices[1] = idx_b;

                int rep_a = spinor_info[idx_a].repno;
                int rep_b = spinor_info[idx_b].repno;
                int rep_ab = mulrep2_abelian(rep_a, rep_b);
                detlist[ndet].sym = rep_ab;
                ms_rep_sizes[rep_ab]++;
                ndet++;
            }
        }
    }
    else if (sect_h == 2 && sect_p == 0) {
        ndet = 0;
        for (int k = 0; k < nacth; k++) {
            for (int m = k + 1; m < nacth; m++) {
                moindex_t idx_k = active_holes_indices[k];
                moindex_t idx_m = active_holes_indices[m];
                detlist[ndet].indices[0] = idx_k;
                detlist[ndet].indices[1] = idx_m;

                int rep_k = spinor_info[idx_k].repno;
                int rep_m = spinor_info[idx_m].repno;
                int rep_km = mulrep2_abelian(rep_k, rep_m);
                detlist[ndet].sym = rep_km;
                ms_rep_sizes[rep_km]++;
                ndet++;
            }
        }
    }
    else if (sect_h == 1 && sect_p == 1) {
        ndet = 0;
        // add the 0h0p vacuum determinant -- if required
        if (cc_opts->mixed) {
            int rep_00 = get_totally_symmetric_irrep(); // is valid only for closed-shell references!
            detlist[0].indices[0] = 0;
            detlist[0].indices[1] = 0;
            detlist[0].sym = rep_00;
            ms_rep_sizes[rep_00]++;
            ndet++;
        }
        // 1h1p determinants
        for (int i = 0; i < nacth; i++) {
            for (int a = 0; a < nactp; a++) {
                moindex_t idx_i = active_holes_indices[i];
                moindex_t idx_a = active_parts_indices[a];
                detlist[ndet].indices[0] = idx_i;
                detlist[ndet].indices[1] = idx_a;

                int rep_i = inverse_irrep_abelian(spinor_info[idx_i].repno);
                //if (rep_i < 4) rep_i += 4;
                //else rep_i -= 4;
                int rep_a = spinor_info[idx_a].repno;
                int rep_ia = mulrep2_abelian(rep_i, rep_a);
                detlist[ndet].sym = rep_ia;
                ms_rep_sizes[rep_ia]++;
                ndet++;
            }
        }
    }
    else if (sect_h == 0 && sect_p == 3) {
        ndet = 0;
        for (int a = 0; a < nactp; a++) {
            for (int b = a + 1; b < nactp; b++) {
                for (int c = b + 1; c < nactp; c++) {
                    moindex_t idx_a = active_parts_indices[a];
                    moindex_t idx_b = active_parts_indices[b];
                    moindex_t idx_c = active_parts_indices[c];
                    detlist[ndet].indices[0] = idx_a;
                    detlist[ndet].indices[1] = idx_b;
                    detlist[ndet].indices[2] = idx_c;

                    int rep_a = spinor_info[idx_a].repno;
                    int rep_b = spinor_info[idx_b].repno;
                    int rep_c = spinor_info[idx_c].repno;
                    int rep_ab = mulrep2_abelian(rep_a, rep_b);
                    int rep_abc = mulrep2_abelian(rep_ab, rep_c);
                    detlist[ndet].sym = rep_abc;

                    ms_rep_sizes[rep_abc]++;
                    ndet++;
                }
            }
        }
    }

    // sort determinants by irrep (in ascending order)
    qsort(detlist, ndet, sizeof(slater_det_t), detcmp);

    cc_free(active_holes_indices);
    cc_free(active_parts_indices);
}


/*******************************************************************************
 * print_slater_det
 *
 * Prints Slater determinant 'det' to the stream 'f'.
 ******************************************************************************/
void print_slater_det(FILE *f, int sect_h, int sect_p, slater_det_t *det)
{
    // здесь нужно сделать исключение -- печать детерминанта из сектора 00, который у нас считается вакуумным
    // детектировать его очень просто -- det->indices[0] == 0 && det->indices[1] == 0
    // нужна отдельная функция, говорящая нам, это вакуумный детерминант или нет
    // и еще отдельная функция, создающая вакуумный детерминант
    if (is_vacuum_det(det)) {
        fprintf(f, "| vac > (%4s)\n", rep_names[det->sym]);
        return;
    }

    fprintf(f, "|");
    for (int ih = 0; ih < sect_h; ih++) {
        int idx = det->indices[ih];
        int rep = spinor_info[idx].repno;
        fprintf(f, " %4s #%4d (%12.6f)", rep_names[rep], idx + 1, spinor_info[idx].eps);
    }
    if (sect_h != 0 && sect_p != 0) {
        fprintf(f, " -> ");
    }
    for (int ip = 0; ip < sect_p; ip++) {
        int idx = det->indices[sect_h + ip];
        int rep = spinor_info[idx].repno;
        fprintf(f, " %4s #%4d (%12.6f)", rep_names[rep], idx + 1, spinor_info[idx].eps);
    }
    fprintf(f, " > (%4s)\n", rep_names[det->sym]);
}


void model_vectors_analysis(int sect_h, int sect_p)
{
    int nrep = 0;
    struct mv_block mv_blocks[64];
    double eff_config[CC_MAX_SPINORS];
    int n_active = 0;
    int active_spinors[CC_MAX_SPINORS];

    get_active_space(sect_h, sect_p, &n_active, active_spinors);

    read_model_vectors_unformatted(sect_h, sect_p, NULL, &nrep, mv_blocks);

    for (size_t irep = 0; irep < nrep; irep++) {
        struct mv_block *block = mv_blocks + irep;

        char *irrep_name = block->rep_name;
        int nroots = block->nroots;
        size_t ms_size = block->ms_size;
        slater_det_t *det_list = block->dets;
        double *energy_cm = block->energy_cm;
        double complex *eigenvals = block->eigval;
        double complex *coef_left = block->vl;
        double complex *coef_right = block->vr;

        int rep0 = 0;

        for (size_t i = 0; i < nroots; i++) {
            print_model_vector(stdout, sect_h, sect_p,
                               irep - rep0 + 1, irrep_name, i, eigenvals[i],
                               ms_size, coef_right + ms_size * i, det_list, 0);

            get_effective_configuration(sect_h, sect_p, ms_size, det_list,
                                        coef_left, coef_right, eff_config);

            printf(" Effective configuration:\n");
            for (size_t j = 0; j < n_active; j++) {
                printf("  %12.6f", eff_config[j]);
                int ispinor = active_spinors[j];
                printf("  %4s #%4d (%12.6f)\n", rep_names[spinor_info[ispinor].repno], ispinor + 1,
                       spinor_info[ispinor].eps);
            }
        }
    }
}


void model_space_props_natorbs(int sect_h, int sect_p)
{
    // calculation of density matrices and natural orbitals
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_denmat) {
        for (int ipair = 0; ipair < cc_opts->n_denmat; ipair++) {
            int sect2_h = cc_opts->denmat_query[ipair].sect2[0];
            int sect2_p = cc_opts->denmat_query[ipair].sect2[1];
            if (!(sect2_h == sect_h && sect2_p == sect_p)) {
                continue;
            }
            int rep1 = get_rep_number(cc_opts->denmat_query[ipair].rep1_name);
            int rep2 = get_rep_number(cc_opts->denmat_query[ipair].rep2_name);
            int state1 = cc_opts->denmat_query[ipair].state1;
            int state2 = cc_opts->denmat_query[ipair].state2;
            // single code for NOs and NTOs:
            // rep1 == rep2 && state1 == state2  =>  single state (NO)
            // else: pair of states (NTO)
            density_matrix(sect_h, sect_p, rep1, state1, rep2, state2);
            if (rep1 == rep2 && state1 == state2) {
                quasi_natural_orbitals_driver(sect_h, sect_p, rep1, state1);
            }
        }
    }

    // transition dipole moments via the DL-TDM techniques
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->do_diplen_tdm) {
        dipole_length_tdms(sect_h, sect_p);
    }

    // estimate properties: only for the target sector
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_ms_props > 0) {
        for (int i = 0; i < cc_opts->n_ms_props; i++) {
            model_space_property(cc_opts->prop_queries + i);
        }
    }
}

void construct_heff(int sect_h, int sect_p, slater_det_t *det_list, size_t *block_sizes, double complex **heff_blocks,
                    int n_diagrams, char **diagram_names);

void diagonalize_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff,
                      double complex **eigvalues, double complex **coef_left, double complex **coef_right);

void print_eigenvalues_table(int sect_h, int sect_p, size_t *block_dims, double complex **eigvalues, double degen_thresh);

double max_energy_of_required_roots(size_t *block_dims, double complex **eigvalues);

void write_formatted_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff);

void restore_intermediate_normalization(size_t total_ms_dim, slater_det_t *det_list, size_t *block_dims,
                                        double complex **heff, double complex **eigvalues);

void print_model_vectors_stdout(int sector_h, int sector_p,
                                size_t total_ms_size, slater_det_t *det_list, size_t *block_dims,
                                double complex **eigvalues, double complex **coef_left, double complex **coef_right);

void write_model_vectors_unformatted(int sector_h, int sector_p,
                                     size_t total_ms_size, slater_det_t *det_list, size_t *block_dims,
                                     double complex **eigvalues, double complex **coef_left, double complex **coef_right);

double get_lowest_eigenvalue(size_t *block_dims, double complex **eigvalues);

size_t first_nonzero_irrep(size_t *block_dims);


/**
 * heff_analysis
 *
 * Constructs and diagonalizes effective Hamiltonian.
 * Performs analysis of model vectors.
 * Writes to disk:
 *   HEFF      formatted file with the effective Hamiltonian
 *   MVCOEF**  unformatted files with model vectors and eigenvalues
 *
 * Agruments:
 *   sect_h, sect_p   number of the Fock-space sector
 *   ...              list of names of the diagrams containing 1-particle,
 *                    2-particle etc parts of the effective interaction operator
 */
void heff_analysis(int sect_h, int sect_p, ...)
{
    double complex *heff[CC_MAX_NUM_IRREPS];
    double complex *eigvalues[CC_MAX_NUM_IRREPS];
    double complex *coef_left[CC_MAX_NUM_IRREPS];
    double complex *coef_right[CC_MAX_NUM_IRREPS];
    size_t block_dims[CC_MAX_NUM_IRREPS];

    printf("\n Effective Hamiltonian analysis\n");

    /*
     * Construct model space
     */
    int nacth, nactp;
    int nspinors = get_num_spinors();
    get_active_space_size(&nacth, &nactp);
    size_t ms_size = get_model_space_size(sect_h, sect_p, nspinors, spinor_info);
    slater_det_t *det_list = (slater_det_t *) cc_malloc(ms_size * sizeof(slater_det_t));
    create_model_dets(sect_h, sect_p, block_dims, det_list);

    printf(" Active space: %d holes, %d particles\n", nacth, nactp);
    printf(" Model space size: %d determinants (total)\n", ms_size);
    if (cc_opts->print_level >= CC_PRINT_HIGH) {
        printf("\n List of model space determinants:\n");
        for (size_t i = 0; i < ms_size; i++) {
            printf(" %3d  ", i);
            print_slater_det(stdout, sect_h, sect_p, det_list+i);
        }
    }
    printf(" Dimensions of symmetry blocks of Heff:\n");
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_dims[irep] == 0) {
            continue;
        }
        printf("  %4s [%d]", rep_names[irep], block_dims[irep]);
    }
    size_t max_heff_size = size_t_max(get_num_irreps(), block_dims);
    printf(" (max %d)", max_heff_size);
    printf("\n");

    va_list args;
    va_start(args, sect_p);
    int n_diagrams = (sect_h + 1) * (sect_p + 1) - 1; // rectangle - (0,0) sector
    char ** diagram_names = (char **) cc_malloc(sizeof(char *) * n_diagrams);
    for (int i = 0; i < n_diagrams; i++) {
        char *dg_name = va_arg(args, char *);
        diagram_names[i] = cc_strdup(dg_name);
    }
    va_end(args);

    construct_heff(sect_h, sect_p, det_list, block_dims, heff, n_diagrams, diagram_names);

    for (int i = 0; i < n_diagrams; i++) {
        cc_free(diagram_names[i]);
    }
    cc_free(diagram_names);

    /*
     * Write Heff to the formatted file
     */
    write_formatted_heff(sect_h, sect_p, block_dims, heff);

    /*
     * Diagonalization
     */
    diagonalize_heff(sect_h, sect_p, block_dims, heff, eigvalues, coef_left, coef_right);

    /*
     * For the 1h1p sector only: restoration of the intermediate normalization
     */
    if (sect_h == 1 && sect_p == 1) {
        restore_intermediate_normalization(ms_size, det_list, block_dims, heff, eigvalues);
    }

    /*
     * Save model vectors and eigenvalues to the formatted files MVCOEF
     */
    write_model_vectors_unformatted(sect_h, sect_p, ms_size, det_list, block_dims,
                                    eigvalues, coef_left, coef_right);

    /*
     * Print model vectors
     */
    print_model_vectors_stdout(sect_h, sect_p, ms_size, det_list, block_dims, eigvalues, coef_left, coef_right);

    /*
     * Print table with energy levels
     */
    print_eigenvalues_table(sect_h, sect_p, block_dims, eigvalues, cc_opts->degen_thresh);

    /*
     * Cleanup
     */
    cc_free(det_list);
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        cc_free(heff[irrep]);
        cc_free(eigvalues[irrep]);
        cc_free(coef_left[irrep]);
        cc_free(coef_right[irrep]);
    }
}


void construct_heff(int sect_h, int sect_p, slater_det_t *det_list, size_t *block_sizes, double complex **heff_blocks,
                    int n_diagrams, char **diagram_names)
{
    /*
     * Allocate memory for the effective Hamiltonian (Heff)
     */
    for (int irep = 0; irep < CC_MAX_NUM_IRREPS; irep++) {
        heff_blocks[irep] = NULL;
    }
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        size_t dim = block_sizes[irep];
        if (dim != 0) {
            heff_blocks[irep] = zzeros(dim, dim);
        }
    }

    /*
     * Construct symmetry blocks of the Heff matrix
     */
    size_t rep_offset = 0;
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_sizes[irep] == 0) {
            continue;
        }

        // get pointer to the current block of det-s with symmetry 'irep'
        slater_det_t *rep_dets = det_list + rep_offset;
        rep_offset += block_sizes[irep];  // go to the next block of dets

        // size of the Heff subblock
        size_t ms_size = block_sizes[irep];

        // construct zero-order Hamiltonian
        zero_order_heff(sect_h, sect_p, ms_size, rep_dets, heff_blocks[irep]);

        // loop over n-particle contibutions (n=1,2,...) to the Veff operator
        for (int ioper = 0; ioper < n_diagrams; ioper++) {
            char *dg_name = diagram_names[ioper];
            diagram_t *dg_veff = diagram_stack_find(dg_name);
            if (dg_veff == NULL) {
                errquit("in construct_heff(): diagram '%s' not found", dg_name);
            }
            setup_slater(dg_veff, (matrix_getter_fun) diagram_get, sect_h, sect_p, sect_h, sect_p, rank(dg_name) / 2);

            // construct matrix of the effective interaction (Veff) in the basis
            // of Slater determinants
            for (int i = 0; i < ms_size; i++) {
                for (int j = 0; j < ms_size; j++) {
                    slater_det_t *bra = rep_dets + i;
                    slater_det_t *ket = rep_dets + j;
                    heff_blocks[irep][i * ms_size + j] += slater(bra, ket);
                }
            }
        }
    }
}


/**
 * Writes blocks of the effective Hamiltonian to the formatted file 'HEFF'
 *
 * @param sect_h
 * @param sect_p
 * @param block_dims
 * @param heff
 */
void write_formatted_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff)
{
    int first_irrep = first_nonzero_irrep(block_dims);

    // open file with formatted Heff
    FILE *hefff = hefff_open(sect_h, sect_p, NULL);

    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_dims[irep] == 0) {
            continue;
        }
        double complex *heff_block = heff[irep];
        size_t ms_size = block_dims[irep];
        hefff_write_block(hefff, carith, irep - first_irrep + 1, ms_size, heff_block);
    }

    hefff_close(hefff);
}


void diagonalize_heff(int sect_h, int sect_p, size_t *block_dims, double complex **heff,
                      double complex **eigvalues, double complex **coef_left, double complex **coef_right)
{
    /*
     * Allocate memory for left, right eigenvectors and eigenvalues
     */
    for (int irep = 0; irep < CC_MAX_NUM_IRREPS; irep++) {
        eigvalues[irep] = NULL;
        coef_left[irep] = NULL;
        coef_right[irep] = NULL;
    }
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        size_t dim = block_dims[irep];
        if (dim != 0) {
            eigvalues[irep] = zzeros(dim, 1);
            coef_left[irep] = zzeros(dim, dim);
            coef_right[irep] = zzeros(dim, dim);
        }
    }

    /*
     * Diagonalize each block of the Heff matrix
     * 1. eigenvalues are sorted in ascending order
     * 2. eigenvectors are resorted too (and biorthonormalized)
     */
    size_t max_dim = size_t_max(get_num_irreps(), block_dims);
    double complex *heff_work = zzeros(max_dim, max_dim);
    for (int irep = 0; irep < get_num_irreps(); irep++) {
        if (block_dims[irep] == 0) {
            continue;
        }
        size_t ms_size = block_dims[irep];
        double complex *heff_block = heff[irep];
        double complex *eigvalues_block = eigvalues[irep];
        double complex *coef_left_block = coef_left[irep];
        double complex *coef_right_block = coef_right[irep];

        memcpy(heff_work, heff_block, sizeof(double complex) * ms_size * ms_size);
        eig(ms_size, heff_work, eigvalues_block, coef_left_block, coef_right_block);
    }
    cc_free(heff_work);

    /*
     * Loewdin orthogonalization of model vectors (if required)
     * NOTE: if no orthogonalization is performed, property matrices will be non-Hermitian
     */
    if (cc_opts->do_hermit == 1) {
        for (int irep = 0; irep < get_num_irreps(); irep++) {
            if (block_dims[irep] == 0) {
                continue;
            }
            size_t ms_size = block_dims[irep];
            double complex *coef_left_block = coef_left[irep];
            double complex *coef_right_block = coef_right[irep];

            loewdin_orth(ms_size, coef_right_block, coef_right_block, cc_opts->print_level == CC_PRINT_DEBUG ? 5 : 0);
            memcpy(coef_left_block, coef_right_block, ms_size * ms_size * sizeof(double complex));
        }
    }
}


double max_energy_of_required_roots(size_t *block_dims, double complex **eigvalues)
{
    double max_energy = -1.0e12;

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {

        // get max number of roots for this irrep
        char *irrep_name = get_irrep_name(irrep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : block_dims[irrep];

        double complex *ev = eigvalues[irrep];
        if (nroots > 0) {
            for (size_t i = 0; i < nroots; i++) {
                if (creal(ev[i]) > max_energy) {
                    max_energy = creal(ev[i]);
                }
            }
        }
    }

    return max_energy;
}


slater_det_t *find_slater_dets_for_irrep(size_t total_ms_dim, slater_det_t *det_list, int irrep)
{
    for (size_t i = 0; i < total_ms_dim; i++) {
        if (det_list[i].sym == irrep) {
            return det_list + i;
        }
    }

    return NULL;
}


void restore_intermediate_normalization(size_t total_ms_dim, slater_det_t *det_list, size_t *block_dims,
                                        double complex **heff, double complex **eigvalues)
{
    if (!(cc_opts->sector_h == 1 && cc_opts->sector_p == 1)) {
        return;
    }

    int first_irrep = first_nonzero_irrep(block_dims);
    int vacuum_irrep = get_vacuum_irrep();
    char *vacuum_irrep_name = get_irrep_name(vacuum_irrep);
    slater_det_t *det_basis = find_slater_dets_for_irrep(total_ms_dim, det_list, vacuum_irrep);
    size_t dim = block_dims[vacuum_irrep];
    double complex *heff_block = heff[vacuum_irrep];

    size_t dim_prime = dim + 1;  // dimension of the Heff' matrix
    double complex *heff_0h0p_1h1p_prime = zzeros(dim_prime, dim_prime);
    double complex *p_omega_p = zzeros(dim_prime, dim_prime);

    // 0h0p-1h1p: transformation of the Heff to with the (P\OmegaP)^-1 matrix.
    // + store this 0h0p+1h1p block; it will be written to the formatted file
    // AFTER all blocks belonging to the "pure" 1h1p sector
    renormalize_wave_operator_0h0p_0h1p(dim, det_basis, heff_block, heff_0h0p_1h1p_prime, p_omega_p);

    // write the 0h0p+1h1p transformed Heff block
    // together with the transformation P \Omega P matrix
    FILE *hefff = hefff_open(1, 1, "0011");
    // heff block
    hefff_write_block(hefff, carith, vacuum_irrep-first_irrep+1, dim_prime, heff_0h0p_1h1p_prime);
    // P \omega P
    for (int i = 0; i < dim_prime * dim_prime; i++) {
        if (carith) {
            fprintf(hefff, "%21.12E%21.12E", creal(p_omega_p[i]), cimag(p_omega_p[i]));
            if (i > 0 && i % 2 != 0) { fprintf(hefff, "\n"); }
        }
        else {
            fprintf(hefff, "%21.12E", creal(p_omega_p[i]));
            if (i > 0 && (i + 1) % 4 == 0) { fprintf(hefff, "\n"); }
        }
    }
    if (dim * dim % 2 != 0) { fprintf(hefff, "\n"); }
    hefff_close(hefff);

    /*
     * The 0h0p-1h1p Heff block should be diagonalized since its eigenvectors
     * are required to construct properties and quasi-natural orbitals.
     * These eigenvectors are to be stored in the MVCOEF0011 unformatted file.
     */
    double complex *vl_0011 = zzeros(dim_prime, dim_prime);
    double complex *vr_0011 = zzeros(dim_prime, dim_prime);
    double complex *ev_0011 = zzeros(dim_prime, 1);
    int nroots_irep = get_nroots_for_irrep(vacuum_irrep_name);
    size_t nroots = cc_opts->nroots_specified ? nroots_irep : dim_prime;

    eig(dim_prime, heff_0h0p_1h1p_prime, ev_0011, vl_0011, vr_0011);
    if (cc_opts->do_hermit == 1) {
        loewdin_orth(dim_prime, vr_0011, vr_0011, cc_opts->print_level == CC_PRINT_DEBUG ? 5 : 0);
        memcpy(vl_0011, vr_0011, dim_prime * dim_prime * sizeof(double complex));
    }

    // print model vectors to the unformatted file MVCOEF**
    slater_det_t *det_list_0011 = (slater_det_t *) cc_malloc(sizeof(slater_det_t) * dim_prime);
    det_list_0011[0].indices[0] = 0;
    det_list_0011[0].indices[1] = 0;
    det_list_0011[0].sym = vacuum_irrep;
    memcpy(det_list_0011 + 1, det_basis, sizeof(slater_det_t) * (dim_prime - 1));
    ev_0011[0] = 0.0 + 0.0 * I;

    int f_mvcoef_0011 = io_open("MVCOEF0011", "w");
    mvcoef_write_vectors_unformatted(f_mvcoef_0011, vacuum_irrep_name, nroots + 1, dim_prime, det_list_0011, ev_0011,
                                     vl_0011, vr_0011);

    //double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);
    mvcoef_close(f_mvcoef_0011, 0.0);  // lowest root == 0.0
    cc_free(vl_0011);
    cc_free(vr_0011);
    cc_free(ev_0011);

    cc_free(heff_0h0p_1h1p_prime);
    cc_free(p_omega_p);
}


/*
 * Prints model vectors to stdout (formatted output).
 * If the 'nroots' directive is given in the input file, only the
 * vectors specified by this directive will be printed and analyzed.
 */
void print_model_vectors_stdout(int sector_h, int sector_p,
                                size_t total_ms_size, slater_det_t *det_list, size_t *block_dims,
                                double complex **eigvalues, double complex **coef_left, double complex **coef_right)
{
    double const COEF_THRESH = 1e-4;   // threshold for printing model vec-s coeff-s
    int first_irrep = first_nonzero_irrep(block_dims);

    printf("\n Sector (%dh,%dp) -- analysis of model vectors (right vectors)\n", sector_h, sector_p);
    printf(" first line : irrep, state number, total energy, eigenvalue\n");
    printf(" other lines: coefficients of contributing determinants (above a threshold of %.1e)\n",
           COEF_THRESH);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t dim = block_dims[irrep];
        slater_det_t *det_basis = find_slater_dets_for_irrep(total_ms_size, det_list, irrep);

        // get max number of roots for this irrep
        char *irrep_name = get_irrep_name(irrep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : block_dims[irrep];
        if (nroots == 0) {
            continue;
        }

        for (size_t i = 0; i < nroots; i++) {
            print_model_vector(stdout, sector_h, sector_p,
                               irrep - first_irrep + 1, irrep_name, i, eigvalues[irrep][i],
                               dim, coef_right[irrep] + dim * i, det_basis, 0);
        }
    }
}


double get_lowest_eigenvalue(size_t *block_dims, double complex **eigvalues)
{
    double min_eig_value = 0.0;

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        if (block_dims[irrep] == 0) {
            continue;
        }
        if (creal(eigvalues[irrep][0]) < min_eig_value) {
            min_eig_value = creal(eigvalues[irrep][0]);
        }
    }

    return min_eig_value;
}


void write_model_vectors_unformatted(int sector_h, int sector_p,
                                     size_t total_ms_size, slater_det_t *det_list, size_t *block_dims,
                                     double complex **eigvalues, double complex **coef_left, double complex **coef_right)
{
    // open unformatted file with model vectors
    int f_mvcoef = mvcoef_open(sector_h, sector_p);

    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t dim = block_dims[irrep];
        slater_det_t *det_basis = find_slater_dets_for_irrep(total_ms_size, det_list, irrep);

        // get max number of roots for this irrep
        char *irrep_name = get_irrep_name(irrep);
        int nroots_irep = get_nroots_for_irrep(irrep_name);
        size_t nroots = cc_opts->nroots_specified ? nroots_irep : block_dims[irrep];
        if (nroots == 0) {
            continue;
        }

        // print model vectors to the unformatted file MVCOEF**
        mvcoef_write_vectors_unformatted(f_mvcoef, irrep_name, nroots, block_dims[irrep], det_basis,
                                         eigvalues[irrep], coef_left[irrep], coef_right[irrep]);
    }

    double lowest_root = get_lowest_eigenvalue(block_dims, eigvalues);
    mvcoef_close(f_mvcoef, lowest_root);
}


void print_eigenvalues_table(int sect_h, int sect_p, size_t *block_dims,
                             double complex **eigvalues, double degen_thresh)
{
    /*
     * Allocate table for eigenvalues
     */
    size_t total_heff_dim = 0;
    for (int i = 0; i < get_num_irreps(); i++) {
        total_heff_dim += block_dims[i];
    }
    if (sect_h == 1 && sect_p == 1) {
        total_heff_dim += 1;  // for the reference energy
    }
    eigval_t *eigenvalues = (eigval_t *) cc_malloc(total_heff_dim * sizeof(eigval_t));
    size_t n_eigenvalues = 0;

    /*
     * Copy eigenvalues to the table
     */
    if (sect_h == 1 && sect_p == 1) {
        eigenvalues[0].eigval = 0.0;
        eigenvalues[0].repno = get_vacuum_irrep();
        n_eigenvalues++;
    }
    for (int irrep = 0; irrep < get_num_irreps(); irrep++) {
        size_t block_dim = block_dims[irrep];
        for (size_t i = 0; i < block_dim; i++) {
            eigenvalues[n_eigenvalues + i].eigval = eigvalues[irrep][i];
            eigenvalues[n_eigenvalues + i].repno = irrep;
        }
        n_eigenvalues += block_dim;
    }

    /*
     * Sort table in ascending order
     */
    qsort(eigenvalues, n_eigenvalues, sizeof(eigval_t), eigval_cmp);

    /*
     * Print table with energy levels
     */
    double complex e0 = eigenvalues[0].eigval;
    double reference_energy = cc_opts->eref;   // (see options.h)
    double max_energy = max_energy_of_required_roots(block_dims, eigvalues);

    printf("\n Heff eigenvalues:\n (degeneracy threshold = %.1e a.u.)\n\n", degen_thresh);
    printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  deg  symmetry\n");
    printf(" ------------------------------------------------------------------------------------------------------------------------\n");

    int ilevel = 1;
    for (size_t i = 0; i < n_eigenvalues; i++) {

        // only energy levels with eigenvalue Ei <= max_energy will be printed
        if (creal(eigenvalues[i].eigval) > max_energy + degen_thresh) {
            break;
        }

        // calculate degeneracy
        int deg = 0;
        for (size_t j = i; j < n_eigenvalues; j++) {
            if (cabs(eigenvalues[i].eigval - eigenvalues[j].eigval) > degen_thresh) {
                break;
            }
            else {
                deg++;
            }
        }

        // print info about (maybe degenerate) level
        double complex ei = eigenvalues[i].eigval;
        double abs_energy = reference_energy + creal(ei);
        double rel_energy = creal(ei) - creal(e0);
        double rel_energy_cm = rel_energy * CODATA_AU_TO_CM;
        double rel_energy_ev = rel_energy * CODATA_AU_TO_EV;
        printf("@%5d%16.10f%10.2e%25.17f%16.10f%16.10f%16.6f  %2d  ",
               ilevel, creal(ei), cimag(ei), abs_energy, rel_energy,
               rel_energy_ev, rel_energy_cm, deg);
        for (size_t irrep = 0; irrep < get_num_irreps(); irrep++) {
            int rep_deg = 0;
            for (size_t j = 0; j < deg; j++) {
                if (eigenvalues[i + j].repno == irrep) {
                    rep_deg++;
                }
            }
            if (rep_deg == 0) { continue; }
            printf(" %s", get_irrep_name(irrep));
            if (rep_deg > 1) {
                printf("(%d)", rep_deg);
            }
        }
        printf("\n");

        // go to the next (maybe degenerate) energy level
        i += deg - 1;
        ilevel++;
    }

    /*
     * Ionization potential
     */
    double lowest_root = creal(eigenvalues[0].eigval);
    if (sect_h != sect_p) {
        printf("\n Ionization potential wrt reference state = %18.12f a.u. = %8.4f eV = %10.2f cm^-1\n",
               -lowest_root, -lowest_root * CODATA_AU_TO_EV, -lowest_root * CODATA_AU_TO_CM);
    }
}


void model_space_properties_and_natural_orbitals(int sect_h, int sect_p)
{
    // calculation of density matrices and natural orbitals
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_denmat) {
        for (int ipair = 0; ipair < cc_opts->n_denmat; ipair++) {
            int sect2_h = cc_opts->denmat_query[ipair].sect2[0];
            int sect2_p = cc_opts->denmat_query[ipair].sect2[1];
            if (!(sect2_h == sect_h && sect2_p == sect_p)) {
                continue;
            }
            int rep1 = get_rep_number(cc_opts->denmat_query[ipair].rep1_name);
            int rep2 = get_rep_number(cc_opts->denmat_query[ipair].rep2_name);
            int state1 = cc_opts->denmat_query[ipair].state1;
            int state2 = cc_opts->denmat_query[ipair].state2;
            // single code for NOs and NTOs:
            // rep1 == rep2 && state1 == state2  =>  single state (NO)
            // else: pair of states (NTO)
            density_matrix(sect_h, sect_p, rep1, state1, rep2, state2);
        }
    }

    // transition dipole moments via the DL-TDM techniques
    // for the target sector only
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->do_diplen_tdm) {
        dipole_length_tdms(sect_h, sect_p);
    }

    // estimate properties: only for the target sector
    if (cc_opts->sector_h == sect_h && cc_opts->sector_p == sect_p &&
        cc_opts->n_ms_props > 0) {
        for (int i = 0; i < cc_opts->n_ms_props; i++) {
            model_space_property(cc_opts->prop_queries + i);
        }
    }
}


size_t first_nonzero_irrep(size_t *block_dims)
{
    for (size_t i = 0; i < get_num_irreps(); i++) {
        if (block_dims[i] != 0) {
            return i;
        }
    }

    return 0;
}


/**
 * construct zero-order effective Hamiltonian (H0) matrix:
 *  - only the block belonging to one irrep;
 *  - only one-electron spinor energies are used.
 * It is a diagonal matrix: heff[i,i] = sum{particles} eps - sum{holes} eps
 *
 * @param sect_h, sect_p  Fock space sector signature
 * @param dim dimension of this Heff block (number of model det-s of one symmetry)
 * @param det_list list of model space determinants (with one symmetry)
 * @param heff Heff matrix, of size dim x dim
 */
void zero_order_heff(int sect_h, int sect_p, size_t dim, slater_det_t *det_list, double complex *heff)
{
    for (size_t i = 0; i < dim; i++) {
        slater_det_t *d = det_list + i;

        double H0 = 0.0;
        for (int h = 0; h < sect_h; h++) {
            H0 -= spinor_info[d->indices[h]].eps;
        }
        for (int p = 0; p < sect_p; p++) {
            H0 += spinor_info[d->indices[sect_h + p]].eps;
        }
        heff[i * dim + i] = H0 + 0.0 * I;
    }
}


/*******************************************************************************
 * eigval_cmp
 *
 * Comparator for eigenvalues (for qsort)
 ******************************************************************************/
int eigval_cmp(const void *aa, const void *bb)
{
    const eigval_t *a = aa, *b = bb;
    return (creal(a->eigval) < creal(b->eigval)) ? -1 : (creal(a->eigval) > creal(b->eigval));
}


/*******************************************************************************
 * eigenvalues_table_old
 *
 * Prints nice table with all information about eigenvalues and symmetries of
 * solutions
 ******************************************************************************/
void eigenvalues_table_old(size_t n_eigenvalues, eigval_t *eigenvalues, double degen_thresh, double max_energy)
{
    size_t i, j;
    size_t ilevel;
    int irep;
    int rep_deg;
    int n_nz_reps;
    int nz_reps[64]; // irreps with non-zero number of electronic states
    double complex e0 = eigenvalues[0].eigval;
    double complex ei;
    double rel_energy;
    double abs_energy;
    double rel_energy_ev;
    double rel_energy_cm;
    double eref = cc_opts->eref;   // (see options.h)

    // find reps with nonzero Heff and their number
    n_nz_reps = 0;
    for (i = 0; i < get_num_irreps(); i++) {
        nz_reps[i] = 0;
    }
    for (i = 0; i < n_eigenvalues; i++) {
        nz_reps[eigenvalues[i].repno]++;
    }
    for (i = 0; i < get_num_irreps(); i++) {
        if (nz_reps[i] != 0) {
            nz_reps[n_nz_reps++] = i;
        }
    }

    printf("\n Heff eigenvalues:\n (degeneracy threshold = %.1e a.u.)\n\n", degen_thresh);

    // print table header (see symmetry.h for 'rep_names[]')
    printf(" Level  Re(eigenvalue)  Im(eigv)               Abs energy  Rel eigenvalue    Rel eigv, eV  Rel eigv, cm-1  deg  symmetry\n");
    /*for (irep = 0; irep < n_nz_reps; irep++) {
        printf("%4s|", rep_names[nz_reps[irep]]);
    }*/
    printf(" ------------------------------------------------------------------------------------------------------------------------\n");

    // print table with energy levels
    // only energy levels with eigenvalue Ei <= max_energy+1e-5 will be printed
    ilevel = 1;
    for (i = 0; i < n_eigenvalues; i++) {
        if (creal(eigenvalues[i].eigval) > max_energy + 1.0e-5) {
            break;
        }

        // calculate degeneracy
        int deg = 0;
        for (j = i; j < n_eigenvalues; j++) {
            if (cabs(eigenvalues[i].eigval - eigenvalues[j].eigval) > degen_thresh) {
                break;
            }
            else {
                deg++;
            }
        }

        // print info about (maybe degenerate) level
        ei = eigenvalues[i].eigval;
        abs_energy = eref + creal(ei);
        rel_energy = creal(ei) - creal(e0);
        rel_energy_cm = rel_energy * CODATA_AU_TO_CM;
        rel_energy_ev = rel_energy * CODATA_AU_TO_EV;
        printf("@%5d%16.10f%10.2e%25.17f%16.10f%16.10f%16.6f  %2d  ",
               ilevel, creal(ei), cimag(ei), abs_energy, rel_energy,
               rel_energy_ev, rel_energy_cm, deg);
        for (irep = 0; irep < n_nz_reps; irep++) {
            rep_deg = 0;
            for (j = 0; j < deg; j++) {
                if (eigenvalues[i + j].repno == nz_reps[irep]) {
                    rep_deg++;
                }
            }
            if (rep_deg == 0) { continue; }
            printf(" %s", rep_names[nz_reps[irep]]);
            if (rep_deg > 1) {
                printf("(%d)", rep_deg);
            }
        }
        printf("\n");

        // go to the next (maybe degenerate) energy level
        i += deg - 1;
        ilevel++;
    }
}


/**
 * Prints model vectors.
 * coef_thresh -- threshold for printing model vec-s coeff-s
 */
void print_model_vector(FILE *f_out, int sect_h, int sect_p,
                        int rep_no, char *rep_name, int state_no, double complex energy,
                        size_t len, double complex *coeffs, slater_det_t *det_list, double coef_thresh)
{
    fprintf(f_out, "\n");
    fprintf(f_out, " Irrep %d (%4s) State %d Energy %23.15f Eigenvalue %14.8f%14.4E\n",
            rep_no, rep_name, state_no + 1, cc_opts->eref + creal(energy), creal(energy), cimag(energy));

    for (size_t j = 0; j < len; j++) {
        double complex coef = coeffs[j];
        if (cabs(coef) < coef_thresh) {
            continue;
        }
        fprintf(f_out, " %10.5f%10.5f ", creal(coef), cimag(coef));
        print_slater_det(f_out, sect_h, sect_p, det_list + j);
    }
}


int get_nroots_for_irrep(char *irrep_name)
{
    int nroots_irep = 0;

    for (int ii = 0; ii < get_num_irreps(); ii++) {
        if (strcmp(irrep_name, cc_opts->nroots_specs.rep_names[ii]) == 0) {
            nroots_irep = cc_opts->nroots_specs.dim[ii];
            break;
        }
    }

    return nroots_irep;
}


/*
 * For transition moments between different Fock space sectors.
 * When using this code, please, cite the following paper:
 *
 * A. Zaitsevskii, A. V. Oleynichenko, E. Eliav,
 * Finite-field calculations of transition properties by the fock space
 * relativistic coupled cluster method: transitions between different
 * Fock space sectors.
 * Symmetry, 2020 (submitted).
 */


/**
 * constructs matrix of the model-space projection of the wave operator \Omega
 * in the basis of 0h0p (vacuum) and 1h1p (singly excited) determinants:
 * P \Omega P, P=0h0p+1h1p
 * note that the vacuum and 1h1p determinants must be of the same symmetry.
 *
 * T1 = exc cluster oper (0h0p)
 * S1 = de-exc cluster oper (1h1p)
 */
void omega_0h0p_0h1p(size_t dim, slater_det_t *det_list_1h1p, double complex *omega,
                     char *exc_oper_name, char *deexc_oper_name)
{
    diagram_t *dg_exc = diagram_stack_find(exc_oper_name);
    diagram_t *dg_deexc = diagram_stack_find(deexc_oper_name);
    assert(dg_exc != NULL && dg_deexc != NULL);

    // < 0h0p | Omega | 0h0p > = 1
    omega[0] = 1.0 + 0.0 * I;

    // < 0h0p | Omega | 1h1p > = S1[a,i]
    for (size_t iket = 0; iket < dim; iket++) {
        slater_det_t *det_ia = det_list_1h1p + iket;
        int i = det_ia->indices[0];
        int a = det_ia->indices[1];
        //omega[iket+1] = diagram_get_2(dg_deexc, a, i);
        omega[(dim + 1) * (iket + 1)] = diagram_get_2(dg_deexc, a, i);
    }

    // < 1h1p | Omega | 0h0p > = T1[i,a]
    for (size_t ibra = 0; ibra < dim; ibra++) {
        slater_det_t *det_ia = det_list_1h1p + ibra;
        int i = det_ia->indices[0];
        int a = det_ia->indices[1];
        //omega[(dim+1)*(ibra+1)] = diagram_get_2(dg_exc, i, a);
        omega[ibra + 1] = diagram_get_2(dg_exc, i, a);
    }

    // < 1h1p i,a | Omega | 1h1p j,b > = \delta_ij \delta_ab + T1[i,a] * S1[j,b]
    for (size_t ibra = 0; ibra < dim; ibra++) {
        for (size_t iket = 0; iket < dim; iket++) {
            slater_det_t *det_ia = det_list_1h1p + ibra;
            slater_det_t *det_jb = det_list_1h1p + iket;
            int i = det_ia->indices[0];
            int a = det_ia->indices[1];
            int j = det_jb->indices[0];
            int b = det_jb->indices[1];
            //omega[(dim+1) * (ibra+1) + iket + 1] = (i==j)*(a==b) + diagram_get_2(dg_exc, i, a) * diagram_get_2(dg_deexc, b, j);
            omega[(dim + 1) * (iket + 1) + ibra + 1] =
                    (i == j) * (a == b) + diagram_get_2(dg_exc, i, a) * diagram_get_2(dg_deexc, b, j);
        }
    }
}


/**
 * Performs renormalization of the wave operator \Omega in order to force
 * intermediate normalization P\OmegaP=P in the P=0h0p+1h1p model space.
 *
 * @param dim dimension of the Heff block in the 1h1p subspace
 * @param det_list list of Slater determinants; their symmetry
 *      must coincide with the symmetry of the vacuum determinant
 * @param heff Heff block in the 1h1p subspace
 * @param heff_prime Heff block in the 0h0p+1h1p subspace, transformed
 *      accordingly with renormalized wave operator \Omega
 *      dim(heff') = dim(heff)+1
 * @param omega transformation matrix (model-space part of the wave operator)
 *      dim(omega) = dim(heff)+1
 */
void renormalize_wave_operator_0h0p_0h1p(size_t dim, slater_det_t *det_list,
                                         double complex *heff, double complex *heff_prime, double complex *omega)
{
    size_t dimx = dim + 1;  // dimension of eXtended (0h0p+1h1p) matrices
    double complex alpha = 1.0 + 0.0 * I;
    double complex beta = 0.0 + 0.0 * I;

    // allocate working arrays, init with zeros
    memset(omega, 0, sizeof(double complex) * dimx * dimx);
    double complex *omega_inv = zzeros(dimx, dimx);
    double complex *heff_0h0p_1h1p = zzeros(dimx, dimx);
    double complex *buf = zzeros(dimx, dimx);

    // construct 0h0p+1h1p (extended) block-diagonal Heff matrix
    // Heff 1h1p block will be placed into the right lower angle.
    // Ecorr is subtracted from the whole diagonal.
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            heff_0h0p_1h1p[dimx * (i + 1) + j + 1] = heff[dim * i + j];
        }
    }

    // construct P \Omega P matrix (P=0h0p+1h1p)
    omega_0h0p_0h1p(dim, det_list, omega, "t1c", "e1c");

    // construct inverse matrix (P \Omega P)^{-1}
    inv(dimx, omega, omega_inv);

    // Heff' = { (P\OmegaP) * { Heff * (P\OmegaP)^-1 } }
    // two matrix multiplications
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                dimx, dimx, dimx, &alpha, heff_0h0p_1h1p, dimx, omega_inv, dimx, &beta, buf, dimx);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                dimx, dimx, dimx, &alpha, omega, dimx, buf, dimx, &beta, heff_prime, dimx);

    // deallocate working arrays
    cc_free(omega_inv);
    cc_free(heff_0h0p_1h1p);
    cc_free(buf);
}


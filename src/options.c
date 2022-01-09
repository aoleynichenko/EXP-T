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
 * options.c
 * =========
 *
 * Functions for object-oriented-like operations with 'struct options'.
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codata.h"
#include "platform.h"
#include "linalg.h"
#include "error.h"
#include "intham1.h"
#include "options.h"
#include "memory.h"

/******************************************************************************/

// arithmetic: CC_ARITH_REAL or CC_ARITH_COMPLEX
int arith = -1;

// complex arithmetic?
int carith = -1;

// sizeof(double) or sizeof(double complex)
int SIZEOF_WORKING_TYPE = sizeof(double complex);

data_type_t WORKING_TYPE = CC_DOUBLE_COMPLEX;

// pointer to options
cc_options_t *cc_opts = NULL;

void print_space(cc_space_t *space);

/******************************************************************************/

// creates new object of type cc_options_t and initialize it with default values
// returns pointer to the created object
cc_options_t *new_options()
{
    int i;

    cc_options_t *opts = (cc_options_t *) cc_malloc(sizeof(cc_options_t));
    if (opts == NULL) {
        errquit("Unable to allocate memory for options");
    }

    // set default values to all fields
    strcpy(opts->title, "no title");
    strcpy(opts->scratch_dir, "scratch");
    opts->clean_scratch = 1;
    opts->print_level = CC_PRINT_MEDIUM;
    opts->recommended_arith = CC_ARITH_REAL;
    opts->max_memory_size = 1024u * 1024u * 1024u;  // 1 Gb
    opts->compress = CC_COMPRESS_NONE;

    // compression of arrays with triples amplitudes
    opts->do_compress_triples = 0;
    opts->compress_triples_thresh = 1e-10;
    opts->compress_triples_data_type = CC_DOUBLE;

    opts->tile_size = 100;
    opts->disk_usage_level = CC_DISK_USAGE_LEVEL_2;  // rank-6+ and pppp on disk
    opts->nthreads = 1;
    opts->cuda_enabled = 0;
    opts->maxiter = 50;
    opts->conv = 1e-9;
    opts->div_thresh = 1e3;
    opts->int_source = CC_INTEGRALS_DIRAC;
    strcpy(opts->integral_file_1, "MRCONEE");
    strcpy(opts->integral_file_2, "MDCINT");
    strcpy(opts->integral_file_prop, "MDPROP");
    opts->x2cmmf = 0;
    // for Daniel Maison integral program
    opts->gaunt_defined = 0;
    opts->breit_defined = 0;
    strcpy(opts->integral_file_gaunt, "MGINT");
    strcpy(opts->integral_file_breit, "MBINT");
    opts->escf = 0.0;
    opts->enuc = 0.0;
    opts->eref = 0.0;
    opts->ground_energy_0h1p = 0.0;
    opts->cc_model = CC_MODEL_CCSD;
    strcpy(opts->cc_model_str, "CCSD");
    opts->sector_h = 0;
    opts->sector_p = 0;
    opts->mixed = 0;
    opts->curr_sector_h = 0;
    opts->curr_sector_p = 0;

    // active space specification
    // [inact H] actsp_min [act H] [act P] actsp_max [inact P]
    // default: no active space
    opts->actsp_defined = CC_ACT_SPACE_UNDEFINED;  // undefined
    opts->actsp_min = +1e10;
    opts->actsp_max = -1e10;
    opts->nacth = 0;
    opts->nactp = 0;
    opts->main_defined = 0;
    opts->main_h = 0;
    opts->main_p = 0;
    memset(opts->active_specs, 0, sizeof(opts->active_specs));
    memset(opts->active_each, 0, sizeof(opts->active_each));

    // number of electrons in each irrep (for open-shell reference states)
    opts->nelec_defined = 0;
    memset(&opts->irrep_occ_numbers, 0, sizeof(opts->irrep_occ_numbers));

    // occupation numbers for each spinor (for open-shell reference states)
    opts->occ_defined = 0;
    for (i = 0; i < CC_MAX_SPINORS; i++) {
        opts->occ[i] = 0;
    }

    // resort integrals or not
    opts->reuse_integrals_1 = 0;
    opts->reuse_integrals_2 = 0;

    // read precomputed amplitudes from disk (as initial guess) or not
    opts->reuse_0h0p = 0;
    opts->reuse_0h1p = 0;
    opts->reuse_1h0p = 0;
    opts->reuse_1h1p = 0;
    opts->reuse_0h2p = 0;
    opts->reuse_2h0p = 0;
    opts->reuse_0h3p = 0;
    opts->reuse_1h2p = 0;

    // flush non-converged amplitudes to disk
    opts->do_flush_iter = 0;

    // denominator shifts
    memset(opts->shifts, 0, sizeof(opts->shifts));
    memset(opts->orbshifts, 0, sizeof(opts->orbshifts));
    opts->shift_type = CC_SHIFT_NONE;

    opts->do_orbshift = 0;
    opts->orbshift = 0.0;
    opts->orbshift_power = 0;

    opts->do_orbshift_0h0p = 0;
    memset(opts->orbshift_0h0p, 0, sizeof(opts->orbshift_0h0p));
    opts->orbshift_0h0p_power = 0;

    // damping
    memset(opts->damping, 0, sizeof(opts->damping));

    // skip some FS sectors
    memset(opts->skip_sector, 0, sizeof(opts->skip_sector));

    // DIIS -- convergence acceleration
    opts->diis_enabled = 1;
    opts->diis_dim = 10;       // dimension of the DIIS subspace
    opts->diis_triples = 0;   // do DIIS for triples amplitudes or not

    // how many roots should be processed
    opts->nroots_specified = 0;     // all roots by default
    // number of lowest roots in each symmetry
    memset(&opts->nroots_specs, 0, sizeof(opts->nroots_specs));
    // roots of interest can be limited by some energy bound
    opts->roots_cutoff_specified = 0;
    opts->roots_cutoff = 0.0;

    // degeneracy threshold (for printing roots)
    opts->degen_thresh = 1e-8;

    // interface to the OneProp code by L. V. Skripnikov (one-elec properties)
    opts->oneprop_on = 0;
    opts->n_oneprop = 0;
    for (int i = 0; i < CC_MAX_NPROP; i++) {
        strcpy(opts->oneprop_file_re[i], "PropInts");
        strcpy(opts->oneprop_file_im[i], "PropIntsIm");
        opts->oneprop_lambda[i] = 0.0 + 0.0 * I;
    }

    // interface to two-electron properties by D. E. Maison and L. V. Skripnikov
    opts->twoprop_on = 0;
    opts->n_twoprop = 0;
    for (int i = 0; i < CC_MAX_NPROP; i++) {
        strcpy(opts->twoprop_file[i], "TWOPROP");
        opts->twoprop_lambda[i] = 0.0 + 0.0 * I;
    }

    // interface to the MDPROP file
    opts->n_mdprop = 0;
    for (int i = 0; i < CC_MAX_NPROP; i++) {
        strcpy(opts->mdprop_file[i], "\0");
        opts->mdprop_lambda[i] = 0.0 + 0.0 * I;
    }

    // perform "hermitization" of Heff or not (default:hermitization is enabled)
    // Bloch (non-Hermitian) or des Cloizeaux (Hermitian) Hamiltonian?
    // Bloch => biorthogonal model vectors => TDM(i->f) != TDM(f->i) =>
    //          not very physical
    // des Cloizeaux => (Loewdin) orthogonal model vectors =>
    //                  much more meaningful TDMs and properties
    opts->do_hermit = 0;

    // perform Dipole-Length (DL) estimation of TDMs or not
    opts->do_diplen_tdm = 0;

    // calculation of density matrices and natural orbitals
    opts->n_denmat = 0;

    // calculation of properties (model-space estimation)
    opts->n_ms_props = 0;

    // calculation of "effective configurations" for the target states
    opts->do_eff_conf = 0;

    // selection of amplitudes
    opts->n_select = 0;

    // (simple) intermediate hamiltonian ("IH-1") parameters
    opts->do_intham1 = 0;
    memset(&opts->ih1_opts, 0, sizeof(opts->ih1_opts));

    // (sophisticated) intermediate hamiltonian parameters
    opts->do_intham = 0;
    memset(&opts->intham_params.main_space, 0, sizeof(opts->intham_params.main_space));

    // restriction of the spinor space for triples
    opts->do_restrict_t3 = 0;
    opts->restrict_t3_bounds[0] = 0.0;
    opts->restrict_t3_bounds[1] = 0.0;

    return opts;
}


// deallocates memory owned by the object of type 'cc_options_t'
void delete_options(cc_options_t *opts)
{
    if (opts == NULL) {
        return;
    }
    cc_free(opts);
}


int triples_enabled()
{
    int model = cc_opts->cc_model;

    if (model == CC_MODEL_CCSDT_1A ||
        model == CC_MODEL_CCSDT_1B ||
        model == CC_MODEL_CCSDT_1B_PRIME ||
        model == CC_MODEL_CCSDT_2 ||
        model == CC_MODEL_CCSDT_3 ||
        model == CC_MODEL_CCSDT) {
        return 1;
    }

    return 0;
}


// print options
void print_options(cc_options_t *opts)
{
    printf("\n");
    printf("\t\t\t\t*************\n");
    printf("\t\t\t\t** Options **\n");
    printf("\t\t\t\t*************\n\n");
    printf(" %-15s  %-40s  %s\n", "title", "title string - comment", opts->title);
    printf(" %-15s  %-40s  %s\n", "scratch_dir", "scratch directory for tmp files", opts->scratch_dir);
    printf(" %-15s  %-40s  %s\n", "--no-clean", "retain scratch directory on exit", opts->clean_scratch ? "yes" : "no");
    printf(" %-15s  %-40s  ", "print", "print level");
    if (opts->print_level == CC_PRINT_LOW) {
        printf("low\n");
    }
    else if (opts->print_level == CC_PRINT_MEDIUM) {
        printf("medium\n");
    }
    else if (opts->print_level == CC_PRINT_HIGH) {
        printf("high\n");
    }
    else if (opts->print_level == CC_PRINT_DEBUG) {
        printf("debug\n");
    }
    else {
        printf("\n");
    }
    printf(" %-15s  %-40s  %s\n", "arith", "recommended arithmetic",
           opts->recommended_arith == CC_ARITH_REAL ? "real" : "complex");
    printf(" %-15s  %-40s  %.1f Mb\n", "memory", "max allowed RAM usage", opts->max_memory_size / (1024.0 * 1024.0));
    printf(" %-15s  %-40s  %s\n", "compress", "compression of integrals on disk", opts->compress ? "LZ4" : "disabled");
    if (opts->do_compress_triples) {
        printf(" %-15s %-40s  yes, thresh=%g, datatype=%s\n", "compress_triples", "compression of triples tensors in RAM",
               opts->compress_triples_thresh, opts->compress_triples_data_type == CC_FLOAT ? "float" : "double");
    }
    else {
        printf(" %-15s %-40s  %s\n", "compress_triples", "compression of triples tensors in RAM", "disabled");
    }
    printf(" %-15s  %-40s  ", "disk_usage", "disk usage level");
    switch (opts->disk_usage_level) {
        case CC_DISK_USAGE_LEVEL_0:
            printf("0 (all data in RAM)\n");
            break;
        case CC_DISK_USAGE_LEVEL_1:
            printf("1 (rank-6+ on disk)\n");
            break;
        case CC_DISK_USAGE_LEVEL_2:
            printf("2 (rank-6+ and pppp on disk)\n");
            break;
        case CC_DISK_USAGE_LEVEL_3:
            printf("3 (rank-6+, pppp and *ppp on disk)\n");
            break;
        case CC_DISK_USAGE_LEVEL_4:
            printf("4 (rank-6+, pppp and *ppp on disk + compression enabled)\n");
            break;
        default:
            printf("\n");
            break;
    }
    printf(" %-15s  %-40s  %d\n", "tilesize", "max dimension of formal blocks (tiles)", opts->tile_size);
    printf(" %-15s  %-40s  %d\n", "nthreads", "number of OpenMP parallel threads", opts->nthreads);
    printf(" %-15s  %-40s  %s\n", "cuda", "calculations on GPU (CUDA)", opts->cuda_enabled ? "enabled" : "disabled");
    printf(" %-15s  %-40s  %d\n", "maxiter", "maximum number of CC iterations", opts->maxiter);
    printf(" %-15s  %-40s  %g\n", "conv", "convergence threshold (by amplitudes)", opts->conv);
    printf(" %-15s  %-40s  %g\n", "div_thresh", "divergence threshold (by amplitudes)", opts->div_thresh);
    printf(" %-15s  %-40s  ", "reuse", "reuse amplitudes and/or integrals");
    if (opts->reuse_0h0p == 0 &&
        opts->reuse_0h1p == 0 &&
        opts->reuse_1h0p == 0 &&
        opts->reuse_1h1p == 0 &&
        opts->reuse_0h2p == 0 &&
        opts->reuse_2h0p == 0 &&
        opts->reuse_0h3p == 0 &&
        opts->reuse_1h2p == 0 &&
        opts->reuse_integrals_1 == 0 &&
        opts->reuse_integrals_2 == 0) {
        printf("nothing ");
    }
    if (opts->reuse_integrals_1) {
        printf("1-integrals ");
    }
    if (opts->reuse_integrals_2) {
        printf("2-integrals ");
    }
    if (opts->reuse_0h0p) {
        printf("0h0p ");
    }
    if (opts->reuse_0h1p) {
        printf("0h1p ");
    }
    if (opts->reuse_1h0p) {
        printf("0h1p ");
    }
    if (opts->reuse_1h1p) {
        printf("1h1p ");
    }
    if (opts->reuse_0h2p) {
        printf("0h2p ");
    }
    if (opts->reuse_2h0p) {
        printf("2h0p ");
    }
    if (opts->reuse_0h3p) {
        printf("0h3p ");
    }
    if (opts->reuse_1h2p) {
        printf("1h2p ");
    }
    printf("\n");

    printf(" %-15s  %-40s  ", "skip", "skip computations in sectors:");
    for (int i = 0; i < MAX_SECTOR_RANK; i++) {
        for (int j = 0; j < MAX_SECTOR_RANK; j++) {
            if (opts->skip_sector[i][j] != 0) {
                printf("%dh%dp ", i, j);
            }
        }
    }
    printf("\n");

    printf(" %-15s  %-40s  ", "flush", "flush amplitudes");
    if (opts->do_flush_iter > 0) {
        printf("each %d iterations\n", opts->do_flush_iter);
    }
    else {
        printf("no\n");
    }

    printf(" %-15s  %-40s  %s\n", "interface", "source of transformed molecular integrals",
           opts->int_source == CC_INTEGRALS_DIRAC ? "DIRAC" : "tm2c");

    printf(" %-15s  %-40s  %s\n", "integrals", "one-electron Hamiltonian integrals file", opts->integral_file_1);
    printf(" %-15s  %-40s  %s\n", "", "two-electron (Coulomb) integrals file", opts->integral_file_2);
    printf(" %-15s  %-40s  %s\n", "", "one-electron property integrals file", opts->integral_file_prop);
    printf(" %-15s  %-40s  ", "x2cmmf", "X2Cmmf Hamiltonian");
    if (opts->x2cmmf) {
        printf("yes\n");
    }
    else {
        printf("no\n");
    }
    printf(" %-15s  %-40s  ", "gaunt", "two-electron (Gaunt) integrals file");
    if (opts->gaunt_defined) {
        printf("%s\n", opts->integral_file_gaunt);
    }
    else {
        printf("not used\n");
    }

    printf(" %-15s  %-40s  %1dh%1dp\n", "sector", "target Fock space sector", opts->sector_h, opts->sector_p);
    /*if (opts->mixed) {
        printf("  Mixed-sector model 0h0p-1h1p\n");
    }*/
    printf(" %-15s  %-40s  ", "model", "coupled cluster model (level of theory)");
    switch (opts->cc_model) {
        case CC_MODEL_CCS:
            printf("CCS");
            break;
        case CC_MODEL_CCD:
            printf("CCD");
            break;
        case CC_MODEL_CCSD:
            printf("CCSD");
            break;
        case CC_MODEL_CCSD_T3:
            printf("CCSD+T(3)");
            break;
        case CC_MODEL_CCSDT_1A:
            printf("CCSDT-1a (dressed Veff)");
            break;
        case CC_MODEL_CCSDT_1B:
            printf("CCSDT-1b (dressed Veff)");
            break;
        case CC_MODEL_CCSDT_1B_PRIME:
            printf("CCSDT-1b\' (undressed Veff)");
            break;
        case CC_MODEL_CCSDT_2:
            printf("CCSDT-2");
            break;
        case CC_MODEL_CCSDT_3:
            printf("CCSDT-3");
            break;
        case CC_MODEL_CCSDT:
            printf("CCSDT");
            break;
        default:
            break;
    }
    printf("\n");
    if (opts->actsp_defined == 1 && opts->nacth == 0 && opts->nactp == 0) {
        printf(" %-15s  %-40s  ", "active", "active space (defined by energy range)");
        printf("%g <= eps <= %g a.u.\n", opts->actsp_min, opts->actsp_max);
    }
    else if (opts->actsp_defined == 2 && (opts->nacth != 0 || opts->nactp != 0)) {
        printf(" %-15s  %-40s  %d spinors\n", "nacth", "number of active hole spinors", opts->nacth);
        printf(" %-15s  %-40s  %d spinors\n", "nactp", "number of active particle spinors", opts->nactp);
    }
    else if (opts->actsp_defined == 3) {
        printf(" %-15s  %-40s  ", "active", "active holes space (defined by irreps)");
        for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
            cc_active_spec_t *sp = &opts->active_specs[i];
            if (sp->nacth > 0) {
                printf("%s:%d ", sp->rep_name, sp->nacth);
            }
        }
        printf("\n");
        printf(" %-15s  %-40s  ", "active", "active particles space (defined by irreps)");
        for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
            cc_active_spec_t *sp = &opts->active_specs[i];
            if (sp->nactp > 0) {
                printf("%s:%d ", sp->rep_name, sp->nactp);
            }
        }
        printf("\n");
    }
    /*else if (opts->actsp_defined == 0) {
        printf("  active space      for each spinor\n");
    }
    else{
        printf("  active space      undefined\n");
    }*/

    if (opts->main_defined) {
        printf(" %-15s  %-40s  %d act holes, %d act particles\n",
               "main", "main active space", opts->main_h, opts->main_p);
    }

    printf(" %-15s  %-40s  ", "shifttype", "formula for denominator shifts");
    if (opts->shift_type == CC_SHIFT_NONE) {
        printf("shifts are disabled\n");
    }
    else if (opts->shift_type == CC_SHIFT_REAL) {
        printf("real shifts\n");
    }
    else if (opts->shift_type == CC_SHIFT_REALIMAG) {
        printf("real simulation of imaginary shifts\n");
    }
    else if (opts->shift_type == CC_SHIFT_IMAG) {
        printf("imaginary shifts\n");
    }
    else if (opts->shift_type == CC_SHIFT_TAYLOR) {
        printf("taylor\n");
    }
    if (opts->shift_type != CC_SHIFT_NONE) {
        printf(" %-15s  %-40s  %s\n", "orbshift", "\"orbital\" shifts (non-trivial sectors)",
               opts->do_orbshift ? "yes" : "no");
        if (opts->do_orbshift) {
            //printf(" %-15s  %-40s  %d\n", "", "\"orbital\" shift attenuation parameter", opts->orbshift_power);
            //printf(" %-15s  %-40s  %f\n", "", "\"orbital\" shift amplitude", opts->orbshift);
            for (int h = 0; h < MAX_SECTOR_RANK; h++) {
                for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                    if (!opts->orbshifts[h][p].enabled) {
                        continue;
                    }
                    printf(" orbshift %dh%dp    %-40s  ", h, p, "orbital shift parameters");
                    cc_shifttype_t type = opts->orbshifts[h][p].type;
                    int power = opts->orbshifts[h][p].power;
                    double *shifts = opts->orbshifts[h][p].shifts;
                    printf("power=%d ", power);
                    printf("shift=%.3f\n", shifts[0]);
                }
            }
        }
        printf(" %-15s  %-40s  %s\n", "orbshift00", "\"orbital\" shifts (0h0p sector)",
               opts->do_orbshift_0h0p ? "yes" : "no");
        if (opts->do_orbshift_0h0p == 1) {
            printf(" %-15s  %-40s  %d\n", "", "\"orbital\" shift (0h0p) attenuation parameter",
                   opts->orbshift_0h0p_power);
            printf(" %-15s  %-40s  %f\n", "", "\"orbital\" shift (0h0p) amplitude", opts->orbshift_0h0p[0]);
        }
        else if (opts->do_orbshift_0h0p == 2) {
            printf(" %-15s  %-40s  %d\n", "", "\"orbital\" shift (0h0p) attenuation parameter",
                   opts->orbshift_0h0p_power);
            printf(" %-15s  \"orbital\" shift (0h0p) amplitudes\n", "");
            printf(" %-15s  singles     1 valence line   %f\n", "", opts->orbshift_0h0p[0]);
            printf(" %-15s              2 valence lines  %f\n", "", opts->orbshift_0h0p[1]);
            printf(" %-15s              1 valence line   %f\n", "", opts->orbshift_0h0p[2]);
            printf(" %-15s              2 valence lines  %f\n", "", opts->orbshift_0h0p[3]);
            printf(" %-15s              3 valence lines  %f\n", "", opts->orbshift_0h0p[4]);
            printf(" %-15s              4 valence lines  %f\n", "", opts->orbshift_0h0p[5]);
        }
        if (!opts->do_orbshift && !opts->do_orbshift_0h0p) {
            for (int h = 0; h < MAX_SECTOR_RANK; h++) {
                for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                    if (!opts->shifts[h][p].enabled) {
                        continue;
                    }
                    printf(" shift %dh%dp       %-40s  ", h, p, "denominator shift parameters");
                    cc_shifttype_t type = opts->shifts[h][p].type;
                    int power = opts->shifts[h][p].power;
                    double *shifts = opts->shifts[h][p].shifts;
                    printf("power=%d ", power);
                    printf("shifts=");
                    for (int i = 0; i < 3; i++) {
                        printf("%.3f ", shifts[i]);
                    }
                    printf("\n");
                }
            }
        }
    }

    printf(" %-15s  %-40s  ", "nroots", "number of roots to be processed");
    if (opts->nroots_specified == 0) {
        printf("all\n");
    }
    else { // number of roots was specified
        printf("lowest ");
        print_space(&opts->nroots_specs);
        printf("\n");
    }
    printf(" %-15s  %-40s  ", "roots_cutoff", "energy cutoff for roots to be processed");
    if (opts->roots_cutoff_specified == 0) {
        printf("all\n");
    }
    else {
        printf("%.0f cm^-1\n", opts->roots_cutoff * CODATA_AU_TO_CM);
    }

    printf(" %-15s  %-40s  %.1e\n", "degen_thresh", "degeneracy threshold (a.u.)", opts->degen_thresh);
    printf(" %-15s  %-40s  ", "occ_irreps", "occupation numbers of spinors");
    if (opts->occ_defined != 0) {
        printf("not used\n");
    }
    else if (opts->nelec_defined == 0) {
        printf("default (from integral interface)\n");
    }
    else {
        printf("(by irreps) ");
        for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
            int occ_number = opts->irrep_occ_numbers.dim[i];
            char *irrep_name = opts->irrep_occ_numbers.rep_names[i];
            if (occ_number != 0) {
                printf("[%s]:%d ", irrep_name, occ_number);
            }
        }
        printf("\n");
    }
    printf(" %-15s  %-40s  ", "occ", "occupation numbers of spinors");
    if (opts->occ_defined == 0) {
        printf("default\n");
    }
    else {
        printf("defined for each spinor\n");
    }

    if (opts->diis_enabled) {
        printf(" %-15s  %-40s  %s\n", "diis", "DIIS technique for convergence", "enabled");
        printf(" %-15s  %-40s  %d\n", "diis <n>", "DIIS subspace dimension", opts->diis_dim);
        printf(" %-15s  %-40s  %s\n", "diis triples", "DIIS for triples amplitudes",
               opts->diis_triples ? "enabled" : "disabled");
    }
    else {
        printf(" %-15s  %-40s  %s\n", "diis", "DIIS technique for convergence", "disabled");
    }
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            if (opts->damping[h][p].enabled) {
                printf(" %-7s %dh%dp     %-30s%dh%dp     ", "damping", h, p,
                       "damping for amplitudes, sector ", h, p);
                printf("  stop at iter %d, old ampl factor %.2f\n",
                       opts->damping[h][p].stop, opts->damping[h][p].factor);
            }
        }
    }

    // interface to one-electron properties by L. V. Skripnikov
    printf(" %-15s  %-40s  %s\n", "oneprop", "interface to OneProp (by L.V.Skripnikov)",
           opts->oneprop_on ? "enabled" : "disabled");
    if (opts->oneprop_on) {
        for (int i = 0; i < opts->n_oneprop; i++) {
            printf(" %-15s  [%1d] %-36s  re=%s im=%s lambda=(%12.4e,%12.4e)\n",
                   "", i, "one-electron property", opts->oneprop_file_re[i], opts->oneprop_file_im[i],
                   creal(opts->oneprop_lambda[i]), cimag(opts->oneprop_lambda[i]));
        }
    }
    // interface to two-electron properties by D. E. Maison and L. V. Skripnikov
    printf(" %-15s  %-40s  %s\n", "twoprop", "interface to TwoProp (by D.E.Maison)",
           opts->twoprop_on ? "enabled" : "disabled");
    if (opts->twoprop_on) {
        for (int i = 0; i < opts->n_twoprop; i++) {
            printf(" %-15s  [%1d] %-36s  file=%s lambda=(%12.4e,%12.4e)\n",
                   "", i, "two-electron property", opts->twoprop_file[i],
                   creal(opts->twoprop_lambda[i]), cimag(opts->twoprop_lambda[i]));
        }
    }

    // interface to the MDPROP files
    printf(" %-15s  %-40s  %s\n", "oneprop", "interface to the MDPROP file",
           (opts->n_mdprop > 0) ? "enabled" : "disabled");
    if (opts->n_mdprop > 0) {
        for (int i = 0; i < opts->n_mdprop; i++) {
            printf(" %-15s  [%1d] %-36s  property=%s lambda=(%12.4e,%12.4e)\n",
                   "", i, "one-electron property from MDPROP", opts->mdprop_file[i],
                   creal(opts->mdprop_lambda[i]), cimag(opts->mdprop_lambda[i]));
        }
    }

    // perform hermitization or not?
    printf(" %-15s  %-40s  %s\n", "nohermit", "hermitization of effective Hamiltonians",
           opts->do_hermit ? "enabled" : "disabled");

    // perform Dipole-Length (DL) estimation of TDMs or not
    printf(" %-15s  %-40s  %s\n", "dltdm", "model-space estimates of tran dipoles",
           opts->do_diplen_tdm ? "enabled" : "disabled");

    // calculation of density matrices and natural orbitals
    if (opts->n_denmat > 0) {
        printf(" %-15s  %-40s  %s\n", "natorb", "model-space natural orbitals", "enabled");
        printf(" %-15s  %-40s\n", "", "target states (sector:rep:state):");

        for (int i = 0; i < opts->n_denmat; i++) {
            cc_denmat_query_t *q = opts->denmat_query + i;
            int *sect1 = q->sect1;
            int *sect2 = q->sect2;
            char *rep1_name = q->rep1_name;
            char *rep2_name = q->rep2_name;
            int state1 = q->state1;
            int state2 = q->state2;
            if (strcmp(rep1_name, rep2_name) == 0 && state1 == state2) {
                printf(" %-15s  %dh%dp:%s:%d (NO)\n", "", sect1[0], sect1[1], rep1_name, state1);
            }
            else {
                printf(" %-15s  %dh%dp:%s:%d-%dh%dp:%s:%d (NTO)\n", "", sect1[0], sect1[1],
                       rep1_name, state1, sect2[0], sect2[1], rep2_name, state2);
            }
        }
    }
    else {
        printf(" %-15s  %-40s  %s\n", "natorb", "model-space natural orbitals", "disabled");
    }

    // calculation of "effective configurations" for the target states
    printf(" %-15s  %-40s  %s\n", "eff_conf", "eff configurations for the target states",
           opts->do_eff_conf ? "enabled" : "disabled");

    if (opts->n_ms_props > 0) {
        printf(" %-15s  %-40s  %s", "prop", "model-space estimates of properties", "enabled,");
        for (int i = 0; i < opts->n_ms_props; i++) {
            if (opts->prop_queries[i].source == CC_PROP_FROM_MDPROP) {
                printf(" %s (md)", opts->prop_queries[i].prop_name);
            }
            else {
                printf(" %s %s (txt)", opts->prop_queries[i].file_real, opts->prop_queries[i].file_imag);
            }
            if (opts->prop_queries[i].do_transpose) {
                printf("^T");
            }
        }
        printf("\n");
    }
    else {
        printf(" %-15s  %-40s  %s\n", "prop", "model-space estimates of properties", "disabled");
    }

    if (opts->n_select == 0) {
        printf(" %-15s  %-40s  %s\n", "select", "selection of cluster amplitudes", "disabled");
    }
    else {
        printf(" %-15s  %-40s  %s\n", "select", "selection of cluster amplitudes", "enabled");
        for (int i = 0; i < opts->n_select; i++) {
            printf(" %-15s  %dh%dp t%d ", "", opts->selects[i].sect_h, opts->selects[i].sect_p,
                   opts->selects[i].rank / 2);
            switch (opts->selects[i].rule) {
                case CC_SELECTION_ALL:
                    printf("all ");
                    break;
                case CC_SELECTION_SPECTATOR:
                    printf("spectator ");
                    break;
                case CC_SELECTION_ACT_TO_ACT:
                    printf("act_to_act ");
                    break;
                case CC_SELECTION_MAX_2_INACT:
                    printf("max2inact ");
                    break;
                case CC_SELECTION_EXC_WINDOW:
                    printf("exc_window [%g;%g] ", opts->selects[i].e1, opts->selects[i].e2);
                    break;
                case CC_SELECTION_EPS_WINDOW:
                    printf("eps_window [%g;%g] ", opts->selects[i].e1, opts->selects[i].e2);
                    break;
                default:
                    break;
            }
            printf("%s 0\n", opts->selects[i].task == CC_SELECTION_SET_ZERO ? "=" : "!=");
        }
    }

    if (opts->do_restrict_t3 == 1) {
        printf(" %-15s  %-40s  enabled, %g < eps < %g\n", "restrict_t3", "restriction of triples", opts->restrict_t3_bounds[0], opts->restrict_t3_bounds[1]);
    }
    else {
        printf(" %-15s  %-40s  %s\n", "restrict_t3", "restriction of triples", "disabled");
    }

    // (simple) intermediate hamiltonian ("IH-1") parameters
    printf(" %-15s  %-40s  ", "intham1", "simple intermediate Hamiltonian");
    if (opts->do_intham1) {
        printf("enabled in sectors: ");
        for (int h = 0; h < MAX_SECTOR_RANK; h++) {
            for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                if (opts->ih1_opts.sectors[h][p]) {
                    printf("%dh%dp ", h, p);
                }
            }
        }
        printf("\n");
    }
    else {
        printf("disabled\n");
    }

    // (sophisticated) intermediate hamiltonian parameters
    printf(" %-15s  %-40s  ", "intham", "intermediate Hamiltonian");
    if (opts->do_intham == 1) {
        printf("enabled\n");
        printf(" %-15s  %-40s  ", "  main", "subspace of target states");
        printf("lowest ");
        print_space(&opts->intham_params.main_space);
        printf("\n");
    }
    else {
        printf("disabled\n");
    }

    printf("\n\n");
}


void print_space(cc_space_t *space)
{
    assert(space->deftype == CC_SPACE_BY_ENERGY ||
           space->deftype == CC_SPACE_BY_TOTAL ||
           space->deftype == CC_SPACE_BY_IRREPS);

    if (space->deftype == CC_SPACE_BY_ENERGY) {
        printf("energy %.f %.f", space->emin, space->emax);
    }
    else if (space->deftype == CC_SPACE_BY_TOTAL) {
        printf("total %d", space->total);
    }
    else { // by irreps
        for (int i = 0; i < CC_MAX_NUM_IRREPS; i++) {
            if (space->dim[i] != 0) {
                printf("%s:%d ", space->rep_names[i], space->dim[i]);
            }
        }
    }
}


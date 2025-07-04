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

/*
 * Functions for object-oriented-like operations with 'struct options'.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codata.h"
#include "linalg.h"
#include "error.h"
#include "intham_imms.h"
#include "options.h"
#include "memory.h"


// arithmetic: CC_ARITH_REAL or CC_ARITH_COMPLEX
int arith = -1;

// sizeof(double) or sizeof(double complex)
int SIZEOF_WORKING_TYPE = sizeof(double complex);
data_type_t WORKING_TYPE = CC_DOUBLE_COMPLEX;

// pointer to options
cc_options_t *cc_opts = NULL;

void print_space(cc_space_t *space);


/**
 * Creates new object of type cc_options_t and initialize it with default values.
 * returns pointer to the created object.
 */
cc_options_t *new_options()
{
    cc_options_t *opts = (cc_options_t *) cc_calloc(1, sizeof(cc_options_t));
    if (opts == NULL) {
        errquit("Unable to allocate memory for options");
    }

    // set default values to all fields
    strcpy(opts->title, "no title");
    strcpy(opts->scratch_dir, "scratch");
    opts->clean_scratch = 1;
    opts->print_level = CC_PRINT_MEDIUM;
    opts->print_model_space = 0;
    opts->print_model_vectors = 0;
    opts->print_eff_config = 0;
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
    opts->openmp_algorithm = CC_OPENMP_ALGORITHM_EXTERNAL;
    opts->cuda_enabled = 0;
    opts->maxiter = 50;
    opts->conv_thresh = 1e-9;
    opts->div_thresh = 1e3;
    opts->int_source = CC_INTEGRALS_DIRAC;
    strcpy(opts->integral_file_1, "MRCONEE");
    strcpy(opts->integral_file_2, "MDCINT");
    strcpy(opts->integral_file_prop, "MDPROP");
    opts->x2cmmf = 0;
    opts->new_sorting = 0;
    // for Daniel Maison integral program
    opts->gaunt_defined = 0;
    opts->breit_defined = 0;
    strcpy(opts->integral_file_gaunt, "MGINT");
    strcpy(opts->integral_file_breit, "MBINT");
    opts->escf = 0.0;
    opts->enuc = 0.0;
    opts->eref = 0.0;
    opts->ground_energy_0h1p = 0.0;
    opts->ground_energy_0h2p = 0.0;
    opts->ground_energy_0h3p = 0.0;
    opts->ground_energy_1h0p = 0.0;
    opts->ground_energy_2h0p = 0.0;
    opts->cc_model = CC_MODEL_CCSD;
    opts->hughes_kaldor_1h2p = 0;
    opts->hughes_kaldor_2h1p = 0;
    opts->sector_h = 0;
    opts->sector_p = 0;
    opts->curr_sector_h = 0;
    opts->curr_sector_p = 0;

    /*
     * flag indicating if we are currently solving lambda equations or not
     */
    opts->curr_in_lambda_equations = 0;

    // active space specification
    // [inact H] act_energy_min [act H] [act P] act_energy_max [inact P]
    // default: no active space
    opts->actsp_defined = CC_ACT_SPACE_UNDEFINED;  // undefined
    opts->act_energy_min = +1e10;
    opts->act_energy_max = -1e10;
    opts->nacth = 0;
    opts->nactp = 0;
    memset(opts->active_specs, 0, sizeof(opts->active_specs));
    memset(opts->active_each, 0, sizeof(opts->active_each));

    // number of electrons in each irrep (for open-shell reference states)
    opts->nelec_defined = 0;
    memset(&opts->irrep_occ_numbers, 0, sizeof(opts->irrep_occ_numbers));

    // occupation numbers for each spinor (for open-shell reference states)
    opts->occ_defined = 0;
    memset(opts->occ, 0, sizeof(opts->occ));

    /*
     * symbolic names of spinors
     */
    memset(opts->spinor_labels, 0, sizeof(char *) * CC_MAX_SPINORS);

    // resort integrals or not
    opts->reuse_integrals_1 = 0;
    opts->reuse_integrals_2 = 0;

    // read precomputed amplitudes from disk (as initial guess) or not
    memset(opts->reuse_amplitudes, 0, sizeof(opts->reuse_amplitudes));

    // flush non-converged amplitudes to disk
    opts->do_flush_iter = 0;

    // denominator shifts
    memset(opts->shifts, 0, sizeof(opts->shifts));
    opts->shift_type = CC_SHIFT_NONE;

    // damping
    memset(opts->damping, 0, sizeof(opts->damping));

    // skip some FS sectors
    memset(opts->skip_sector, 0, sizeof(opts->skip_sector));

    // DIIS -- convergence acceleration
    opts->diis_enabled = 1;
    opts->diis_dim = 10;       // dimension of the DIIS subspace
    opts->diis_triples = 0;   // do DIIS for triples amplitudes or not

    // CROP -- convergence acceleration
    opts->crop_enabled = 0;
    opts->crop_dim = 10;       // dimension of the CROP subspace
    opts->crop_triples = 0;   // do CROP for triples amplitudes or not

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

    // for analytic calculation of properties
    opts->n_analyt_prop = 0;
    for (int i = 0; i < CC_MAX_NPROP; i++) {
        strcpy(opts->analyt_prop_files[i], "\0");
    }

    memset(opts->calc_density, 0, sizeof(opts->calc_density));
    memset(opts->density_num_states, 0, sizeof(opts->density_num_states));
    memset(opts->density_target_states, 0, sizeof(opts->density_target_states));
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            opts->calc_density[h][p] = CC_DENSITY_MATRIX_DISABLED;
        }
    }

    opts->calc_lambda_0h1p = 0;
    opts->lambda_0h1p_num_states = 0;
    memset(opts->lambda_0h1p_states, 0, sizeof(opts->lambda_0h1p_num_states));

    // perform "hermitization" of Heff or not (default:hermitization is enabled)
    // Bloch (non-Hermitian) or des Cloizeaux (Hermitian) Hamiltonian?
    // Bloch => biorthogonal model vectors => TDM(i->f) != TDM(f->i) =>
    //          not very physical
    // des Cloizeaux => (Loewdin) orthogonal model vectors =>
    //                  much more meaningful TDMs and properties
    opts->do_hermit = 0;

    // perform Dipole-Length (DL) estimation of TDMs or not
    opts->calc_model_space_tdms = 0;

    // calculation of density matrices and natural orbitals
    opts->n_denmat = 0;

    // calculation of properties (model-space estimation)
    opts->n_model_space_props = 0;

    // selection of amplitudes
    opts->n_select = 0;

    // (simple) intermediate hamiltonian ("IH-1") parameters
    opts->do_intham_imms = 0;
    memset(&opts->intham_imms_opts, 0, sizeof(opts->intham_imms_opts));
    opts->intham_imms_opts.scale_shift = 1.0;
    opts->intham_imms_opts.shift_type = CC_SHIFT_REALIMAG;
    opts->intham_imms_opts.npower = 5;
    for (int i = 0; i < MAX_SECTOR_RANK; i++) {
        for (int j = 0; j < MAX_SECTOR_RANK; j++) {
            opts->intham_imms_opts.det_shift_auto[i][j] = IH_IMMS_FRONTIER_ENERGY_UPPER_BOUND;
        }
    }

    // restriction of the spinor space for triples
    opts->do_restrict_t3 = 0;
    opts->restrict_t3_bounds[0] = 0.0;
    opts->restrict_t3_bounds[1] = 0.0;

    /*
     * a flag detecting whether cluster amplitudes should be flushed
     * into formatted txt files or not
     */
    opts->do_flush_amplitudes_txt = 0;

    /*
     * calculate overlap integrals for target wavefunctions
     * (in the given sector)
     */
    memset(opts->calc_overlap, 0, sizeof(opts->calc_overlap));

    /*
     * data from the MRCONEE (DIRAC) file
     */
    opts->mrconee_data = NULL;

    /*
     * do not recalculate orbital energies
     */
    opts->use_oe = 0;

    /*
     * use tensor trains in different situations
     */
    opts->tt_options.tt_module_enabled = 0;
    opts->tt_options.use_tensor_trains = 1;
    opts->tt_options.svd_thresh = 1e-9;
    opts->tt_options.cholesky_thresh = 1e-6;
    opts->tt_options.use_pyscf_integrals = 0;
    opts->tt_options.pyscf_integrals_path[0] = '\0';

    /*
     * use Goldstone formalism
     */
    opts->use_goldstone = 0;

    return opts;
}


/*
 * deallocates memory owned by the object of type 'cc_options_t'
 */
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


/*
 * print options
 */
void print_options(cc_options_t *opts)
{
    printf("\n");
    printf("\t\t\t\t*************\n");
    printf("\t\t\t\t** Options **\n");
    printf("\t\t\t\t*************\n\n");
    printf(" %-15s  %-40s  %s\n", "title", "title string - comment", opts->title);
    printf(" %-15s  %-40s  %s\n", "scratch_dir", "scratch directory for tmp files", opts->scratch_dir);
    printf(" %-15s  %-40s  %s\n", "--no-clean", "retain scratch directory on exit", opts->clean_scratch ? "no" : "yes");

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

    printf(" %-15s  %-40s  %s\n", "flush_amplitude", "write formatted files with cluser ampl-s",
           opts->do_flush_amplitudes_txt ? "yes" : "no");

    printf(" %-15s  %-40s  %s\n", "arith", "recommended arithmetic",
           opts->recommended_arith == CC_ARITH_REAL ? "real" : "complex");
    printf(" %-15s  %-40s  %.1f Mb\n", "memory", "max allowed RAM usage", opts->max_memory_size / (1024.0 * 1024.0));
    printf(" %-15s  %-40s  %s\n", "compress", "compression of integrals on disk", opts->compress ? "LZ4" : "disabled");
    if (opts->do_compress_triples) {
        printf(" %-15s %-40s  yes, thresh=%g, datatype=%s\n", "compress_triples",
               "compression of triples tensors in RAM",
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
    printf(" %-15s  %-40s  %s\n", "openmp_algorithm", "parallelization algorithm for mult",
           opts->openmp_algorithm == CC_OPENMP_ALGORITHM_EXTERNAL ? "external" : "internal");
    printf(" %-15s  %-40s  %s\n", "cuda", "calculations on GPU (CUDA)", opts->cuda_enabled ? "enabled" : "disabled");
    printf(" %-15s  %-40s  %d\n", "maxiter", "maximum number of CC iterations", opts->maxiter);
    printf(" %-15s  %-40s  %g\n", "conv_thresh", "convergence threshold (by amplitudes)", opts->conv_thresh);
    printf(" %-15s  %-40s  %g\n", "div_thresh", "divergence threshold (by amplitudes)", opts->div_thresh);

    printf(" %-15s  %-40s  ", "reuse", "reuse amplitudes and/or integrals");
    int num_reused = 0;
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            if (opts->reuse_amplitudes[h][p]) {
                printf("%dh%dp ", h, p);
                num_reused++;
            }
        }
    }
    if (opts->reuse_integrals_1) {
        printf("1-integrals ");
        num_reused++;
    }
    if (opts->reuse_integrals_2) {
        printf("2-integrals ");
        num_reused++;
    }
    if (num_reused == 0) {
        printf("nothing\n");
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
           opts->int_source == CC_INTEGRALS_DIRAC ? "DIRAC" : "PySCF");

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
    printf(" %-15s  %-40s  ", "new_sorting", "new implementation of integral sorting");
    if (opts->new_sorting) {
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
        printf("%g <= eps <= %g a.u.\n", opts->act_energy_min, opts->act_energy_max);
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

    printf(" %-15s  %-40s  %s\n", "goldstone", "use Goldstone formalism", opts->use_goldstone ? "yes" : "no");

    printf(" %-15s  %-40s  ", "shift_type", "formula for denominator shifts");
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

    printf(" %-15s  %-40s  %.5e\n", "degen_thresh", "degeneracy threshold (a.u.)", opts->degen_thresh);
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

    if (opts->crop_enabled) {
        printf(" %-15s  %-40s  %s\n", "crop", "CROP technique for convergence", "enabled");
        printf(" %-15s  %-40s  %d\n", "crop <n>", "CROP subspace dimension", opts->crop_dim);
        printf(" %-15s  %-40s  %s\n", "crop triples", "CROP for triples amplitudes",
               opts->crop_triples ? "enabled" : "disabled");
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

    // density matrix in the 0h0p sector
    printf(" %-15s  %-40s  ", "density 0h0p", "construct density matrix in 0h0p");

    if (opts->calc_density[0][0] == CC_DENSITY_MATRIX_DISABLED) {
        printf("disabled\n");
    }
    else if (opts->calc_density[0][0] == CC_DENSITY_MATRIX_LAMBDA) {
        printf("lambda equations\n");
    }
    else if (opts->calc_density[0][0] == CC_DENSITY_MATRIX_EXPECTATION) {
        printf("expectation value\n");
    }

    // density matrix calculations
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {

            char sector_label[16];
            char keyword[16];
            char explanation[64];
            sprintf(sector_label, "%dh%dp", h, p);
            sprintf(keyword, "density %s", sector_label);
            sprintf(explanation, "construct density matrix in %s", sector_label);

            printf(" %-15s  %-40s  ", keyword, explanation);
            if (opts->calc_density[h][p] == CC_DENSITY_MATRIX_DISABLED) {
                printf("disabled\n");
            }
            else if (opts->calc_density[h][p] == CC_DENSITY_MATRIX_LAMBDA) {
                if (opts->density_num_states[h][p] == -1) {
                    printf("for all target states\n");
                }
                else {
                    printf("for states");
                    for (int i = 0; i < opts->density_num_states[h][p]; i++) {
                        cc_denmat_query_t *dm_query = opts->density_target_states[h][p] + i;
                        if (dm_query->state1 == dm_query->state2 &&
                            strcmp(dm_query->rep1_name, dm_query->rep2_name) == 0) {
                            printf(" %s:%d", dm_query->rep1_name, dm_query->state1 + 1);
                        }
                        else {
                            printf(" %s:%d-%s:%d", dm_query->rep1_name, dm_query->state1 + 1,
                                   dm_query->rep2_name, dm_query->state2 + 1);
                        }
                    }
                    printf("\n");
                }
            }
            else if (opts->calc_density[h][p] == CC_DENSITY_MATRIX_EXPECTATION) {
                printf("expectation value\n");
            }

        }
    }

    // lambda equations & exact analytic density matrix in the 0h1p sector
    printf(" %-15s  %-40s  ", "lambda 0h1p", "solve lambda equations in the 0h1p sector");
    if (opts->calc_lambda_0h1p == 0) {
        printf("disabled\n");
    }
    else {
        printf("for states");
        for (int i = 0; i < opts->lambda_0h1p_num_states; i++) {
            printf(" %s:%d", opts->lambda_0h1p_states[i].rep1_name, opts->lambda_0h1p_states[i].state1 + 1);
        }
        printf("\n");
    }

    // overlap integrals in non-trivial sectors
    printf(" %-15s  %-40s  ", "overlap", "calculate overlap int-s for target wfns");
    int n_overlap_sectors = 0;
    for (int h = 0; h < MAX_SECTOR_RANK; h++) {
        for (int p = 0; p < MAX_SECTOR_RANK; p++) {
            if (opts->calc_overlap[h][p] != 0) {
                n_overlap_sectors += 1;
                printf("%dh%dp ", h, p);
            }
        }
    }
    if (n_overlap_sectors == 0) {
        printf("disabled");
    }
    printf("\n");

    // perform hermitization or not?
    printf(" %-15s  %-40s  %s\n", "nohermit", "hermitization of effective Hamiltonians",
           opts->do_hermit ? "enabled" : "disabled");

    // perform Dipole-Length (DL) estimation of TDMs or not
    printf(" %-15s  %-40s  %s\n", "dltdm", "model-space estimates of tran dipoles",
           opts->calc_model_space_tdms ? "enabled" : "disabled");

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

    if (opts->n_model_space_props > 0) {
        printf(" %-15s  %-40s  %s", "prop", "model-space estimates of properties", "enabled:");
        for (int i = 0; i < opts->n_model_space_props; i++) {
            cc_ms_prop_query_t *query = opts->prop_queries + i;

            if (query->source == CC_PROP_FROM_MDPROP) {
                printf(" %s (md)", query->prop_name);
            }
            else {
                printf(" %s %s (txt)", query->file_real, query->file_imag);
            }
            if (query->do_transpose) {
                printf("^T");
            }
            if (!(strcmp(query->irrep_name, "") == 0)) {
                printf(" %s", query->irrep_name);
            }
            if (query->approx_numerator > 0 || query->approx_denominator > 0) {
                printf(" %d/%d", query->approx_numerator, query->approx_denominator);
            }
            printf(";");
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
        printf(" %-15s  %-40s  enabled, %g < eps < %g\n", "restrict_t3", "restriction of triples",
               opts->restrict_t3_bounds[0], opts->restrict_t3_bounds[1]);
    }
    else {
        printf(" %-15s  %-40s  %s\n", "restrict_t3", "restriction of triples", "disabled");
    }

    /*
     * do not recalculate orbital energies
     */
    printf(" %-15s  %-40s  %s\n", "use_oe", "use orbital energies from DIRAC", opts->use_oe ? "enabled" : "disabled");

    // (simple) intermediate hamiltonian ("IH-1") parameters
    printf(" %-15s  %-40s  ", "ih_imms", "simple intermediate Hamiltonian");
    if (opts->do_intham_imms) {
        printf("enabled in sectors: ");
        for (int h = 0; h < MAX_SECTOR_RANK; h++) {
            for (int p = 0; p < MAX_SECTOR_RANK; p++) {
                if (opts->intham_imms_opts.sectors[h][p]) {
                    printf("%dh%dp ", h, p);
                }
            }
        }
        printf("\n");
    }
    else {
        printf("disabled\n");
    }

    /*
     * use tensor trains in different situations
     * if the TT library was linked to EXP-T
     */
#ifdef TENSOR_TRAIN
    printf(" %-15s  %-40s  %s\n", "tensor_trains", "tensor-train-based CCSD in the 0h0p sector", opts->tt_options.tt_module_enabled ? "enabled" : "disabled");
#endif // TENSOR_TRAIN

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


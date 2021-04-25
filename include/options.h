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
* options.h
* =========
*
* Options: all data read from the input file.
* See the code for detailed description of every option.
*
* 2018-2021 Alexander Oleynichenko
******************************************************************************/

#ifndef CC_OPTIONS_H_INCLUDED
#define CC_OPTIONS_H_INCLUDED

#include <stdbool.h>
#include <complex.h>

#include "comdef.h"
#include "intham.h"
#include "linalg.h"
#include "selection.h"

#define CC_ARITH_REAL    0
#define CC_ARITH_COMPLEX 1

// maximum number of queues for properties calculations
#define CC_MAX_NPROP 256
#define CC_MAX_SELECTION 64

// boolean constants
enum {
    CC_DISABLED = 0,
    CC_ENABLED = 1
};

// print level
enum {
    CC_PRINT_LOW,
    CC_PRINT_MEDIUM,
    CC_PRINT_HIGH,
    CC_PRINT_DEBUG
};

// integral source
typedef enum {
    CC_INTEGRALS_DIRAC,
    CC_INTEGRALS_TM2C
} cc_interface_t;

// source of property matrices
enum {
    CC_PROP_FROM_MDPROP,
    CC_PROP_FROM_TXTPROP
};

// level of disk usage
enum {
    CC_DISK_USAGE_LEVEL_0,  // all data in RAM
    CC_DISK_USAGE_LEVEL_1,  // rank 6+ diagrams on disk (default)
    CC_DISK_USAGE_LEVEL_2,  // rank 6+, pppp on disk
    CC_DISK_USAGE_LEVEL_3,  // rank 6+, pppp, *ppp on disk
    CC_DISK_USAGE_LEVEL_4   // rank 6+, pppp, *ppp on disk + compression enabled
};

typedef enum {
    CC_MODEL_CCS,
    CC_MODEL_CCD,
    CC_MODEL_CCSD,
    CC_MODEL_CCSD_T3,        // CCSD+T(3)
    CC_MODEL_CCSD_T3_STAR,   // CCSD+T*(3)
    CC_MODEL_CCSD_T4,        // CCSD+T(4)
    CC_MODEL_CCSD_T4_STAR,   // CCSD+T*(4)
    CC_MODEL_CCSDT_1A,       // CCSDT-1a
    CC_MODEL_CCSDT_1B,       // CCSDT-1b
    CC_MODEL_CCSDT_1B_PRIME, // CCSDT-1b'
    CC_MODEL_CCSDT_2,        // CCSDT-2
    CC_MODEL_CCSDT_3,        // CCSDT-3
    CC_MODEL_CCSDT
} cc_model_t;

typedef enum {
    CC_SHIFT_NONE,
    CC_SHIFT_REAL,
    CC_SHIFT_REALIMAG,
    CC_SHIFT_IMAG,
    CC_SHIFT_TAYLOR
} cc_shifttype_t;

// damping parameters
typedef struct {
    bool enabled;   // flag
    int stop;      // stop damping at this step
    double factor;  // damping factor
} cc_damping_params_t;

// dynamic denominator shifts
typedef struct {
    bool enabled;         // flag
    cc_shifttype_t type;  // formula of the shift: real, taylor, etc
    int power;           // "attenuation parameter"
    double shifts[CC_DIAGRAM_MAX_RANK];  // shifts for singles, doubles, triples, ...
} cc_shift_params_t;

// for active space specification
typedef struct {
    int nacth;
    int nactp;
    char rep_name[32];
} cc_active_spec_t;

enum {
    CC_SPACE_BY_IRREPS,
    CC_SPACE_BY_ENERGY,
    CC_SPACE_BY_TOTAL
};

enum {
    CC_ACT_SPACE_UNDEFINED = 0,
    CC_ACT_SPACE_SPEC_ENERGY,
    CC_ACT_SPACE_SPEC_TOTAL,
    CC_ACT_SPACE_SPEC_IRREPS,
    CC_ACT_SPACE_SPEC_BINARY
};

// for definition of subspaces (of orbitals or electronic states)
typedef struct {
    int deftype;
    int total;
    double emin;
    double emax;
    int dim[CC_MAX_NUM_IRREPS];
    char rep_names[CC_MAX_NUM_IRREPS][32];
} cc_space_t;

// for calculations of density matrices and natural orbitals (including transition ones)
typedef struct {
    int sect1[2];
    int sect2[2];
    char rep1_name[32];
    char rep2_name[32];
    int state1;
    int state2;
} cc_denmat_query_t;

// for calculation of properties
typedef struct {
    int  source;
    char prop_name[32];
    char file_real[CC_MAX_FILE_NAME_LENGTH];
    char file_imag[CC_MAX_FILE_NAME_LENGTH];
    int  do_transpose;     // optional transposition
} cc_ms_prop_query_t;

// data compression algorithm
typedef enum {
    CC_COMPRESS_NONE,
    CC_COMPRESS_LZ4
} cc_compression_type_t;

typedef struct {
    cc_space_t main_space;
} intham_param_t;

// sizeof(double) or sizeof(double complex)
extern int SIZEOF_WORKING_TYPE;

extern data_type_t WORKING_TYPE;

struct cc_options {

    // title (suggested by user)
    char title[CC_MAX_TITLE];

    // scratch directory
    char scratch_dir[CC_MAX_PATH_LENGTH];

    // remove tmp files from scratch directory on exit
    int clean_scratch;

    // print level
    int print_level;

    // recommended arithmetic
    int recommended_arith;

    // max memory usage (bytes)
    size_t max_memory_size;

    // data compression
    cc_compression_type_t compress;

    // orbital range tiling
    int tile_size;

    // level of disk usage (all in RAM, rank6 on disk, pppp on disk, etc)
    int disk_usage_level;

    // number of lightweight (OpenMP) threads
    int nthreads;

    // is CUDA enabled or not
    int cuda_enabled;

    // convergence options
    int maxiter;    // max number of iterations
    double conv;    // convergence threshold (by amplitudes)
    double div_thresh; // divergence threshold (by amplitudes)

    // CC model: CCSD, CCSD-T(3), CCSDT-1, etc
    cc_model_t cc_model;
    char cc_model_str[32];

    // target Fock-space sector
    int sector_h;
    int sector_p;
    int mixed;

    // current Fock-space sector
    int curr_sector_h;
    int curr_sector_p;

    // active space specification
    // [inact H] actsp_min [act H] [act P] actsp_max [inact P]
    // default: no active space
    int actsp_defined;  // 1 if via min/max, 2 if via nacth/p, 3 if via active_specs, 4 via 1/0 numbers
    double actsp_min;
    double actsp_max;
    int nacth;
    int nactp;
    int main_defined;
    int main_h;
    int main_p;

    cc_active_spec_t active_specs[CC_MAX_NUM_IRREPS];
    int active_each[CC_MAX_SPINORS];

    // number of electrons in each irrep -- will be useful for open-shell references
    int nelec_defined;
    int nelec[CC_MAX_NUM_IRREPS];
    // occupation numbers of each spinor -- for open-shell references
    int occ_defined;
    char occ[CC_MAX_SPINORS];

    // interface to molecular integrals
    cc_interface_t int_source;   // source of MO integrals: DIRAC, etc
    char integral_file_1[CC_MAX_PATH_LENGTH];  // name of the file with 1-el MO integrals and basic info
    char integral_file_2[CC_MAX_PATH_LENGTH];  // name of the file with 2-el MO integrals
    char integral_file_prop[CC_MAX_PATH_LENGTH];  // name of the file with properties MO integrals
    int x2cmmf;    // x2c molecular mean field Hamiltonian enabled in DIRAC
    // for Daniel Maison integral program:
    int gaunt_defined;
    int breit_defined;
    char integral_file_gaunt[CC_MAX_PATH_LENGTH];
    char integral_file_breit[CC_MAX_PATH_LENGTH];

    // to avoid resorting of integrals
    int reuse_integrals_1;
    int reuse_integrals_2;

    // precomputed cluster amplitudes as initial guess
    int reuse_0h0p;
    int reuse_0h1p;
    int reuse_1h0p;
    int reuse_1h1p;
    int reuse_0h2p;
    int reuse_2h0p;
    int reuse_0h3p;

    // flush non-converged amplitudes to disk
    int do_flush_iter;

    // SCF energy
    double escf;

    // nuclear repulsion energy
    double enuc;

    // CCSD(0,0) energy
    double eref;

    // denominator shifts
    cc_shifttype_t shift_type;
    cc_shift_params_t shifts[MAX_SECTOR_RANK][MAX_SECTOR_RANK];

    int do_orbshift;
    int orbshift_power;
    double orbshift;

    int do_orbshift_0h0p;
    int orbshift_0h0p_power;
    double orbshift_0h0p[32];

    // DIIS -- convergence acceleration
    int diis_enabled;
    int diis_dim;       // dimension of the DIIS subspace
    int diis_triples;   // do DIIS for triples amplitudes or not

    // damping
    cc_damping_params_t damping[MAX_SECTOR_RANK][MAX_SECTOR_RANK];

    // skip calculations in some sectors
    int skip_sector[MAX_SECTOR_RANK][MAX_SECTOR_RANK];

    // how many roots should be processed
    int nroots_specified;  // is 'nroots' option was specified in input file
    cc_space_t nroots_specs;  // number of lowest roots in each symmetry

    // degeneracy threshold (used for representing results)
    double degen_thresh;

    // interface to the OneProp code by L. V. Skripnikov (one-elec properties)
    int oneprop_on;
    int n_oneprop;
    char oneprop_file_re[CC_MAX_NPROP][CC_MAX_PATH_LENGTH];
    char oneprop_file_im[CC_MAX_NPROP][CC_MAX_PATH_LENGTH];
    double complex oneprop_lambda[CC_MAX_NPROP];

    // interface to two-electron properties by D. E. Maison and L. V. Skripnikov
    int twoprop_on;
    int n_twoprop;
    char twoprop_file[CC_MAX_NPROP][CC_MAX_PATH_LENGTH];
    double complex twoprop_lambda[CC_MAX_NPROP];

    // properties from the MDPROP file
    int n_mdprop;
    char mdprop_file[CC_MAX_NPROP][CC_MAX_PATH_LENGTH];
    double complex mdprop_lambda[CC_MAX_NPROP];

    // perform Dipole-Length (DL) estimation of TDMs or not
    int do_diplen_tdm;

    // for calculation of density matrices and natural orbitals
    int do_hermit;      // perform "hermitization" of Heff or not
    int n_denmat;
    cc_denmat_query_t denmat_query[CC_MAX_NPROP];

    // for model-space estimation of properties
    int n_ms_props;
    cc_ms_prop_query_t prop_queries[CC_MAX_NPROP];

    // selection of amplitudes
    // (does not affect performance, just for fast implementation of new approximate schemes)
    int n_select;
    ampl_selection_t selects[CC_MAX_SELECTION];

    // intermediate hamiltonian parameters
    int do_intham;
    intham_param_t intham_params;
};

typedef struct cc_options cc_options_t;

// and general options:

// arithmetic: CC_ARITH_REAL or CC_ARITH_COMPLEX
extern int arith;

// complex arithmetic?
extern int carith;

extern cc_options_t *cc_opts;

// helper functions
cc_options_t *new_options();

void delete_options(cc_options_t *);

void print_options(cc_options_t *);

#endif /* CC_OPTIONS_H_INCLUDED */

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
 * dirac_interface.c
 * =================
 *
 * Implementation of the interface to the DIRAC program package.
 * Reads info about spinors and transformed MO integrals from the unformatted
 * file produced by DIRAC (with default names MRCONEE, MDCINT, MDPROP)
 *
 * 2018-2021 Alexander Oleynichenko
 ******************************************************************************/

#include "interfaces.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "engine.h"
#include "error.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"

#define MAX_INTS_FILE_NAME_LENGTH 1024

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// arguments to be passed to the Fortran code
// (via the 'mrconee_mdcint' common block, see dirac_binary.f90)
extern struct {
    char oneel_file[MAX_INTS_FILE_NAME_LENGTH];
    char twoel_file[MAX_INTS_FILE_NAME_LENGTH];
    char prop_file[MAX_INTS_FILE_NAME_LENGTH];
    char gaunt_file[MAX_INTS_FILE_NAME_LENGTH];
    int64_t gaunt_enabled;
    int64_t n_twoprop;
    char twoprop_files[CC_MAX_NPROP][MAX_INTS_FILE_NAME_LENGTH];
} mrconee_mdcint;

// data which are returned by the Fortran code;
// see common /dirac_data/ in dirac_binary.f90
extern struct {
    int32_t nspinors;        // total number of spinors
    int32_t breit;           // breit active in DHF
    double enuc;             // core energy (inactive energy + nuclear repulsion)
    int32_t invsym;          // inversion symmetry (yes : 2; no : 1)
    int32_t nz_arith;        // group type (1 real, 2 complex, 4 quaternion)
    int32_t is_spinfree;     // spinfree formalism
    int32_t norb_total;      // total number of orbitals (so including frozen or deleted orbitals)
    double escf;             // total SCF energy (from DIRAC)
    int32_t nsymrp;          // number of fermion irreps in parent group
    char repnames[112];      // their names (type: string)
    int32_t nactive[8];      // number of occupied spinors in each irrep active in transf-n
    int32_t nstr[2];         // total number of orbitals of this ircop
    int32_t nfrozen[6];      // number of occupied frozen (core) spinors
    int32_t ndelete[8];      // number of deleted spinors
    int32_t nsymrpa;         // number of fermion irreps in the Abelian subgroup
    char repanames[256];     // names of these irreps
    int32_t multb[64][64];   // multiplication table for direct products in the Abelian subgroup
    int32_t nbsymrp;         // (not used) number of boson symmetry reps (for LUCITA)
    int32_t irpmo[CC_MAX_SPINORS];   // Irrep in parent group (1:gerade, 2:ungerade)
    int32_t irpamo[CC_MAX_SPINORS];  // Irrep in Abelian subgroup
    int32_t ibspi[CC_MAX_SPINORS];   // Approximate boson irrep identification (needed in LUCIAREL and LUCITA)
    int32_t iocc[CC_MAX_SPINORS];    // spinor occupation numbers
    double eorbmo[CC_MAX_SPINORS];   // Orbital energy taken from DHF
    int32_t place4;                  // placeholder for alignment
} dirac_data;

extern struct {
    int32_t tile_size;
} tilesize;

extern struct {
    int32_t recommended_carith;
} recommendedcarith;

extern struct {
    int32_t skip_mdcint;
} skip_sorting;

void invoke_dirac_interface_binary(cc_options_t *opts);
void dirac_interface_binary(int *err);
void dirac_interface_debug_print();
char **dirac_interface_extract_rep_names(int nsym);
void detect_dirac_point_group(int nsym, char **rep_names, char *group_name, int *fully_sym_irrep);
void rename_irreps_dirac_to_expt(int nsym, char **rep_names);
void extract_direct_product_table_dirac(int nsym, int *mult_table);
void strcpy_to_fortran(size_t dest_len, char *dest, char *src);
void create_spinor_info_dirac(int32_t n_spinors, int32_t *irpamo_array, double *spinor_energies, int32_t *occ_numbers);
void check_irrep_names();
int rep_name_exists(int nrep, char **list, char *query);


/**
 * Reads info about spinors and transformed MO integrals from formatted file.
 * NOTE: no error checking at all!
 */
void dirac_interface(cc_options_t *opts)
{
    int point_group_type;
    int nsym;
    int irrep_a1;
    char point_group_name[MAX_IRREP_NAME];
    int *mult_table = NULL;

    // call the Fortran subroutine
    invoke_dirac_interface_binary(opts);

    // extract data from the "Fortran common block"
    // and reinterpret them

    opts->enuc = dirac_data.enuc;
    opts->escf = dirac_data.escf;

    // nz = arithmetic = 1:real, 2:complex, 4:quaternion
    int nz = dirac_data.nz_arith;
    if (nz == 1) {
        point_group_type = CC_GROUP_REAL;
    }
    else if (nz == 2) {
        point_group_type = CC_GROUP_COMPLEX;
    }
    else {
        point_group_type = CC_GROUP_QUATERNION;
    }

    if (nz == 1 || dirac_data.is_spinfree) {
        arith = CC_ARITH_REAL;
        carith = 0;
    }
    else {
        arith = CC_ARITH_COMPLEX;
        carith = 1;
    }
    // what about arithmetic recommended by the user ?
    if (opts->recommended_arith == CC_ARITH_COMPLEX) {
        carith = 1;
    }
    WORKING_TYPE = (carith == 1) ? CC_DOUBLE_COMPLEX : CC_DOUBLE;
    SIZEOF_WORKING_TYPE = (carith == 1) ? sizeof(double complex) : sizeof(double);

    if (opts->print_level >= CC_PRINT_DEBUG) {
        if (dirac_data.is_spinfree) { printf(" Spin-orbit: off (spinfree)\n"); }
        else { printf(" Spin-orbit: on\n"); }
    }

    // Number of double group irreps
    nsym = dirac_data.nsymrpa * 2;

    // Representation names
    char **rep_names = dirac_interface_extract_rep_names(nsym);
    detect_dirac_point_group(nsym, rep_names, point_group_name, &irrep_a1);
    rename_irreps_dirac_to_expt(nsym, rep_names);

    // "Multiplication table" ( = dir prod table for Abelian groups).
    // for (quasi) infinite groups Cinfv and Dinfh multiplication tables provided by DIRAC
    // are not sufficient for CCSDT and/or 0h3p calculations (too few irreps). So for these
    // groups additional irreps will be generated

    if (strcmp(point_group_name, "Cinfv") == 0 || strcmp(point_group_name, "Dinfh") == 0) {
        int num_irreps;
        int max_omega_x2;
        char **new_irrep_names = NULL;

        if (strcmp(point_group_name, "Cinfv") == 0) { // Cinfv
            max_omega_x2 = CC_MAX_NUM_IRREPS / 2 - 2;
            new_irrep_names = generate_irreps_Cinfv(max_omega_x2, &num_irreps);
            irrep_a1 = search_string("0", new_irrep_names, num_irreps);
            mult_table = construct_direct_product_table(num_irreps, new_irrep_names, multiply_irreps_Cinfv);
        }
        else { // Dinfh
            max_omega_x2 = CC_MAX_NUM_IRREPS / 4 - 2;
            new_irrep_names = generate_irreps_Dinfh(max_omega_x2, &num_irreps);
            irrep_a1 = search_string("0g", new_irrep_names, num_irreps);
            mult_table = construct_direct_product_table(num_irreps, new_irrep_names, multiply_irreps_Dinfh);
        }

        printf("\n");
        printf(" %s group is to be extended:\n", point_group_name);
        printf(" max 2x|Omega| = %d\n", max_omega_x2);
        printf(" number of irreps = %d\n", num_irreps);
        printf("\n");

        // for each spinor: update irrep enumeration
        for (int i = 0; i < dirac_data.nspinors; i++) {
            int old_irrep_no = dirac_data.irpamo[i];
            char *old_irrep_name = rep_names[old_irrep_no];
            int new_irrep_no = search_string(old_irrep_name, new_irrep_names, num_irreps);
            dirac_data.irpamo[i] = new_irrep_no;
        }

        rep_names = new_irrep_names;
        nsym = num_irreps;
    }
    else { // all groups except for Cinfv and Dinfh
        mult_table = (int *) cc_malloc(nsym * nsym * sizeof(int));
        extract_direct_product_table_dirac(nsym, mult_table);
    }

    setup_symmetry(point_group_type, point_group_name, nsym, rep_names, irrep_a1, mult_table);

    // information about spinors and active space
    create_spinor_info_dirac(dirac_data.nspinors, dirac_data.irpamo, dirac_data.eorbmo, dirac_data.iocc);
    create_spinor_blocks(opts->tile_size);
    setup_occupation_numbers(opts, dirac_data.nspinors, spinor_info);
    setup_active_space(opts);
    setup_fast_access_spinor_lists();

    print_symmetry_info();
    print_spinor_info_table();

    // additional checks:
    // 1. names of irreps in some directives of the input file can be incorrect
    check_irrep_names();
}


void invoke_dirac_interface_binary(cc_options_t *opts)
{
    int err;

    // prepare Fortran-compatible names of integral files
    strcpy_to_fortran(MAX_INTS_FILE_NAME_LENGTH, mrconee_mdcint.oneel_file, opts->integral_file_1);
    strcpy_to_fortran(MAX_INTS_FILE_NAME_LENGTH, mrconee_mdcint.twoel_file, opts->integral_file_2);
    strcpy_to_fortran(MAX_INTS_FILE_NAME_LENGTH, mrconee_mdcint.prop_file, opts->integral_file_prop);
    strcpy_to_fortran(MAX_INTS_FILE_NAME_LENGTH, mrconee_mdcint.gaunt_file, opts->integral_file_gaunt);
    for (int i = 0; i < CC_MAX_NPROP; i++) {
        strcpy_to_fortran(MAX_INTS_FILE_NAME_LENGTH, mrconee_mdcint.twoprop_files[i], opts->twoprop_file[i]);
    }
    mrconee_mdcint.gaunt_enabled = opts->gaunt_defined ? 1 : 0;
    mrconee_mdcint.n_twoprop = opts->n_twoprop;

    // get tile size and pass it to the Fortran side
    tilesize.tile_size = opts->tile_size;
    skip_sorting.skip_mdcint = opts->reuse_integrals_2;
    recommendedcarith.recommended_carith = (opts->recommended_arith == CC_ARITH_COMPLEX) ? 1 : 0;

    fflush(stdout);
    timer_new_entry("dirac", "DIRAC interface (MRCONEE/MDCINT)");
    timer_start("dirac");
    dirac_interface_binary(&err);
    timer_stop("dirac");
    if (err == 1) {
        errquit("ERROR IN THE BINARY INTERFACE\n");
    }

    if (opts->print_level >= CC_PRINT_DEBUG) {
        dirac_interface_debug_print();
    }
}


void strcpy_to_fortran(size_t dest_len, char *dest, char *src)
{
    for (size_t i = 0; i < dest_len; i++) {
        dest[i] = ' ';
    }
    for (size_t i = 0; i < strlen(src); i++) {
        dest[i] = src[i];
    }
}


void detect_dirac_point_group(int nsym, char **rep_names, char *group_name, int *fully_sym_irrep)
{
    if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "A  b") == 0) {
        *fully_sym_irrep = 4;
        strcpy(group_name, "C1");
    }
    else if (strcmp(rep_names[0], "Ag a") == 0 && strcmp(rep_names[1], "Au a") == 0) {
        *fully_sym_irrep = 8;
        strcpy(group_name, "Ci");
    }
    else if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "B  a") == 0) {
        *fully_sym_irrep = 8;
        strcpy(group_name, "C2");
    }
    else if (strcmp(rep_names[0], "A' a") == 0 && strcmp(rep_names[1], "A\" a") == 0) {
        *fully_sym_irrep = 8;
        strcpy(group_name, "Cs");
    }
    else if (strcmp(rep_names[0], "A1 a") == 0) {
        *fully_sym_irrep = 16;
        strcpy(group_name, "C2v");
    }
    else if (strcmp(rep_names[0], "A  a") == 0) {
        *fully_sym_irrep = 16;
        strcpy(group_name, "D2");
    }
    else if (strcmp(rep_names[0], "Ag a") == 0 && strcmp(rep_names[1], "Bg a") == 0) {
        *fully_sym_irrep = 16;
        strcpy(group_name, "C2h");
    }
    else if (strcmp(rep_names[0], "Ag a") == 0) {
        *fully_sym_irrep = 32;
        strcpy(group_name, "D2h");
    }
        // double groups
    else if (strcmp(rep_names[0], "   A") == 0 && strcmp(rep_names[1], "   a") == 0) {
        *fully_sym_irrep = 1;
        strcpy(group_name, "C1");
    }
    else if (strcmp(rep_names[0], "  AG") == 0 && strcmp(rep_names[1], "  AU") == 0) {
        *fully_sym_irrep = 2;
        strcpy(group_name, "Ci");
    }
    else if (strcmp(rep_names[0], "  1E") == 0 && strcmp(rep_names[1], "  2E") == 0) {
        *fully_sym_irrep = 2;
        strcpy(group_name, "C2, Cs, C2v or D2");
    }
    else if (strcmp(rep_names[0], " 1Eg") == 0 && strcmp(rep_names[1], " 2Eg") == 0) {
        *fully_sym_irrep = 4;
        strcpy(group_name, "C2h or D2h");
    }
    else if (strcmp(rep_names[0], "   1") == 0 && strcmp(rep_names[1], "  -1") == 0) {
        *fully_sym_irrep = 32;
        strcpy(group_name, "Cinfv");
    }
    else if (strcmp(rep_names[0], "  1g") == 0 && strcmp(rep_names[1], " -1g") == 0) {
        *fully_sym_irrep = 32;
        strcpy(group_name, "Dinfh");
    }
    else {
        *fully_sym_irrep = 0;
        strcpy(group_name, "undetected");
    }
}


/**
 * Changes names of irreps for DIRAC to EXP-T notation
 * in order to make them more readable.
 *
 * rules for nonrel groups:
 *  left -- L. Visscher's notation, right -- A. Oleynichenko's notation
 * (Ms projection)
 * a -> a
 * b -> b
 * 3 -> -3/2
 * 3 -> +3/2
 * 0 -> 0
 * 4 -> 2
 * 2 -> +1
 * 2 -> -1
 */
void rename_irreps_dirac_to_expt(int nsym, char **rep_names)
{
    // C1 nonrel
    if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "A  b") == 0) {
        char *translation[] = {
                "A_a", "A_b", "A_-3/2", "A_+3/2", "A_0", "A_2", "A_+1", "A_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // C2 nonrel
    if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "B  a") == 0) {
        char *translation[] = {
                "A_a", "B_a",
                "A_b", "B_b",
                "A_-3/2", "B_-3/2",
                "A_+3/2", "B_+3/2",
                "A_0", "B_0",
                "A_2", "B_2",
                "A_+1", "B_+1",
                "A_-1", "B_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // Cs nonrel
    if (strcmp(rep_names[0], "A' a") == 0 && strcmp(rep_names[1], "A\" a") == 0) {
        char *translation[] = {
                "A'_a", "A\"_a",
                "A'_b", "A\"_b",
                "A'_-3/2", "A\"_-3/2",
                "A'_+3/2", "A\"_+3/2",
                "A'_0", "A\"_0",
                "A'_2", "A\"_2",
                "A'_+1", "A\"_+1",
                "A'_-1", "A\"_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // Ci nonrel
    if (rep_name_exists(nsym, rep_names, "Ag a") &&
        rep_name_exists(nsym, rep_names, "Au a") &&
        rep_name_exists(nsym, rep_names, "Ag b") &&
        rep_name_exists(nsym, rep_names, "Au b") &&
        !rep_name_exists(nsym, rep_names, "Bg a") &&   // not C2h
        !rep_name_exists(nsym, rep_names, "B3ua")      // not D2h
            ) {
        char *translation[] = {
                "Ag_a", "Au_a",
                "Ag_b", "Au_b",
                "Ag_-3/2", "Au_-3/2",
                "Ag_+3/2", "Au_+3/2",
                "Ag_0", "Au_0",
                "Ag_2", "Au_2",
                "Ag_+1", "Au_+1",
                "Ag_-1", "Au_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // C2v nonrel
    if (strcmp(rep_names[0], "A1 a") == 0 && strcmp(rep_names[1], "B2 a") == 0 ||
        strcmp(rep_names[0], "A1 a") == 0 && strcmp(rep_names[1], "B1 a") == 0) {
        char *translation[] = {
                "A1_a", "B2_a", "B1_a", "A2_a",
                "A1_b", "B2_b", "B1_b", "A2_b",
                "A1_-3/2", "B2_-3/2", "B1_-3/2", "A2_-3/2",
                "A1_+3/2", "B2_+3/2", "B1_+3/2", "A2_+3/2",
                "A1_0", "B2_0", "B1_0", "A2_0",
                "A1_2", "B2_2", "B1_2", "A2_2",
                "A1_+1", "B2_+1", "B1_+1", "A2_+1",
                "A1_-1", "B2_-1", "B1_-1", "A2_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // C2h nonrel
    if (rep_name_exists(nsym, rep_names, "Ag a") &&
        rep_name_exists(nsym, rep_names, "Au a") &&
        rep_name_exists(nsym, rep_names, "Ag b") &&
        rep_name_exists(nsym, rep_names, "Au b") &&
        rep_name_exists(nsym, rep_names, "Bg a") &&
        rep_name_exists(nsym, rep_names, "Bu a") &&
        rep_name_exists(nsym, rep_names, "Bg b") &&
        rep_name_exists(nsym, rep_names, "Bu b")) {
        char *translation[] = {
                "Ag_a", "Bg_a", "Bu_a", "Au_a",
                "Ag_b", "Bg_b", "Bu_b", "Au_b",
                "Ag_-3/2", "Bg_-3/2", "Bu_-3/2", "Au_-3/2",
                "Ag_+3/2", "Bg_+3/2", "Bu_+3/2", "Au_+3/2",
                "Ag_0", "Bg_0", "Bu_0", "Au_0",
                "Ag_2", "Bg_2", "Bu_2", "Au_2",
                "Ag_+1", "Bg_+1", "Bu_+1", "Au_+1",
                "Ag_-1", "Bg_-1", "Bu_-1", "Au_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // D2 nonrel
    if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "B3 a") == 0) {
        char *translation[] = {
                "A_a", "B3_a", "B1_a", "B2_a",
                "A_b", "B3_b", "B1_b", "B2_b",
                "A_-3/2", "B3_-3/2", "B1_-3/2", "B2_-3/2",
                "A_+3/2", "B3_+3/2", "B1_+3/2", "B2_+3/2",
                "A_0", "B3_0", "B1_0", "B2_0",
                "A_2", "B3_2", "B1_2", "B2_2",
                "A_+1", "B3_+1", "B1_+1", "B2_+1",
                "A_-1", "B3_-1", "B1_-1", "B2_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // D2h nonrel
    if (strcmp(rep_names[0], "Ag a") == 0 && (
            strcmp(rep_names[1], "B1ua") == 0 ||
            strcmp(rep_names[1], "B2ua") == 0 ||
            strcmp(rep_names[1], "B3ua") == 0 ||
            strcmp(rep_names[1], "B1ga") == 0 ||
            strcmp(rep_names[1], "B2ga") == 0 ||
            strcmp(rep_names[1], "B3ga") == 0)) {
        char *translation[] = {
                "Ag_a", "B1u_a", "B2u_a", "B3g_a", "B3u_a", "B2g_a", "B1g_a", "Au_a",
                "Ag_b", "B1u_b", "B2u_b", "B3g_b", "B3u_b", "B2g_b", "B1g_b", "Au_b",
                "Ag_-3/2", "B1u_-3/2", "B2u_-3/2", "B3g_-3/2", "B3u_-3/2", "B2g_-3/2", "B1g_-3/2", "Au_-3/2",
                "Ag_+3/2", "B1u_+3/2", "B2u_+3/2", "B3g_+3/2", "B3u_+3/2", "B2g_+3/2", "B1g_+3/2", "Au_+3/2",
                "Ag_0", "B1u_0", "B2u_0", "B3g_0", "B3u_0", "B2g_0", "B1g_0", "Au_0",
                "Ag_2", "B1u_2", "B2u_2", "B3g_2", "B3u_2", "B2g_2", "B1g_2", "Au_2",
                "Ag_+1", "B1u_+1", "B2u_+1", "B3g_+1", "B3u_+1", "B2g_+1", "B1g_+1", "Au_+1",
                "Ag_-1", "B1u_-1", "B2u_-1", "B3g_-1", "B3u_-1", "B2g_-1", "B1g_-1", "Au_-1"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // C1 rel
    if (strcmp(rep_names[0], "   A") == 0 && strcmp(rep_names[1], "   a") == 0) {
        strcpy(rep_names[0], "A");
        strcpy(rep_names[1], "a");
    }
    // Ci rel
    if (strcmp(rep_names[0], "  AG") == 0 && strcmp(rep_names[1], "  AU") == 0) {
        strcpy(rep_names[0], "AG");
        strcpy(rep_names[1], "AU");
        strcpy(rep_names[2], "ag");
        strcpy(rep_names[3], "au");
    }
    // C2, Cs, C2v, D2 rel
    if (strcmp(rep_names[0], "  1E") == 0 && strcmp(rep_names[1], "  2E") == 0) {
        strcpy(rep_names[0], "1E");
        strcpy(rep_names[1], "2E");
        strcpy(rep_names[2], "a");
        strcpy(rep_names[3], "b");
    }
    // C2h, D2h rel
    if (strcmp(rep_names[0], " 1Eg") == 0 && strcmp(rep_names[1], " 2Eg") == 0) {
        strcpy(rep_names[0], "1Eg");
        strcpy(rep_names[1], "2Eg");
        strcpy(rep_names[2], "1Eu");
        strcpy(rep_names[3], "2Eu");
        strcpy(rep_names[4], "ag");
        strcpy(rep_names[5], "bg");
        strcpy(rep_names[6], "au");
        strcpy(rep_names[7], "bu");
    }
    // nonrel Cinfv = C2v, nonrel Dinfh = D2h
    // subgroups of the D2h point group, relativistic case -- nothing to do
    // Cinfv rel
    if (strcmp(rep_names[0], "   1") == 0 && strcmp(rep_names[1], "  -1") == 0) {
        char *translation[] = {
                "1/2+", "1/2-", "3/2+", "3/2-", "5/2+", "5/2-", "7/2+", "7/2-",
                "9/2+", "9/2-", "11/2+", "11/2-", "13/2+", "13/2-", "15/2+", "15/2-",
                "17/2+", "17/2-", "19/2+", "19/2-", "21/2+", "21/2-", "23/2+", "23/2-",
                "25/2+", "25/2-", "27/2+", "27/2-", "29/2+", "29/2-", "31/2+", "31/2-",
                "0", "1+", "1-", "2+", "2-", "3+", "3-", "4+",
                "4-", "5+", "5-", "6+", "6-", "7+", "7-", "8+",
                "8-", "9+", "9-", "10+", "10-", "11+", "11-", "12+",
                "12-", "13+", "13-", "14+", "14-", "15+", "15-", "16+"
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // Dinfh rel
    if (strcmp(rep_names[0], "  1g") == 0 && strcmp(rep_names[1], " -1g") == 0) {
        char *translation[] = {
                "1/2g+", "1/2g-", "3/2g+", "3/2g-", "5/2g+", "5/2g-", "7/2g+", "7/2g-",
                "9/2g+", "9/2g-", "11/2g+", "11/2g-", "13/2g+", "13/2g-", "15/2g+", "15/2g-",
                "1/2u+", "1/2u-", "3/2u+", "3/2u-", "5/2u+", "5/2u-", "7/2u+", "7/2u-",
                "9/2u+", "9/2u-", "11/2u+", "11/2u-", "13/2u+", "13/2u-", "15/2u+", "15/2u-",
                "0g", "1g+", "1g-", "2g+", "2g-", "3g+", "3g-", "4g+",
                "4g-", "5g+", "5g-", "6g+", "6g-", "7g+", "7g-", "8g+",
                "0u", "1u+", "1u-", "2u+", "2u-", "3u+", "3u-", "4u+",
                "4u-", "5u+", "5u-", "6u+", "6u-", "7u+", "7u-", "8u+",
        };
        for (int irep = 0; irep < nsym; irep++) {
            strcpy(rep_names[irep], translation[irep]);
        }
    }
}


char **dirac_interface_extract_rep_names(int nsym)
{
    char **rep_names = (char **) cc_malloc(nsym * sizeof(char *));

    for (int i = 0; i < nsym; i++) {
        rep_names[i] = (char *) cc_malloc(sizeof(char) * 8);
        for (int j = 0; j < 4; j++) {
            rep_names[i][j] = dirac_data.repanames[4 * i + j];
        }
        rep_names[i][4] = '\0';  // cut end of line (only 4 char-s are significant)
    }

    return rep_names;
}


void extract_direct_product_table_dirac(int nsym, int *mult_table)
{
    memset(mult_table, 0, sizeof(int) * nsym * nsym);

    for (int i = 0; i < nsym; i++) {
        for (int j = 0; j < nsym; j++) {
            mult_table[i * nsym + j] = dirac_data.multb[i][j] - 1; // -1 since C arrays start from zero
        }
    }
}


/**
 * adapter to the 'create_spinor_info' subroutine,
 * is needed to encapsulate operations with fixed-width integers.
 */
void create_spinor_info_dirac(int32_t n_spinors, int32_t *irpamo_array, double *spinor_energies, int32_t *occ_numbers)
{
    int *irrep_numbers_buf = cc_malloc(n_spinors * sizeof(int));
    int *occ_numbers_buf = cc_malloc(n_spinors * sizeof(int));

    for (int i = 0; i < n_spinors; i++) {
        irrep_numbers_buf[i] = irpamo_array[i] - 1;
        occ_numbers_buf[i] = occ_numbers[i];
    }

    create_spinor_info(n_spinors, irrep_numbers_buf, spinor_energies, occ_numbers_buf);

    cc_free(irrep_numbers_buf);
    cc_free(occ_numbers_buf);
}


void check_irrep_names()
{
    // now we can check all symbols of irreps in the input data
    // (1) properties calculations
    for (size_t i = 0; i < cc_opts->n_denmat; i++) {
        cc_denmat_query_t *query = cc_opts->denmat_query + i;
        if (get_rep_number(query->rep1_name) == -1) {
            errquit("Error in the input file: wrong irrep name '%s'", query->rep1_name);
        }
        if (get_rep_number(query->rep2_name) == -1) {
            errquit("Error in the input file: wrong irrep name '%s'", query->rep2_name);
        }
    }
    // (2) active space specification
    for (size_t i = 0; i < CC_MAX_NUM_IRREPS; i++) {
        cc_active_spec_t *sp = cc_opts->active_specs + i;
        if (*(sp->rep_name)) {
            if (get_rep_number(sp->rep_name) == -1) {
                errquit("Error in the input file: wrong irrep name in the active space specification: '%s'",
                        sp->rep_name);
            }
        }
    }
    // (3) number of roots to be analyzed
    cc_space_t *nroots_specs = &cc_opts->nroots_specs;
    for (size_t i = 0; i < CC_MAX_NUM_IRREPS; i++) {
        char *rep_name = nroots_specs->rep_names[i];
        if (*(rep_name)) {
            if (get_rep_number(rep_name) == -1) {
                errquit("Error in the input file: wrong irrep name in the nroots specification: '%s'",
                        rep_name);
            }
        }
    }
}


// prints string with current time
void c_asctime()
{
    time_t tm = time(0);
    printf(" %s", asctime(localtime(&tm)));
}


int rep_name_exists(int nrep, char **list, char *query)
{
    for (int i = 0; i < nrep; i++) {
        if (strcmp(list[i], query) == 0) {
            return 1;
        }
    }
    return 0;
}


void dirac_interface_debug_print()
{
    printf("nspinors = %d\n", dirac_data.nspinors);
    printf("breit = %d\n", dirac_data.breit);
    printf("enuc = %.8f\n", dirac_data.enuc);
    printf("invsym = %d\n", dirac_data.invsym);
    printf("nz_arith = %d\n", dirac_data.nz_arith);
    printf("is_spinfree = %d\n", dirac_data.is_spinfree);
    printf("norb_total = %d\n", dirac_data.norb_total);
    printf("escf = %.8f\n", dirac_data.escf);
    printf("enuc = %.8f\n", dirac_data.enuc);
    printf("nsymrp = %d\n", dirac_data.nsymrp);
    printf("repnames = \n");
    for (int i = 0; i < dirac_data.nsymrp; i++) {
        printf("  [%d] ", i);
        for (int j = 0; j < 14; j++) {
            printf("%c", dirac_data.repnames[14 * i + j]);
        }
        printf("\n");
    }
    printf("nactive = ");
    for (int i = 0; i < 8; i++) {
        printf("%d ", dirac_data.nactive[i]);
    }
    printf("\n");
    printf("nstr = %d %d\n", dirac_data.nstr[0], dirac_data.nstr[1]);
    printf("nfrozen = ");
    for (int i = 0; i < 6; i++) {
        printf("%d ", dirac_data.nfrozen[i]);
    }
    printf("\n");
    printf("ndelete = ");
    for (int i = 0; i < 8; i++) {
        printf("%d ", dirac_data.ndelete[i]);
    }
    printf("\n");
    printf("nsymrpa = %d\n", dirac_data.nsymrpa);
    printf("repanames = \n");
    for (int i = 0; i < 2 * dirac_data.nsymrpa; i++) {
        printf("  [%d] ", i);
        for (int j = 0; j < 4; j++) {
            printf("%c", dirac_data.repanames[4 * i + j]);
        }
        printf("\n");
    }
    printf("irrep direct product table =\n");
    for (int i = 0; i < 2 * dirac_data.nsymrpa; i++) {
        for (int j = 0; j < 2 * dirac_data.nsymrpa; j++) {
            printf("%2d ", dirac_data.multb[i][j]);
        }
        printf("\n");
    }
    printf("nbsymrp = %d\n", dirac_data.nbsymrp);
    printf("spinor data: ");
    printf("no, irpmo, irpamo, iocc, eorbmo, ibspi\n");
    for (int i = 0; i < dirac_data.nspinors; i++) {
        printf("%4d%4d%4d%4d%15.8f%4d\n", i + 1, dirac_data.irpmo[i], dirac_data.irpamo[i],
               dirac_data.iocc[i], dirac_data.eorbmo[i], dirac_data.ibspi[i]);
    }
}

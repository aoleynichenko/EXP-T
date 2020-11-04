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
 * dirac_interface.c
 * =================
 *
 * Implementation of the interface to the DIRAC program package.
 * Reads info about spinors and transformed MO integrals from the unformatted
 * file produced by DIRAC (with default names MRCONEE, MDCINT, MDPROP)
 *
 * 2018 Alexander Oleynichenko
 ******************************************************************************/

#include "interfaces.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "engine.h"
#include "platform.h"
#include "error.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "symmetry.h"
#include "timer.h"


void dirac_interface_binary(int *err);

// arguments to be passed to the Fortran code
// (via the 'mrconee_mdcint' common block, see dirac_binary.f90)
extern struct {
    char oneel_file[1024];
    char twoel_file[1024];
    char prop_file[1024];
    char gaunt_file[1024];
    int64_t gaunt_enabled;
    int64_t n_twoprop;
    char twoprop_files[64][1024];
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
    int32_t multb[CC_MAX_NREP][CC_MAX_NREP];   // multiplication table for direct products in the Abelian subgroup
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


// prints string with current time
void c_asctime()
{
    time_t tm = time(0);
    printf(" %s", asctime(localtime(&tm)));
}


// reads info about spinors and transformed MO integrals from formatted file
// NOTE: no error checking at all!
void dirac_interface(char *moints_file_1, char *moints_file_2, char *moints_file_prop, cc_options_t *opts)
{
    int i, j, k;
    int err;
    int nz;           // arithmetic, 1:real, 2:complex, 4:quaternion

    // get tile size and pass it to the Fortran side
    tilesize.tile_size = opts->tile_size;
    skip_sorting.skip_mdcint = opts->reuse_integrals_2;
    recommendedcarith.recommended_carith = (opts->recommended_arith == CC_ARITH_COMPLEX) ? 1 : 0;

    // prepare Fortran-compatible names of integral files
    for (i = 0; i < 1024; i++){
        mrconee_mdcint.oneel_file[i] = ' ';
        mrconee_mdcint.twoel_file[i] = ' ';
        mrconee_mdcint.prop_file [i] = ' ';
        mrconee_mdcint.gaunt_file[i] = ' ';
        for (j = 0; j < CC_MAX_NPROP; j++) {
            mrconee_mdcint.twoprop_files[j][i] = ' ';
        }
    }

    for (i = 0; i < strlen(moints_file_1); i++){
        mrconee_mdcint.oneel_file[i] = moints_file_1[i];
    }
    for (i = 0; i < strlen(moints_file_2); i++){
        mrconee_mdcint.twoel_file[i] = moints_file_2[i];
    }
    for (i = 0; i < strlen(moints_file_prop); i++){
        mrconee_mdcint.prop_file[i] = moints_file_prop[i];
    }
    for (i = 0; i < strlen(opts->integral_file_gaunt); i++){
        mrconee_mdcint.gaunt_file[i] = opts->integral_file_gaunt[i];
    }
    for (j = 0; j < CC_MAX_NPROP; j++) {
        for (i = 0; i < strlen(opts->twoprop_file[j]); i++) {
            mrconee_mdcint.twoprop_files[j][i] = opts->twoprop_file[j][i];
        }
    }
    mrconee_mdcint.gaunt_enabled = opts->gaunt_defined ? 1 : 0;
    mrconee_mdcint.n_twoprop = opts->n_twoprop;

    // call the Fortran subroutine
    fflush(stdout);
    timer_new_entry("dirac", "DIRAC interface (MRCONEE/MDCINT)");
    timer_start("dirac");
    dirac_interface_binary(&err);
    timer_stop("dirac");
    if (err == 1) {
        errquit("ERROR IN THE BINARY INTERFACE\n");
    }

    if (opts->print_level >= CC_PRINT_DEBUG) {
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
        for (i = 0; i < dirac_data.nsymrp; i++){
            printf("  [%d] ", i);
            for (j = 0; j < 14; j++){
                printf("%c", dirac_data.repnames[14 * i + j]);
            }
            printf("\n");
        }
        printf("nactive = ");
        for (i = 0; i < 8; i++){
            printf("%d ", dirac_data.nactive[i]);
        }
        printf("\n");
        printf("nstr = %d %d\n", dirac_data.nstr[0], dirac_data.nstr[1]);
        printf("nfrozen = ");
        for (i = 0; i < 6; i++){
            printf("%d ", dirac_data.nfrozen[i]);
        }
        printf("\n");
        printf("ndelete = ");
        for (i = 0; i < 8; i++){
            printf("%d ", dirac_data.ndelete[i]);
        }
        printf("\n");
        printf("nsymrpa = %d\n", dirac_data.nsymrpa);
        printf("repanames = \n");
        for (i = 0; i < 2 * dirac_data.nsymrpa; i++){
            printf("  [%d] ", i);
            for (j = 0; j < 4; j++){
                printf("%c", dirac_data.repanames[4 * i + j]);
            }
            printf("\n");
        }
        printf("irrep direct product table =\n");
        for (i = 0; i < 2 * dirac_data.nsymrpa; i++){
            for (j = 0; j < 2 * dirac_data.nsymrpa; j++){
                printf("%2d ", dirac_data.multb[i][j]);//64*i+j]);
            }
            printf("\n");
        }
        printf("nbsymrp = %d\n", dirac_data.nbsymrp);
        printf("spinor data: ");
        printf("no, irpmo, irpamo, iocc, eorbmo, ibspi\n");
        for (i = 0; i < dirac_data.nspinors; i++){
            printf("%4d%4d%4d%4d%15.8f%4d\n", i + 1, dirac_data.irpmo[i], dirac_data.irpamo[i],
                   dirac_data.iocc[i], dirac_data.eorbmo[i], dirac_data.ibspi[i]);
        }
    } // end of debug print

    // extract data from the "Fortran common block"
    // and reinterpret them

    opts->enuc = dirac_data.enuc;
    opts->escf = dirac_data.escf;

    nz = dirac_data.nz_arith;
    if (nz == 1) {
        point_group_nz = CC_GROUP_REAL;
    }
    else if (nz == 2) {
        point_group_nz = CC_GROUP_COMPLEX;
    }
    else {
        point_group_nz = CC_GROUP_QUATERNION;
    }

    if (nz == 1 || dirac_data.is_spinfree) {
        arith = CC_ARITH_REAL;
        carith = 0;
    }
    else{
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
        if (dirac_data.is_spinfree) printf(" Spin-orbit: off (spinfree)\n");
        else printf(" Spin-orbit: on\n");
    }

    // Number of double group irreps
    nsym = dirac_data.nsymrpa * 2;

    // Representation names
    rep_names = (char **) cc_malloc(nsym * sizeof(char *));
    for (i = 0; i < nsym; i++){
        rep_names[i] = (char *) cc_malloc(sizeof(char) * 8);
        for (j = 0; j < 4; j++){
            rep_names[i][j] = dirac_data.repanames[4 * i + j];
        }
        rep_names[i][4] = '\0';  // cut end of line (only 4 char-s are significant)
    }

    // detect fully symmetric irrep
    if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "A  b") == 0) {
        irrep_a1 = 4;
        set_point_group_name("C1");
    }
    else if (strcmp(rep_names[0], "Ag a") == 0 && strcmp(rep_names[1], "Au a") == 0) {
        irrep_a1 = 8;
        set_point_group_name("Ci");
    }
    else if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "B  a") == 0) {
        irrep_a1 = 8;
        set_point_group_name("C2");
    }
    else if (strcmp(rep_names[0], "A' a") == 0 && strcmp(rep_names[1], "A\" a") == 0) {
        irrep_a1 = 8;
        set_point_group_name("Cs");
    }
    else if (strcmp(rep_names[0], "A1 a") == 0) {
        irrep_a1 = 16;
        set_point_group_name("C2v");
    }
    else if (strcmp(rep_names[0], "A  a") == 0) {
        irrep_a1 = 16;
        set_point_group_name("D2");
    }
    else if (strcmp(rep_names[0], "Ag a") == 0 && strcmp(rep_names[1], "Bg a") == 0) {
        irrep_a1 = 16;
        set_point_group_name("C2h");
    }
    else if (strcmp(rep_names[0], "Ag a") == 0) {
        irrep_a1 = 32;
        set_point_group_name("D2h");
    }
        // double groups
    else if (strcmp(rep_names[0], "   A") == 0 && strcmp(rep_names[1], "   a") == 0) {
        irrep_a1 = 1;
        set_point_group_name("C1");
    }
    else if (strcmp(rep_names[0], "  AG") == 0 && strcmp(rep_names[1], "  AU") == 0) {
        irrep_a1 = 2;
        set_point_group_name("Ci");
    }
    else if (strcmp(rep_names[0], "  1E") == 0 && strcmp(rep_names[1], "  2E") == 0) {
        irrep_a1 = 2;
        set_point_group_name("C2, Cs, C2v or D2");
    }
    else if (strcmp(rep_names[0], " 1Eg") == 0 && strcmp(rep_names[1], " 2Eg") == 0) {
        irrep_a1 = 4;
        set_point_group_name("C2h or D2h");
    }
    else if (strcmp(rep_names[0], "   1") == 0 && strcmp(rep_names[1], "  -1") == 0) {
        irrep_a1 = 32;
        set_point_group_name("C32 aka Cinfv");
    }
    else if (strcmp(rep_names[0], "  1g") == 0 && strcmp(rep_names[1], " -1g") == 0) {
        irrep_a1 = 32;
        set_point_group_name("C16h aka Dinfh");
    }
    else{
        irrep_a1 = 0;
        set_point_group_name("undetected");
    }

    // let us use new irrep symbols
    /* for nonrel groups
       left -- L. Visscher's notation, right -- A. Oleynichenko's notation
     (Ms projection)
     a -> a
     b -> b
     3 -> -3/2
     3 -> +3/2
     0 -> 0
     4 -> 2
     2 -> +1
     2 -> -1
     */
    // C1 nonrel
    if (strcmp(rep_names[0], "A  a") == 0 && strcmp(rep_names[1], "A  b") == 0) {
        char *translation[] = {
                "A_a", "A_b", "A_-3/2", "A_+3/2", "A_0", "A_2", "A_+1", "A_-1"
        };
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // Ci nonrel
    if (strcmp(rep_names[0], "Ag a") == 0 && strcmp(rep_names[1], "Au a") == 0) {
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
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
            strcpy(rep_names[irep], translation[irep]);
        }
    }
    // C2h nonrel
    if (strcmp(rep_names[0], "Ag a") == 0 && strcmp(rep_names[1], "Bg a") == 0) {
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
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
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
        for (int irep = 0; irep < nsym; irep++){
            strcpy(rep_names[irep], translation[irep]);
        }
    }

    // "Multiplication table" ( = dir prod table for Abelian groups)
    // in fact, 3d array

    // allocate direct product table
    memset(dir_prod_table_abelian, 0, sizeof(int)*CC_MAX_NREP*CC_MAX_NREP);
    dir_prod_table = (int ***) cc_malloc(nsym * sizeof(int **));
    for (i = 0; i < nsym; i++){
        dir_prod_table[i] = (int **) cc_malloc(nsym * sizeof(int *));
        for (j = 0; j < nsym; j++){
            dir_prod_table[i][j] = (int *) cc_malloc((nsym + 1) * sizeof(int));
            for (k = 0; k < nsym + 1; k++){
                dir_prod_table[i][j][k] = -1;
            }
            // get values from the dirac_data structure (from MRCONEE files)
            dir_prod_table[i][j][0] = 1;
            dir_prod_table[i][j][1] = dirac_data.multb[i][j] - 1; // -1 since C arrays start from zero
            dir_prod_table_abelian[i][j] = dirac_data.multb[i][j] - 1;
        }
    }
    is_abelian = 1;

    /*printf("Mult table:\n");
    for (i = 0; i < nsym; i++) {
        for (j = 0; j < nsym; j++)
            printf("%3d", dir_prod_table[i][j][1]);
        printf("\n");
    }*/

    nspinors = dirac_data.nspinors;

    // fill spinor info array
    //  *** Spinor, irrep, occupation, energy ***
    spinor_info = (spinor_attr_t *) cc_malloc(sizeof(spinor_attr_t) * nspinors);
    for (i = 0; i < nspinors; i++){
        spinor_info[i].seqno = i;    // numeration will start from 0
        spinor_info[i].repno = dirac_data.irpamo[i] - 1; // and for rep-s too
        spinor_info[i].occ = dirac_data.iocc[i];
        spinor_info[i].eps = dirac_data.eorbmo[i];
        spinor_info[i].active = 0;
        if (opts->print_level >= CC_PRINT_DEBUG)
            printf("%4d%4d%4d%4d%8.4f\n", spinor_info[i].seqno, spinor_info[i].repno,
                   spinor_info[i].occ, spinor_info[i].active, spinor_info[i].eps);
    }

    create_spinor_blocks(opts->tile_size);
    classify_spinors(opts->actsp_min, opts->actsp_max,
                     opts->nacth, opts->nactp);  // => quasiparticles: IH, AH, AP, IP

    print_symmetry_info();
    print_spinor_info();

    // now we can check all symbols of irreps in the input data
    // (1) properties calculations
    for (size_t i = 0; i < cc_opts->n_denmat; i++){
        cc_denmat_query_t *query = cc_opts->denmat_query + i;
        if (get_rep_number(query->rep1_name) == -1) {
            errquit("Error in the input file: wrong irrep name '%s'", query->rep1_name);
        }
        if (get_rep_number(query->rep2_name) == -1) {
            errquit("Error in the input file: wrong irrep name '%s'", query->rep2_name);
        }
    }
    // (2) active space specification
    for (size_t i = 0; i < CC_MAX_NREP; i++) {
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
    for (size_t i = 0; i < CC_MAX_NREP; i++) {
        char *rep_name = nroots_specs->rep_names[i];
        if (*(rep_name)) {
            if (get_rep_number(rep_name) == -1) {
                errquit("Error in the input file: wrong irrep name in the nroots specification: '%s'",
                        rep_name);
            }
        }
    }
}

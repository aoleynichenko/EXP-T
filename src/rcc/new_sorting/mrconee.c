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


#include "mrconee.h"

#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "libunf.h"
#include "../../transform_natorb/errquit.h"

enum {
    DIRAC_INT_4 = 4,
    DIRAC_INT_8 = 8
};

static void new_detect_dirac_point_group(
    char **rep_names, char *group_name, int *fully_sym_irrep);

int new_rep_name_exists(int nrep, char **list, char *query);

void new_rename_irreps_dirac_to_expt(int nsym, char **rep_names);

int test_dirac_integer_size(char *path);

int mrconee_read_header(unf_file_t *file, mrconee_data_t *data);

int mrconee_read_fermion_irrep_occs(unf_file_t *file, mrconee_data_t *data, int *fermion_irrep_occs);

int mrconee_read_abelian_irreps(unf_file_t *file, mrconee_data_t *data);

int mrconee_read_multiplication_table(unf_file_t *file, mrconee_data_t *data);

int mrconee_read_spinor_info(unf_file_t *file, mrconee_data_t *data, int *fermion_irrep_occs);

int mrconee_read_fock(unf_file_t *file, mrconee_data_t *data);


mrconee_data_t *read_mrconee(char *path)
{
    printf("\n");
    printf("                ***********************************\n");
    printf("                ** DIRAC integrals: MRCONEE file **\n");
    printf("                ***********************************\n");
    printf("\n");

    unf_file_t *file = unf_open(path, "r", UNF_ACCESS_SEQUENTIAL);
    if (file == NULL) {
        errquit("cannot open MRCONEE unformatted file\n");
        return NULL;
    }

    mrconee_data_t *data = (mrconee_data_t *) calloc(1, sizeof(mrconee_data_t));

    /*
     * check size of integer numbers used by DIRAC: 4-byte or 8-byte
     */
    int dirac_int_size = test_dirac_integer_size(path);
    data->dirac_int_size = dirac_int_size;

    /*
     * record 1
     * header: number of spinors, SCF and nuclear repulsion energy,
     * inversion symmetry, point group type, spinfree
     */
    int error_code = mrconee_read_header(file, data);
    if (error_code == EXIT_FAILURE) {
        errquit("while reading MRCONEE file (record 1)");
        free_mrconee_data(data);
        return NULL;
    }

    /*
     * record 2
     * number of electrons in each fermion irrep
     */
    int fermion_irrep_occs[8];
    error_code = mrconee_read_fermion_irrep_occs(file, data, fermion_irrep_occs);
    if (error_code == EXIT_FAILURE) {
        errquit("while reading MRCONEE file (record 2)");
        free_mrconee_data(data);
        return NULL;
    }

    /*
     * record 3
     * irreps of abelian subgroup
     */
    error_code = mrconee_read_abelian_irreps(file, data);
    if (error_code == EXIT_FAILURE) {
        errquit("while reading MRCONEE file (record 3)");
        free_mrconee_data(data);
        return NULL;
    }

    /*
     * record 4
     * irrep multiplication table
     */
    error_code = mrconee_read_multiplication_table(file, data);
    if (error_code == EXIT_FAILURE) {
        errquit("while reading MRCONEE file (record 4)");
        free_mrconee_data(data);
        return NULL;
    }

    /*
     * record 5
     * information about spinors
     */
    error_code = mrconee_read_spinor_info(file, data, fermion_irrep_occs);
    if (error_code == EXIT_FAILURE) {
        errquit("while reading MRCONEE file (record 5)");
        free_mrconee_data(data);
        return NULL;
    }

    /*
     * record 6
     * Fock matrix
     */
    error_code = mrconee_read_fock(file, data);
    if (error_code == EXIT_FAILURE) {
        errquit("while reading MRCONEE file (record 6)");
        free_mrconee_data(data);
        return NULL;
    }

    return data;
}


int mrconee_read_header(unf_file_t *file, mrconee_data_t *data)
{
    // total number of spinors
    int32_t num_spinors;
    // Breit interaction active in SCF
    int32_t breit;
    // core energy (inactive energy + nuclear repulsion)
    double enuc;
    // inversion symmetry (yes : 2; no : 1)
    int32_t invsym;
    // group type (1 real, 2 complex, 4 quaternion)
    int32_t nz_arith;
    // spinfree formalism
    int32_t is_spinfree;
    // total number of orbitals (so including frozen or deleted orbitals)
    int32_t norb_total;
    // total SCF energy
    double scf_energy;

    int nread;
    if (data->dirac_int_size == DIRAC_INT_4) {
        nread = unf_read(file, "2i4,r8,4i4,r8", &num_spinors, &breit, &enuc, &invsym, &nz_arith, &is_spinfree,
                                 &norb_total, &scf_energy);
    }
    else {
        int64_t num_spinors_8;
        int64_t breit_8;
        int64_t invsym_8;
        int64_t nz_arith_8;
        int64_t is_spinfree_8;
        int64_t norb_total_8;

        nread = unf_read(file, "2i8,r8,4i8,r8", &num_spinors_8, &breit_8, &enuc, &invsym_8, &nz_arith_8,
                                 &is_spinfree_8,
                                 &norb_total_8, &scf_energy);

        num_spinors = (int32_t) num_spinors_8;
        breit = (int32_t) breit_8;
        invsym = (int32_t) invsym_8;
        nz_arith = (int32_t) nz_arith_8;
        is_spinfree = (int32_t) is_spinfree_8;
        norb_total = (int32_t) norb_total_8;
    }

    if (nread != 8) {
        return EXIT_FAILURE;
    }

    data->num_spinors = num_spinors;
    data->nuc_rep_energy = enuc;
    data->group_arith = nz_arith;
    data->is_spinfree = is_spinfree;
    data->scf_energy = scf_energy;
    data->invsym = invsym;

    return EXIT_SUCCESS;
}


int mrconee_read_fermion_irrep_occs(unf_file_t *file, mrconee_data_t *data, int *fermion_irrep_occs)
{
    int invsym = data->invsym;

    // number of fermion irreps in parent group
    int32_t nsymrp;
    // names of these irreps (gerade, ungerade)
    char repnames[14 * 8];
    // number of spinors active in the transf-n (not valence!)
    int32_t nactive[8];
    // total number of orbitals of this ircop
    int32_t nstr[2];
    // number of occupied frozen (core) spinors, [0]: total, [1]: positive energy, [2]: negative energy
    int32_t nfrozen[3][2];
    // number of deleted spinors
    int32_t ndelete[8];

    int nread;
    if (data->dirac_int_size == DIRAC_INT_4) {
        nread = unf_read(file, "i4,c14[i4],6i4[i4]",
                                 &nsymrp, repnames, &nsymrp, nactive, &nsymrp,
                                 nstr, &invsym, nfrozen[0], &invsym, nfrozen[1], &invsym, nfrozen[2], &invsym, ndelete,
                                 &invsym);
    }
    else {
        int64_t nsymrp_8;
        int64_t nactive_8[8];
        int64_t nstr_8[2];
        int64_t nfrozen_8[3][2];
        int64_t ndelete_8[8];

        nread = unf_read(file, "i8,c14[i4],6i8[i4]",
                                 &nsymrp_8, repnames, &nsymrp_8, nactive_8, &nsymrp_8,
                                 nstr_8, &invsym, nfrozen_8[0], &invsym, nfrozen_8[1], &invsym, nfrozen_8[2], &invsym,
                                 ndelete_8, &invsym);

        for (int i = 0; i < nsymrp_8; i++) {
            nactive[i] = (int32_t) nactive_8[i];
        }
    }

    if (nread != 8) {
        return EXIT_FAILURE;
    }

    for (int i = 0; i < nsymrp; i++) {
        fermion_irrep_occs[i] = nactive[i];
    }

    return EXIT_SUCCESS;
}


int mrconee_read_abelian_irreps(unf_file_t *file, mrconee_data_t *data)
{
    // number of fermion irreps in the Abelian subgroup
    int32_t nsymrpa;
    // names of these irreps
    char repanames[4 * 4 * 64];

    int nread;
    if (data->dirac_int_size == DIRAC_INT_4) {
        unf_read(file, "i4", &nsymrpa);
        int32_t size_repanames = 2 * nsymrpa;
        unf_backspace(file);
        nread = unf_read(file, "i4,c4[i4]", &nsymrpa, repanames, &size_repanames);
    }
    else {
        int64_t nsymrpa_8;

        unf_read(file, "i8", &nsymrpa_8);
        int32_t size_repanames = (int32_t) (2 * nsymrpa_8);
        unf_backspace(file);
        nread = unf_read(file, "i8,c4[i4]", &nsymrpa_8, repanames, &size_repanames);

        nsymrpa = (int32_t) nsymrpa_8;
    }
    if (nread != 2) {
        return EXIT_FAILURE;
    }

    data->num_irreps = 2 * nsymrpa;
    data->irrep_names = (char **) calloc(data->num_irreps, sizeof(char *));
    for (int i = 0; i < data->num_irreps; i++) {
        char *name = (char *) calloc(32, sizeof(char));
        name[0] = repanames[4 * i];
        name[1] = repanames[4 * i + 1];
        name[2] = repanames[4 * i + 2];
        name[3] = repanames[4 * i + 3];
        name[4] = '\0';
        data->irrep_names[i] = name;
    }

    data->point_group = (char *) calloc(16, sizeof(char));
    new_detect_dirac_point_group(data->irrep_names, data->point_group, &data->totally_sym_irrep);
    new_rename_irreps_dirac_to_expt(data->num_irreps, data->irrep_names);

    return EXIT_SUCCESS;
}


int mrconee_read_multiplication_table(unf_file_t *file, mrconee_data_t *data)
{
    int nread;
    int32_t multb[64 * 64];
    int multb_size = data->num_irreps * data->num_irreps;

    if (data->dirac_int_size == DIRAC_INT_4) {
        nread = unf_read(file, "i4[i4]", multb, &multb_size);
    }
    else {
        int64_t multb_8[64 * 64];
        nread = unf_read(file, "i8[i4]", multb_8, &multb_size);
        for (int i = 0; i < data->num_irreps * data->num_irreps; i++) {
            multb[i] = (int32_t) multb_8[i];
        }
    }
    if (nread != 1) {
        return EXIT_FAILURE;
    }

    data->mult_table = (int *) calloc(data->num_irreps * data->num_irreps, sizeof(int *));
    for (int i = 0; i < data->num_irreps * data->num_irreps; i++) {
        data->mult_table[i] = multb[i] - 1; // Fortran to C conversion
    }

    return EXIT_SUCCESS;
}


int mrconee_read_spinor_info(unf_file_t *file, mrconee_data_t *data, int *fermion_irrep_occs)
{
    int invsym = data->invsym;
    int num_spinors = data->num_spinors;

    // approximate boson irrep identification (needed in LUCIAREL and LUCITA)
    int32_t *ibspi = (int32_t *) calloc(num_spinors, sizeof(int32_t));
    // number of boson symmetry reps (for LUCITA)
    int32_t nbsymrp;
    // (?) number of Kramers pairs in each u/g irrep
    int32_t norb[2];

    struct irp_irpa_eorb {
        int32_t irp; // irrep in parent group (1:gerade, 2:ungerade)
        int32_t irpa; // irrep in Abelian subgroup
        double eorb; // one-electron energies taken from SCF
    };

    struct irp_irpa_eorb *buf_irp_irpa_eorb = (struct irp_irpa_eorb *) calloc(
        num_spinors, sizeof(struct irp_irpa_eorb));
    int size_irp_irpa_eorb = num_spinors * sizeof(struct irp_irpa_eorb);

    int nread;
    if (data->dirac_int_size == DIRAC_INT_4) {
        nread = unf_read(file, "c[i4],2i4[i4],i4", buf_irp_irpa_eorb, &size_irp_irpa_eorb, ibspi,
                                 &num_spinors, norb, &invsym, &nbsymrp);
        if (nread != 4) {
            free(ibspi);
            free(buf_irp_irpa_eorb);
            return EXIT_FAILURE;
        }
    }
    else {
        int64_t *ibspi_8 = (int64_t *) calloc(data->num_spinors, sizeof(int64_t));
        int64_t nbsymrp_8;
        int64_t norb_8[2];

        struct irp_irpa_eorb_8 {
            int64_t irp;
            int64_t irpa;
            double eorb;
        };

        struct irp_irpa_eorb_8 *buf_irp_irpa_eorb_8 = (struct irp_irpa_eorb_8 *) calloc(
            data->num_spinors, sizeof(struct irp_irpa_eorb_8));
        int size_irp_irpa_eorb_8 = data->num_spinors * sizeof(struct irp_irpa_eorb_8);

        nread = unf_read(file, "c[i4],2i8[i4],i8", buf_irp_irpa_eorb_8, &size_irp_irpa_eorb_8, ibspi_8,
                                 &num_spinors, norb_8, &invsym, &nbsymrp_8);

        if (nread != 4) {
            free(ibspi);
            free(ibspi_8);
            free(buf_irp_irpa_eorb);
            free(buf_irp_irpa_eorb_8);
            return EXIT_FAILURE;
        }

        for (int i = 0; i < data->num_spinors; i++) {
            buf_irp_irpa_eorb[i].irp = (int32_t) buf_irp_irpa_eorb_8[i].irp;
            buf_irp_irpa_eorb[i].irpa = (int32_t) buf_irp_irpa_eorb_8[i].irpa;
            buf_irp_irpa_eorb[i].eorb = buf_irp_irpa_eorb_8[i].eorb;
        }

        free(ibspi_8);
        free(buf_irp_irpa_eorb_8);
    }

    data->occ_numbers = (int *) calloc(data->num_spinors, sizeof(int));
    data->spinor_irreps = (int *) calloc(data->num_spinors, sizeof(int));
    data->spinor_energies = (double *) calloc(data->num_spinors, sizeof(double));

    for (int i = 0; i < data->num_spinors; i++) {
        data->spinor_irreps[i] = buf_irp_irpa_eorb[i].irpa - 1;
        data->spinor_energies[i] = buf_irp_irpa_eorb[i].eorb;
        int irp = buf_irp_irpa_eorb[i].irp;
        if (fermion_irrep_occs[irp - 1] > 0) {
            fermion_irrep_occs[irp - 1] -= 1;
            data->occ_numbers[i] = 1;
        }
    }

    free(ibspi);
    free(buf_irp_irpa_eorb);

    return EXIT_SUCCESS;
}


int mrconee_read_fock(unf_file_t *file, mrconee_data_t *data)
{
    data->fock = (double _Complex *) calloc(data->num_spinors * data->num_spinors, sizeof(double _Complex));
    int32_t fock_size = data->num_spinors * data->num_spinors;

    int nread = unf_read(file, "z8[i4]", data->fock, &fock_size);
    if (nread != 1) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


void free_mrconee_data(mrconee_data_t *data)
{
    if (data->occ_numbers) {
        free(data->occ_numbers);
    }

    if (data->spinor_energies) {
        free(data->spinor_energies);
    }

    if (data->spinor_irreps) {
        free(data->spinor_irreps);
    }

    if (data->irrep_names) {
        for (int i = 0; i < data->num_irreps; i++) {
            char *name = data->irrep_names[i];
            if (name) {
                free(name);
            }
        }
        free(data->irrep_names);
    }

    if (data->point_group) {
        free(data->point_group);
    }

    if (data->mult_table) {
        free(data->mult_table);
    }

    if (data->fock) {
        free(data->fock);
    }

    free(data);
}


int test_dirac_integer_size(char *path)
{
    unf_file_t *file = unf_open(path, "r", UNF_ACCESS_SEQUENTIAL);
    if (file == NULL) {
        return -1;
    }

    double enuc;
    double scf_energy;
    int32_t num_spinors;
    int32_t breit;
    int32_t invsym;
    int32_t nz_arith;
    int32_t is_spinfree;
    int32_t norb_total;

    unf_read(file, "2i4,r8,4i4,r8", &num_spinors, &breit, &enuc, &invsym, &nz_arith, &is_spinfree,
                     &norb_total, &scf_energy);

    if (!(invsym == 1 || invsym == 2)) {
        return DIRAC_INT_8;
    }

    if (!(nz_arith == 1 || nz_arith == 2 || nz_arith == 4)) {
        return DIRAC_INT_8;
    }

    if (!(is_spinfree == -1 || is_spinfree == 0 || is_spinfree == 1)) {
        return DIRAC_INT_8;
    }

    if (norb_total < 0) {
        return DIRAC_INT_8;
    }

    return DIRAC_INT_4;
}


void print_mrconee_data(FILE *out, mrconee_data_t *data)
{
    fprintf(out, " size of integer numbers in DIRAC                   %d bytes\n", data->dirac_int_size);
    fprintf(out, " number of spinors                                  %d\n", data->num_spinors);
    fprintf(out, " core energy (inactive energy + nuclear repulsion)  %.12f a.u.\n", data->nuc_rep_energy);
    fprintf(out, " total SCF energy                                   %.12f a.u.\n", data->scf_energy);
    fprintf(out, " double group type                                  ");
    if (data->group_arith == 1) {
        fprintf(out, "real\n");
    }
    else if (data->group_arith == 2) {
        fprintf(out, "complex\n");
    }
    else if (data->group_arith == 4) {
        fprintf(out, "quaternion\n");
    }
    else {
        fprintf(out, "unknown\n");
    }
    fprintf(out, " spin-free                                          %s\n", data->is_spinfree ? "yes" : "no");
    fprintf(out, " Abelian subgroup                                   %s\n",
            data->point_group ? data->point_group : "n/a");
    fprintf(out, " totally symmetric irrep                            %s\n",
            data->irrep_names ? data->irrep_names[data->totally_sym_irrep] : "n/a");
    fprintf(out, " number of irreps in the Abelian subgroup           %d\n", data->num_irreps);
    fprintf(out, "\n");

    /*fprintf(out, " spinors info:\n");
    fprintf(out, " -----------------------------------------------------\n");
    fprintf(out, "   no       irrep     occ      one-electron energy    \n");
    fprintf(out, " -----------------------------------------------------\n");
    for (int i = 0; i < data->num_spinors; i++) {
        char *irrep_name = data->irrep_names[data->spinor_irreps[i]];
        fprintf(out, " %4d%12s%8d%25.8f\n", i + 1, irrep_name, data->occ_numbers[i], data->spinor_energies[i]);
    }
    fprintf(out, " -----------------------------------------------------\n");

    fprintf(out, "\n");*/
}


static void new_detect_dirac_point_group(char **rep_names, char *group_name, int *fully_sym_irrep)
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
void new_rename_irreps_dirac_to_expt(int nsym, char **rep_names)
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
    if (new_rep_name_exists(nsym, rep_names, "Ag a") &&
        new_rep_name_exists(nsym, rep_names, "Au a") &&
        new_rep_name_exists(nsym, rep_names, "Ag b") &&
        new_rep_name_exists(nsym, rep_names, "Au b") &&
        !new_rep_name_exists(nsym, rep_names, "Bg a") && // not C2h
        !new_rep_name_exists(nsym, rep_names, "B3ua") // not D2h
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
    if (new_rep_name_exists(nsym, rep_names, "Ag a") &&
        new_rep_name_exists(nsym, rep_names, "Au a") &&
        new_rep_name_exists(nsym, rep_names, "Ag b") &&
        new_rep_name_exists(nsym, rep_names, "Au b") &&
        new_rep_name_exists(nsym, rep_names, "Bg a") &&
        new_rep_name_exists(nsym, rep_names, "Bu a") &&
        new_rep_name_exists(nsym, rep_names, "Bg b") &&
        new_rep_name_exists(nsym, rep_names, "Bu b")) {
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


int new_rep_name_exists(int nrep, char **list, char *query)
{
    for (int i = 0; i < nrep; i++) {
        if (strcmp(list[i], query) == 0) {
            return 1;
        }
    }
    return 0;
}

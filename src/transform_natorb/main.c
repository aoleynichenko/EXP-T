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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "read_natural_spinors.h"
#include "read_dirac_output.h"
#include "print_vectors.h"
#include "mapping.h"
#include "transform.h"


void print_header();

void print_footer();

void parse_arguments(int argc, char **argv, char **path_natorb_file, char **path_dirac_output);

void print_spinor_energies(int n_spinors, double *spinor_energies, char *comment);

void print_natural_spinors_info(int n_nat_spinors, int n_mol_spinors, double *occ_numbers, double complex *nto_coeffs);


int main(int argc, char **argv)
{
    /*
     * parse command line arguments
     */

    char *path_natorb_file = NULL;
    char *path_dirac_output = NULL;
    parse_arguments(argc, argv, &path_natorb_file, &path_dirac_output);

    print_header();

    printf(" > path to the file containing natural transition spinors:\n");
    printf("   %s\n\n", path_natorb_file);
    printf(" > path to the DIRAC output file containing molecular spinors in the AO basis:\n");
    printf("   %s\n", path_dirac_output);
    printf("\n");

    /*
     * read natural spinors expanded in the basis of molecular spinors
     */

    int n_act_spinors = 0;
    int n_nat_spinors = 0;
    double *act_spinor_energies = NULL;
    double *occ_numbers = NULL;
    double complex *nat_spinors_coeffs = NULL;

    read_natural_spinors(path_natorb_file,
                         &n_act_spinors, &act_spinor_energies,
                         &n_nat_spinors, &occ_numbers, &nat_spinors_coeffs
    );

    print_spinor_energies(n_act_spinors, act_spinor_energies, "energies of active molecular spinors");
    print_natural_spinors_info(n_nat_spinors, n_act_spinors, occ_numbers, nat_spinors_coeffs);

    /*
     * read molecular spinors from the DIRAC output file
     */

    int n_mol_spinors = 0;
    int n_basis_fun = 0;
    double complex *cmo_alpha = NULL;
    double complex *cmo_beta = NULL;
    double *mol_spinor_energies = NULL;
    char **basis_fun_labels = NULL;

    read_dirac_output(path_dirac_output, &n_mol_spinors, &n_basis_fun, &mol_spinor_energies,
                      &cmo_alpha, &cmo_beta, &basis_fun_labels);

    print_spinor_energies(n_mol_spinors, mol_spinor_energies, "energies of molecular spinors");

    /*
     * construct mapping: active spinor index -> molecular spinor index
     */

    int *act_to_all = NULL;
    int *kramers_conj = NULL;
    construct_spinor_mapping(n_act_spinors, act_spinor_energies,
                             n_mol_spinors, mol_spinor_energies, &act_to_all, &kramers_conj);

    /*
     * transform natural spinors to the AO basis
     * and print quaternion coefficients of their expansions to the standard output.
     * DIRAC output format will be simulated.
     */

    double complex *cno_alpha = (double complex *) calloc(n_nat_spinors * n_basis_fun, sizeof(double complex));
    double complex *cno_beta = (double complex *) calloc(n_nat_spinors * n_basis_fun, sizeof(double complex));

    transform(n_nat_spinors, n_act_spinors, n_basis_fun,
              act_to_all, kramers_conj, nat_spinors_coeffs,
              cmo_alpha, cmo_beta, cno_alpha, cno_beta);

    print_vectors(stdout, n_nat_spinors, n_basis_fun, occ_numbers, cno_alpha, cno_beta, basis_fun_labels);

    /*
     * deallocate memory and exit
     */

    free(act_spinor_energies);
    free(occ_numbers);
    free(nat_spinors_coeffs);
    free(mol_spinor_energies);
    for (int i = 0; i < n_basis_fun; i++) {
        free(basis_fun_labels[i]);
    }
    free(basis_fun_labels);
    free(cmo_alpha);
    free(cmo_beta);
    free(cno_alpha);
    free(cno_beta);
    free(act_to_all);
    free(kramers_conj);

    print_footer();
}


void print_header()
{
    char buf[100];

    sprintf(buf, "%s %s", __DATE__, __TIME__);

    printf("\n");
    printf("\t**********************************************************************************\n"
           "\t**                                                                              **\n"
           "\t**           program for transformation of natural (transition) spinors         **\n"
           "\t**                        to the atomic orbital basis                           **\n"
           "\t**                                                                              **\n"
           "\t**                         2021-2022 A. Oleynichenko                            **\n"
           "\t**                       compiled %-30s                **\n"
           "\t**                                                                              **\n"
           "\t**********************************************************************************\n\n",
           buf);

    time_t lt = time(NULL);
    struct tm *ptr = localtime(&lt);
    printf(" > started at %s\n", asctime(ptr));
}


void print_footer()
{
    time_t lt = time(NULL);
    struct tm *ptr = localtime(&lt);
    printf(" > finished at %s\n", asctime(ptr));
}


void parse_arguments(int argc, char **argv, char **path_natorb_file, char **path_dirac_output)
{
    if (argc != 3) {
        printf("\n");
        printf(" usage: %s <path-to-nts-file> <path-to-dirac-output>\n", argv[0]);
        printf("\n");
        printf("  where:\n");
        printf("  <path-to-nts-file>     path to the file containing natural transition spinors\n");
        printf("                         expanded in the basis of molecular spinors\n");
        printf("  <path-to-dirac-output> path to the DIRAC output file containing molecular spinors\n");
        printf("                         expanded in the basis of atomic orbitals\n");
        printf("\n");
        exit(1);
    }

    *path_natorb_file = argv[1];
    *path_dirac_output = argv[2];
}


void print_spinor_energies(int n_spinors, double *spinor_energies, char *comment)
{
    printf(" > %s\n", comment);

    for (int ispinor = 0; ispinor < n_spinors; ispinor++) {
        printf("%6d%20.12f\n", ispinor + 1, spinor_energies[ispinor]);
    }

    printf("\n");
}


void print_natural_spinors_info(int n_nat_spinors, int n_mol_spinors, double *occ_numbers, double complex *nto_coeffs)
{
    printf(" > (quasi) natural orbitals\n");
    printf("\n");

    for (int i_nts = 0; i_nts < n_nat_spinors; i_nts++) {
        printf(" [%d] occ number = %.8f\n\n", i_nts + 1, occ_numbers[i_nts]);
        for (int ispinor = 0; ispinor < n_mol_spinors; ispinor++) {
            double coef = nto_coeffs[n_mol_spinors * i_nts + ispinor];
            printf("%6d%24.12e%24.12e\n", ispinor + 1, creal(coef), cimag(coef));
        }
        printf("\n");
    }
}


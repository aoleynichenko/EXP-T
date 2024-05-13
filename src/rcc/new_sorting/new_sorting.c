/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2024 The EXP-T developers.
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


#include "new_sorting.h"

#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "mdprop.h"
#include "mrconee.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"


void check_irrep_names();


void new_sorting(cc_options_t *options)
{
    /*
     * basic information + Fock matrix
     */
    mrconee_data_t *mrconee_data = read_mrconee(options->integral_file_1);
    print_mrconee_data(stdout, mrconee_data);

    options->enuc = mrconee_data->nuc_rep_energy;
    options->escf = mrconee_data->scf_energy;

    /*
     * one-electron property matrices (calculated by DIRAC)
     */
    read_mdprop(options->integral_file_prop, mrconee_data);

    /*
     * setup arithmetic: real or complex
     */
    if (mrconee_data->group_arith == 1 || mrconee_data->is_spinfree) {
        arith = CC_ARITH_REAL;
    }
    else {
        arith = CC_ARITH_COMPLEX;
    }

    // what about arithmetic recommended by the user ?
    if (options->recommended_arith == CC_ARITH_COMPLEX) {
        arith = CC_ARITH_COMPLEX;
    }

    if (arith == CC_ARITH_COMPLEX) {
        WORKING_TYPE = CC_DOUBLE_COMPLEX;
        SIZEOF_WORKING_TYPE = sizeof(double complex);
    }
    else {
        WORKING_TYPE = CC_DOUBLE;
        SIZEOF_WORKING_TYPE = sizeof(double);
    }

    /*
     * setup symmetry.
     * "Multiplication table" ( = direct product table for Abelian groups).
     * for (quasi) infinite groups Cinfv and Dinfh multiplication tables provided by DIRAC
     * are not sufficient for CCSDT and/or 0h3p calculations (too few irreps). So for these
     * groups additional irreps will be generated.
     */

    if (strcmp(mrconee_data->point_group, "Cinfv") == 0 || strcmp(mrconee_data->point_group, "Dinfh") == 0) {
        int new_num_irreps;
        int max_omega_x2;
        char **new_irrep_names = NULL;
        int new_irrep_a1;
        int *new_mult_table = NULL;

        if (strcmp(mrconee_data->point_group, "Cinfv") == 0) {
            // Cinfv
            max_omega_x2 = CC_MAX_NUM_IRREPS / 2 - 2;
            new_irrep_names = generate_irreps_Cinfv(max_omega_x2, &new_num_irreps);
            new_irrep_a1 = search_string("0", new_irrep_names, new_num_irreps);
            new_mult_table = construct_direct_product_table(new_num_irreps, new_irrep_names, multiply_irreps_Cinfv);
        }
        else {
            // Dinfh
            max_omega_x2 = CC_MAX_NUM_IRREPS / 4 - 2;
            new_irrep_names = generate_irreps_Dinfh(max_omega_x2, &new_num_irreps);
            new_irrep_a1 = search_string("0g", new_irrep_names, new_num_irreps);
            new_mult_table = construct_direct_product_table(new_num_irreps, new_irrep_names, multiply_irreps_Dinfh);
        }

        printf("\n");
        printf(" %s group is to be extended:\n", mrconee_data->point_group);
        printf(" max 2x|Omega| = %d\n", max_omega_x2);
        printf(" number of irreps = %d\n", new_num_irreps);
        printf("\n");

        // for each spinor: update irrep enumeration
        for (int i = 0; i < mrconee_data->num_spinors; i++) {
            int old_irrep_no = mrconee_data->spinor_irreps[i];
            char *old_irrep_name = mrconee_data->irrep_names[old_irrep_no];
            int new_irrep_no = search_string(old_irrep_name, new_irrep_names, new_num_irreps);
            mrconee_data->spinor_irreps[i] = new_irrep_no;
        }

        // update list of irreps and multiplication table
        for (int i = 0; i < mrconee_data->num_irreps; i++) {
            free(mrconee_data->irrep_names[i]);
        }
        free(mrconee_data->irrep_names);
        free(mrconee_data->mult_table);

        mrconee_data->num_irreps = new_num_irreps;
        mrconee_data->irrep_names = new_irrep_names;
        mrconee_data->mult_table = new_mult_table;
        mrconee_data->totally_sym_irrep = new_irrep_a1;
    }

    setup_symmetry(mrconee_data->group_arith, mrconee_data->point_group, mrconee_data->num_irreps,
                   mrconee_data->irrep_names, mrconee_data->totally_sym_irrep, mrconee_data->mult_table);

    /*
     * setup information about spinors and active space
     */
    //create_spinor_info_dirac(dirac_data.nspinors, dirac_data.irpamo, dirac_data.eorbmo, dirac_data.iocc);
    create_spinor_info(mrconee_data->num_spinors, mrconee_data->spinor_irreps,
        mrconee_data->spinor_energies, mrconee_data->occ_numbers);
    create_spinor_blocks(options->tile_size);
    setup_occupation_numbers(options, mrconee_data->num_spinors, spinor_info);
    setup_active_space(options);
    setup_fast_access_spinor_lists();

    print_symmetry_info();
    print_spinor_info_table();

    /*
     * additional checks:
     * names of irreps in some directives of the input file can be incorrect
     */
    check_irrep_names();

    /*
     * save data read from mrconee
     */
    options->mrconee_data = mrconee_data;
}

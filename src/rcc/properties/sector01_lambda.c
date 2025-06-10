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
 * Lambda equations in the 0h1p Fock space sector.
 *
 * Basic theory can be found in:
 * K. R. Shamasundar, S. Asokan, S. Pal
 * A constrained variational approach for energy derivatives in
 * Fock-space multireference coupled-cluster theory.
 * J. Chem. Phys. 120, 6381 (2004)
 * doi: 10.1063/1.1652436
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../models/ccutils.h"
#include "engine.h"
#include "../models/diis.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/mvcoef.h"


void construct_lambda_singles_0h1p();

void construct_lambda_doubles_0h1p();

void construct_lambda_triples_0h1p();

void sector_0h1p_lambda_equations_guess();

void construct_disconnected_rank2_rank2(char *src1_name, char *src2_name, char *dst);

void sector_0h0p_analytic_density_matrix_lambda();

void calculate_natural_spinors_and_occ_numbers(char *path_to_dm_file, char *path_to_nat_spinors_file);

void calc_L00_perturbative_triples_correction();

void diagram_clear_off_diagonal(char *diagram_name);

void construct_one_part_active_density_matrix_diagram_0h1p(
    char *diagram_name,
    char *bra_irrep_name, int bra_state_number,
    char *ket_irrep_name, int ket_state_number
);


/**
 * Iterative solution of lambda equations.
 * Finally density matrix is constructed and written to disk.
 */
int sector_0h1p_lambda_equations(char *irrep_name, int state_number)
{
    int triples = triples_enabled();
    cc_opts->curr_in_lambda_equations = 1;

    printf("\n");
    printf("\n");
    printf(" Lambda equations in the 0h1p sector for the state %d (%s)\n", state_number + 1, irrep_name);
    printf("\n");

    construct_one_part_active_density_matrix_diagram_0h1p("dm1", irrep_name, state_number, 0, 0);

    /*
     * initialize L1 and L2 amplitudes (diagrams L00_1c and L00_2c)
     */
    sector_0h1p_lambda_equations_guess();

    /*
     * for the CCSD(T) model additional one- and two-particle contributions
     * to lambda equations must be calculated and then added at each iteration
     */
    /*if (cc_opts->cc_model == CC_MODEL_CCSD_T3) {
        calc_L00_perturbative_triples_correction();
    }*/

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
        0, 1,
        "L01_1c", "L01_1nw",
        "L01_2c", "L01_2nw",
        NULL, NULL, NULL,
        construct_lambda_singles_0h1p,
        construct_lambda_doubles_0h1p,
        construct_lambda_triples_0h1p,
        NULL
    );
    cc_opts->curr_in_lambda_equations = 0;

    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze lambda amplitudes and save them to disk
     */
    print_cluster_operator_analysis(0, 1, "L01_1c", "L01_2c", NULL);
    save_cluster_amplitudes(0, 1, "L01_1c", "L01_2c", NULL, NULL);
    if (cc_opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(0, 1, "L01_1c", "L01_2c", NULL);
    }

    return EXIT_SUCCESS;
}


/**
 * Constructs initial guess for Lambda amplitudes
 * (or read them from disk)
 */
void sector_0h1p_lambda_equations_guess()
{
    int calc_L1 = 1;
    int calc_L2 = 1;
    int calc_L3 = 1;

    if (cc_opts->reuse_amplitudes[0][1]) {
        printf(" Trying to read amplitudes from disk ...\n");
        if (diagram_read_binary("L01_1c.dg") != NULL) {
            printf(" L1 amplitudes successfully read from disk\n");
            calc_L1 = 0;
        }
        if (diagram_read_binary("L01_2c.dg") != NULL) {
            printf(" L2 amplitudes successfully read from disk\n");
            calc_L2 = 0;
        }
        printf("\n");
        /*if (triples) {
            if (diagram_read_binary("t3c.dg") != NULL) {
                printf(" T3 amplitudes successfully read from disk\n");
                calc_t3 = 0;
            }
            else {
                printf(" T3 amplitudes will be calculated\n");
            }
        }*/
    }

    if (calc_L1) {
        tmplt("L01_1c", "pp", "01", "12", NOT_PERM_UNIQUE);
//        copy("pv", "L01_1c");
//        diveps("L01_1c");
    }

    if (calc_L2) {
        tmplt("L01_2c", "ppph", "0010", "1234", NOT_PERM_UNIQUE);
//        copy("ppvh", "L01_2c");
//        diveps("L01_2c");
    }
}


void construct_lambda_singles_0h1p()
{
    tmplt("L01_1nw", "pp", "01", "12", NOT_PERM_UNIQUE);
}


void construct_lambda_doubles_0h1p()
{
    tmplt("L01_2nw", "ppph", "0010", "1234", NOT_PERM_UNIQUE);
    tmplt("cross", "pp", "11", "12", NOT_PERM_UNIQUE); // intermediate diagram

    dg_stack_pos_t pos = get_stack_pos();

    // L01_D8
    /*reorder("L01_2c", "r1", "3412");
    mult("r1", "pp", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("L01_2nw", 1.0, "r3");
    restore_stack_pos(pos);

    // L01_D9
    reorder("hh", "r1", "21");
    mult("L01_2c", "r1", "r2", 1);
    update("L01_2nw", -1.0, "r2");
    restore_stack_pos(pos);

    // L01_D10
    reorder("L01_2c", "r1", "1342");
    reorder("hpph", "r2", "2413");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    update("L01_2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L01_D11
    reorder("L01_2c", "r1", "3412");
    mult("pppp", "r1", "r2", 2);
    update("L01_2nw", 0.5, "r2");
    restore_stack_pos(pos);

    // L01_D21
    reorder("pphh", "r1", "3412");
    mult("r1", "t2c", "r2", 3);
    mult("L01_2c", "r2", "r3", 1);
    update("L01_2nw", -0.5, "r3");
    restore_stack_pos(pos);

    // L01_D22
    reorder("t2c", "r1", "3412");
    reorder("L01_2c", "r2", "3412");
    mult("pphh", "r1", "r3", 3);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    perm("r5", "(12)");
    update("L01_2nw", -0.5, "r5");
    restore_stack_pos(pos);

    // L01_D24
    reorder("L01_2c", "r1", "3412");
    mult("r1", "t2c", "r2", 3);
    mult("pphh", "r2", "r3", 1);
    reorder("r3", "r4", "2143"); // interchange electrons
    update("L01_2nw", -0.5, "r4");
    restore_stack_pos(pos);

    // L01_D25
    reorder("L01_2c", "r1", "1342");
    reorder("t2c", "r2", "1324");
    reorder("pphh", "r3", "2431");
    mult("r3", "r2", "r4", 2);
    mult("r1", "r4", "r5", 2);
    reorder("r5", "r6", "1324");
    perm("r6", "(12)");
    update("L01_2nw", 1.0, "r6");
    restore_stack_pos(pos);

    // L01_D27
    reorder("L01_2c", "r1", "3412");
    mult("r1", "t2c", "r2", 2);
    mult("pphh", "r2", "r3", 2);
    update("L01_2nw", 0.25, "r3");
    restore_stack_pos(pos);*/

    /*
     * contributions from the energy term
     */

    //clear("L01_1c");
    //clear("L01_2c");

    // intermediate diagram - "cross"
    copy("dm1", "cross");

    /*reorder("L01_2c", "r1", "3412");
    mult("s2c", "r1", "r2", 3);
    update("cross", -0.5, "r2");
    restore_stack_pos(pos);*/

    // L01_DE1
    construct_disconnected_rank2_rank2("cross", "ph", "r1");
    tmplt("r2", "ppph", "0010", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r1", "r2");
    perm("r2", "(12)");
    //update("L01_2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // L01_DE2
    reorder("cross", "r1", "21");
    reorder("ppvh", "r2", "1243");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    //update("L01_2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L01_DE3 ("folded")
    /*reorder("veff01", "r1", "21");
    reorder("L01_2c", "r2", "1243");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    update("L01_2nw", -1.0, "r4");
    restore_stack_pos(pos);*/

    /*clear("L01_1nw");
    clear("L01_2nw");*/
}


void construct_lambda_triples_0h1p()
{

}
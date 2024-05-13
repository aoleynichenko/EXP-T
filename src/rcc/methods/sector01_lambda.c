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

#include "ccutils.h"
#include "engine.h"
#include "diis.h"
#include "options.h"
#include "spinors.h"
#include "symmetry.h"
#include "utils.h"
#include "../heff/mvcoef.h"


extern int diveps_lambda;

void construct_lambda_singles_0h1p();

void construct_lambda_doubles_0h1p();

void construct_lambda_triples_0h1p();

void sector_0h1p_lambda_equations_guess();

void construct_disconnected_rank2_rank2(char *src1_name, char *src2_name, char *dst);

void sector_0h0p_analytic_density_matrix_lambda();

void sector_0h0p_calculate_natural_spinors();

void calc_L00_perturbative_triples_correction();

void diagram_clear_off_diagonal(char *diagram_name);

void construct_active_space_density_matrix_diagram_0h1p(char *diagram_name, char *irrep_name, int state_number);


/**
 * Iterative solution of lambda equations.
 * Finally density matrix is constructed and written to disk.
 */
int sector_0h1p_lambda_equations(char *irrep_name, int state_number)
{
    int triples = triples_enabled();

    printf("\n");
    printf("\n");
    printf(" Lambda equations in the 0h1p sector for the state %d (%s)\n", state_number + 1, irrep_name);
    printf("\n");

    construct_active_space_density_matrix_diagram_0h1p("dm1", irrep_name, state_number);

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
    diveps_lambda = 1;
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
    diveps_lambda = 0;
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze lambda amplitudes and save them to disk
     */
    print_cluster_operator_analysis(0, 0, "L01_1c", "L01_2c", NULL);
    save_cluster_amplitudes(0, 0, "L01_1c", "L01_2c", NULL, NULL);
    if (cc_opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(0, 0, "L01_1c", "L01_2c", NULL);
    }

    return EXIT_SUCCESS;
}


void construct_active_space_density_matrix_diagram_0h1p(char *diagram_name, char *irrep_name, int state_number)
{
    int ms_dim = 0;
    slater_det_t *det_list = NULL;
    double eigenvalue = 0.0;
    double exc_energy_cm = 0.0;
    double complex *coef_left = NULL;
    double complex *coef_right = NULL;

    mvcoef_read_vectors_unformatted_for_state(
        0, 1, "MVCOEF01", irrep_name, state_number,
        &ms_dim, &det_list, &eigenvalue, &exc_energy_cm, &coef_left, &coef_right
    );

    printf("mod sp dim = %d\n", ms_dim);
    printf("eigenvalue = %.8f\n", eigenvalue);
    printf("exc ene cm = %.2f\n", exc_energy_cm);

    double left_norm = 0.0;
    printf("left vector:\n");
    for (int i = 0; i < ms_dim; i++) {
        printf("%20.12f%20.12f", creal(coef_left[i]), cimag(coef_left[i]));
        print_slater_det(stdout, 0, 1, det_list + i);
        left_norm += cabs(coef_left[i]) * cabs(coef_left[i]);
    }
    printf("<left|left> = %20.12f\n", left_norm);

    double right_norm = 0.0;
    printf("right vector:\n");
    for (int i = 0; i < ms_dim; i++) {
        printf("%20.12f%20.12f", creal(coef_right[i]), cimag(coef_right[i]));
        print_slater_det(stdout, 0, 1, det_list + i);
        right_norm += cabs(coef_right[i]) * cabs(coef_right[i]);
    }
    printf("<right|right> = %20.12f\n", right_norm);

    double complex dot_prod = 0.0;
    for (int i = 0; i < ms_dim; i++) {
        dot_prod += conj(coef_left[i]) * coef_right[i];
    }
    printf("dot = %20.12f%20.12f\n", creal(dot_prod), cimag(dot_prod));

    tmplt(diagram_name, "pp", "11", "12", NOT_PERM_UNIQUE);
    diagram_t *diag_dm = diagram_stack_find(diagram_name);

    for (int i = 0; i < ms_dim; i++) {
        for (int j = 0; j < ms_dim; j++) {
            slater_det_t *det_i = det_list + i;
            slater_det_t *det_j = det_list + j;
            moindex_t p = det_i->indices[0];
            moindex_t q = det_j->indices[0];
            double complex C_p = coef_left[i];
            double complex C_q = coef_right[j];
            double complex gamma_qp = conj(C_p) * C_q;
            //printf("gamma[%d,%d] = %20.12f%20.12f\n", q+1, p+1, creal(gamma_qp), cimag(gamma_qp));

            int idx[2];
            idx[0] = q;
            idx[1] = p;
            diagram_set(diag_dm, gamma_qp, idx);
        }
    }

    prt(diagram_name);

    tmplt("heff01", "pp", "11", "12", NOT_PERM_UNIQUE);
    diagram_t *dg_h0 = diagram_stack_find("heff01");
    for (int i = 0; i < ms_dim; i++) {
        slater_det_t *det_i = det_list + i;
        moindex_t index = det_i->indices[0];

        int repno, occ, active;
        double eps = 0.0;
        get_spinor_info(index, &repno, &occ, &active, &eps);

        int idx[] = {index, index};
        diagram_set(dg_h0, eps, idx);
    }

    update("heff01", 1.0, "veff01");

    //prt(diagram_name);

    double E = scalar_product("N", "N", diagram_name, "heff01");

    printf("Re E = %24.16f\n", creal(E));
    printf("Im E = %24.16f\n", cimag(E));
    printf("eigv = %24.16f\n", eigenvalue);

    /*int num_irreps = 0;
    struct mv_block mv_blocks[CC_MAX_NUM_IRREPS];
    mvcoef_read_vectors_unformatted(0, 1, NULL, &num_irreps, mv_blocks);

    for (size_t irep = 0; irep < num_irreps; irep++) {
        struct mv_block *model_vectors = mv_blocks + irep;
        printf("[%d] %d\n", irep, model_vectors->ms_size);
    }

    for (size_t irep = 0; irep < num_irreps; irep++) {
        mvblock_free(mv_blocks + irep);
    }*/
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
    reorder("L01_2c", "r1", "3412");
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
    restore_stack_pos(pos);

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
    update("L01_2nw", 1.0, "r2");
    restore_stack_pos(pos);

    // L01_DE2
    reorder("cross", "r1", "21");
    reorder("ppvh", "r2", "1243");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    update("L01_2nw", 1.0, "r4");
    restore_stack_pos(pos);

    // L01_DE3 ("folded")
    /*reorder("veff01", "r1", "21");
    reorder("L01_2c", "r2", "1243");
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    update("L01_2nw", 1.0, "r4");
    restore_stack_pos(pos);*/

    /*clear("L01_1nw");
    clear("L01_2nw");*/
}


void construct_lambda_triples_0h1p()
{

}
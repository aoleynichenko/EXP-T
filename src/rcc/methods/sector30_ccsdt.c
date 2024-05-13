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
 * Multireference Fock-space coupled cluster method for model spaces
 * with three holes (Fock space sector 3h0p)
 *
 * Symbolic names used for diagrams:
 * k3c   current approximation to T^{3h0p}_3 amplitudes
 * k3nw  next (new) approximation to T^{3h0p}_3 amplitudes
 *
 * The CCSD model does not require solution of any amplitude equations.
 * Heff is constructed in one step from molecular integrals and Singles and
 * Doubles amplitudes from lower sectors:
 *   Heff(3h0p,CCSD) = f(T1,T2)
 *
 * This code was obtained by modification of the CC(0,0) program:
 * all we need is to flip particle creation lines to turn them to valence hole
 * annihilation lines. See Kaldor, J. Comp. Chem. V. 8, P.448 (1987) for details.
 */

#include "methods.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ccutils.h"
#include "engine.h"
#include "diis.h"
#include "heff.h"
#include "options.h"
#include "sort.h"

void sort_integrals_3h0p();

void init_amplitudes_3h0p();

void const_terms_3h0p();

void construct_triples_3h0p();

void construct_folded_3h0p();

void heff30_ccsd();

void heff30_triples();

void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence);

void diag_heff_3h0p_ccsd_k_ij_a_bc(int pt_order);
void diag_heff_3h0p_ccsd_i_jk_c_ab(int pt_order);
void diag_heff_3h0p_ccsd_ijk_c_ab(int pt_order);
void diag_heff_3h0p_ccsd_k_ij_abc(int pt_order);

void diag_T3_3h0p_i_jk_c_ab(int pt_order);
void diag_T3_3h0p_ijk_c_ab(int pt_order);
void diag_T3_3h0p_k_ij_a_bc(int pt_order);
void diag_T3_3h0p_k_ij_abc(int pt_order);
void diag_T3_3h0p_c_ab(int pt_order);
void diag_T3_3h0p_a_bc(int pt_order);
void diag_T3_3h0p_abc(int pt_order);
void diag_T3_3h0p_k_ij_c_ab(int pt_order);
void diag_T3_3h0p_ijk(int pt_order);
void diag_T3_3h0p_i_jk(int pt_order);
void diag_T3_3h0p_k_ij(int pt_order);

/**
 * Solves amplitude equations and find correlation energy for ground state
 * (Fock space sector 3h0p).
 */
int sector_3h0p(cc_options_t *opts)
{
    print_sector_banner(3, 0);
    opts->curr_sector_h = 3;
    opts->curr_sector_p = 0;

    sort_integrals_3h0p();

    /* CCSD, CCSDT-1 */
    if (cc_opts->cc_model < CC_MODEL_CCSDT_2) {
        heff30_ccsd();
        if (cc_opts->cc_model >= CC_MODEL_CCSD_T3) {
            heff30_triples();
        }
    }
        /* CCSDT-2, CCSDT-3, CCSDT */
    else {

        // setup shifts for the IH-like technique by A. V. Zaitsevskii
        if (intham_imms_in_sector(3, 0)) {
            intham_imms_setup(3, 0);
        }

        const_terms_3h0p();
        init_amplitudes_3h0p();
        intruder_state_analysis(NULL, NULL, "k3c");

        /*
         * main loop
         */
        int success = solve_amplitude_equations(
                3, 0,
                NULL, NULL,
                NULL, NULL,
                "k3c", "k3nw",
                "veff30",
                NULL,
                NULL,
                construct_triples_3h0p,
                construct_folded_3h0p
        );
        if (success == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }

        /*
         * analyze converged amplitudes and save them to disk
         */
        print_cluster_operator_analysis(3, 0, NULL, NULL, "k3c");
        save_cluster_amplitudes(3, 0, NULL, NULL, "k3c", "veff30");
        if (opts->do_flush_amplitudes_txt) {
            save_cluster_amplitudes_formatted(3, 0, NULL, NULL, "k3c");
        }
    }

    // следующий код загружает всю диаграмму в память для построения эфф гам-на:
    diagram_t *dg = diagram_stack_find("veff30");
    for (size_t isb = 0; isb < dg->n_blocks; isb++) {
        block_t *sb = dg->blocks[isb];
        block_load(sb);
        sb->storage_type = CC_DIAGRAM_IN_MEM;
    }

    /*
     * construct and diagonalize effective Hamiltonian matrix,
     * analyze its eigenvectors & eigenvalues
     */
    heff_analysis(3, 0, "veff10", "veff20", "veff30");

    /*
     * Ionization potentials
     */
    double ion_pot_1 = cc_opts->ground_energy_1h0p;
    double ion_pot_2 = cc_opts->ground_energy_2h0p - cc_opts->ground_energy_1h0p;
    double ion_pot_3 = cc_opts->ground_energy_3h0p - cc_opts->ground_energy_2h0p;

    print_ionization_potential("Ionization potential 0h0p -> 1h0p", ion_pot_1);
    print_ionization_potential("Ionization potential 1h0p -> 2h0p", ion_pot_2);
    print_ionization_potential("Ionization potential 2h0p -> 3h0p", ion_pot_3);

    /*
     * calculate quasi-natural orbitals and model-space estimates of properties
     */
    calculate_model_space_properties(3, 0);
    calculate_model_space_natural_spinors(3, 0);

    return EXIT_SUCCESS;
}


/**
 * Sorting of integrals used in the CC(0h3p) amplitude equations.
 */
void sort_integrals_3h0p()
{
    /*request_sorting("vphp", "pphp", "1000", "1234");
    request_sorting("ppph", "ppph", "0000", "1234");
    request_sorting("vphh", "pphh", "1000", "1234");
    request_sorting("pvvv", "pppp", "0111", "1234");
    request_sorting("vvhv", "pphp", "1101", "1234");
    request_sorting("pphv", "pphp", "0001", "1234");
    request_sorting("ppvvr", "pppp", "0011", "3412");
    request_sorting("vphv", "pphp", "1001", "1234");
    request_sorting("ppvh", "ppph", "0010", "1234");
    request_sorting("vvvv", "pppp", "1111", "1234");
    request_sorting("pvph", "ppph", "0100", "1234");
    request_sorting("pvhv", "pphp", "0101", "1234");*/
    request_sorting("pggg", "phhh", "0111", "1234");
    request_sorting("pghh", "phhh", "0100", "1234");
    request_sorting("gghg", "hhhh", "1101", "1234");
    request_sorting("gghh", "hhhh", "1100", "1234");
    request_sorting("gphh", "hphh", "1000", "1234");
    request_sorting("pghg", "phhh", "0101", "1234");
    perform_sorting();

    if (cc_opts->print_level >= CC_PRINT_DEBUG) {
        diagram_stack_print();
    }
}


void const_terms_3h0p()
{
    double t1, t2;

    printf(" Construction of S^(3,0)-independent contributions to the FSCC-equations ...\n");
    timer_new_entry("const_s30", "Constant part of 3h0p amplitudes");
    timer_start("const_s30");

    tmplt("k3_0", "hhhhhh", "000111", "123456", IS_PERM_UNIQUE);

    t1 = abs_time();

    diag_T3_3h0p_i_jk_c_ab(PT_INF);
    diag_T3_3h0p_ijk_c_ab(PT_INF);
    diag_T3_3h0p_k_ij_a_bc(PT_INF);
    diag_T3_3h0p_k_ij_abc(PT_INF);
    diag_T3_3h0p_c_ab(PT_INF);
    diag_T3_3h0p_a_bc(PT_INF);
    diag_T3_3h0p_abc(PT_INF);
    diag_T3_3h0p_k_ij_c_ab(PT_INF);

    if (cc_opts->cc_model >= CC_MODEL_CCSDT_2) {
        /*diag_T3_0h3p_k_ij_a_bc();
        diag_T3_0h3p_k_ij_abc();
        diag_T3_0h3p_i_jk_c_ab();
        diag_T3_0h3p_ijk_c_ab();*/
    }
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        /*diag_T3_0h3p_k_ij_c_ab();
        diag_T3_0h3p_i_jk();
        diag_T3_0h3p_ijk();
        diag_T3_0h3p_k_ij();*/
    }
    t2 = abs_time();
    printf(" done in %.2f sec\n", t2 - t1);

    timer_stop("const_s30");
}


void init_amplitudes_3h0p()
{
    int calc_t3 = 1, calc_veff = 1;

    printf("\n Initialization of T{3h0p}_3 amplitudes ...\n");

    if (cc_opts->reuse_amplitudes[3][0]) {
        printf(" Trying to read amplitudes (sector 3h0p) from disk ...\n");
        if (diagram_read_binary("k3c.dg") != NULL) {
            printf(" T{3h0p}_3 amplitudes successfully read from disk\n");
            calc_t3 = 0;
        }
        else {
            printf(" T{3h0p}_3 amplitudes will be calculated\n");
        }
        if (diagram_read_binary("veff30.dg") != NULL) {
            printf(" Heff{3h0p} diagram successfully read from disk\n");
            calc_veff = 0;
        }
        else {
            printf(" Heff{3h0p} diagram will be calculated\n");
        }
    }

    if (calc_t3) {
        copy("k3_0", "k3c");
        closed("k3c", "veff30");
        diveps("k3c");
    }

    printf(" done\n\n");
}


// non-iterative construction of Heff at the CCSD level of theory
// diagrams are based on CCSDT-3 diagrams (lines were simply "reversed")
void heff30_ccsd()
{
    dg_stack_pos_t pos;

    timer_new_entry("30-CCSD-Heff", "3h0p -- CCSD part of eff Hamiltonian");
    timer_start("30-CCSD-Heff");

    printf("\n restrict valence ... ");
    restrict_valence("t1c", "t1c-g", "10", 0);
    restrict_valence("t2c", "t2c-g12", "1100", 0);
    restrict_valence("h2c", "h2c-g12", "1110", 0);
    restrict_valence("h2c", "h2c-g1", "1010", 0);
    restrict_valence("g2c", "g2c-g1", "1011", 0);
    printf("done\n");

    printf(" T1 and T2 contributions to Heff(3h0p)\n");

    tmplt("veff30", "hhhhhh", "111111", "123456", NOT_PERM_UNIQUE);

    diag_heff_3h0p_ccsd_k_ij_a_bc(PT_INF);
    diag_heff_3h0p_ccsd_i_jk_c_ab(PT_INF);
    diag_heff_3h0p_ccsd_ijk_c_ab(PT_INF);
    diag_heff_3h0p_ccsd_k_ij_abc(PT_INF);

    timer_stop("30-CCSD-Heff");
}


void diag_heff_3h0p_ccsd_k_ij_a_bc(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pggg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3e
        reorder("g2c", "r1", "3412");
        mult("r1", "pghh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        reorder("ppgg", "ppggr", "3412");
        mult("ppggr", "t1c-g", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("h1c", "r1", "21");
        mult("r1", "pghh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("g2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c-g", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("h1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c-g12", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|4/56)");
    update("veff30", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_heff_3h0p_ccsd_i_jk_c_ab(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("gghg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME && pt_order == PT_INF) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3a
        reorder("ph", "phr_", "21");
        reorder("h2c-g12", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("t2c-g12", "r1", "3412");
        reorder("pphg", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
            goto finish;
        }

        // T4b
        reorder("h1c", "r2", "21");
        mult("gghh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c-g", "r1", "r2", 1);
        mult("t1c-g", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("h2c-g12", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("t2c-g12", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c-g", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g2c-g1", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(1/23|6/45)");
    update("veff30", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_3h0p_ccsd_ijk_c_ab(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("h2c-g1", "r1", "2413");
        reorder("gphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        mult("t1c-g", "pghg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("h1c", "r1", "21");
        reorder("gphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("h2c-g1", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c-g", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g2c-g1", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(123|6/45)");
    update("veff30", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_3h0p_ccsd_k_ij_abc(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("h2c-g1", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("pghg", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c-g", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("h2c-g1", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c-g12", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|456)");
    update("veff30", 1.0, "r3");
    restore_stack_pos(pos);
}


void heff30_triples()
{

}


/**
 * Triples amplitude equations.
 * This code is invoked only if cc_model >= CCSDT
 */
void construct_triples_3h0p()
{
    timer_new_entry("30-T3", "3h0p -- Triples equations (T3)");
    timer_start("30-T3");

    copy("k3_0", "k3nw");

    // constant terms
    //diag_T3_3h0p_i_jk_c_ab(PT_INF);
    //diag_T3_3h0p_ijk_c_ab(PT_INF);
    //diag_T3_3h0p_k_ij_a_bc(PT_INF);
    //diag_T3_3h0p_k_ij_abc(PT_INF);
    //diag_T3_3h0p_c_ab(PT_INF);
    //diag_T3_3h0p_a_bc(PT_INF);
    //diag_T3_3h0p_abc(PT_INF);
    //diag_T3_3h0p_k_ij_c_ab(PT_INF);

    diag_T3_3h0p_ijk(PT_INF);
    diag_T3_3h0p_i_jk(PT_INF);
    diag_T3_3h0p_k_ij(PT_INF);

    // constant terms
    //diag_T3_0h3p_k_ij_a_bc();
    //diag_T3_0h3p_k_ij_abc();
    //diag_T3_0h3p_i_jk_c_ab();
    //diag_T3_0h3p_ijk_c_ab();

    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        /*diag_T3_0h3p_c_ab();
        diag_T3_0h3p_a_bc();
        diag_T3_0h3p_abc();*/
    }

    // constant terms
    //diag_T3_0h3p_k_ij_c_ab();
    //diag_T3_0h3p_i_jk();
    //diag_T3_0h3p_ijk();
    //diag_T3_0h3p_k_ij();

    timer_stop("30-T3");
}


void construct_folded_3h0p()
{
    dg_stack_pos_t pos = get_stack_pos();
    timer_new_entry("folded_K3", "Folded contribution to T{3h0p}_3");
    timer_start("folded_K3");

    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        // F1
        reorder("veff10", "r1", "21");
        mult("k3c", "r1", "r2", 1);
        update("k3nw", 1.0, "r2");
        restore_stack_pos(pos);
    }

    // F2
    reorder("veff20", "r1", "2341");
    mult("g2c", "r1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    tmplt("r3_ext", "hhhhhh", "000111", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r3_ext");
    diagram_stack_erase("r3");
    perm("r3_ext", "(3/12|4/56)");
    update("k3nw", 1.0, "r3_ext");
    restore_stack_pos(pos);

    // F3
    reorder("veff20", "r1", "3412");
    mult("r1", "h1c", "r2", 1);
    reorder("r2", "r3", "1243");
    mult("g2c", "r3", "r4", 1);
    reorder("r4", "r5", "126345");
    diagram_stack_erase("r4");
    perm("r5", "(3/12|4/56)");
    update("k3nw", -1.0, "r5");
    restore_stack_pos(pos);

    // F4
    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        reorder("veff20", "r1", "3412");
        mult("k3c", "r1", "r2", 2);
        perm("r2", "(4/56)");
        update("k3nw", -0.5, "r2");
        restore_stack_pos(pos);
    }

    // F5
    reorder("veff30", "r1", "456123");
    mult("r1", "h1c", "r2", 1);
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    tmplt("r3_ext", "hhhhhh", "000111", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r3_ext");
    diagram_stack_erase("r3");
    perm("r3_ext", "(3/12)");
    update("k3nw", 1.0, "r3_ext");
    restore_stack_pos(pos);

    // F6
    reorder("veff30", "r1", "456123");
    mult("r1", "g2c", "r2", 2);
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    tmplt("r3_ext", "hhhhhh", "000111", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r3_ext");
    diagram_stack_erase("r3");
    perm("r3_ext", "(1/23)");
    update("k3nw", -0.5, "r3_ext");
    restore_stack_pos(pos);

    // F7
    reorder("veff30", "r1", "345612");
    mult("h1c", "r1", "r2", 1);
    diagram_stack_erase("r1");
    mult("h1c", "r2", "r3", 1);
    diagram_stack_erase("r2");
    tmplt("r3_ext", "hhhhhh", "000111", "123456", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r3_ext");
    diagram_stack_erase("r3");
    perm("r3_ext", "(3/12)");
    update("k3nw", -1.0, "r3_ext");
    restore_stack_pos(pos);

    // F8
    reorder("veff30", "r1", "456123");
    mult("h1c", "r1", "r2", 1);
    diagram_stack_erase("r1");
    mult("h1c", "r2", "r3", 1);
    diagram_stack_erase("r2");
    mult("h1c", "r3", "r4", 1);
    update("k3nw", 1.0, "r4");
    restore_stack_pos(pos);

    // F9
    reorder("veff30", "r1", "456123");
    mult("g2c", "r1", "r2", 2);
    diagram_stack_erase("r1");
    mult("h1c", "r2", "r3", 1);
    diagram_stack_erase("r2");
    perm("r3", "(1/23)");
    update("k3nw", -0.5, "r3");
    restore_stack_pos(pos);

    if (cc_opts->cc_model >= CC_MODEL_CCSDT) {
        // F10
        reorder("veff30", "r1", "456123");
        mult("k3c", "r1", "r2", 3);
        update("k3nw", 0.125, "r2");
        restore_stack_pos(pos);
    }

    timer_stop("folded_K3");
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_T3_3h0p_i_jk_c_ab(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hhhg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME && pt_order == PT_INF) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3a
        reorder("h2c", "h2c_21", "2143");
        reorder("ph", "phr_", "21");
        reorder("h2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphg", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        if (cc_opts->cc_model <= CC_MODEL_CCSDT_2 && pt_order == PT_INF) {
            goto finish;
        }

        // T4b
        reorder("h1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("t2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5d
        reorder("g3c", "g3c_132", "132465");  // TODO: merge two transpositions
        reorder("g3c_132", "t5d_r1", "146235");
        reorder("pphh", "t5d_r2", "2341");
        mult("t5d_r1", "t5d_r2", "t5d_r3", 3);
        mult("t5d_r3", "h2c", "t5d_r4", 1);
        reorder("t5d_r4", "t5d_r5", "154236");
        update("r4", -0.5, "t5d_r5");
        restore_stack_pos(pos2);
    }

    perm("r4", "(1/23|6/45)");
    update("k3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_T3_3h0p_ijk_c_ab(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4c
        reorder("phhg", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7b
        reorder("h1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("h2c", "h2c_21", "2143");
        reorder("pphh", "r1", "3142");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(123|6/45)");
    update("k3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_T3_3h0p_k_ij_a_bc(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhp", "0110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phgg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }
        if (pt_order <= PT_2) {
            goto finish;
        }

        // T3e
        reorder("g2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4a
        reorder("ppgg", "ppggr", "3412");
        mult("ppggr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7d
        reorder("h1c", "r1", "21");
        mult("r1", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("g2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T10b
        reorder("h1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");

    if ((pt_order == PT_INF && cc_opts->cc_model > CC_MODEL_CCSDT_3) ||
        (pt_order != PT_INF && pt_order >= PT_4)) {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T5e
        reorder("h3c", "t5e_r1", "134562");
        reorder("pphh", "t5e_r2", "4123");
        reorder("g2c", "t5e_r3", "2341");
        mult("t5e_r1", "t5e_r2", "t5e_r4", 3);
        diagram_stack_erase("t5e_r1");
        mult("t5e_r4", "t5e_r3", "t5e_r5", 1);
        diagram_stack_erase("t5e_r3");
        reorder("t5e_r5", "t5e_r6", "124356");
        diagram_stack_erase("t5e_r5");
        update("r3", -0.5, "t5e_r6");
        restore_stack_pos(pos2);
    }

    perm("r3", "(3/12|4/56)");
    update("k3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_T3_3h0p_k_ij_abc(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhp", "0110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("phhg", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(3/12|456)");
    update("k3_0", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (6/45) aka (c/ab)
 * Diagrams included:
 * T2a T5c T6b T6e T9b
 */
void diag_T3_3h0p_c_ab(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hp", "10", "21", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2a
        reorder("pg", "pgr_", "21");
        update("i1", 1.0, "pgr_");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5c
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 3);
        update("i1", -0.5, "r2");
        restore_stack_pos(pos2);

        // T6b
        reorder("h1c", "r1", "21");
        mult("r1", "ph", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T6e
        reorder("ppgh", "r1", "3142");
        mult("r1", "t1c", "r2", 2);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9b
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        mult("r1", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("g3c", "i1", "r1", 1);
    perm("r1", "(6/45)");
    update("k3_0", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (4/56) aka (a/bc)
 * Diagrams included:
 * (a/bc): T2c T5g T9d
 */
void diag_T3_3h0p_a_bc(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhpp", "1100", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2c
        reorder("ppgg", "ppggr", "3412");
        update("i1", 0.5, "ppggr");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5g
        reorder("g2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        update("i1", 0.25, "r2");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9d
        reorder("h1c", "r1", "21");
        reorder("h1c", "r2", "21");
        mult("r1", "pphh", "r3", 1);
        mult("r2", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h3c", "i1", "r1", 2);
    perm("r1", "(4/56)");
    update("k3_0", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (456) aka (abc)
 * Diagrams included:
 * (abc):  T6g
 */
void diag_T3_3h0p_abc(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
        return;
    }

    tmplt("i1", "hhpp", "1100", "3412", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6g
        reorder("h1c", "r1", "21");
        mult("ppgh", "r1", "r2", 1);
        reorder("r2", "r3", "3412");
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);
    }

    mult("h3c", "i1", "r1", 2);
    perm("r1", "(456)");
    update("k3_0", 1.0, "r1");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|6/45) aka (k/ij|c/ab)
 * Diagrams included:
 * T2e T5a T6c T6d T9e
 */
void diag_T3_3h0p_k_ij_c_ab(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhp", "0100", "2431", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2e
        reorder("phhg", "r0", "2431");
        update("i1", 1.0, "r0");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5a
        reorder("pphh", "r2", "3142");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r3", "2413");
        mult("r3", "r2", "r4", 2);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T6c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T6d
        reorder("phhh", "r1", "2314");
        reorder("h1c", "r2", "21");
        mult("r2", "r1", "r3", 1);
        reorder("r3", "r4", "2134");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9e
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3124");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g3c", "r1", "124536");
    mult("r1", "i1", "r3", 2);
    diagram_stack_erase("r1");
    reorder("r3", "r4", "125346");
    diagram_stack_erase("r3");
    perm("r4", "(3/12|6/45)");
    update("k3_0", 1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123) aka (ijk)
 * Diagrams included:
 * (ijk):  T6h
 */
void diag_T3_3h0p_ijk(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order < PT_4) {
        return;
    }

    tmplt("i1", "hhhh", "0000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T6h
        reorder("hphh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    reorder("k3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(123)");
    update("k3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23) aka (i/jk)
 * Diagrams included:
 * (i/jk): T2d T5f T9c
 */
void diag_T3_3h0p_i_jk(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hhhh", "0000", "1234", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2d
        update("i1", 0.5, "hhhh");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5f
        reorder("pphh", "r2", "3412");
        mult("t2c", "r2", "r3", 2);
        update("i1", 0.25, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9c
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 0.5, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("k3c", "r1", "456123");
    mult("r1", "i1", "r2", 2);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(1/23)");
    update("k3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12) aka (k/ij)
 * Diagrams included:
 * T2b T5b T6a T6f T9a
 */
void diag_T3_3h0p_k_ij(int pt_order)
{
    dg_stack_pos_t pos, pos2;
    pos = get_stack_pos();

    if (cc_opts->cc_model < CC_MODEL_CCSDT && pt_order == PT_INF) {
        return;
    }
    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hh", "00", "12", NOT_PERM_UNIQUE);

    {
        pos2 = get_stack_pos();

        // T2b
        update("i1", -1.0, "hh");

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T5b
        reorder("pphh", "r2", "3412");
        mult("t2c", "r2", "r3", 3);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T6a
        reorder("ph", "r2", "21");
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T6f
        reorder("hphh", "r2", "1342");
        mult("r2", "t1c", "r3", 2);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        if (pt_order <= PT_4) {
            goto finish;
        }

        // T9a
        reorder("pphh", "r2", "3142");
        mult("r2", "t1c", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("k3c", "r1", "456123");
    mult("r1", "i1", "r2", 1);
    diagram_stack_erase("r1");
    reorder("r2", "r3", "456123");
    diagram_stack_erase("r2");
    perm("r3", "(3/12)");
    update("k3nw", 1.0, "r3");
    restore_stack_pos(pos);
}


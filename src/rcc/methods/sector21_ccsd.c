/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2023 The EXP-T developers.
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
 * Fock space sector 2h1p. CCSD model
 *
 * The theory for the 2h1p sector is analogous to that for the 1h2p sector.
 * For the latter, see:
 * S. R. Hughes, U. Kaldor. Fock-space coupled-cluster method: The (1,2) sector.
 * Phys. Rev. A 47, 4705 (1993). doi: 10.1103/PhysRevA.47.4705
 *
 * FS-CCSD calculations in the 2h1p sector were also reported for the bromine atom,
 * S. R. Hughes, U. Kaldor. The Fock-space coupled cluster method extended
 * to higher sectors. Condensed Matter Theories. PP. 385-394 (1993).
 *
 */

#include "methods.h"

#include <math.h>
#include <stdlib.h>

#include "ccutils.h"
#include "diis.h"
#include "engine.h"
#include "heff.h"
#include "datamodel.h"
#include "options.h"
#include "sort.h"
#include "spinors.h"
#include "utils.h"

void heff21_ccsd_2h1p_dependent();

void heff21_ccsd_const();

void sort_integrals_2h1p();

int solve_amplitude_equations_2h1p();

void linsys_solver();

int iterative_solver();

void add_eps_to_diagonal(char *name);

void calc_2h1p_intermediates();

void init_amplitudes_2h1p();

void construct_doubles_2h1p();

void construct_folded_2h1p();


void ampl_equations_2h1p_D1();

void ampl_equations_2h1p_D2a();

void ampl_equations_2h1p_D2b_1();

void ampl_equations_2h1p_D2b_2();

void ampl_equations_2h1p_D2c();

void ampl_equations_2h1p_D2d();

void ampl_equations_2h1p_D2e_1();

void ampl_equations_2h1p_D2e_2();

void ampl_equations_2h1p_D3a();

void ampl_equations_2h1p_D3b();

void ampl_equations_2h1p_D3c_1();

void ampl_equations_2h1p_D3c_2();

void ampl_equations_2h1p_D3d();

void ampl_equations_2h1p_D4a_1();

void ampl_equations_2h1p_D4a_2();

void ampl_equations_2h1p_D4b();

void ampl_equations_2h1p_D5a_1();

void ampl_equations_2h1p_D5a_2();

void ampl_equations_2h1p_D5b();

void ampl_equations_2h1p_D5c_1();

void ampl_equations_2h1p_D5c_2();

void ampl_equations_2h1p_D5d_1();

void ampl_equations_2h1p_D5d_2();

void ampl_equations_2h1p_D5e();

void ampl_equations_2h1p_D5f_1();

void ampl_equations_2h1p_D5f_2();

void ampl_equations_2h1p_D5g();

void ampl_equations_2h1p_D5h_1();

void ampl_equations_2h1p_D5h_2();

void ampl_equations_2h1p_D6a();

void ampl_equations_2h1p_D6b();

void ampl_equations_2h1p_D6c_1();

void ampl_equations_2h1p_D6c_2();

void ampl_equations_2h1p_D7a();

void ampl_equations_2h1p_D7b();

void ampl_equations_2h1p_D7c_1();

void ampl_equations_2h1p_D7c_2();

void ampl_equations_2h1p_D7d_1();

void ampl_equations_2h1p_D7d_2();

void ampl_equations_2h1p_D7e();

void ampl_equations_2h1p_D8a();

void ampl_equations_2h1p_D8b_1();

void ampl_equations_2h1p_D8b_2();

void ampl_equations_2h1p_D9();

void ampl_equations_2h1p_folded_F1();

void ampl_equations_2h1p_folded_F2();

void ampl_equations_2h1p_folded_F3();

void ampl_equations_2h1p_folded_F4();

void ampl_equations_2h1p_folded_F5();

void ampl_equations_2h1p_folded_F6();

void ampl_equations_2h1p_folded_F7();

void ampl_equations_2h1p_folded_F8();

void ampl_equations_2h1p_folded_F9();

void ampl_equations_2h1p_folded_F10();


void diag_heff_2h1p_i_jk_c_ab__inv_h12_p3_const();

void diag_heff_2h1p_i_jk_c_ab__inv_h23_p1_const();

void diag_heff_2h1p_i_jk_c_ab__inv_h13_p2_const();

void diag_heff_2h1p_i_jk_c_ab__inv_h12_p1_w2c();


void diag_heff_2h1p_ijk_c_ab__inv_h23_p1_const();

void diag_heff_2h1p_ijk_c_ab__inv_h13_p2_const();

void diag_heff_2h1p_ijk_c_ab__inv_h12_p2_const();

void diag_heff_2h1p_ijk_c_ab__inv_h13_p3_const();

void diag_heff_2h1p_ijk_c_ab__inv_h12_p3_const();

void diag_heff_2h1p_ijk_c_ab__inv_h12_p1_w2c();


void diag_heff_2h1p_k_ij_a_bc__inv_h12_p3_const();

void diag_heff_2h1p_k_ij_a_bc__inv_h23_p1_const();

void diag_heff_2h1p_k_ij_a_bc__inv_h13_p2_const();

void diag_heff_2h1p_k_ij_a_bc__inv_h23_p3_const();


void diag_heff_2h1p_k_ij_abc__inv_h12_p2_const();

void diag_heff_2h1p_k_ij_abc__inv_h23_p1_const();

void diag_heff_2h1p_k_ij_abc__inv_h13_p2_const();

void diag_heff_2h1p_k_ij_abc__inv_h12_p3_const();

void diag_heff_2h1p_k_ij_abc__inv_h23_p3_const();

void diag_heff_2h1p_k_ij_abc__inv_h13_p3_const();


void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence);


int sector_2h1p(cc_options_t *opts)
{
    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 2h1p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    if (opts->hughes_kaldor_2h1p) {
        printf("\n");
        printf(" The Hughes-Kaldor version of the FS-CCSD(2h1p) method will be used\n");
        printf(" (folded terms with Heff{2h1p} will be excluded).\n");
        printf("\n");
    }

    opts->curr_sector_h = 2;
    opts->curr_sector_p = 1;

    // setup shifts for the IH-like technique by A. V. Zaitsevskii
    if (intham_imms_in_sector(2, 1)) {
        intham_imms_setup(2, 1);
    }

    sort_integrals_2h1p();
    init_amplitudes_2h1p();
    predict_intruders_for_diagram("w2nw", 10);

    /*
     * initial approximation to Veff{2h1p}.
     * includes diagrams independent on the T{2h1p} cluster amplitudes.
     */
    heff21_ccsd_const();

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            2, 1, NULL, NULL, "w2c", "w2nw", NULL, NULL, NULL,
            NULL, construct_doubles_2h1p, NULL, construct_folded_2h1p
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze cluster amplitudes and write them to disk
     */
    print_cluster_operator_analysis(2, 1, NULL, "w2c", NULL);
    save_cluster_amplitudes(2, 1, NULL, "w2c", NULL, NULL);
    if (cc_opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(2, 1, NULL, "w2c", NULL);
    }

    /*
     * effective Hamiltonian construction and analysis
     */
    copy("veff21_const", "veff21");
    heff21_ccsd_2h1p_dependent();

    heff_analysis(2, 1, "veff01", "veff10", "veff11", "veff20", "veff21");

    /*
     * calculate quasi-natural orbitals and model-space estimates of properties
     */
    calculate_model_space_properties(2, 1);
    calculate_model_space_natural_spinors(2, 1);

    return EXIT_SUCCESS;
}


void sort_integrals_2h1p()
{
    /*
     * two-electron integrals
     */
    request_sorting("vphh", "pphh", "1000", "1234");
    request_sorting("pppg", "ppph", "0001", "1234");
    request_sorting("pvgg", "pphh", "0111", "1234");
    request_sorting("ppgh", "pphh", "0010", "1234");
    request_sorting("ppph", "ppph", "0000", "1234");
    request_sorting("hvgg", "hphh", "0111", "1234");
    request_sorting("hvhg", "hphh", "0101", "1234");
    perform_sorting();

    /*
     * prepare cluster amplitudes from previous sectors
     */
}


void init_amplitudes_2h1p()
{
    tmplt("w2nw", "hphh", "0111", "1234", NOT_PERM_UNIQUE);
    copy("w2nw", "w2c");

    if (cc_opts->reuse_amplitudes[2][1]) {
        printf(" Trying to read amplitudes and effective interaction (sector 2h1p) from disk ...\n");
        // Doubles
        if (diagram_read("w2c.dg") != NULL) {
            printf(" T{2h1p}_1 amplitudes successfully read from disk\n");
        }
        else {
            printf(" T{2h1p}_1 amplitudes will be calculated\n");
        }
        if (diagram_read("veff21.dg") != NULL) {
            printf(" Heff{2h1p} diagram successfully read from disk\n");
        }
        else {
            printf(" Heff{2h1p} diagram will be calculated\n");
        }
    }
    printf("\n");
}


void construct_doubles_2h1p()
{
    /*
     * (V\Omega)_conn
     */
    ampl_equations_2h1p_D1();
    ampl_equations_2h1p_D2a();
    ampl_equations_2h1p_D2b_1();
    ampl_equations_2h1p_D2b_2();
    ampl_equations_2h1p_D2c();
    ampl_equations_2h1p_D2d();
    ampl_equations_2h1p_D2e_1();
    ampl_equations_2h1p_D2e_2();
    ampl_equations_2h1p_D3a();
    ampl_equations_2h1p_D3b();
    ampl_equations_2h1p_D3c_1();
    ampl_equations_2h1p_D3c_2();
    ampl_equations_2h1p_D3d();
    ampl_equations_2h1p_D4a_1();
    ampl_equations_2h1p_D4a_2();
    ampl_equations_2h1p_D4b();
    ampl_equations_2h1p_D5a_1();
    ampl_equations_2h1p_D5a_2();
    ampl_equations_2h1p_D5b();
    ampl_equations_2h1p_D5c_1();
    ampl_equations_2h1p_D5c_2();
    ampl_equations_2h1p_D5d_1();
    ampl_equations_2h1p_D5d_2();
    ampl_equations_2h1p_D5e();
    ampl_equations_2h1p_D5f_1();
    ampl_equations_2h1p_D5f_2();
    ampl_equations_2h1p_D5g();
    ampl_equations_2h1p_D5h_1();
    ampl_equations_2h1p_D5h_2();
    ampl_equations_2h1p_D6a();
    ampl_equations_2h1p_D6b();
    ampl_equations_2h1p_D6c_1();
    ampl_equations_2h1p_D6c_2();
    ampl_equations_2h1p_D7a();
    ampl_equations_2h1p_D7b();
    ampl_equations_2h1p_D7c_1();
    ampl_equations_2h1p_D7c_2();
    ampl_equations_2h1p_D7d_1();
    ampl_equations_2h1p_D7d_2();
    ampl_equations_2h1p_D7e();
    ampl_equations_2h1p_D8a();
    ampl_equations_2h1p_D8b_1();
    ampl_equations_2h1p_D8b_2();
    ampl_equations_2h1p_D9();
}


void construct_folded_2h1p()
{
    timer_new_entry("21-CCSD-Folded", "2h1p -- Folded diagrams");
    timer_start("21-CCSD-Folded");

    /*
     * effective interaction will be used in folded diagrams F9 and F10
     */
    if (cc_opts->hughes_kaldor_2h1p == 0) {
        copy("veff21_const", "veff21");
        heff21_ccsd_2h1p_dependent();
    }

    ampl_equations_2h1p_folded_F1();
    ampl_equations_2h1p_folded_F2();
    ampl_equations_2h1p_folded_F3();
    ampl_equations_2h1p_folded_F4();
    ampl_equations_2h1p_folded_F5();
    ampl_equations_2h1p_folded_F6();
    ampl_equations_2h1p_folded_F7();
    ampl_equations_2h1p_folded_F8();

    /*
     * folded terms involving the Heff{2h1p} operator
     */
    if (cc_opts->hughes_kaldor_2h1p == 0) {
        ampl_equations_2h1p_folded_F9();
        ampl_equations_2h1p_folded_F10();
    }

    timer_stop("21-CCSD-Folded");
}


/**
 * The part of Heff{2h1p} which is constant in the 2h1p sector
 * (no diagrams depending on the T{2h1p} amplitudes).
 * This part is calculated only once at start.
 * Result is stored in the "veff21_const" diagram.
 */
void heff21_ccsd_const()
{
    timer_new_entry("21-CCSD-Heff", "2h1p -- CCSD part of eff Hamiltonian");
    timer_start("21-CCSD-Heff");

    printf("\n Construction of the CCSD eff Hamiltonian in the 2h1p sector (constant part)\n");

    tmplt("veff21_const", "hhphhp", "111111", "123456", NOT_PERM_UNIQUE);

    printf(" [1/4] P(i/jk|c/ab)\n");

    diag_heff_2h1p_i_jk_c_ab__inv_h12_p3_const();
    diag_heff_2h1p_i_jk_c_ab__inv_h23_p1_const();
    diag_heff_2h1p_i_jk_c_ab__inv_h13_p2_const();
    //diag_heff_2h1p_i_jk_c_ab__inv_h12_p1_w2c(); // depends on 2h1p

    printf(" [2/4] P(ijk|c/ab)\n");
    diag_heff_2h1p_ijk_c_ab__inv_h23_p1_const();
    diag_heff_2h1p_ijk_c_ab__inv_h13_p2_const();
    diag_heff_2h1p_ijk_c_ab__inv_h12_p2_const();
    diag_heff_2h1p_ijk_c_ab__inv_h13_p3_const();
    diag_heff_2h1p_ijk_c_ab__inv_h12_p3_const();
    //diag_heff_2h1p_ijk_c_ab__inv_h12_p1_w2c(); // depends on 2h1p

    printf(" [3/4] P(k/ij|a/bc)\n");
    diag_heff_2h1p_k_ij_a_bc__inv_h12_p3_const();
    diag_heff_2h1p_k_ij_a_bc__inv_h23_p1_const();
    diag_heff_2h1p_k_ij_a_bc__inv_h13_p2_const();
    diag_heff_2h1p_k_ij_a_bc__inv_h23_p3_const();

    printf(" [4/4] P(k/ij|abc)\n");
    diag_heff_2h1p_k_ij_abc__inv_h12_p2_const();
    diag_heff_2h1p_k_ij_abc__inv_h23_p1_const();
    diag_heff_2h1p_k_ij_abc__inv_h13_p2_const();
    diag_heff_2h1p_k_ij_abc__inv_h12_p3_const();
    diag_heff_2h1p_k_ij_abc__inv_h23_p3_const();
    diag_heff_2h1p_k_ij_abc__inv_h13_p3_const();

    printf("\n");
    timer_stop("21-CCSD-Heff");
}


/**
 * The part of Heff{2h1p} which is dependent on the T{2h1p} amplitudes.
 * Should be recalculated at each iteration.
 * Result is added to the "veff21" diagram.
 */
void heff21_ccsd_2h1p_dependent()
{
    timer_new_entry("21-CCSD-Heff", "2h1p -- CCSD part of eff Hamiltonian");
    timer_start("21-CCSD-Heff");

    diag_heff_2h1p_i_jk_c_ab__inv_h12_p1_w2c();
    diag_heff_2h1p_ijk_c_ab__inv_h12_p1_w2c();

    timer_stop("21-CCSD-Heff");
}


/*
 * Amplitude equations for doubles: direct terms, (V\Omega)_conn
 */


void ampl_equations_2h1p_D1()
{
    tmplt("w2nw", "hphh", "0111", "1234", NOT_PERM_UNIQUE);
    update("w2nw", 1.0, "hvgg");
}


void ampl_equations_2h1p_D2a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pg", "r0", "21");
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "r0", "r1", 1);
    perm("r1", "(34)");
    update("w2nw", 1.0, "r1");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D2b_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("g2c", "r0", "3412");
    mult("r0", "vh", "r1", 1);
    reorder("r1", "r2", "3412");
    update("w2nw", -1.0, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D2b_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("w2c", "w2c_21", "2143");
    reorder("w2c_21", "r0", "3412");
    mult("r0", "hh", "r1", 1);
    reorder("r1", "r2", "3412");
    reorder("r2", "r3", "2143"); // interchange electrons: 1 <-> 2
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D2c()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ppgg", "ppggr", "3412");
    reorder("s2c", "s2c_21", "2143");
    mult("ppggr", "s2c_21", "r1", 2);
    reorder("r1", "r2", "3412");
    update("w2nw", 0.5, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D2d()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("g2c", "r1", "3412");
    mult("hvhh", "r1", "r2", 2);
    update("w2nw", 0.5, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D2e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "r1", "1324");
    mult("r1", "pvhg", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(34)");
    update("w2nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D2e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("phhg", "r1", "2431");
    reorder("e2c", "r2", "1423");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(34)");
    reorder("r4", "r5", "2143"); // interchange electrons: 1 <-> 2
    update("w2nw", -1.0, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D3a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("g2c", "r0", "3412");
    reorder("s2c", "s2c_21", "2143");
    reorder("pphh", "v1", "3412");
    mult("s2c_21", "v1", "r1", 2);
    mult("r1", "r0", "r2", 2);
    update("w2nw", 0.25, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D3b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "r1", "1324");
    reorder("pphh", "r2", "4231");
    reorder("e2c", "r3", "1423");
    mult("r1", "r2", "r4", 2);
    mult("r4", "r3", "r5", 2);
    reorder("r5", "r6", "1324");
    perm("r6", "(34)");
    update("w2nw", -1.0, "r6");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D3c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v1", "3412");
    mult("t2c", "v1", "r1", 3);
    reorder("w2c", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    update("w2nw", -0.5, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D3c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "3412");
    reorder("g2c", "r2", "1342");
    mult("s2c", "r1", "r3", 3);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1423");
    update("w2nw", -0.5, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D3d()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "r1", "3412");
    reorder("e2c", "e2c_21", "2143");
    mult("r1", "pphh", "r2", 3);
    mult("e2c_21", "r2", "r3", 1);
    perm("r3", "(34)");
    update("w2nw", -0.5, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D4a_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvgg", "r1", "2341");
    mult("t1c", "r1", "r2", 1);
    update("w2nw", 1.0, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D4a_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("phgg", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    reorder("r2", "r3", "2143"); // interchange electrons: 1 <-> 2
    update("w2nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D4b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("hvhg", "r1", "1243");
    reorder("h1c", "r0", "21");
    mult("r1", "r0", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5a_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "r1", "21");
    reorder("w2c", "r2", "2341");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    update("w2nw", -1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5a_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "r1", "21");
    reorder("g2c", "r0", "3412");
    mult("s1c", "r1", "r2", 1);
    mult("r0", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("w2nw", -1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h1c", "r1", "21");
    reorder("e2c", "e2c_21", "2143");
    mult("r1", "ph", "r2", 1);
    mult("e2c_21", "r2", "r3", 1);
    perm("r3", "(34)");
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ppgh", "r1", "3421");
    reorder("e2c", "r2", "1423");
    mult("t1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(34)");
    update("w2nw", -1.0, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "r1", "1324");
    reorder("pphg", "r2", "4312");
    mult("s1c", "r2", "r3", 1);
    mult("r1", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(34)");
    update("w2nw", 1.0, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5d_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("hphh", "r1", "1243");
    reorder("h1c", "r2", "21");
    reorder("e2c", "r3", "1423");
    mult("r1", "r2", "r4", 1);
    reorder("r4", "r5", "1432");
    mult("r5", "r3", "r6", 2);
    reorder("r6", "r7", "1324");
    perm("r7", "(34)");
    update("w2nw", 1.0, "r7");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5d_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "r1", "1324");
    reorder("pvhh", "r2", "2314");
    reorder("h1c", "r3", "21");
    mult("r3", "r2", "r4", 1);
    mult("r1", "r4", "r5", 2);
    reorder("r5", "r6", "1423");
    perm("r6", "(34)");
    update("w2nw", -1.0, "r6");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5e()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4312");
    reorder("h1c", "r0", "21");
    reorder("s2c", "s2c_21", "2143");
    mult("s2c_21", "r1", "i1", 2);
    mult("i1", "r0", "r2", 1);
    reorder("r2", "r3", "1243");
    perm("r3", "(34)");
    update("w2nw", -0.5, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5f_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2341");
    reorder("g2c", "r0", "3412");
    mult("t1c", "r1", "i1", 1);
    mult("i1", "r0", "r2", 2);
    update("w2nw", 0.5, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5f_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("phhh", "r1", "2341");
    reorder("g2c", "r0", "3412");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "r0", "r2", 2);
    reorder("r2", "r3", "2143"); // interchange electrons: 1 <-> 2
    update("w2nw", 0.5, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5g()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ppgh", "r1", "3142");
    reorder("e2c", "e2c_21", "2143");
    mult("r1", "t1c", "r2", 2);
    mult("e2c_21", "r2", "r3", 1);
    perm("r3", "(34)");
    update("w2nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5h_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("phhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("w2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D5h_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("g2c", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    reorder("r3", "r4", "2143");
    update("w2nw", -1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D6a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ppgg", "ppggr", "3412");
    mult("ppggr", "s1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("t1c", "r1", "r2", 1);
    update("w2nw", 1.0, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D6b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h1c", "r0", "21");
    mult("r0", "hvhh", "i1", 1);
    mult("r0", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("w2nw", 1.0, "r1");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D6c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h1c", "r2", "21");
    mult("t1c", "pvhg", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("w2nw", -1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D6c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("phhg", "r1", "2431");
    reorder("h1c", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "2134"); // [1243] + interchange electrons 1 <-> 2
    perm("r5", "(34)");
    update("w2nw", -1.0, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "vr1", "3412");
    reorder("g2c", "r0", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("r2", "r0", "r3", 2);
    update("w2nw", 0.5, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "vr1", "3412");
    reorder("h1c", "r0", "21");
    mult("s2c", "vr1", "r1", 2);
    mult("r0", "r1", "r2", 1);
    mult("r0", "r2", "r3", 1);
    reorder("r3", "r4", "4321"); // [3412] + interchange electrons 1 <-> 2
    update("w2nw", 0.5, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "1342");
    reorder("e2c", "r2", "1423");
    reorder("h1c", "r3", "21");
    mult("r2", "r1", "r4", 2);
    mult("r3", "r4", "r5", 1);
    mult("t1c", "r5", "r6", 1);
    reorder("r6", "r7", "1324");
    perm("r7", "(34)");
    update("w2nw", 1.0, "r7");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "3142");
    reorder("h2c", "h2c_21", "2143");
    reorder("h2c_21", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "2134"); // [1243] + interchange electrons 1 <-> 2
    perm("r7", "(34)");
    update("w2nw", -1.0, "r7");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7d_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r0", "2431");
    mult("r0", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("t1c", "r2", "r3", 1);
    reorder("w2c", "w2r", "2341");
    mult("r3", "w2r", "r4", 1);
    update("w2nw", -1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7d_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r0", "2431");
    mult("r0", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    reorder("g2c", "g2r", "2341");
    mult("r3", "g2r", "r4", 1);
    reorder("r4", "r5", "2143"); // interchange electrons 1 <-> 2
    update("w2nw", -1.0, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D7e()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r0", "2431");
    reorder("h1c", "h1cr", "21");
    mult("r0", "t1c", "r1", 2);
    mult("h1cr", "r1", "r2", 1);
    reorder("e2c", "e2r", "1243");
    mult("e2r", "r2", "r3", 1);
    reorder("r3", "r4", "2134"); // [1243] + interchange electrons 1 <-> 2
    perm("r4", "(34)");
    update("w2nw", -1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D8a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4312");
    reorder("h1c", "h1r", "21");
    mult("s1c", "r1", "r2", 1);
    mult("t1c", "r2", "r3", 1);
    mult("r3", "h1r", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(34)");
    update("w2nw", -1.0, "r5");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D8b_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2341");
    reorder("h1c", "h1r", "21");
    mult("t1c", "r1", "i1", 1);
    mult("h1r", "i1", "i2", 1);
    mult("h1r", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    update("w2nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D8b_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("phhh", "r1", "2341");
    reorder("h1c", "h1r", "21");
    mult("s1c", "r1", "i1", 1);
    mult("h1r", "i1", "i2", 1);
    mult("h1r", "i2", "r2", 1);
    reorder("r2", "r3", "4321"); // [3412] + interchange electrons 1 <-> 2
    update("w2nw", 1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_D9()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h1c", "h1r", "21");
    reorder("pphh", "vr1", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("t1c", "r1", "r2", 1);
    mult("h1r", "r2", "r3", 1);
    mult("h1r", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("w2nw", 1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Amplitude equations for doubles: folded terms
 */


void ampl_equations_2h1p_folded_F1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("w2c", "r1", "1243");
    reorder("veff10", "r2", "21");
    mult("r1", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(34)");
    update("w2nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("w2c", "r1", "3412");
    mult("r1", "veff01", "r2", 1);
    reorder("r2", "r3", "3412");
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F3()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff20", "r1", "3412");
    mult("w2c", "r1", "r2", 2);
    update("w2nw", -0.5, "r2");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F4()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff20", "r1", "3412");
    mult("r1", "e1c", "r2", 1);
    reorder("r2", "r3", "3412");
    tmplt("r4", "hphh", "0111", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r4");
    update("w2nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F5()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff20", "r1", "3412");
    mult("e1c", "r1", "r2", 1);
    mult("h1c", "r2", "r3", 1);
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F6()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("w2c", "r1", "1432");
    reorder("veff11", "veff11_21", "2143");
    reorder("veff11_21", "r2", "2314");
    mult("r1", "r2", "r3", 2);
    reorder("r3", "r4", "1342");
    perm("r4", "(34)");
    update("w2nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F7()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff11", "veff11_21", "2143");
    reorder("e1c", "r1", "21");
    mult("veff11_21", "r1", "r2", 1);
    set_order("r2", "1234");
    tmplt("r3", "hphh", "0111", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r2", "r3");
    perm("r3", "(34)");
    update("w2nw", -1.0, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F8()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff11", "veff11_21", "2143");
    reorder("veff11_21", "r1", "2341");
    reorder("e1c", "r2", "21");
    mult("h1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    perm("r4", "(34)");
    update("w2nw", 1.0, "r4");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F9()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff21", "r1", "345612");
    mult("r1", "w2c", "r2", 3);
    reorder("r2", "r3", "4123");
    update("w2nw", -0.5, "r3");
    restore_stack_pos(pos);
}


void ampl_equations_2h1p_folded_F10()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff21", "r1", "134562");
    mult("r1", "e1c", "r2", 2);
    tmplt("r3", "hphh", "0111", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r2", "r3");
    update("w2nw", 1.0, "r3");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_heff_2h1p_i_jk_c_ab__inv_h12_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hvhp", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("s2c", "s2c_21", "2143");
        reorder("pphp", "r2", "4312");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c", "r2", "21");
        mult("hvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphp", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        reorder("s2c", "s2c_21", "2143");
        mult("s2c_21", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
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
    perm("r4", "(12)");
    closed("r4", "r5");
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_i_jk_c_ab__inv_h23_p1_const()
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

        // T3a
        reorder("ph", "phr_", "21");
        reorder("h2c", "h2c_21", "2143");
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

        // T4b
        reorder("h1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    perm("r4", "(56)");
    closed("r4", "r5");
    reorder("r5", "r6", "321654");
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_i_jk_c_ab__inv_h13_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vhhg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("e2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphg", "r2", "4312");
        mult("s2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("h1c", "r2", "21");
        mult("vhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("e2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("s2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    closed("r4", "r5");
    perm("r5", "(13|46)");
    reorder("r5", "r6", "132465"); // interchange electrons 2 <-> 3
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


// contribution to Heff dependent on w2c (ampl-s in 2h1p)
void diag_heff_2h1p_i_jk_c_ab__inv_h12_p1_w2c()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhph", "0000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("hhhp", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("t2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("t2c", "r1", "3412");
        reorder("pphp", "r2", "4312");
        mult("t2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c", "r2", "21");
        mult("hhhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphp", "r1", "4312");
        mult("t1c", "r1", "r2", 1);
        mult("t1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("t2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c", "r3", "21");
        mult("t2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("w2c", "r1", "4321"); // [3412] + interchange electrons 1 <-> 2
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    diagram_stack_erase("r3");
    closed("r4", "r5");
    reorder("r5", "r6", "321456"); // and the sign must be changed, +1 -> -1
    update("veff21", -1.0, "r6");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_heff_2h1p_ijk_c_ab__inv_h23_p1_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhh", "0010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("phhg", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

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
    reorder("e2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(23|56)");
    reorder("r5", "r6", "321654"); // interchange electrons 1 <-> 3
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_ijk_c_ab__inv_h13_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1010", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("phhg", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(13|46)");
    reorder("r5", "r6", "132465"); // interchange electrons 2 <-> 3
    update("veff21_const", 1.0, "r6");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_ijk_c_ab__inv_h12_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phph", "1000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("s1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("t2c", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(13)");
    reorder("r5", "r6", "132456"); // + sign will be changed!
    update("veff21_const", -1.0, "r6");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_ijk_c_ab__inv_h13_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphh", "0110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("t1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(12|46)");
    reorder("r5", "r6", "123465"); // the sign will be changed!
    update("veff21_const", -1.0, "r6");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_ijk_c_ab__inv_h12_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "0100", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pvhp", "r0", "2431");
        mult("t1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("g2c", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(12)");
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


// TODO: contribution to Heff dependent on w2c (ampl-s in 2h1p)
void diag_heff_2h1p_ijk_c_ab__inv_h12_p1_w2c()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhph", "0000", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("t2c", "r1", "2413");
        reorder("hphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("t1c", "phhp", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c", "r1", "21");
        reorder("hphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("t2c", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("w2c", "r1", "4321"); // [3412] + interchange electrons 1 <-> 2
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    closed("r4", "r5");
    perm("r5", "(23)");
    reorder("r5", "r6", "321456"); // the sign will be changed, +1 -> -1
    update("veff21", -1.0, "r6");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_heff_2h1p_k_ij_a_bc__inv_h12_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvgp", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("h2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        mult("ppgp", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c", "r0", "21");
        reorder("h1c", "r1", "21");
        mult("r0", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c", "r0", "21");
        reorder("h1c", "r1", "21");
        mult("r0", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    perm("r4", "(45)");
    update("veff21_const", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_a_bc__inv_h23_p1_const()
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

        // T3e
        reorder("g2c", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppgg", "ppggr", "3412");
        mult("ppggr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

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

        // T10b
        reorder("h1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    perm("r4", "(23)");
    reorder("r4", "r5", "321654"); // interchange electrons: 1 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_a_bc__inv_h13_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("phpg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "phhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("pppg", "pppgr", "3412");
        mult("pppgr", "t1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "phhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "r0", "21");
        reorder("t1c", "r1", "21");
        mult("r0", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons: 2 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_a_bc__inv_h23_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvgg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("g2c", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppgg", "ppggr", "3412");
        mult("ppggr", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("g2c", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    diagram_stack_erase("r2");
    closed("r3", "r4");
    reorder("r4", "r5", "123654"); // the sign of the prefactor will be changed, +1 -> -1
    update("veff21_const", -1.0, "r5");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_heff_2h1p_k_ij_abc__inv_h12_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhpp", "0100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("t2c", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("phhp", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("t2c", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(13|45)");
    reorder("r4", "r5", "132456"); // the sign will be changed
    update("veff21_const", -1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_abc__inv_h23_p1_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhhp", "0110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("phhg", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

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
    mult("s2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(23|56)");
    reorder("r4", "r5", "321654"); // interchange electrons: 1 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_abc__inv_h13_p2_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "0010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("h2c", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("phhg", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "t1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
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
    reorder("e2c", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons: 2 <-> 3
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_abc__inv_h12_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1100", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhp", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("s2c", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
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
    closed("r3", "r4");
    perm("r4", "(45)");
    update("veff21_const", 1.0, "r4");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_abc__inv_h23_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5");  // the sign was changed, +1 -> -1
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7");  // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    finish:
    mult("t2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(56)");
    reorder("r4", "r5", "123654"); // the sign will be changed, +1 -> -1
    update("veff21_const", -1.0, "r5");
    restore_stack_pos(pos);
}


void diag_heff_2h1p_k_ij_abc__inv_h13_p3_const()
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphp", "1010", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppph", "r1", "1342");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5"); // the sign was changed, +1 -> -1
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("t1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppph", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c", "r1", "21");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7"); // the sign was changed, -1 -> +1
        restore_stack_pos(pos2);
    }

    finish:
    mult("h2c", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    closed("r3", "r4");
    perm("r4", "(46)");
    reorder("r4", "r5", "123465"); // the sign will be changed, +1 -> -1
    update("veff21_const", 1.0, "r5");
    restore_stack_pos(pos);
}

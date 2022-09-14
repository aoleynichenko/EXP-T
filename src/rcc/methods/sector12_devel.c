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

/*
 * Fock space sector 1h2p.
 *
 * For more details, see:
 * S. R. Hughes, U. Kaldor. Fock-space coupled-cluster method: The (1,2) sector.
 * Phys. Rev. A 47, 4705 (1993). doi: 10.1103/PhysRevA.47.4705
 *
 * 2021-2022 Alexander Oleynichenko
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

void heff12_ccsd_1h2p_dependent();

void heff12_ccsd_const();

void sort_integrals_1h2p();

int solve_amplitude_equations_1h2p();

void linsys_solver();

int iterative_solver();

void add_eps_to_diagonal(char *name);

void calc_1h2p_intermediates();

void init_amplitudes_1h2p();

void construct_doubles_1h2p();

void construct_folded_1h2p();

void ampl_equations_1h2p_D1();

void ampl_equations_1h2p_D2a_1();

void ampl_equations_1h2p_D2a_2();

void ampl_equations_1h2p_D2b();

void ampl_equations_1h2p_D2c();

void ampl_equations_1h2p_D2d();

void ampl_equations_1h2p_D2e_1();

void ampl_equations_1h2p_D2e_2();

void ampl_equations_1h2p_D3a();

void ampl_equations_1h2p_D3b();

void ampl_equations_1h2p_D3c();

void ampl_equations_1h2p_D3d_1();

void ampl_equations_1h2p_D3d_2();

void ampl_equations_1h2p_D4a();

void ampl_equations_1h2p_D4b_1();

void ampl_equations_1h2p_D4b_2();

void ampl_equations_1h2p_D5a();

void ampl_equations_1h2p_D5b_1();

void ampl_equations_1h2p_D5b_2();

void ampl_equations_1h2p_D5c_1();

void ampl_equations_1h2p_D5c_2();

void ampl_equations_1h2p_D5d_1();

void ampl_equations_1h2p_D5d_2();

void ampl_equations_1h2p_D5e_1();

void ampl_equations_1h2p_D5e_2();

void ampl_equations_1h2p_D5f();

void ampl_equations_1h2p_D5g_1();

void ampl_equations_1h2p_D5g_2();

void ampl_equations_1h2p_D5h();

void ampl_equations_1h2p_D6a();

void ampl_equations_1h2p_D6b();

void ampl_equations_1h2p_D6c_1();

void ampl_equations_1h2p_D6c_2();

void ampl_equations_1h2p_D7a();

void ampl_equations_1h2p_D7b();

void ampl_equations_1h2p_D7c_1();

void ampl_equations_1h2p_D7c_2();

void ampl_equations_1h2p_D7d();

void ampl_equations_1h2p_D7e_1();

void ampl_equations_1h2p_D7e_2();

void ampl_equations_1h2p_D8a_1();

void ampl_equations_1h2p_D8a_2();

void ampl_equations_1h2p_D8b();

void ampl_equations_1h2p_D9();

void ampl_equations_1h2p_folded_F1();

void ampl_equations_1h2p_folded_F2();

void ampl_equations_1h2p_folded_F3();

void ampl_equations_1h2p_folded_F4();

void ampl_equations_1h2p_folded_F5();

void ampl_equations_1h2p_folded_F6();

void ampl_equations_1h2p_folded_F7();

void ampl_equations_1h2p_folded_F8();

void ampl_equations_1h2p_folded_F9();

void ampl_equations_1h2p_folded_F10();

void diag_heff_1h2p_i_jk_c_ab__inv_h2_p13_const(int pt_order);

void diag_heff_1h2p_i_jk_c_ab__inv_h3_p12_const(int pt_order);

void diag_heff_1h2p_i_jk_c_ab__inv_h1_p23_const(int pt_order);

void diag_heff_1h2p_i_jk_c_ab__inv_h3_p23_m2c(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h2_p13_const(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h3_p12_const(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h1_p23_const(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h3_p13_const(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h2_p12_const(int pt_order);

void diag_heff_1h2p_ijk_c_ab__inv_h3_p23_const(int pt_order);

void diag_heff_1h2p_k_ij_a_bc__inv_h1_p23_const(int pt_order);

void diag_heff_1h2p_k_ij_a_bc__inv_h3_p12_const(int pt_order);

void diag_heff_1h2p_k_ij_a_bc__inv_h2_p13_const(int pt_order);

void diag_heff_1h2p_k_ij_a_bc__inv_h1_p12_m2c(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h1_p23_const(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h2_p13_const(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h3_p12_const(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h1_p12_m2c(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h2_p12_const(int pt_order);

void diag_heff_1h2p_k_ij_abc__inv_h3_p13_const(int pt_order);

void restrict_valence(char *src_name /*large*/, char *tgt_name /*small*/, char *new_valence, int extract_valence);


int sector_1h2p(cc_options_t *opts)
{
    printf("\n");
    printf("\t\t\t\t*****************\n");
    printf("\t\t\t\t** Sector 1h2p **\n");
    printf("\t\t\t\t*****************\n");
    printf("\n");

    if (opts->hughes_kaldor_1h2p) {
        printf("\n");
        printf(" The Hughes-Kaldor version of the FS-CCSD(1h2p) method will be used\n");
        printf(" (folded terms with Heff{1h2p} will be excluded).\n");
        printf("\n");
    }

    opts->curr_sector_h = 1;
    opts->curr_sector_p = 2;

    // setup shifts for the IH-like technique by A. V. Zaitsevskii
    if (intham_imms_in_sector(1, 2)) {
        intham_imms_setup(1, 2);
    }

    sort_integrals_1h2p();
    init_amplitudes_1h2p();
    predict_intruders_for_diagram("m2nw", 10);

    /*
     * initial approximation to Veff{1h2p}.
     * includes diagrams independent on the T{1h2p} cluster amplitudes.
     */
    heff12_ccsd_const();

    /*
     * main loop
     */
    int success = solve_amplitude_equations(
            1, 2, NULL, NULL, "m2c", "m2nw", NULL, NULL, NULL,
            NULL, construct_doubles_1h2p, NULL, construct_folded_1h2p
    );
    if (success == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    /*
     * analyze cluster amplitudes and write them to disk
     */
    print_cluster_operator_analysis(1, 2, NULL, "m2c", NULL);
    save_cluster_amplitudes(1, 2, NULL, "m2c", NULL, NULL);
    if (cc_opts->do_flush_amplitudes_txt) {
        save_cluster_amplitudes_formatted(1, 2, NULL, "m2c", NULL);
    }

    /*
     * effective Hamiltonian construction and analysis
     */
    copy("veff12_const", "veff12");
    heff12_ccsd_1h2p_dependent();

    heff_analysis(1, 2, "veff01", "veff10", "veff11", "veff02", "veff12");

    /*
     * calculate quasi-natural orbitals and model-space estimates of properties
     */
    calculate_model_space_properties(1, 2);
    calculate_model_space_natural_spinors(1, 2);

    return EXIT_SUCCESS;
}


void sort_integrals_1h2p()
{
    /*
     * two-electron integrals
     */
    request_sorting("gvhv", "hphp", "1101", "1234");
    request_sorting("pphv", "pphp", "0001", "1234");
    request_sorting("gvhh", "hphh", "1100", "1234");
    request_sorting("vghg", "phhh", "1101", "1234");
    request_sorting("vghh", "phhh", "1100", "1234");
    request_sorting("vvhv", "pphp", "1101", "1234");
    request_sorting("gphh", "hphh", "1000", "1234");
    request_sorting("pvhv", "pphp", "0101", "1234");
    request_sorting("pghg", "phhh", "0101", "1234");
    request_sorting("pghv", "phhp", "0101", "2431");
    request_sorting("pvvv", "pppp", "0111", "1234");
    request_sorting("ppvv", "pppp", "0011", "1234");
    request_sorting("pgvg", "phph", "0111", "1234");
    request_sorting("pghh", "phhh", "0100", "1234");
    request_sorting("ppvg", "ppph", "0011", "1234");
    request_sorting("pvgv", "pphp", "0111", "1234");
    request_sorting("ppgv", "pphp", "0011", "3412");
    request_sorting("pgvv", "phpp", "0111", "1234");
    request_sorting("ppvh", "ppph", "0010", "1234");

    request_sorting("vphh", "pphh", "1000", "1234");
    request_sorting("phpg", "phph", "0001", "1234");
    request_sorting("pppg", "ppph", "0001", "1234");
    request_sorting("ppph", "ppph", "0000", "1234");
    request_sorting("ppgh", "pphh", "0010", "1234");
    request_sorting("vvhg", "pphh", "1101", "1234");
    request_sorting("vvpg", "ppph", "1101", "1234");
    request_sorting("pvpg", "ppph", "0101", "1234");
    perform_sorting();

    /*
     * prepare cluster amplitudes from previous sectors
     */
    restrict_valence("t1c", "t1c-v", "01", 0);
    restrict_valence("t1c", "t1c-g", "10", 0);
    restrict_valence("t2c", "t2c-v12-g1", "1011", 0);
    restrict_valence("t2c", "t2c-v1-g1", "0101", 0);
    restrict_valence("t2c", "t2c-v12", "0011", 0);

    restrict_valence("s2c", "s2c-v1-g", "1110", 0);
    restrict_valence("s2c", "s2c-v12", "1011", 0);
    restrict_valence("s2c", "s2c-g", "1100", 0);
    restrict_valence("s2c", "s2c-v1", "1010", 0);

    restrict_valence("h2c", "h2c-g1-v", "1011", 0);
    restrict_valence("h2c", "h2c-g1", "1010", 0);
    restrict_valence("h2c", "h2c-v", "0011", 0);

    restrict_valence("e2c", "e2c-v1", "1011", 0);
    restrict_valence("e2c", "e2c-v", "1011", 0);
    restrict_valence("e2c", "e2c-g", "1101", 0);

    restrict_valence("x2c", "x2c-v1", "1110", 0);
    restrict_valence("x2c", "x2c-v2", "1101", 0);
}


void init_amplitudes_1h2p()
{
    tmplt("m2nw", "ppph", "1101", "1234", NOT_PERM_UNIQUE);
    copy("m2nw", "m2c");

    if (cc_opts->reuse_amplitudes[1][2]) {
        printf(" Trying to read amplitudes and effective interaction (sector 1h2p) from disk ...\n");
        // Doubles
        if (diagram_read("m2c.dg") != NULL) {
            printf(" T{12}_1 amplitudes successfully read from disk\n");
        }
        else {
            printf(" T{12}_1 amplitudes will be calculated\n");
        }
    }
    printf("\n");
}


void construct_doubles_1h2p()
{
    /*
     * (V\Omega)_conn
     */
    ampl_equations_1h2p_D1();
    ampl_equations_1h2p_D2a_1();
    ampl_equations_1h2p_D2a_2();
    ampl_equations_1h2p_D2b();
    ampl_equations_1h2p_D2c();
    ampl_equations_1h2p_D2d();
    ampl_equations_1h2p_D2e_1();
    ampl_equations_1h2p_D2e_2();
    ampl_equations_1h2p_D3a();
    ampl_equations_1h2p_D3b();
    ampl_equations_1h2p_D3c();
    ampl_equations_1h2p_D3d_1();
    ampl_equations_1h2p_D3d_2();
    ampl_equations_1h2p_D4a();
    ampl_equations_1h2p_D4b_1();
    ampl_equations_1h2p_D4b_2();
    ampl_equations_1h2p_D5a();
    ampl_equations_1h2p_D5b_1();
    ampl_equations_1h2p_D5b_2();
    ampl_equations_1h2p_D5c_1();
    ampl_equations_1h2p_D5c_2();
    ampl_equations_1h2p_D5d_1();
    ampl_equations_1h2p_D5d_2();
    ampl_equations_1h2p_D5e_1();
    ampl_equations_1h2p_D5e_2();
    ampl_equations_1h2p_D5f();
    ampl_equations_1h2p_D5g_1();
    ampl_equations_1h2p_D5g_2();
    ampl_equations_1h2p_D5h();
    ampl_equations_1h2p_D6a();
    ampl_equations_1h2p_D6b();
    ampl_equations_1h2p_D6c_1();
    ampl_equations_1h2p_D6c_2();
    ampl_equations_1h2p_D7a();
    ampl_equations_1h2p_D7b();
    ampl_equations_1h2p_D7c_1();
    ampl_equations_1h2p_D7c_2();
    ampl_equations_1h2p_D7d();
    ampl_equations_1h2p_D7e_1();
    ampl_equations_1h2p_D7e_2();
    ampl_equations_1h2p_D8a_1();
    ampl_equations_1h2p_D8a_2();
    ampl_equations_1h2p_D8b();
    ampl_equations_1h2p_D9();
}

void ampl_equations_1h2p_D1()
{
    tmplt("m2nw", "ppph", "1101", "1234", NOT_PERM_UNIQUE);
    update("m2nw", 1.0, "vvpg");
}

void ampl_equations_1h2p_D2a_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pg", "pgr", "21");
    mult("x2c", "pgr", "r1", 1);
    update("m2nw", 1.0, "r1");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2a_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pp", "ppr_", "21");
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "ppr_", "r1", 1);
    reorder("r1", "r2", "2143");
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e2c", "e2cr_", "3412");
    mult("e2cr_", "vh", "r1", 1);
    reorder("r1", "r2", "3412");
    perm("r2", "(12)");
    update("m2nw", -1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2c()
{
    dg_stack_pos_t pos = get_stack_pos();
    mult("ppgp", "x2c", "r1", 2);
    reorder("r1", "r2", "4321");
    update("m2nw", 0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2d()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "h2c_21", "2143");
    reorder("h2c_21", "r1", "3412");
    mult("vvhh", "r1", "r2", 2);
    update("m2nw", 0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s2c", "r1", "1324");
    mult("r1", "pvhg", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    update("m2nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D2e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e2c", "e2c_1243", "1243");
    reorder("e2c_1243", "r1", "1324");
    reorder("pvhp", "r0", "2431");
    mult("r1", "r0", "r3", 2);
    reorder("r3", "r4", "1324");
    perm("r4", "(12)");
    reorder("r4", "r5", "2143"); // 1 <-> 2
    update("m2nw", -1.0, "r5");  // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "h2c_21", "2143");
    reorder("h2c_21", "h2cr_", "3412");
    reorder("pphh", "v1", "3412");
    mult("x2c", "v1", "r1", 2);
    mult("r1", "h2cr_", "r2", 2);
    update("m2nw", 0.25, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e2c", "r1", "1423");
    reorder("pphh", "r2", "3142");
    reorder("s2c", "r3", "1324");
    mult("r1", "r2", "r4", 2);
    mult("r3", "r4", "r5", 2);
    reorder("r5", "r6", "1324");
    perm("r6", "(12)");
    update("m2nw", -1.0, "r6"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3c()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v1", "3412");
    mult("s2c", "v1", "r1", 3);
    reorder("e2c", "e2c_2134", "2134");
    reorder("e2c_2134", "r2", "2341");
    mult("r1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("m2nw", -0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3d_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("h2c", "t2cr_", "3412");
    mult("t2cr_", "pphh", "r1", 3);
    mult("x2c", "r1", "r2", 1);
    update("m2nw", -0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D3d_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t2c", "t2cr_", "3412");
    mult("t2cr_", "pphh", "r1", 3);
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "r1", "r2", 1);
    reorder("r2", "r3", "2143");
    update("m2nw", -0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D4a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvpg", "r1", "2341");
    mult("s1c", "r1", "r2", 1);
    perm("r2", "(12)");
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D4b_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("vvhg", "r1", "1243");
    reorder("t1c", "t1cr", "21");
    mult("r1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D4b_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("vvhp", "r1", "1243");
    reorder("h1c", "t1cr", "21");
    mult("r1", "t1cr", "r2", 1);
    reorder("r2", "r3", "2134");   // 1243 + interchange electrons 1 <-> 2
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "r1", "21");
    reorder("e2c", "r2", "1342");
    mult("s1c", "r1", "r3", 1);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1423");
    perm("r5", "(12)");
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5b_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "1243");
    reorder("ph", "r2", "21");
    reorder("t1c", "r3", "21");
    mult("r1", "r2", "r4", 1);
    mult("r4", "r3", "r5", 1);
    reorder("r5", "r6", "1243");
    update("m2nw", -1.0, "r6");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5b_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ph", "r1", "21");
    reorder("h1c", "r2", "21");
    mult("x2c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    update("m2nw", -1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("ppph", "r1", "3421");
    reorder("e2c", "e2c_2134", "2134");  // sign will be changed
    reorder("e2c_2134", "r2", "2413");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    update("m2nw", -1.0, "r5"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s2c", "r0", "2143");
    reorder("r0", "r3", "2413");
    reorder("ppgh", "r1", "3421");
    mult("s1c", "r1", "r2", 1);
    mult("r2", "r3", "r4", 2);
    reorder("r4", "r5", "1324");
    perm("r5", "(12)");
    reorder("r5", "r6", "2143");
    update("m2nw", 1.0, "r6");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5d_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2314");
    reorder("t1c", "t1cr", "21");
    mult("t1cr", "r1", "i1", 1);
    reorder("e2c", "e2c_1243", "1243"); // sign will be changed
    reorder("e2c_1243", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12)");
    reorder("r3", "r4", "2143");
    update("m2nw", 1.0, "r4"); // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5d_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2314");
    reorder("h1c", "t1cr", "21");
    mult("t1cr", "r1", "i1", 1);
    reorder("s2c", "r2", "1324");
    mult("r2", "i1", "i2", 2);
    reorder("i2", "r3", "1423");
    perm("r3", "(12)");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("x2c", "r1", "i1", 2);
    mult("i1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    update("m2nw", -0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphp", "r1", "4312");
    reorder("h1c", "t1cr", "21");
    mult("x2c", "r1", "i1", 2);
    mult("i1", "t1cr", "r2", 1);
    reorder("r2", "r3", "1243");
    reorder("r3", "r4", "2143");
    update("m2nw", -0.5, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5f()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2341");
    reorder("h2c", "h2c_2143", "2143");
    reorder("h2c_2143", "t2cr", "3412");
    mult("s1c", "r1", "i1", 1);
    mult("i1", "t2cr", "r2", 2);
    perm("r2", "(12)");
    update("m2nw", 0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5g_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphp", "r1", "4231");
    reorder("m2c", "r2", "1243");
    mult("r1", "t1c", "r3", 2);
    mult("r2", "r3", "r4", 1);
    reorder("r4", "r5", "1243");
    update("m2nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5g_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4231");
    mult("r1", "t1c", "i1", 2);
    mult("x2c", "i1", "r2", 1);
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D5h()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2431");
    mult("r1", "t1c", "i1", 2);
    reorder("e2c", "e2c_2134", "2134");
    reorder("e2c_2134", "r2", "2341");
    mult("i1", "r2", "r3", 1);
    perm("r3", "(12)");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6a()
{
    dg_stack_pos_t pos = get_stack_pos();
    mult("ppgp", "s1c", "i1", 1);
    reorder("i1", "r1", "4123");
    mult("s1c", "r1", "r2", 1);
    reorder("r2", "r3", "2143"); // interchange electrons 1 <-> 2
    update("m2nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "t1c_21", "21");
    reorder("h1c", "h1c_21", "21");
    mult("h1c_21", "vvhh", "i1", 1);
    mult("t1c_21", "i1", "i2", 1);
    reorder("i2", "r1", "3412");
    update("m2nw", 1.0, "r1");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "r1", "21");
    mult("s1c", "pvhg", "r2", 1);
    mult("r2", "r1", "r3", 1);
    reorder("r3", "r4", "1243");
    perm("r4", "(12)");
    update("m2nw", -1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D6c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhp", "r1", "2431");
    reorder("h1c", "r2", "21");
    mult("s1c", "r1", "r3", 1);
    mult("r3", "r2", "r4", 1);
    reorder("r4", "r5", "1243");
    perm("r5", "(12)");
    reorder("r5", "r6", "2143");  // interchange electrons 1 <-> 2
    update("m2nw", -1.0, "r6");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7a()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "vr1", "3412");
    reorder("h2c", "h2c_2143", "2143");
    reorder("h2c_2143", "t2cr", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("r2", "t2cr", "r3", 2);
    update("m2nw", 0.5, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "vr1", "3412");
    reorder("t1c", "t1cr", "21");
    reorder("h1c", "h1cr", "21");
    mult("x2c", "vr1", "r1", 2);
    mult("h1cr", "r1", "r2", 1);
    mult("t1cr", "r2", "r3", 1);
    reorder("r3", "r4", "3412");
    update("m2nw", 0.5, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7c_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "3142");
    reorder("e2c", "e2c_2134", "2134");  // the sign will be changed
    reorder("e2c_2134", "r2", "2413");
    reorder("t1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    perm("r7", "(12)");
    update("m2nw", 1.0, "r7");  // sign is changed
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7c_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "3142");
    reorder("s2c", "s2c_2143", "2143");
    reorder("s2c_2143", "r2", "2413");
    reorder("h1c", "r5", "21");
    mult("r2", "r1", "r3", 2);
    mult("s1c", "r3", "r4", 1);
    mult("r4", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    perm("r7", "(12)");
    reorder("r7", "r8", "2143");
    update("m2nw", -1.0, "r8");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7d()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v_", "2431");
    mult("v_", "t1c", "r1", 2);
    reorder("r1", "r2", "21");
    mult("s1c", "r2", "r3", 1);
    reorder("e2c", "e2c_2134", "2134");
    reorder("e2c_2134", "t2r", "2341");
    mult("r3", "t2r", "r4", 1);
    perm("r4", "(12)");
    update("m2nw", -1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7e_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "r1", "2431");
    reorder("t1c", "r2", "21");
    reorder("m2c", "r3", "1243");
    mult("r1", "t1c", "r4", 2);
    mult("r2", "r4", "r5", 1);
    mult("r3", "r5", "r6", 1);
    reorder("r6", "r7", "1243");
    update("m2nw", -1.0, "r7");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D7e_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphh", "v_", "2431");
    reorder("h1c", "h1cr", "21");
    mult("v_", "t1c", "r1", 2);
    mult("h1cr", "r1", "r2", 1);
    reorder("x2c", "x2r", "1243");
    mult("x2r", "r2", "r3", 1);
    reorder("r3", "r4", "1243");
    reorder("r4", "r5", "2143");  // interchange electrons 1 <-> 2
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D8a_1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphp", "r1", "4312");
    reorder("h1c", "h1cr", "21");
    mult("s1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "h1cr", "r4", 1);
    reorder("r4", "r5", "2134");  // 1243 + interchange electrons 1<->2
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D8a_2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pphg", "r1", "4312");
    reorder("t1c", "t1cr", "21");
    mult("s1c", "r1", "r2", 1);
    mult("s1c", "r2", "r3", 1);
    mult("r3", "t1cr", "r4", 1);
    reorder("r4", "r5", "1243");
    update("m2nw", -1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D8b()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("pvhh", "r1", "2341");
    reorder("t1c", "t1cr", "21");
    reorder("h1c", "h1cr", "21");
    mult("s1c", "r1", "i1", 1);
    mult("h1cr", "i1", "i2", 1);
    mult("t1cr", "i2", "r2", 1);
    reorder("r2", "r3", "3412");
    perm("r3", "(12)");
    update("m2nw", 1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_D9()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("t1c", "t1cr", "21");
    reorder("h1c", "h1cr", "21");
    reorder("pphh", "vr1", "3412");
    mult("s1c", "vr1", "r1", 1);
    mult("s1c", "r1", "r2", 1);
    mult("h1cr", "r2", "r3", 1);
    mult("t1cr", "r3", "r4", 1);
    reorder("r4", "r5", "3412");
    update("m2nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void construct_folded_1h2p()
{
    /*
     * effective interaction will be used in folded diagrams F9 and F10
     */
    if (cc_opts->hughes_kaldor_1h2p == 0) {
        copy("veff12_const", "veff12");
        heff12_ccsd_1h2p_dependent();
    }

    ampl_equations_1h2p_folded_F1();
    ampl_equations_1h2p_folded_F2();
    ampl_equations_1h2p_folded_F3();
    ampl_equations_1h2p_folded_F4();
    ampl_equations_1h2p_folded_F5();
    ampl_equations_1h2p_folded_F6();
    ampl_equations_1h2p_folded_F7();
    ampl_equations_1h2p_folded_F8();

    /*
     * folded terms involving the Heff{1h2p} operator
     */
    if (cc_opts->hughes_kaldor_1h2p == 0) {
        ampl_equations_1h2p_folded_F9();
        ampl_equations_1h2p_folded_F10();
    }
}

void ampl_equations_1h2p_folded_F1()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "2341");
    mult("veff01", "r1", "r2", 1);
    perm("r2", "(12)");
    update("m2nw", -1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F2()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff10", "r1", "21");
    mult("m2c", "r1", "r2", 1);
    update("m2nw", 1.0, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F3()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "3412");
    mult("veff02", "r1", "r2", 2);
    update("m2nw", -0.5, "r2");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F4()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("e1c", "e1c_21", "21");
    mult("veff02", "e1c_21", "r1", 1);
    tmplt("r1_ext", "ppph", "1101", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r1", "r1_ext");
    update("m2nw", -1.0, "r1_ext");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F5()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("s1c", "s1c_21", "21");
    reorder("e1c", "e1c_21", "21");
    mult("e1c_21", "veff02", "r1", 1);
    mult("s1c_21", "r1", "r2", 1);
    reorder("r2", "r3", "3412");
    update("m2nw", -1.0, "r3");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F6()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("m2c", "r1", "2341");
    reorder("veff11", "r2", "1423");
    mult("r2", "r1", "r3", 2);
    reorder("r3", "r4", "1342");
    perm("r4", "(12)");
    update("m2nw", 1.0, "r4");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F7()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff11", "r1", "3412");
    mult("r1", "e1c", "r2", 1);
    reorder("r2", "r3", "3412");
    tmplt("r3_ext", "ppph", "1101", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r3", "r3_ext");
    perm("r3_ext", "(12)");
    update("m2nw", 1.0, "r3_ext");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F8()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff11", "r1", "1423");
    reorder("s1c", "r2", "21");
    mult("r2", "r1", "r3", 1);
    mult("r3", "e1c", "r4", 1);
    reorder("r4", "r5", "2413");
    perm("r5", "(12)");
    update("m2nw", 1.0, "r5");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F9()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff12", "r1", "124653");
    mult("r1", "e1c", "r2", 2);
    tmplt("r2_ext", "ppph", "1101", "1234", NOT_PERM_UNIQUE);
    expand_diagram("r2", "r2_ext");
    update("m2nw", 1.0, "r2_ext");
    restore_stack_pos(pos);
}

void ampl_equations_1h2p_folded_F10()
{
    dg_stack_pos_t pos = get_stack_pos();
    reorder("veff12", "r1", "126345");
    reorder("m2c", "r2", "3412");
    mult("r1", "r2", "r3", 3);
    reorder("r3", "r4", "1243");
    update("m2nw", 0.5, "r4");
    restore_stack_pos(pos);
}


/**
 * The part of Heff{1h2p} which is constant in the 1h2p sector
 * (no diagrams depending on the T{1h2p} amplitudes).
 * This part is calculated only once at start.
 * Result is stored in the "veff12_const" diagram.
 */
void heff12_ccsd_const()
{
    timer_new_entry("12-CCSD-Heff", "1h2p -- CCSD part of eff Hamiltonian");
    timer_start("12-CCSD-Heff");

    printf("\n Construction of the CCSD eff Hamiltonian in the 1h2p sector (constant part)\n");

    tmplt("veff12_const", "pphpph", "111111", "123456", NOT_PERM_UNIQUE);

    printf(" [1/4] P(i/jk|c/ab)\n");
    diag_heff_1h2p_i_jk_c_ab__inv_h2_p13_const(PT_INF);
    diag_heff_1h2p_i_jk_c_ab__inv_h3_p12_const(PT_INF);
    diag_heff_1h2p_i_jk_c_ab__inv_h1_p23_const(PT_INF);
    //diag_heff_1h2p_i_jk_c_ab__inv_h3_p23_m2c(PT_INF);

    printf(" [2/4] P(ijk|c/ab)\n");
    diag_heff_1h2p_ijk_c_ab__inv_h2_p13_const(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h3_p12_const(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h1_p23_const(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h3_p13_const(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h2_p12_const(PT_INF);
    diag_heff_1h2p_ijk_c_ab__inv_h3_p23_const(PT_INF);

    printf(" [3/4] P(k/ij|a/bc)\n");
    diag_heff_1h2p_k_ij_a_bc__inv_h1_p23_const(PT_INF);
    diag_heff_1h2p_k_ij_a_bc__inv_h3_p12_const(PT_INF);
    diag_heff_1h2p_k_ij_a_bc__inv_h2_p13_const(PT_INF);
    //diag_heff_1h2p_k_ij_a_bc__inv_h1_p12_m2c(PT_INF);   // S2_1h2p is required

    printf(" [4/4] P(k/ij|abc)\n");
    diag_heff_1h2p_k_ij_abc__inv_h1_p23_const(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h2_p13_const(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h3_p12_const(PT_INF);
    //diag_heff_1h2p_k_ij_abc__inv_h1_p12_m2c(PT_INF);   // S2_1h2p is required
    diag_heff_1h2p_k_ij_abc__inv_h2_p12_const(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h3_p13_const(PT_INF);

    printf("\n");

    timer_stop("12-CCSD-Heff");
}


/**
 * The part of Heff{1h2p} which is dependent on the T{1h2p} amplitudes.
 * Should be recalculated at each iteration.
 * Result is added to the "veff12" diagram.
 */
void heff12_ccsd_1h2p_dependent()
{
    timer_new_entry("12-CCSD-Heff", "1h2p -- CCSD part of eff Hamiltonian");
    timer_start("12-CCSD-Heff");

    diag_heff_1h2p_i_jk_c_ab__inv_h3_p23_m2c(PT_INF);
    diag_heff_1h2p_k_ij_a_bc__inv_h1_p12_m2c(PT_INF);
    diag_heff_1h2p_k_ij_abc__inv_h1_p12_m2c(PT_INF);

    timer_stop("12-CCSD-Heff");
}



/*
 * Final T3 permutation: (1/23|6/45) aka (i/jk|c/ab)
 * Diagrams included:
 * (i/jk|c/ab): T1b T3a T3d T4b T5d T7c T8a T8e T10a
 */
void diag_heff_1h2p_i_jk_c_ab__inv_h2_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("gvhv", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("s2c-v1-g", "s2c_21", "2143");
        reorder("s2c_21", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("s2c-g", "s2c_21", "2143");
        reorder("pphv", "r2", "4312");
        mult("s2c_21", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c-v", "r2", "21");
        mult("gvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphv", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("t1c-g", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("s2c-v1-g", "s2c_21", "2143");
        reorder("s2c_21", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c-v", "r3", "21");
        reorder("s2c-g", "s2c_21", "2143");
        mult("s2c_21", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c-v", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c-v1", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons 2 <-> 3
    update("veff12_const", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_i_jk_c_ab__inv_h3_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vghg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("e2c-g", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphg", "r2", "4312");
        mult("s2c-g", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("h1c", "r2", "21");
        mult("vghh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphg", "r1", "4312");
        mult("t1c-g", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("e2c-g", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("s2c-g", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("t1c-g", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c-v12", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12)");
    update("veff12_const", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_i_jk_c_ab__inv_h1_p23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhv", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        reorder("ph", "phr_", "21");
        reorder("x2c-v2", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphv", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("t1c-v", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphv", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        reorder("pphh", "r1", "4231");
        reorder("x2c-v2", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("t1c-v", "r3", "21");
        mult("x2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("t1c-v", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c-g1-v", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(56)");
    reorder("r4", "r5", "321654");
    update("veff12_const", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_i_jk_c_ab__inv_h3_p23_m2c(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1b
        reorder("vvhg", "r2", "1243");
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3a
        // m2c-dependent
        reorder("ph", "phr_", "21");
        reorder("m2c", "r1", "1243");
        mult("r1", "phr_", "r2", 1);
        update("i1", -1.0, "r2");
        restore_stack_pos(pos2);

        // T3d
        reorder("pphg", "r2", "4312");
        mult("x2c", "r2", "r3", 2);
        update("i1", -0.5, "r3");
        restore_stack_pos(pos2);

        // T4b
        reorder("h1c", "r2", "21");
        mult("vvhh", "r2", "r3", 1);
        reorder("r3", "r4", "1243");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T7c
        reorder("pphg", "r1", "4312");
        mult("s1c", "r1", "r2", 1);
        mult("s1c", "r2", "r3", 1);
        update("i1", -1.0, "r3");
        restore_stack_pos(pos2);

        // T8a
        // m2c-dependent
        reorder("pphh", "r1", "4231");
        reorder("m2c", "r3", "1243");
        mult("r1", "t1c", "r4", 2);
        mult("r3", "r4", "r5", 1);
        update("i1", -1.0, "r5");
        restore_stack_pos(pos2);

        // T8e
        reorder("pphh", "r1", "3412");
        reorder("h1c", "r3", "21");
        mult("x2c", "r1", "r2", 2);
        mult("r2", "r3", "r4", 1);
        reorder("r4", "r5", "1243");
        update("i1", 0.5, "r5");
        restore_stack_pos(pos2);

        // T10a
        reorder("h1c", "r1", "21");
        reorder("pphh", "r2", "3412");
        mult("s1c", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "1243");
        update("i1", 1.0, "r6");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c-v12-g1", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    reorder("r4", "r5", "231546");
    update("veff12", -1.0, "r5");  // the sign is changed from +1 to -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (123|6/45) aka (ijk|c/ab)
 * Diagrams included:
 * (ijk|c/ab):  T3c T4c T7b T8b
 */
void diag_heff_1h2p_ijk_c_ab__inv_h2_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hpph", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("s2c-v1", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("gphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pvhv", "r0", "2431");
        mult("t1c-g", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c-v", "r1", "21");
        reorder("gphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c-v1", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c-g", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c-v", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13|46)");
    reorder("r4", "r5", "132465"); // interchange electrons: 2 <-> 3
    update("veff12_const", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h3_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phhh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("h2c-g1", "h2c_21", "2143");
        reorder("h2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pghg", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("h2c-g1", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c-v12", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12)");
    update("veff12_const", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h1_p23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "ppph", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("s2c-v1", "s2c_21", "2143");
        reorder("s2c_21", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        reorder("pvhv", "r0", "2431");
        mult("s1c", "r0", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c-v", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("s2c-v1", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("h2c-g1-v", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23|56)");
    reorder("r4", "r5", "321654");
    update("veff12_const", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h3_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r1", "2413");
        reorder("gphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4"); // change sign: x -1
        restore_stack_pos(pos2);

        // T4c
        mult("t1c-g", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("gphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("e2c", "e2c_2134", "2134");
        reorder("pphh", "r1", "3142");
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("t1c-g", "r3", "r4", 1);
        update("i1", 1.0, "r4"); // change sign: x -1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("s2c-v12", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(13)");
    reorder("r4", "r5", "132456");
    update("veff12_const", 1.0, "r5");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h2_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phph", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("t2c-v1-g1", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T4c
        mult("s1c", "pghv", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("t1c-v", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("pphh", "r1", "3142");
        reorder("t2c-v1-g1", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c-v", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(12|46)");
    reorder("r4", "r5", "123465");     // <ijk|abc> -> <ijk|acb>
    update("veff12_const", -1.0, "r5"); // change sign: x -1
    restore_stack_pos(pos);
}

void diag_heff_1h2p_ijk_c_ab__inv_h3_p23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphh", "1110", "1243", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3c
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r1", "2413");
        reorder("vphh", "r2", "3142");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");  // change sign: x -1
        restore_stack_pos(pos2);

        // T4c
        mult("s1c", "pvhg", "r1", 1);
        update("i1", -1.0, "r1");
        restore_stack_pos(pos2);

        // T7b
        reorder("h1c", "r1", "21");
        reorder("vphh", "r2", "1324");
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        reorder("r4", "r5", "3124");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T8b
        reorder("e2c", "e2c_2134", "2134");
        reorder("pphh", "r1", "3142");
        reorder("e2c_2134", "r2", "2413");
        mult("r2", "r1", "r3", 2);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");  // change sign: x -1
        restore_stack_pos(pos2);
    }

    finish:
    reorder("t2c-v12-g1", "r1", "3412");
    mult("r1", "i1", "r3", 1);
    reorder("r3", "r4", "345126");
    perm("r4", "(23)");
    reorder("r4", "r5", "231546");
    update("veff12_const", -1.0, "r5"); // change sign: x -1
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|4/56) aka (k/ij|a/bc)
 * Diagrams included:
 * (k/ij|a/bc): T1a T3e T4a T5e T7d T8d T10b
 */
void diag_heff_1h2p_k_ij_a_bc__inv_h1_p23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvvv", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("t2c-v12", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppvv", "r0", "3412");
        mult("r0", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c-v", "r1", "21");
        mult("r1", "pvhh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("t2c-v12", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c-v", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("e2c-g", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(23)");
    reorder("r3", "r4", "321654");  // 1 <-> 3
    update("veff12_const", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_a_bc__inv_h3_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pgvg", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("h2c-v", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "pghh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppvg", "r0", "3412");
        mult("r0", "t1c-g", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("h1c", "r0", "21");
        reorder("t1c-v", "r1", "21");
        mult("r0", "pghh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c-v", "h2c_21", "2143");
        reorder("h2c_21", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c-g", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "r0", "21");
        reorder("t1c-v", "r1", "21");
        mult("r0", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c-v1", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(45)");
    update("veff12_const", 1.0, "r3");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_a_bc__inv_h2_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pvgv", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("h2c-v", "r1", "3412");
        mult("r1", "pvhh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        mult("ppgv", "s1c", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c-v", "t1c_21", "21");
        reorder("h1c", "h1c_21", "21");
        mult("t1c_21", "pvhh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("h2c-v", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("s1c", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("h1c", "h1c_21", "21");
        reorder("t1c-v", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("h1c_21", "r2", "r3", 1);
        mult("s1c", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c-v1-g", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|46)");
    reorder("r3", "r4", "132465");
    update("veff12_const", 1.0, "r4");
    restore_stack_pos(pos);
}

// S_2^{1h2p} is required
void diag_heff_1h2p_k_ij_a_bc__inv_h1_p12_m2c(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hppp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T1a
        reorder("pgvv", "r1", "2341");
        update("i1", 1.0, "r1");
        restore_stack_pos(pos2);

        // T3e
        reorder("t2c-v12", "r1", "3412");
        mult("r1", "pghh", "r2", 2);
        reorder("r2", "r3", "4123");
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T4a
        reorder("ppvv", "r0", "3412");
        mult("r0", "t1c-g", "r1", 1);
        reorder("r1", "r2", "4123");
        update("i1", 1.0, "r2");
        restore_stack_pos(pos2);

        // T7d
        reorder("t1c-v", "r1", "21");
        mult("r1", "pghh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "4123");
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);

        // T8d
        reorder("t2c-v12", "r1", "3412");
        mult("r1", "pphh", "r2", 2);
        mult("t1c-g", "r2", "r3", 1);
        update("i1", 0.5, "r3");
        restore_stack_pos(pos2);

        // T10b
        reorder("t1c-v", "r1", "21");
        mult("r1", "pphh", "r2", 1);
        mult("r1", "r2", "r3", 1);
        mult("t1c-g", "r3", "r4", 1);
        update("i1", 1.0, "r4");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    reorder("r3", "r4", "123654");
    update("veff12", -1.0, "r4");
    restore_stack_pos(pos);
}


/*
 * Final T3 permutation: (3/12|456) aka (k/ij|abc)
 * Diagrams included:
 * (k/ij|abc):  T3b T4d T7a T8c
 */
void diag_heff_1h2p_k_ij_abc__inv_h1_p23_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pppp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppvh", "r1", "1342");
        reorder("s2c-v1", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhv", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c-v", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c-v", "r1", "21");
        reorder("ppvh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c-v", "r1", "21");
        reorder("s2c-v1", "s2c_21", "2143");
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
    reorder("e2c-g", "e2c_21", "2143");
    mult("e2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(23|56)");
    reorder("r3", "r4", "321654");  // interchange electrons 1 <-> 3
    update("veff12_const", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h2_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "phpp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("s2c-v1", "s2c_21", "2143");
        reorder("s2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhv", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c-v", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("s2c-v1", "s2c_21", "2143");
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
    mult("s2c-v1-g", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|46)");
    reorder("r3", "r4", "132465"); // interchange electrons 2 <-> 3
    update("veff12_const", 1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h3_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hphp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppvh", "r1", "1342");
        reorder("h2c-g1", "h2c_21", "2143");
        reorder("h2c_21", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pghg", "r0", "2431");
        reorder("r0", "r1", "4123");
        reorder("t1c-v", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppvh", "r2", "1342");
        mult("r2", "t1c-g", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c-v", "r1", "21");
        reorder("h2c-g1", "h2c_21", "2143");
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
    mult("x2c-v1", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(45)");
    update("veff12_const", 1.0, "r3");
    restore_stack_pos(pos);
}

// S_2^{1h2p} is required
void diag_heff_1h2p_k_ij_abc__inv_h1_p12_m2c(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    if (pt_order <= PT_2) {
        return;
    }

    tmplt("i1", "hppp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_1B_PRIME) {
            goto finish;
        }

        // T3b
        reorder("ppvh", "r1", "1342");
        reorder("t2c-v1-g1", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        if (pt_order == PT_INF && cc_opts->cc_model <= CC_MODEL_CCSDT_2) {
            goto finish;
        }

        // T4d
        reorder("pghv", "r1", "4123");
        reorder("t1c-v", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        if (pt_order <= PT_3) {
            goto finish;
        }

        // T7a
        reorder("t1c-v", "r1", "21");
        reorder("ppvh", "r2", "1342");
        mult("r2", "t1c-g", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c-v", "r1", "21");
        reorder("t2c-v1-g1", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    reorder("m2c", "m2c_21", "2143");
    mult("m2c_21", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(56)");
    reorder("r3", "r4", "123654");
    update("veff12", -1.0, "r4");
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h2_p12_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "hhpp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppgh", "r1", "1342");
        reorder("t2c-v1-g1", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", 1.0, "r5");
        restore_stack_pos(pos2);

        // T4d
        reorder("pghv", "r1", "4123");
        reorder("h1c", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("t1c-v", "r1", "21");
        reorder("ppgh", "r2", "1342");
        mult("r2", "t1c-g", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("h1c", "r1", "21");
        reorder("t2c-v1-g1", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", -1.0, "r7");
        restore_stack_pos(pos2);
    }

    finish:
    mult("x2c-v1", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(46)");
    reorder("r3", "r4", "123465");
    update("veff12_const", -1.0, "r4"); // change sign
    restore_stack_pos(pos);
}

void diag_heff_1h2p_k_ij_abc__inv_h3_p13_const(int pt_order)
{
    dg_stack_pos_t pos;
    pos = get_stack_pos();

    tmplt("i1", "pphp", "1110", "2341", NOT_PERM_UNIQUE);

    {
        dg_stack_pos_t pos2;
        pos2 = get_stack_pos();

        // T3b
        reorder("ppvh", "r1", "1342");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        mult("r1", "r2", "r3", 2);
        reorder("r3", "r4", "1324");
        reorder("r4", "r5", "2341");
        update("i1", -1.0, "r5");  // change sign: x -1
        restore_stack_pos(pos2);

        // T4d
        reorder("pvhg", "r1", "4123");
        reorder("t1c-v", "r2", "21");
        mult("r1", "r2", "r3", 1);
        reorder("r3", "r4", "2431");
        update("i1", -1.0, "r4");
        restore_stack_pos(pos2);

        // T7a
        reorder("h1c", "r1", "21");
        reorder("ppvh", "r2", "1342");
        mult("r2", "s1c", "r3", 1);
        reorder("r3", "r4", "1423");
        mult("r4", "r1", "r5", 1);
        reorder("r5", "r6", "2341");
        update("i1", -1.0, "r6");
        restore_stack_pos(pos2);

        // T8c
        reorder("t1c-v", "r1", "21");
        reorder("e2c", "e2c_2134", "2134");
        reorder("e2c_2134", "r2", "2413");
        reorder("pphh", "r3", "1342");
        mult("r3", "r2", "r4", 2);
        reorder("r4", "r5", "1342");
        mult("r5", "r1", "r6", 1);
        reorder("r6", "r7", "2431");
        update("i1", 1.0, "r7");  // change sign: x -1
        restore_stack_pos(pos2);
    }

    finish:
    mult("s2c-v1-g", "i1", "r2", 1);
    reorder("r2", "r3", "124356");
    perm("r3", "(13|45)");
    reorder("r3", "r4", "132456");
    update("veff12_const", -1.0, "r4");  // change sign: x -1
    restore_stack_pos(pos);
}

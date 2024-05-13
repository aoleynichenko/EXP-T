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
 * Utility functions for FSCC models.
 */

#ifndef CC_CCUTILS_H_INCLUDED
#define CC_CCUTILS_H_INCLUDED

#include "comdef.h"
#include "options.h"

void print_sector_banner(int sect_h, int sect_p);

void damping(int h, int p, char *old_ampl, char *new_ampl, int iter);

void print_max_ampl(char *label, char *dg_name);

double t1_diagnostic(char *dg_name, int nelec);

int check_divergence(char *dg_name);

void print_ampl_norm(char *label, char *dg_name);

void get_cc_model_name(int sect_h, int sect_p, cc_model_t cc_model, char *cc_model_name);

void print_ionization_potential(char *label, double ion_pot);

void print_sector_header(int sect_h, int sect_p);

void intruder_state_analysis(char *oper_t1, char *oper_t2, char *oper_t3);

void save_cluster_amplitudes(int sect_h, int sect_p, char *oper_t1, char *oper_t2, char *oper_t3, char *oper_veff);

void save_cluster_amplitudes_formatted(int sect_h, int sect_p, char *oper_t1, char *oper_t2, char *oper_t3);

void print_cluster_operator_analysis(int sect_h, int sect_p, char *oper_t1, char *oper_t2, char *oper_t3);

void print_iterations_table_header(int sect_h, int sect_p, int singles_column, int doubles_column, int triples_column);

void print_iterations_table_footer(int singles_column, int doubles_column, int triples_column, int num_iterations,
                                   int converged);

int solve_amplitude_equations(
        int sector_h, int sector_p,
        char *singles, char *singles_buf,
        char *doubles, char *doubles_buf,
        char *triples, char *triples_buf,
        char *veff,
        void (*construct_singles)(),
        void (*construct_doubles)(),
        void (*construct_triples)(),
        void (*construct_folded)()
);

#endif /* CC_CCUTILS_H_INCLUDED */

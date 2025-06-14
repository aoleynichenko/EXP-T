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

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "ecp.h"
#include "elements.h"
#include "molecule.h"

#define MAX_ANG_MOM    10
#define MAX_CONTRACTED 100
#define MAX_PRIMITIVES 30
#define MAX(a, b) (((a)>(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))

void gen_dirac_mol_file(FILE *mol_file, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);

void gen_dirac_inp_file(FILE *inp_file, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib);

void gen_dirac_mol_file_symmetry_label(char *buf, molecule_t *mol);

void basis_get_blocks_for_dirac(basis_t *bas, int *n_blocks, int *n_block_sizes, int len_block_sizes);

char angular_momentum_to_char(int l);

int basis_get_max_ang_mom(basis_t *bas);

void print_l_blocks_dirac(FILE *out, basis_t *bas, int L);

int basis_get_number_l_blocks(basis_t *bas, int L);

int has_pseudopotential(ecp_lib_t *ecp_lib);

/**
 * generates DIRAC 'mol' file.
 * about mol format:
 * http://www.diracprogram.org/doc/release-19/molecule_and_basis/molecule_using_mol.html
 */
void gen_pam(
    FILE *out, char *mol_file_name, char *inp_file_name,
    molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib
)
{
    FILE *f_mol = fopen(mol_file_name, "w");
    gen_dirac_mol_file(f_mol, mol, bas_lib, ecp_lib);
    fprintf(out, "geometry, basis set, ecp -> %s\n", mol_file_name);
    fclose(f_mol);

    FILE *f_inp = fopen(inp_file_name, "w");
    gen_dirac_inp_file(f_inp, mol, bas_lib, ecp_lib);
    fprintf(out, "template for settings -> %s\n", inp_file_name);
    fclose(f_inp);
}


/**
 * generates DIRAC 'mol' file.
 * about mol format:
 * http://www.diracprogram.org/doc/release-19/molecule_and_basis/molecule_using_mol.html
 */
void gen_dirac_mol_file(FILE *mol_file, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    // first three lines -- comment
    fprintf(mol_file, "\n");
    fprintf(mol_file, "\n");
    fprintf(mol_file, "\n");

    // for the Cinfv symmetry the ghost atom is added to ensure the axial symmetry
    if (mol->sym_group.group == SYMMETRY_Cinfv) {
        // add a ghost atom
        molecule_add_atom(mol, 0, mol->sym_group.xyz[0] * 100,
                          mol->sym_group.xyz[1] * 100, mol->sym_group.xyz[2] * 100);
        // and its absent (NULL) basis set
    }

    // 'C' for Cartesian GTO basis;
    // number of atom types
    // 'A' label for angstroms if needed
    int n_atom_types = molecule_n_atom_types(mol);
    int n_point_charge_types = molecule_n_point_charge_types(mol);
    fprintf(mol_file, "C  %2d", n_atom_types + n_point_charge_types);

    // symmetry data
    char sym_buf[16];
    gen_dirac_mol_file_symmetry_label(sym_buf, mol);
    fprintf(mol_file, "%s\n", sym_buf);

    // loop over atom types
    // heavier elements first
    for (int ielem = N_CHEM_ELEMENTS - 1; ielem >= 0; ielem--) {
        char elem_sym[MAX_ELEMENT_SYMBOL];

        int n_atoms = molecule_n_atoms_of(mol, ielem);
        if (n_atoms == 0) {
            continue;
        }

        get_element_symbol(ielem, elem_sym);

        // generate XYZ coordinates for atoms of the given type
        fprintf(mol_file, "      %3d.   %2d\n", ielem, n_atoms);
        for (int j = 0; j < mol->n_atoms; j++) {
            if (mol->charges[j] == ielem) {
                fprintf(mol_file, "%s%16.10f%16.10f%16.10f\n", elem_sym, mol->x[j], mol->y[j], mol->z[j]);
            }
        }

        // generate basis
        basis_t *bas = basis_lib_get(bas_lib, elem_sym);
        if (bas == NULL) {
            fprintf(mol_file, "LARGE 0\n");
            continue;
        }

        int max_ang_mom = basis_get_max_ang_mom(bas);

        fprintf(mol_file, "LARGE EXPLICIT%5d", max_ang_mom + 1);
        for (int l = 0; l <= max_ang_mom; l++) {
            int n_blocks_dirac = basis_get_number_l_blocks(bas, l);
            fprintf(mol_file, "%5d", n_blocks_dirac);
        }
        fprintf(mol_file, "\n");

        for (int l = 0; l <= max_ang_mom; l++) {
            print_l_blocks_dirac(mol_file, bas, l);
        }

        // print ECP if needed
        ecp_t *ecp = ecp_lib_get(ecp_lib, elem_sym);
        if (ecp == NULL) {
            continue;
        }
        int len_arep, len_esop;
        ecp_get_len(ecp, &len_arep, &len_esop);
        fprintf(mol_file, "ECP %d %d %d\n", ecp->n_core_elec, len_arep, len_esop);

        for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
            ecp_expansion_t *f = ecp->arep[i];
            if (f == NULL && i < len_arep) { // no expansion, but it is needed
                fprintf(mol_file, "1\n");
                fprintf(mol_file, "%4d%24.12e%24.12e\n", 0, 1.0, 0.0);
                continue;
            }
            else if (f == NULL) {
                continue;
            }
            if (i == 0) {
                fprintf(mol_file, "# ul\n");
            }
            else {
                fprintf(mol_file, "# %c\n", tolower(angular_momentum_to_char(i - 1)));
            }
            fprintf(mol_file, "%d\n", f->nprim);
            for (int j = 0; j < f->nprim; j++) {
                fprintf(mol_file, "%4d%24.12e%24.12e\n", f->powers[j], f->e[j], f->c[j]);
            }
        }
        for (int i = 0; i <= ECP_MAX_ANG_MOM; i++) {
            ecp_expansion_t *f = ecp->esop[i];
            if (f == NULL && i < len_esop && i >= 2) { // no expansion, but it is needed
                fprintf(mol_file, "1\n");
                fprintf(mol_file, "%4d%24.12e%24.12e\n", 0, 1.0, 0.0);
                continue;
            }
            else if (f == NULL) {
                continue;
            }
            else {
                fprintf(mol_file, "# %c (spin-orbit)\n", tolower(angular_momentum_to_char(i - 1)));
            }
            fprintf(mol_file, "%d\n", f->nprim);
            for (int j = 0; j < f->nprim; j++) {
                fprintf(mol_file, "%4d%24.12e%24.12e\n", f->powers[j], f->e[j], f->c[j]);
            }
        }
    }

    // generate blocks of point charges
    if (mol->n_point_charges > 0) {

        double point_charge_types[1000];
        int n_point_charges_each_type[1000];
        int n_point_charge_types = 0;

        memset(point_charge_types, 0, sizeof(point_charge_types));
        memset(n_point_charges_each_type, 0, sizeof(n_point_charges_each_type));

        for (int i = 0; i < mol->n_point_charges; i++) {
            double q = mol->point_charges[i];

            int charge_found = 0;
            for (int j = 0; j < n_point_charge_types; j++) {
                if (fabs(point_charge_types[j] - q) < 1e-6) {
                    charge_found = 1;
                    n_point_charges_each_type[j]++;
                    break;
                }
            }

            // new type of point charge
            if (!charge_found) {
                point_charge_types[n_point_charge_types] = q;
                n_point_charges_each_type[n_point_charge_types] = 1;
                n_point_charge_types++;
            }
        }

        for (int qtype = 0; qtype < n_point_charge_types; qtype++) {

            double q = point_charge_types[qtype];

            fprintf(mol_file, "%10.6f%5d\n", q, n_point_charges_each_type[qtype]);

            for (int i = 0; i < mol->n_point_charges; i++) {
                if (fabs(mol->point_charges[i] - q) < 1e-6) {
                    fprintf(mol_file, "%s%16.10f%16.10f%16.10f\n", "Gh", mol->qx[i], mol->qy[i], mol->qz[i]);
                }
            }

            fprintf(mol_file, "LARGE 0\n");
        }

    }

    // generate charged cage (geodesic polyhedron) if needed
    if (mol->cage.n_points > 0) {
        fprintf(mol_file, "%10.6f%5d\n", mol->cage.total_charge / mol->cage.n_points, mol->cage.n_points);
        for (int i = 0; i < mol->cage.n_points; i++) {
            double x = mol->cage.points[3 * i];
            double y = mol->cage.points[3 * i + 1];
            double z = mol->cage.points[3 * i + 2];
            fprintf(mol_file, "%s%16.10f%16.10f%16.10f\n", "Z", x, y, z);
        }
        fprintf(mol_file, "LARGE POINTCHARGE\n");
    }

    fprintf(mol_file, "FINISH\n");
}


/**
 * Writes a string containing basic symmetry operations into a buffer
 * (DIRAC format, see https://www.diracprogram.org/doc/release-24/molecule_and_basis/molecule_using_mol.html#how-to-force-specific-symmetry )
 */
void gen_dirac_mol_file_symmetry_label(char *buf, molecule_t *mol)
{
    sprintf(buf, ""); // default - empty string

    if (mol->sym_group.group == SYMMETRY_C1) {
        sprintf(buf, "    0");
    }
    else if (mol->sym_group.group == SYMMETRY_Ci) {
        sprintf(buf, "    1XYZ");
    }
    else if (mol->sym_group.group == SYMMETRY_Cs) {
        if (mol->sym_group.xyz[0] == 1 && mol->sym_group.xyz[1] == 1) {
            sprintf(buf, "    1  Z");
        }
        else if (mol->sym_group.xyz[0] == 1 && mol->sym_group.xyz[2] == 1) {
            sprintf(buf, "    1  Y");
        }
        else if (mol->sym_group.xyz[1] == 1 && mol->sym_group.xyz[2] == 1) {
            sprintf(buf, "    1  X");
        }
    }
    else if (mol->sym_group.group == SYMMETRY_C2) {
        if (mol->sym_group.xyz[0] == 1) {
            sprintf(buf, "    1 YZ");
        }
        else if (mol->sym_group.xyz[1] == 1) {
            sprintf(buf, "    1X Z");
        }
        else if (mol->sym_group.xyz[2] == 1) {
            sprintf(buf, "    1XY");
        }
    }
    else if (mol->sym_group.group == SYMMETRY_C2v) {
        if (mol->sym_group.xyz[0] == 1) {
            sprintf(buf, "    2 Y   Z");
        }
        else if (mol->sym_group.xyz[1] == 1) {
            sprintf(buf, "    2  ZX");
        }
        else if (mol->sym_group.xyz[2] == 1) {
            sprintf(buf, "    2 Y X");
        }
    }
    else if (mol->sym_group.group == SYMMETRY_C2h) {
        if (mol->sym_group.xyz[0] == 1) {
            sprintf(buf, "    2X  XYZ");
        }
        else if (mol->sym_group.xyz[1] == 1) {
            sprintf(buf, "    2 Y XYZ");
        }
        else if (mol->sym_group.xyz[2] == 1) {
            sprintf(buf, "    2  ZXYZ");
        }
    }
    else if (mol->sym_group.group == SYMMETRY_D2) {
        sprintf(buf, "    2XY  YZ");
    }
    else if (mol->sym_group.group == SYMMETRY_D2h) {
        sprintf(buf, "    3  Z  Y  X");
    }
    else {
        // pass
    }
}


/**
 * generates DIRAC 'inp' file.
 * about dirac input files:
 * https://www.diracprogram.org/doc/release-24/manual/index.html
 */
void gen_dirac_inp_file(FILE *inp_file, molecule_t *mol, basis_lib_t *bas_lib, ecp_lib_t *ecp_lib)
{
    // if the symmetry is not specified explicitly, has_inversion == 0,
    // exception: atoms (always have inversion)
    int has_inversion = molecule_has_inversion_center(mol);
    int has_rpp = has_pseudopotential(ecp_lib);

    fprintf(inp_file, "!\n");
    fprintf(inp_file, "**GENERAL\n");
    fprintf(inp_file, ".PCMOUT\n");
    fprintf(inp_file, "**DIRAC\n");
    fprintf(inp_file, ".WAVE F\n");
    fprintf(inp_file, ".4INDEX\n");
    fprintf(inp_file, "**HAMILTONIAN\n");
    if (has_rpp) {
        fprintf(inp_file, ".ECP\n");
    }
    else {
        fprintf(inp_file, "!.GAUNT\n");
        fprintf(inp_file, ".X2Cmmf\n");
    }
    fprintf(inp_file, "**WAVE FUNCTIONS\n");
    fprintf(inp_file, ".SCF\n");
    fprintf(inp_file, "*SCF\n");
    fprintf(inp_file, ".CLOSED SHELL\n");
    fprintf(inp_file, " %s\n", has_inversion ? "<?> <?>" : "<?>");
    fprintf(inp_file, "**MOLTRA\n");
    fprintf(inp_file, ".PRPTRA\n");
    fprintf(inp_file, ".ACTIVE\n");
    if (has_inversion) {
        fprintf(inp_file, "all\n");
        fprintf(inp_file, "all\n");
    }
    else {
        fprintf(inp_file, "all\n");
    }
    fprintf(inp_file, "*END OF DIRAC\n");
}


/**
 * returns number of L-blocks for given angular momentum L
 */
int basis_get_number_l_blocks(basis_t *bas, int L)
{
    bfn_t *first_basis_functions[MAX_CONTRACTED];   // functions which are first in each block
    int block_sizes[MAX_CONTRACTED];
    int n_blocks = 0;
    int has_primitive_functions = 0;
    int n_primitive_functions = 0;

    memset(first_basis_functions, 0, sizeof(first_basis_functions));
    memset(block_sizes, 0, sizeof(block_sizes));

    for (int ibf = 0; ibf < bas->nfun; ibf++) {
        bfn_t *bf = bas->functions[ibf];
        if (bf->l != L) {
            continue;
        }

        // find block for this function
        int ib;
        if (bf->nprim == 1) {
            has_primitive_functions = 1;
            n_primitive_functions++;
            continue;
        }
        for (ib = 0; ib < n_blocks; ib++) {
            if (bfn_same_exponents(first_basis_functions[ib], bf)) {
                block_sizes[ib]++;
                break;
            }
        }
        if (ib == n_blocks) { // block for this function was not found
            first_basis_functions[n_blocks] = bf;
            block_sizes[n_blocks] = 1;
            n_blocks++;
        }
    }

    int n_blocks_dirac = 0;
    for (int i = 0; i < n_blocks; i++) {
        n_blocks_dirac += block_sizes[i] / 4 + (block_sizes[i] % 4 != 0);
    }
    // DIRAC limitation: max 30 gaussian primitives in the block
    if (has_primitive_functions) {
        n_blocks_dirac += (n_primitive_functions / MAX_PRIMITIVES) +
                          (n_primitive_functions % MAX_PRIMITIVES > 0);
    }

    return n_blocks_dirac;
}


/**
 * prints basis functions with given angular momentum L
 */
void print_l_blocks_dirac(FILE *out, basis_t *bas, int L)
{
    bfn_t *basis_functions[MAX_CONTRACTED][MAX_CONTRACTED];   // functions which are first in each block
    int block_sizes[MAX_CONTRACTED];
    int n_blocks = 0;
    int has_primitive_functions = 0;
    int n_primitive_functions = 0;

    memset(basis_functions, 0, sizeof(basis_functions));
    memset(block_sizes, 0, sizeof(block_sizes));

    for (int ibf = 0; ibf < bas->nfun; ibf++) {
        bfn_t *bf = bas->functions[ibf];
        if (bf->l != L) {
            continue;
        }

        // find block for this function
        int ib;
        if (bf->nprim == 1) {
            has_primitive_functions = 1;
            n_primitive_functions++;
            continue;
        }
        for (ib = 0; ib < n_blocks; ib++) {
            if (bfn_same_exponents(basis_functions[ib][0], bf)) { // add function to the existing block
                int bs = block_sizes[ib];
                basis_functions[ib][bs] = bf;
                block_sizes[ib]++;
                break;
            }
        }
        if (ib == n_blocks) { // block for this function was not found; create new block
            basis_functions[n_blocks][0] = bf;
            block_sizes[n_blocks] = 1;
            n_blocks++;
        }
    }

    if (n_blocks == 0 && has_primitive_functions == 0) {
        return;
        // nothing to print
    }

    fprintf(out, "# %c functions\n", tolower(angular_momentum_to_char(L)));

    for (int ibl = 0; ibl < n_blocks; ibl++) {
        int block_size = block_sizes[ibl];
        bfn_t **funcs = basis_functions[ibl];
        int nprim = funcs[0]->nprim;
        int nprint = 0;

        // max 4 contracted functions in block
        for (int icontr = 0; icontr < block_size; icontr += 4) {
            nprint = (block_size - icontr > 4) ? 4 : (block_size - icontr);
            fprintf(out, "F%3d%3d\n", nprim, nprint);
            for (int i = 0; i < nprim; i++) {
                fprintf(out, "%20.8f", funcs[0]->e[i]);
                for (int j = 0; j < nprint; j++) {
                    fprintf(out, "%14.8f", funcs[icontr + j]->c[i]);
                }
                fprintf(out, "\n");
            }
        }
    }

    // print primitive Gaussians
    if (has_primitive_functions) {
        int n_prim_blocks = (n_primitive_functions / MAX_PRIMITIVES) +
                            (n_primitive_functions % MAX_PRIMITIVES > 0);
        int ibf = 0;
        for (int i = 0; i < n_prim_blocks; i++) {
            int ibf_begin = i * MAX_PRIMITIVES;
            int ibf_end = MIN((i + 1) * MAX_PRIMITIVES, n_primitive_functions);
            int n_prim_i = ibf_end - ibf_begin;
            fprintf(out, "F%3d%3d\n", n_prim_i, 0);

            for (; ibf < bas->nfun; ibf++) {
                bfn_t *bf = bas->functions[ibf];
                if (bf->l != L || bf->nprim != 1) {
                    continue;
                }
                if (ibf_begin == ibf_end) {
                    break;  // go to the next block with primitives
                }
                ibf_begin++;
                fprintf(out, "%20.8f\n", bf->e[0]);
            }
        }
    }
}


/**
 * counts number of L-blocks for each L (in DIRAC format)
 * note that in mol format max length of line == 80 =>
 * only 4 contracted functions can be printed simultaneously
 */
void basis_get_blocks_for_dirac(basis_t *bas, int *n_blocks, int *n_block_sizes, int len_block_sizes)
{
    int n_cont_fun_by_l[MAX_ANG_MOM];
    int n_prim_fun_by_l[MAX_ANG_MOM];

    memset(n_cont_fun_by_l, 0, sizeof(n_cont_fun_by_l));
    memset(n_prim_fun_by_l, 0, sizeof(n_prim_fun_by_l));

    for (int i = 0; i < len_block_sizes; i++) {
        n_block_sizes[i] = 0;
    }

    // count functions in each irrep of SO(3) (L)
    for (int ifun = 0; ifun < bas->nfun; ifun++) {
        bfn_t *bf = bas->functions[ifun];
        int ang_mom = bf->l;
        if (bf->nprim == 1) {
            n_prim_fun_by_l[ang_mom]++;
        }
        else {
            n_cont_fun_by_l[ang_mom]++;
        }
    }

    // for each angulam momentum: max 4 contracted functions in block
    // + 1 separate block for uncontracted Gaussians
    for (int i = 0; i < len_block_sizes; i++) {
        n_block_sizes[i] = n_cont_fun_by_l[i] / 4 + (n_cont_fun_by_l[i] % 4 != 0) + (n_prim_fun_by_l[i] > 0);
    }

    *n_blocks = 0;
    for (int i = 0; i < len_block_sizes; i++) {
        if (n_block_sizes[i] != 0) {
            *n_blocks = i + 1;
        }
    }
}


int has_pseudopotential(ecp_lib_t *ecp_lib)
{
    for (int i = N_CHEM_ELEMENTS - 1; i >= 0; i--) {
        char elem_sym[MAX_ELEMENT_SYMBOL];

        get_element_symbol(i, elem_sym);
        ecp_t *ecp = ecp_lib_get(ecp_lib, elem_sym);
        if (ecp != NULL) {
            return 1;
        }
    }

    return 0;
}

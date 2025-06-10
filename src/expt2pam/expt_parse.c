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

/**
 * parser
 *
 * Format of EXP-T input file for molecules is very similar to NWChem's one:
 * http://www.nwchem-sw.org/index.php/Main_Page
 *
 * Example (for UO2 molecule, adapted from http://www.nwchem-sw.org/index.php/SODFT):
 *
 * geometry units angstrom
 *   U     0.00000      0.00000     0.00000
 *   O     0.00000      0.00000     1.68000
 *   O     0.00000      0.00000    -1.68000
 * end
 *
 * basis
 * U    S
 *       12.12525300         0.02192100
 *        7.16154500        -0.22516000
 *        4.77483600         0.56029900
 *        2.01169300        -1.07120900
 * . . .
 * U    P
 *       17.25477000         0.00139800
 *        7.73535600        -0.03334600
 *        5.15587800         0.11057800
 *        2.24167000        -0.31726800
 * . . .
 *
 * O    S
 *       47.10551800        -0.01440800
 *        5.91134600         0.12956800
 *        0.97648300        -0.56311800
 * . . .
 * O    P
 *       16.69221900         0.04485600
 *        3.90070200         0.22261300
 *        1.07825300         0.50018800
 * . . .
 * end
 *
 * ecp
 * U nelec 78
 * U s
 *  2          4.06365300        112.92010300
 *  2          1.88399500         15.64750000
 *  2          0.88656700         -3.68997100
 * U p
 *  2          3.98618100        118.75801600
 *  2          2.00016000         15.07722800
 *  2          0.96084100          0.55672000
 * U d
 *  2          4.14797200         60.85589200
 *  2          2.23456300         29.28004700
 *  2          0.91369500          4.99802900
 * U f
 *  2          3.99893800         49.92403500
 *  2          1.99884000        -24.67404200
 *  2          0.99564100          1.38948000
 * O nelec 2
 * O s
 *  2         10.44567000         50.77106900
 * O p
 *  2         18.04517400         -4.90355100
 * O d
 *  2          8.16479800         -3.31212400
 * end
 *
 * #spin-orbit part
 * so
 * U p
 *  2    3.986181      1.816350
 *  2    2.000160     11.543940
 *  2    0.960841      0.794644
 * U d
 *  2    4.147972      0.353683
 *  2    2.234563      3.499282
 *  2    0.913695      0.514635
 * U f
 *  2    3.998938      4.744214
 *  2    1.998840     -5.211731
 *  2    0.995641      1.867860
 * end
 *
 */

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "geodesic.h"
#include "ecp.h"
#include "error.h"
#include "expt_lexer.h"
#include "molecule.h"
#include "zmatrix.h"

#define MAX_LINE_LEN   1024
#define MAX_FILE_NAME  1024
#define MAX_TAG_LEN    64
#define MAX_PRIMITIVES 100
#define MAX_CONTRACTED 100

void directive_geometry(molecule_t *mol);
void directive_cage(molecule_t *mol);
void directive_zmatrix(molecule_t *mol);
void directive_basis(basis_lib_t *basis);
void directive_basis_block(basis_lib_t *basis);
void directive_ecp(ecp_lib_t *ecp_lib, int pot_type);

void yyerror(char *s);
int next_token();
void put_back(int token_type);
int match(int required_type);
void str_tolower(char *s);
int parse_state_spec(char *buf, char *rep_name, int *state_no);
void set_input_file_name(char *file_name);
void parse_sector_label(char *s, int *h, int *p);
int angular_momentum_to_int(char *lstr);


/**
 * main subroutine of the parser
 */
void expt_parse(char *file_name, molecule_t *mol, basis_lib_t *basis, ecp_lib_t *ecp_lib)
{
    int token_type;

    yyin = fopen(file_name, "r");
    if (yyin == NULL) {
        errquit("Input file \"%s\" not found", file_name);
    }

    set_input_file_name(file_name);

    token_type = next_token();
    while (token_type != END_OF_FILE) {
        switch (token_type) {
            case KEYWORD_GEOMETRY:
                directive_geometry(mol);
                break;
            case KEYWORD_BASIS:
                directive_basis(basis);
                break;
            case KEYWORD_ECP:
                directive_ecp(ecp_lib, ECP_AVERAGED);
                break;
            case KEYWORD_SO:
                directive_ecp(ecp_lib, ECP_SPIN_ORBIT);
                break;
            case END_OF_LINE:
                // nothing to do
                break;
            default:
                yyerror("unknown keyword");
                break;
        }

        /* each "non-newline" directive must end with EOL or EOF */
        if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
            token_type = next_token();
            if (token_type != END_OF_LINE && token_type != END_OF_FILE) {
                yyerror("end of line is expected");
            }
        }

        /* go to the next directive */
        token_type = next_token();
    }

    fclose(yyin);
}


void directive_geometry(molecule_t *mol)
{
    static char *err_no_newline = "expected end of line";
    static char *err_no_elem = "wrong element symbol";
    static char *err_no_float = "float number is expected";
    static char *err_no_sym = "Schoenflies symbol is expected";
    static char *err_wrong_sym = "unknown symmetry group";
    static char *err_eol_or_xyz = "end of line or orientation specification are expected";
    static char *err_wrong_xyz = "wrong orientation specification";
    static char *err_xyz_used = "orientation specification is unnecessary";
    static char *err_wrong_units = "wrong units (allowed: au, atomic, bohr, angstrom)";
    int token_type;
    double factor = 1.0;

    token_type = next_token();

    if (token_type == TT_WORD && strcmp(yytext, "units") == 0) {
        token_type = next_token();
        if (token_type == TT_WORD && strcmp(yytext, "au") == 0) {
            factor = 1.0;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "atomic") == 0) {
            factor = 1.0;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "bohr") == 0) {
            factor = 1.0;
        }
        else if (token_type == TT_WORD && strcmp(yytext, "angstrom") == 0) {
            factor = 1.0 / 0.529177210903;
        }
        else {
            yyerror(err_wrong_units);
        }
        token_type = next_token();
    }
    else if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    while (1) {
        token_type = next_token();
        str_tolower(yytext);

        // symmetry specification
        if (token_type == TT_WORD && strcmp(yytext, "symmetry") == 0) {
            token_type = next_token();
            if (token_type != TT_WORD) {
                yyerror(err_no_sym);
            }
            str_tolower(yytext);
            if (strcmp(yytext, "c1") == 0) {
                mol->sym_group.group = SYMMETRY_C1;
            }
            else if (strcmp(yytext, "c2") == 0) {
                mol->sym_group.group = SYMMETRY_C2;
                mol->sym_group.xyz[2] = 1;
            }
            else if (strcmp(yytext, "cs") == 0) {
                mol->sym_group.group = SYMMETRY_Cs;
                mol->sym_group.xyz[0] = 1;
                mol->sym_group.xyz[1] = 1;
            }
            else if (strcmp(yytext, "ci") == 0) {
                mol->sym_group.group = SYMMETRY_Ci;
            }
            else if (strcmp(yytext, "c2v") == 0) {
                mol->sym_group.group = SYMMETRY_C2v;
                mol->sym_group.xyz[2] = 1;
            }
            else if (strcmp(yytext, "c2h") == 0) {
                mol->sym_group.group = SYMMETRY_C2h;
                mol->sym_group.xyz[2] = 1;
            }
            else if (strcmp(yytext, "d2") == 0) {
                mol->sym_group.group = SYMMETRY_D2;
            }
            else if (strcmp(yytext, "d2h") == 0) {
                mol->sym_group.group = SYMMETRY_D2h;
            }
            else if (strcmp(yytext, "cinfv") == 0) {
                mol->sym_group.group = SYMMETRY_Cinfv;
                mol->sym_group.xyz[2] = 1;
            }
            else if (strcmp(yytext, "dinfh") == 0) {
                mol->sym_group.group = SYMMETRY_Dinfh;
                mol->sym_group.xyz[2] = 1;
            }
            else {
                yyerror(err_wrong_sym);
            }

            token_type = next_token();
            if (token_type == END_OF_LINE) {  // default orientation will be used
                continue;
            }
            else if (token_type == TT_WORD) {
                str_tolower(yytext);
                if (mol->sym_group.group == SYMMETRY_C2 ||
                    mol->sym_group.group == SYMMETRY_C2v ||
                    mol->sym_group.group == SYMMETRY_C2h ||
                    mol->sym_group.group == SYMMETRY_Cinfv ||
                    mol->sym_group.group == SYMMETRY_Dinfh) {
                    if (strcmp(yytext, "x") == 0) {
                        mol->sym_group.xyz[0] = 1;
                        mol->sym_group.xyz[1] = 0;
                        mol->sym_group.xyz[2] = 0;
                    }
                    else if (strcmp(yytext, "y") == 0) {
                        mol->sym_group.xyz[0] = 0;
                        mol->sym_group.xyz[1] = 1;
                        mol->sym_group.xyz[2] = 0;
                    }
                    else if (strcmp(yytext, "z") == 0) {
                        mol->sym_group.xyz[0] = 0;
                        mol->sym_group.xyz[1] = 0;
                        mol->sym_group.xyz[2] = 1;
                    }
                    else {
                        yyerror(err_wrong_xyz);
                    }
                }
                else if (mol->sym_group.group == SYMMETRY_Cs) {
                    if (strcmp(yytext, "xy") == 0 || strcmp(yytext, "yx") == 0) {
                        mol->sym_group.xyz[0] = 1;
                        mol->sym_group.xyz[1] = 1;
                        mol->sym_group.xyz[2] = 0;
                    }
                    else if (strcmp(yytext, "yz") == 0 || strcmp(yytext, "zy") == 0) {
                        mol->sym_group.xyz[0] = 0;
                        mol->sym_group.xyz[1] = 1;
                        mol->sym_group.xyz[2] = 1;
                    }
                    else if (strcmp(yytext, "xz") == 0 || strcmp(yytext, "zx") == 0) {
                        mol->sym_group.xyz[0] = 1;
                        mol->sym_group.xyz[1] = 0;
                        mol->sym_group.xyz[2] = 1;
                    }
                    else {
                        yyerror(err_wrong_xyz);
                    }
                }
                else {
                    yyerror(err_xyz_used);
                }
            }
            else {
                yyerror(err_eol_or_xyz);
            }

            token_type = next_token();
            if (token_type != END_OF_LINE) {  // default orientation will be used
                yyerror(err_no_newline);
            }
        }
            // z-matrix specification
        else if (token_type == TT_WORD &&
                 (strcmp(yytext, "zmatrix") == 0 ||
                  strcmp(yytext, "zmat") == 0 ||
                  strcmp(yytext, "zmt") == 0)) {
            directive_zmatrix(mol);
        }
            // cage specification
        else if (token_type == TT_WORD && strcmp(yytext, "cage") == 0) {
            directive_cage(mol);
        }
            // xyz geometry
        else if (token_type == TT_WORD) { // add new atom
            int nuc_charge;
            double r[] = {0.0, 0.0, 0.0};

            // read element symbol -> nuclear charge
            nuc_charge = get_element_nuc_charge(yytext);
            if (nuc_charge == -1) {
                yyerror(err_no_elem);
            }

            // read nuclear coordinates x,y,z
            for (int i = 0; i < 3; i++) {
                token_type = next_token();
                if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                    yyerror(err_no_float);
                }
                r[i] = atof(yytext);
            }

            // optional keyword: 'charge'
            token_type = next_token();
            int do_add_atom = 1;
            if (token_type == TT_WORD && strcmp(yytext, "charge") == 0) {
                token_type = next_token();
                if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                    yyerror("float number (charge) is expected");
                }
                double q = atof(yytext);
                if (nuc_charge == 0) {
                    do_add_atom = 0;
                }

                molecule_add_point_charge(mol, q, factor * r[0], factor * r[1], factor * r[2]);
            }
            else {
                put_back(token_type);
            }

            // check for the end of line
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }

            if (do_add_atom) {
                molecule_add_atom(mol, nuc_charge, factor * r[0], factor * r[1], factor * r[2]);
            }
        }
        else if (token_type == END_OF_LINE) {
            /* empty line, nothing to do */
            continue;
        }
        else if (token_type == KEYWORD_END) {
            return;
        }
        else {
            yyerror("unexpected token");
        }
    }
}


enum {
    TYPE_INT = 0,
    TYPE_DOUBLE
};

typedef struct {
    char *name;
    int type;
    union {
        int ival;
        double dval;
    };
} var_t;

typedef struct {
    var_t *vars;
    size_t len;
    size_t cap;
} var_list_t;


var_list_t *var_list_new()
{
    var_list_t *li = (var_list_t *) malloc(1 * sizeof(var_list_t));
    if (li != NULL) {
        li->len = 0;
        li->cap = 1;
        li->vars = (var_t *) malloc(sizeof(var_t) * li->cap);
    }
    return li;
}


void var_list_delete(var_list_t *li)
{
    for (size_t i = 0; i < li->len; i++) {
        if (li->vars[i].name != NULL) {
            free(li->vars[i].name);
        }
    }
    free(li->vars);
    free(li);
}


void var_list_print(var_list_t *li)
{
    printf("List len=%d cap=%d\n", li->len, li->cap);
    for (size_t i = 0; i < li->len; i++) {
        printf(" [%d] %-10s %s ", i, li->vars[i].name, (li->vars[i].type == TYPE_INT) ? "int" : "dbl");
        if (li->vars[i].type == TYPE_INT) {
            printf("%d\n", li->vars[i].ival);
        }
        else {
            printf("%f\n", li->vars[i].dval);
        }
    }
}


void var_list_grow(var_list_t *li)
{
    li->cap *= 2;
    li->vars = realloc(li->vars, sizeof(var_t) * li->cap);
}


void var_list_set(var_list_t *li, char *key, int type, double val)
{
    // try to find the key
    for (size_t i = 0; i < li->len; i++) {
        if (strcmp(li->vars[i].name, key) == 0) {
            if (type == TYPE_INT) {
                li->vars[i].ival = (int) val;
                li->vars[i].type = TYPE_INT;
            }
            else { // TYPE_DOUBLE
                li->vars[i].dval = val;
                li->vars[i].type = TYPE_DOUBLE;
            }
            return;
        }
    }

    // add the (key,val) pair
    if (li->len + 1 == li->cap) {
        var_list_grow(li);
    }

    li->vars[li->len].name = (char *) malloc(sizeof(char) * strlen(key));
    strcpy(li->vars[li->len].name, key);
    if (type == TYPE_INT) {
        li->vars[li->len].ival = (int) val;
        li->vars[li->len].type = TYPE_INT;
    }
    else { // TYPE_DOUBLE
        li->vars[li->len].dval = val;
        li->vars[li->len].type = TYPE_DOUBLE;
    }
    li->len++;
}


int var_list_get(var_list_t *li, char *key, int *type, double *val)
{
    assert(li != NULL);
    assert(key != NULL);
    assert(type != NULL);
    assert(val != NULL);

    var_list_print(li);

    // try to find the key
    for (size_t i = 0; i < li->len; i++) {
        if (strcmp(li->vars[i].name, key) == 0) {
            if (li->vars[i].type == TYPE_INT) {
                *val = (double) li->vars[i].ival;
                *type = li->vars[i].type;
            }
            else { // TYPE_DOUBLE
                *val = li->vars[i].dval;
                *type = li->vars[i].type;
            }
            return 1;
        }
    }

    return 0;
}


int arr_to_pos_integer(char *s)
{
    for (size_t i = 0; i < strlen(s); i++) {
        if (!isdigit(s[i])) {
            return -1;
        }
    }
    return atoi(s);
}


int tag_to_element_nuc_charge(char *s)
{
    char *buf;
    size_t i, j;
    int nuc_charge;

    // if tag contains only one integer == nuc_charge
    nuc_charge = arr_to_pos_integer(s);
    if (nuc_charge >= 0) {
        return nuc_charge;
    }

    // else: extract only first substring of alphabetic symbols
    buf = (char *) malloc(sizeof(char) * strlen(s));
    for (i = 0, j = 0; i < strlen(s); i++) {
        if (isalpha(s[i])) {
            buf[j++] = s[i];
        }
    }
    buf[j] = '\0';

    nuc_charge = get_element_nuc_charge(buf);
    free(buf);

    return nuc_charge;
}


void directive_zmatrix(molecule_t *mol)
{
    static char *err_no_newline = "end of line is expected";
    static char *err_end_of_file = "reached end of file";
    static char *err_no_elem = "wrong element symbol";
    static char *err_no_elem_tag = "wrong (unknown) element tag";
    static char *err_wrong_end = "unexpected end of line";

    int atomlist[MAX_N_ATOMS];
    int rconnect[MAX_N_ATOMS];
    int aconnect[MAX_N_ATOMS];
    int dconnect[MAX_N_ATOMS];
    double rlist[MAX_N_ATOMS];
    double alist[MAX_N_ATOMS];
    double dlist[MAX_N_ATOMS];

    memset(atomlist, 0, sizeof(atomlist));
    memset(rconnect, 0, sizeof(rconnect));
    memset(aconnect, 0, sizeof(aconnect));
    memset(dconnect, 0, sizeof(dconnect));
    memset(rlist, 0, sizeof(rlist));
    memset(alist, 0, sizeof(alist));
    memset(dlist, 0, sizeof(dlist));

    if(!match(END_OF_LINE)) {
        yyerror(err_no_newline);
    }

    int count_atoms = 0;

    int token_type = next_token();
    while (token_type != KEYWORD_END) {

        if (token_type == TT_WORD) {

            // read element symbol -> nuclear charge
            int nuc_charge = get_element_nuc_charge(yytext);
            if (nuc_charge == -1) {
                yyerror(err_no_elem);
            }
            atomlist[count_atoms] = nuc_charge;

            for (int icoord = 0; icoord < 3; icoord++) {
                int atom_number = 0;
                double coord_value = 0.0;

                if (icoord > count_atoms - 1) {
                    break;
                }

                token_type = next_token();
                if (token_type == TT_INTEGER) {
                    atom_number = atoi(yytext);
                }
                else {
                    yyerror("integer positive number is expected");
                }

                token_type = next_token();
                if (token_type == TT_INTEGER || token_type == TT_FLOAT) {
                    coord_value = atof(yytext);
                }
                else {
                    yyerror("integer positive number is expected");
                }

                if (icoord == 0) {
                    rconnect[count_atoms] = atom_number;
                    rlist[count_atoms] = coord_value;
                }
                else if (icoord == 1) {
                    aconnect[count_atoms] = atom_number;
                    alist[count_atoms] = coord_value;
                }
                else {
                    dconnect[count_atoms] = atom_number;
                    dlist[count_atoms] = coord_value;
                }
            }

            //printf("%4d  %4d%14.8f  %4d%14.8f  %4d%14.8f\n", nuc_charge, rconnect[count_atoms], rlist[count_atoms],
            //       aconnect[count_atoms], alist[count_atoms], dconnect[count_atoms], dlist[count_atoms]);

            if(!match(END_OF_LINE)) {
                yyerror(err_no_newline);
            }

            count_atoms += 1;
        }
        else if (token_type == END_OF_LINE) {
            // nothing to do; skip empty line
        }
        else if (token_type == END_OF_FILE) {
            yyerror(err_end_of_file);
        }

        token_type = next_token();
    }

    /*
     * check numbers of atoms
     */
    for (int iatom = 0; iatom < count_atoms; iatom++) {
        if (iatom >= 1 && (rconnect[iatom] <= 0 || rconnect[iatom] > iatom)) {
            errquit("in zmatrix: wrong connectivity (atom number must be in range [1,%d])", iatom);
        }
        if (iatom >= 2 && (aconnect[iatom] <= 0 || aconnect[iatom] > iatom)) {
            errquit("in zmatrix: wrong connectivity (atom number must be in range [1,%d])", iatom);
        }
        if (iatom >= 3 && (dconnect[iatom] <= 0 || dconnect[iatom] > iatom)) {
            errquit("in zmatrix: wrong connectivity (atom number must be in range [1,%d])", iatom);
        }
    }

    /*
     * transform coordinates: zmatrix -> cartesian
     */
    convert_zmatrix_to_xyz(count_atoms, atomlist, rconnect, rlist, aconnect, alist, dconnect, dlist, mol);
}


/**
 * constructs cage of positively charged points around the molecule.
 * geodesic polyhedra are used as a cage.
 * following options are allowed:
 *   icosahedron + level of tesselation
 *   tetrahedron + level of tesselation
 *   octahedron + level of tesselation
 *   dodecahedron (no tesselation)
 * sphere can be scaled to the given radius.
 *
 * syntax:
 * cage
 *   origin <x> <y> <z>
 *   charge <z>
 *   radius <r>
 *   [icosahedron|tetrahedron|octahedron|dodecahedron] <int level (optional)>
 * end
 */
void directive_cage(molecule_t *mol)
{
    static char *err_no_newline = "end of line is expected";
    static char *err_end_of_file = "reached end of file";
    static char *err_wrong_end = "unexpected end of line";
    static char *err_no_float = "float number is expected";
    static char *err_positive_f = "positive float number is expected";
    static char *err_positive_i = "positive integer number is expected";

    int token_type;
    int polyhedron_type = CAGE_ICOSAHEDRON;
    int tessel_level = 1;
    double radius = 1.0;
    double charge = 0.0;
    double origin[] = {0.0, 0.0, 0.0};

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    while (1) {
        token_type = next_token();
        str_tolower(yytext);

        // symmetry specification
        if (token_type == TT_WORD && strcmp(yytext, "origin") == 0) {
            // read origin coordinates x,y,z
            for (int i = 0; i < 3; i++) {
                token_type = next_token();
                if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                    yyerror(err_no_float);
                }
                origin[i] = atof(yytext);
            }
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == TT_WORD && strcmp(yytext, "radius") == 0) {
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror(err_no_float);
            }
            radius = atof(yytext);
            if (radius <= 0) {
                yyerror(err_positive_f);
            }
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == TT_WORD && strcmp(yytext, "charge") == 0) {
            token_type = next_token();
            if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                yyerror(err_no_float);
            }
            charge = atof(yytext);
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == TT_WORD && strcmp(yytext, "icosahedron") == 0) {
            polyhedron_type = CAGE_ICOSAHEDRON;
            token_type = next_token();
            if (token_type == TT_INTEGER) {
                tessel_level = atoi(yytext);
                token_type = next_token();
            }
            if (tessel_level <= 0) {
                yyerror(err_positive_i);
            }
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == TT_WORD && strcmp(yytext, "tetrahedron") == 0) {
            polyhedron_type = CAGE_TETRAHEDRON;
            token_type = next_token();
            if (token_type == TT_INTEGER) {
                tessel_level = atoi(yytext);
                token_type = next_token();
            }
            if (tessel_level <= 0) {
                yyerror(err_positive_i);
            }
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == TT_WORD && strcmp(yytext, "octahedron") == 0) {
            polyhedron_type = CAGE_OCTAHEDRON;
            token_type = next_token();
            if (token_type == TT_INTEGER) {
                tessel_level = atoi(yytext);
                token_type = next_token();
            }
            if (tessel_level <= 0) {
                yyerror(err_positive_i);
            }
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == TT_WORD && strcmp(yytext, "dodecahedron") == 0) {
            polyhedron_type = CAGE_DODECAHEDRON;
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
        }
        else if (token_type == END_OF_LINE) {
            /* empty line, nothing to do */
            continue;
        }
        else if (token_type == KEYWORD_END) {
            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_no_newline);
            }
            break;
        }
        else if (token_type == END_OF_FILE) {
            yyerror(err_end_of_file);
        }
        else {
            yyerror("unexpected token");
        }
    }

    // generate polyhedron and transform it wrt origin and radius

    if (polyhedron_type == CAGE_ICOSAHEDRON ||
        polyhedron_type == CAGE_TETRAHEDRON ||
        polyhedron_type == CAGE_OCTAHEDRON) {
        geodesicSphere poly;
        if (polyhedron_type == CAGE_ICOSAHEDRON) {
            poly = icosahedronSphere(tessel_level);
        }
        else if (polyhedron_type == CAGE_TETRAHEDRON) {
            poly = tetrahedronSphere(tessel_level);
        }
        else {   // CAGE_OCTAHEDRON
            poly = octahedronSphere(tessel_level);
        }
        mol->cage.total_charge = charge;
        mol->cage.n_points = poly.numPoints;
        mol->cage.points = (double *) malloc(3 * sizeof(double) * poly.numPoints);
        for (int i = 0; i < poly.numPoints; i++) {
            for (int j = 0; j < 3; j++) {
                mol->cage.points[3 * i + j] = origin[j] + radius * poly.points[3 * i + j];
            }
        }
        deleteGeodesicSphere(&poly);
    }
    else {  // DODECAHEDRON
        mol->cage.total_charge = charge;
        mol->cage.n_points = DODECAHEDRON_POINT_COUNT;
        mol->cage.points = (double *) malloc(3 * sizeof(double) * DODECAHEDRON_POINT_COUNT);
        for (int i = 0; i < DODECAHEDRON_POINT_COUNT; i++) {
            for (int j = 0; j < 3; j++) {
                mol->cage.points[3 * i + j] = origin[j] + radius * _dodecahedron_points[3 * i + j];
            }
        }
    }
}


void directive_ecp(ecp_lib_t *ecp_lib, int pot_type)
{
    static char *err_no_newline = "end of line is expected";
    static char *err_ang_mom = "angular momentum symbol (one of SPDFGHIKLM) or keyword 'nelec' is expected";
    static char *err_end_line = "end of line is expected";
    static char *err_end_float = "float number or end of line are expected";
    static char *err_float = "float number is expected";
    static char *err_integer = "non-negative integer number is expected";
    static char *err_no_fun = "no functions found in this block";
    static char *err_no_elem = "wrong element symbol";
    int token_type;

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    // read functions (ECP expansions for each L) until the 'end' keyword

    while (1) {
        int ang_mom, nelec, read_expansion = 0;
        int n_primitives = 0;
        int f_n[MAX_PRIMITIVES];
        double f_e[MAX_PRIMITIVES];
        double f_c[MAX_PRIMITIVES];
        token_type = next_token();

        if (token_type == TT_WORD) {

            // get ECP for this element
            str_tolower(yytext);
            ecp_t *ecp = ecp_lib_get(ecp_lib, yytext);
            if (ecp == NULL) {
                int z = get_element_nuc_charge(yytext);
                if (z == -1) {
                    yyerror(err_no_elem);
                }
                ecp = ecp_new(z);
            }

            // get angular momentum or number of core electrons
            token_type = next_token();
            if (token_type != TT_WORD) {
                yyerror(err_ang_mom);
            }
            str_tolower(yytext);
            if (strcmp(yytext, "nelec") == 0) {  // number of core electrons
                token_type = next_token();
                if (token_type != TT_INTEGER) {
                    yyerror(err_integer);
                }
                nelec = atoi(yytext);
                if (nelec < 0) {
                    yyerror(err_integer);
                }
                ecp->n_core_elec = nelec;
            }
                // angular momentum: U_L or S, P, D, ...
            else if (strcmp(yytext, "ul") == 0) {
                ang_mom = ECP_UL;
                read_expansion = 1;
            }
            else {
                ang_mom = angular_momentum_to_int(yytext);
                if (ang_mom == -1) {
                    yyerror(err_ang_mom);
                }
                read_expansion = 1;
            }

            token_type = next_token();
            if (token_type != END_OF_LINE) {
                yyerror(err_end_line);
            }

            // read Gaussian expansion if required
            if (!read_expansion) {
                ecp_lib_set(ecp_lib, ecp);
                continue;
            }

            n_primitives = 0;
            token_type = next_token();
            while (token_type == TT_INTEGER) {
                int n = atoi(yytext);
                token_type = next_token();
                if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                    yyerror(err_float);
                }
                double exponent = atof(yytext);
                token_type = next_token();
                if (token_type != TT_INTEGER && token_type != TT_FLOAT) {
                    yyerror(err_float);
                }
                double coef = atof(yytext);

                token_type = next_token();
                if (token_type != END_OF_LINE) {
                    yyerror(err_end_line);
                }

                f_n[n_primitives] = n;
                f_e[n_primitives] = exponent;
                f_c[n_primitives] = coef;
                n_primitives++;

                token_type = next_token();
            }
            put_back(token_type);

            // add expansion to ECP
            ecp_add_function(ecp, pot_type, ang_mom, n_primitives, f_n, f_e, f_c);
            ecp_lib_set(ecp_lib, ecp);
        }
        else if (token_type == END_OF_LINE) {
            /* empty line, nothing to do */
            continue;
        }
        else if (token_type == KEYWORD_END) {
            return;
        }
        else {
            yyerror("unexpected token");
        }
    }
}


void directive_basis(basis_lib_t *basis)
{
    static char *err_no_newline = "expected end of line";
    int token_type;

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_no_newline);
    }

    // read blocks of functions until the 'end' keyword

    while (1) {
        token_type = next_token();

        if (token_type == TT_WORD) {
            put_back(token_type);
            directive_basis_block(basis);
        }
        else if (token_type == END_OF_LINE) {
            /* empty line, nothing to do */
            continue;
        }
        else if (token_type == KEYWORD_END) {
            return;
        }
        else {
            yyerror("unexpected token");
        }
    }
}


void directive_basis_block(basis_lib_t *basis_lib)
{
    static char *err_ang_mom = "angular momentum symbol (one of SPDFGHIKLM) is expected";
    static char *err_end_line = "end of line is expected";
    static char *err_end_float = "float number or end of line are expected";
    static char *err_float = "float number is expected";
    static char *err_no_fun = "no functions found in this block";
    static char *err_no_elem = "wrong element symbol";

    int token_type;
    static double matrix[MAX_CONTRACTED][MAX_PRIMITIVES];
    int n_primitives = 0;
    int n_contracted = 0;
    basis_t *basis = NULL;

    memset(matrix, 0, sizeof(matrix));

    token_type = next_token();
    str_tolower(yytext);
    basis = basis_lib_get(basis_lib, yytext);
    if (basis == NULL) {
        int z = get_element_nuc_charge(yytext);
        if (z == -1) {
            yyerror(err_no_elem);
        }
        basis = basis_new(z);
    }

    token_type = next_token();
    if (token_type != TT_WORD) {
        yyerror(err_ang_mom);
    }
    int ang_mom = angular_momentum_to_int(yytext);
    if (ang_mom == -1) {
        yyerror(err_ang_mom);
    }

    token_type = next_token();
    if (token_type != END_OF_LINE) {
        yyerror(err_end_line);
    }

    while (1) { // loop over rows
        token_type = next_token();
        if (token_type == END_OF_LINE) { // just empty line
            continue;
        }
        else if (token_type == TT_WORD || token_type == KEYWORD_END) { // end of block
            put_back(token_type);
            break;
        }
        else if (token_type == TT_INTEGER || token_type == TT_FLOAT) {
            // loop over columns
            double exponent = atof(yytext);
            matrix[0][n_primitives] = exponent;
            n_contracted = 0;
            while ((token_type = next_token()) == TT_INTEGER || token_type == TT_FLOAT) {
                double coef = atof(yytext);
                matrix[1 + n_contracted][n_primitives] = coef;
                n_contracted++;
            }
            if (token_type != END_OF_LINE) {
                yyerror(err_end_float);
            }
            n_primitives++;
        }
        else {
            yyerror(err_float);
        }
    }

    if (n_primitives == 0) {
        yyerror(err_no_fun);
    }

    if (n_contracted == 0) { // special case of uncontracted functions
        for (int i = 0; i < n_primitives; i++) {
            double coef_1 = 1.0;
            basis_add_function(basis, ang_mom, 1, &matrix[0][i], &coef_1);
        }
    }
    else {
        for (int i = 0; i < n_contracted; i++) {
            basis_add_function(basis, ang_mom, n_primitives, matrix[0], matrix[1 + i]);
        }
    }

    basis_lib_set(basis_lib, basis);
}


/**
 * Converts string containing angular momentum symbol (S,P,D,...) to its integer
 * value: S -> 0, P -> 1, ...
 * returns -1 in case of error.
 */
int angular_momentum_to_int(char *lstr)
{
    if (strlen(lstr) != 1) {
        return -1;
    }

    switch (lstr[0]) {
        case 'S':
        case 's':
            return 0;
        case 'P':
        case 'p':
            return 1;
        case 'D':
        case 'd':
            return 2;
        case 'F':
        case 'f':
            return 3;
        case 'G':
        case 'g':
            return 4;
        case 'H':
        case 'h':
            return 5;
        case 'I':
        case 'i':
            return 6;
        case 'K':
        case 'k':
            return 7;
        case 'L':
        case 'l':
            return 8;
        case 'M':
        case 'm':
            return 9;
        default:
            return -1;
    }
}


/*******************************************************************************
                     PARSER OF INPUT FILES - GLOBAL VARIABLES
 ******************************************************************************/

// flag: was the last token returned to the lexer?
static int pushed_back = 0;

// name of the input file with language to be processed
static char input_file_name[MAX_FILE_NAME];


/*******************************************************************************
                     PARSER OF INPUT FILES - HELPER FUNCTIONS
 ******************************************************************************/


void set_input_file_name(char *file_name)
{
    strcpy(input_file_name, file_name);
}


void print_line_by_num(char *file_name, int required_lineno)
{
    FILE *f;
    int c;
    int lineno;

    f = fopen(file_name, "r");
    if (f == NULL) {
        errquit("print_line_by_num(): cannot open file '%s'", file_name);
    }

    lineno = 1;
    while ((c = fgetc(f)) != EOF) {
        if (lineno == required_lineno) {
            printf(" ");
            do {
                putchar(c);
                c = fgetc(f);
            } while (c != '\n' && c != EOF);
            printf("\n");
            fclose(f);
            return;
        }
        if (c == '\n') {
            lineno++;
        }
    }

    errquit("print_line_by_num(): no line number %d in '%s'", required_lineno, file_name);
    fclose(f);
}


void yyerror(char *s)
{
    int lineno = yylineno;

    // lineno = yylineno-1 if last token was END_OF_LINE
    if (strcmp(yytext, "\n") == 0) {
        lineno = yylineno - 1;
    }
    fclose(yyin);

    printf("%s:%d:%d: error: %s\n\n", input_file_name, lineno, yycol, s);
    print_line_by_num(input_file_name, lineno);

    // print pointer to the erroneous token
    printf(" ");
    for (int i = 0; i < yycol; i++) {
        printf(" ");
    }
    printf("^\n\n");
    printf("execution terminated.\n");

    exit(1);
}


int next_token()
{
    if (pushed_back) {
        int token_type = pushed_back;
        pushed_back = 0;
        return token_type;
    }
    return yylex();
}


void put_back(int token_type)
{
    pushed_back = token_type;
}


int match(int required_type)
{
    int token_type = next_token();
    if (token_type != required_type) {
        return 0;
    }
    return 1;
}


void str_tolower(char *s)
{
    for (; *s; s++) {
        *s = tolower(*s);
    }
}


void str_toupper(char *s)
{
    for (; *s; s++) {
        *s = toupper(*s);
    }
}

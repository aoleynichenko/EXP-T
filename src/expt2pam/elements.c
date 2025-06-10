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

#include "elements.h"

#include <assert.h>
#include <ctype.h>
#include <string.h>

element_t periodic_table[] = {
        {0,   "Gh"},
        {1,   "H"},
        {2,   "He"},
        {3,   "Li"},
        {4,   "Be"},
        {5,   "B"},
        {6,   "C"},
        {7,   "N"},
        {8,   "O"},
        {9,   "F"},
        {10,  "Ne"},
        {11,  "Na"},
        {12,  "Mg"},
        {13,  "Al"},
        {14,  "Si"},
        {15,  "P"},
        {16,  "S"},
        {17,  "Cl"},
        {18,  "Ar"},
        {19,  "K"},
        {20,  "Ca"},
        {21,  "Sc"},
        {22,  "Ti"},
        {23,  "V"},
        {24,  "Cr"},
        {25,  "Mn"},
        {26,  "Fe"},
        {27,  "Co"},
        {28,  "Ni"},
        {29,  "Cu"},
        {30,  "Zn"},
        {31,  "Ga"},
        {32,  "Ge"},
        {33,  "As"},
        {34,  "Se"},
        {35,  "Br"},
        {36,  "Kr"},
        {37,  "Rb"},
        {38,  "Sr"},
        {39,  "Y"},
        {40,  "Zr"},
        {41,  "Nb"},
        {42,  "Mo"},
        {43,  "Tc"},
        {44,  "Ru"},
        {45,  "Rh"},
        {46,  "Pd"},
        {47,  "Ag"},
        {48,  "Cd"},
        {49,  "In"},
        {50,  "Sn"},
        {51,  "Sb"},
        {52,  "Te"},
        {53,  "I"},
        {54,  "Xe"},
        {55,  "Cs"},
        {56,  "Ba"},
        {57,  "La"},
        {58,  "Ce"},
        {59,  "Pr"},
        {60,  "Nd"},
        {61,  "Pm"},
        {62,  "Sm"},
        {63,  "Eu"},
        {64,  "Gd"},
        {65,  "Tb"},
        {66,  "Dy"},
        {67,  "Er"},
        {68,  "Ho"},
        {69,  "Tm"},
        {70,  "Yb"},
        {71,  "Lu"},
        {72,  "Hf"},
        {73,  "Ta"},
        {74,  "W"},
        {75,  "Re"},
        {76,  "Os"},
        {77,  "Ir"},
        {78,  "Pt"},
        {79,  "Au"},
        {80,  "Hg"},
        {81,  "Tl"},
        {82,  "Pb"},
        {83,  "Bi"},
        {84,  "Po"},
        {85,  "At"},
        {86,  "Rn"},
        {87,  "Fr"},
        {88,  "Ra"},
        {89,  "Ac"},
        {90,  "Th"},
        {91,  "Pa"},
        {92,  "U"},
        {93,  "Np"},
        {94,  "Pu"},
        {95,  "Am"},
        {96,  "Cm"},
        {97,  "Bk"},
        {98,  "Cf"},
        {99,  "Es"},
        {100, "Fm"},
        {101, "Md"},
        {102, "Nb"},
        {103, "Lr"},
        {104, "Rf"},
        {105, "Db"},
        {106, "Sg"},
        {107, "Bh"},
        {108, "Hs"},
        {109, "Mt"},
        {110, "Ds"},
        {111, "Rg"},
        {112, "Cn"},
        {113, "Nh"},
        {114, "Fl"},
        {115, "Mc"},
        {116, "Lv"},
        {117, "Ts"},
        {118, "Og"},
        {119, "E119"},
        {120, "E120"},
        {121, "E121"},
        {122, "E122"},
        {123, "E123"},
        {124, "E124"},
        {125, "E125"},
        {126, "E126"},
        {127, "E127"},
        {128, "E128"},
        {129, "E129"},
        {130, "E130"}
};


void get_element_symbol(int z, char *sym)
{
    assert(z >= 0 && z <= N_CHEM_ELEMENTS);

    strcpy(sym, periodic_table[z].sym);
}


int get_element_nuc_charge(char *sym)
{
    const int n_elements = sizeof(periodic_table) / sizeof(element_t);
    char buf[MAX_ELEMENT_SYMBOL];

    // cast 'sym' to the 'standard form' (only the first letter is capitalized)
    strncpy(buf, sym, MAX_ELEMENT_SYMBOL);
    buf[MAX_ELEMENT_SYMBOL - 1] = '\0';
    buf[0] = toupper(buf[0]);
    for (int i = 1; i < MAX_ELEMENT_SYMBOL - 1; i++) {
        buf[i] = tolower(buf[i]);
    }

    for (int i = 0; i < n_elements; i++) {
        if (strcmp(periodic_table[i].sym, buf) == 0) {
            return i;
        }
    }

    return -1;
}

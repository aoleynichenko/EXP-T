/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2022 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

%{

#include "lexer.h"

#define MAX_TOKEN_LEN 1024

void inc_col();
int yycol = 0;
int prev_len = -1;

%}

%option case-insensitive

%%
"#"[^\n]*       { inc_col(); } /* single-line comments beginning with '#' */
[ \t\r]         { inc_col(); } /* ignore whitespace */
title           { inc_col(); return KEYWORD_TITLE;        }
print           { inc_col(); return KEYWORD_PRINT;        }
scratch_dir     { inc_col(); return KEYWORD_SCRATCH_DIR;  }
degen_thresh    { inc_col(); return KEYWORD_DEGEN_THRESH; }
sector          { inc_col(); return KEYWORD_SECTOR;       }
model           { inc_col(); return KEYWORD_MODEL;        }
occ             { inc_col(); return KEYWORD_OCC;          }
occ_irreps      { inc_col(); return KEYWORD_OCC_IRREPS;   }
active          { inc_col(); return KEYWORD_ACTIVE;       }
nactp           { inc_col(); return KEYWORD_NACTP;        }
nacth           { inc_col(); return KEYWORD_NACTH;        }
hermit          { inc_col(); return KEYWORD_HERMIT;       }
mstdm           { inc_col(); return KEYWORD_MSTDM;        }
natorb          { inc_col(); return KEYWORD_NATORB;       }
nroots          { inc_col(); return KEYWORD_NROOTS;       }
roots_cutoff    { inc_col(); return KEYWORD_ROOTS_CUTOFF; }
maxiter         { inc_col(); return KEYWORD_MAXITER;      }
conv            { inc_col(); return KEYWORD_CONV;         }
div_thresh      { inc_col(); return KEYWORD_DIV_THRESH;   }
damping         { inc_col(); return KEYWORD_DAMPING;      }
diis            { inc_col(); return KEYWORD_DIIS;         }
shifttype       { inc_col(); return KEYWORD_SHIFTTYPE;    }
shift           { inc_col(); return KEYWORD_SHIFT;        }
reuse           { inc_col(); return KEYWORD_REUSE;        }
flush           { inc_col(); return KEYWORD_FLUSH;        }
integrals       { inc_col(); return KEYWORD_INTEGRALS;    }
oneprop         { inc_col(); return KEYWORD_ONEPROP;      }
twoprop         { inc_col(); return KEYWORD_TWOPROP;      }
memory          { inc_col(); return KEYWORD_MEMORY;       }
tilesize        { inc_col(); return KEYWORD_TILESIZE;     }
disk_usage      { inc_col(); return KEYWORD_DISK_USAGE;   }
compress        { inc_col(); return KEYWORD_COMPRESS;     }
compress_triples { inc_col(); return KEYWORD_COMPRESS_TRIPLES; }
nthreads        { inc_col(); return KEYWORD_NTHREADS;     }
openmp          { inc_col(); return KEYWORD_OPENMP;       }
cuda            { inc_col(); return KEYWORD_CUDA;         }
arith           { inc_col(); return KEYWORD_ARITH;        }
mdprop          { inc_col(); return KEYWORD_MDPROP;       }
txtprop         { inc_col(); return KEYWORD_TXTPROP;      }
end             { inc_col(); return KEYWORD_END;          }
select          { inc_col(); return KEYWORD_SELECT;       }
intham1         { inc_col(); return KEYWORD_IH_IMMS;      }
ih_imms         { inc_col(); return KEYWORD_IH_IMMS;      }
gaunt           { inc_col(); return KEYWORD_GAUNT;        }
skip            { inc_col(); return KEYWORD_SKIP;         }
interface       { inc_col(); return KEYWORD_INTERFACE;    }
breit           { inc_col(); return KEYWORD_BREIT;        }
x2cmmf          { inc_col(); return KEYWORD_X2CMMF;       }
restrict_t3     { inc_col(); return KEYWORD_RESTRICT_T3;  }
spinor_labels   { inc_col(); return KEYWORD_SPINOR_LABELS;}
flush_amplitudes_txt { inc_col(); return KEYWORD_FLUSH_AMPLITUDES_TXT; }
analyt_prop     { inc_col(); return KEYWORD_ANALYT_PROP;  }
density         { inc_col(); return KEYWORD_DENSITY;      }
overlap         { inc_col(); return KEYWORD_OVERLAP;      }
hughes_kaldor_1h2p { inc_col(); return KEYWORD_HUGHES_KALDOR; }
!=              { inc_col(); return TT_NEQ;               }
=               { inc_col(); return TT_EQ;                }
\".*\"          { inc_col(); return TT_QUOTE;             }
[0-9]h[0-9]p    { inc_col(); return TT_SECTOR;            }
[-]             { inc_col(); return TT_HYPHEN;            }
[*]             { inc_col(); return TT_STAR;              }
\[[^ \t\n\[\]]+\]:[0-9]+   { inc_col(); return TT_ELEC_STATE; }
[0-9]*          { inc_col(); return TT_INTEGER;           }
[+-]?[0-9]*(\.[0-9]*)?([eE][+-]?[0-9]+)?   { inc_col(); return TT_FLOAT; }
([^ *\t\n\[\]])+ { inc_col(); return TT_WORD;     }
\n              { inc_col(); yylineno++; return END_OF_LINE; }
<<EOF>>         { return END_OF_FILE; }
%%

/* just dummy */
int yywrap() {
   return 1;
}

/**
 * updates the number of the current column
 */
void inc_col()
{
    // update yycol
    if (prev_len == -1) {
        yycol = 0;
    }
    else {
        yycol += prev_len;
    }

    // update prev_len
    if (strcmp(yytext, "\n") == 0) {
        prev_len = -1;
    }
    else {
        prev_len = strlen(yytext);
    }
}


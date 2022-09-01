#!/usr/bin/env python

# Test: finite-field calculation of transition dipole moments
# in the rubidium atom (the 0h1p Fock space sectors)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# obtain transformed integrals
dirac_inp = 'TRA.inp'
dirac_mol = "Rb.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

# Field -F
execute("expt.x --no-clean input_F- > ccsd_F-.out")
execute("mv scratch/HEFF HEFFM")

# Field +F
execute("expt.x --no-clean input_F+ > ccsd_F+.out")
execute("mv scratch/HEFF HEFFP")

# TDM calculation
filter_tdm12 = Filter("1  ->  2", [12555.138, None, 0.000000, 2.627905, 1.621082], [1e-1, None, 1e-4, 1e-4, 1e-4])
filter_tdm13 = Filter("1  ->  3", [12786.559, None, 0.000000, 5.173064, 2.274437], [1e-1, None, 1e-4, 1e-4, 1e-4])
filter_tdm23 = Filter("2  ->  3", [  231.421, None, 0.000000, 0.000000, 0.000012], [1e-1, None, 1e-4, 1e-4, 1e-4])

ret = Test("", "ff.inp", filters=[filter_tdm12,filter_tdm13,filter_tdm23], binary="heffman.x < ").run()

execute("mv ff.inp.test.out ff.out")
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(ret)

#!/usr/bin/env python

# Test: high-spin openshell CC calculation (sector 0h0p)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

dirac_inp = "TRA.inp"
dirac_mol = "O2.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")

t1_scf_read   = Filter("Total SCF energy = ",       -149.686661451845, 1e-7)
t1_scf_recalc = Filter("SCF reference energy = ",   -149.718144633814, 1e-7)
t1_mp2_corr   = Filter("MP2 correlation energy = ",   -0.383061310628, 1e-7)
t1_mp2_total  = Filter("Total MP2 energy = ",       -150.101205944442, 1e-7)
t1_ccsd_corr  = Filter("CCSD correlation energy = ",  -0.366958682628, 1e-7)
t1_ccsd_total = Filter("Total CCSD energy = ",      -150.085103316442, 1e-7)
ret_code = Test('O2 high spin CCSD (triplet)', "ccsd.inp", filters=[t1_scf_read,t1_mp2_total,t1_mp2_corr,t1_ccsd_corr,t1_ccsd_total,t1_scf_recalc]).run()

execute("mv ccsd.inp.test.out ccsd.out")

# cleanup
execute("rm -rf HINT VINT*")
execute("rm -rf MRCONEE* MDCINT* scratch")

sys.exit(ret_code)

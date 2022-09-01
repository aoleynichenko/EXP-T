#!/usr/bin/env python

# Test:
# calculation of the dipole moment of the LiF molecule using exact analytic
# density matrix in the CCSD(T) approximation (in the 0h0p sector).
# Lambda equations are solved to construct the CCSD(T) density matrix.
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
#symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']
symmetries = ['Cs', 'Cinfv']

for sym in symmetries:
    dirac_inp = "TRA.inp"
    dirac_mol = "LiF-%s.mol" % (sym)
    execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
    
    t1_scf  = Filter("Total SCF energy = ",       -107.057122389928, 1e-7)
    t1_mp2  = Filter("Total MP2 energy = ",       -107.284230678125, 1e-7)
    t1_ccsd = Filter("CCSD correlation energy = ",  -0.229801015210, 1e-7)
    t1_ccsdt= Filter("Total CCSD(T) energy = ",   -107.290154157247, 1e-7)
    t1_norm = Filter("norm                       =", 1.105784970064, 1e-7)
    t1_prop = Filter("Correlation contribution",    -0.069153508909, 1e-7)
    
    ret = Test(sym, "ccsd_t.inp", filters=[t1_scf,t1_ccsd,t1_ccsdt,t1_mp2,t1_norm,t1_prop]).run()
    ret_codes.append(ret)

    execute("mv ccsd_t.inp.test.out ccsd_t_%s.out" % (sym))
    execute("rm -rf MRCONEE* MDCINT* MR* MD*")
    execute("rm -rf scratch")

sys.exit(1 if any(ret_codes) else 0)

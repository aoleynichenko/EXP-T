#!/usr/bin/env python

# Test:
# calculation of the dipole moment of NO using exact analytic density
# matrix (in the 0h0p sector). Lambda equations are solved to construct
# the density matrix. Note that non-canonical orbitals are used in order
# to verify the correctness of diagrams involving the Fock operator.
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
#symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']
symmetries = ['C1', 'C2', 'C2v', 'Cinfv']

# todo: here we have troubles with Cs: wrong correlation energy & norm
# todo: wrong MP2 energies in C1 and Cs
# bugs in DIRAC? (imaginary matrix elements in HINT)
# 

for sym in symmetries:
    dirac_inp = "TRA.inp"
    dirac_mol = "NO-%s.mol" % (sym)
    execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
    t1_scf  = Filter("Total SCF energy = ",       -129.33752817437897, 1e-7)
    #t1_mp2  = Filter("Total MP2 energy = ",       -129.680552037539, 1e-7)
    t1_ccsd = Filter("CCSD correlation energy = ",  -0.338389211976, 1e-7)
    t1_norm = Filter("norm                       =", 1.192409560987, 1e-7)
    t1_prop = Filter("Correlation contribution",    -0.192944953356, 1e-7)
    
    ret = Test(sym, "ccsd.inp", filters=[t1_scf,t1_ccsd,t1_norm,t1_prop]).run()
    ret_codes.append(ret)
    
    execute("mv ccsd.inp.test.out ccsd_%s.out" % (sym))
    execute("rm -rf MRCONEE* MDCINT* MR* MD*")
    execute("rm -rf scratch")

sys.exit(1 if any(ret_codes) else 0)


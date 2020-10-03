#!/usr/bin/env python

# Test: external field at CC level (large T1 amplitudes)
# TODO: problem with the Cinfv test -- wrong SCF energy (both in DIRAC and EXPT)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute

print('>>> efield/LiNa/relativistic(4c-DC)')

symmetries = ['C1', 'C2', 'Cs', 'C2v']

for sym in symmetries:
    dirac_scf = "scf.inp"
    dirac_tra = "TRA.inp"
    dirac_mol = "LiNa-%s.mol" % (sym)
    execute("pam --nobackup --noarch --inp=" + dirac_scf + " --mol=" + dirac_mol + " --outcmo")
    execute("pam --nobackup --noarch --inp=" + dirac_tra + " --mol=" + dirac_mol + " --incmo --get=\"MRCONEE MDCINT\"")
    execute("mv MRCONEE MRCONEE-%s" % (sym))
    execute("mv MDCINT MDCINT-%s" % (sym))
    t1_scf  = Filter("SCF reference energy = ",   -169.073977054460   , 1e-8)
    t1_mp2c = Filter("MP2 correlation energy = ",   -0.039137634196386, 1e-8)
    t1_mp2  = Filter("Total MP2 energy = ",       -169.113114688656   , 1e-8)
    t1_ccsd = Filter("CCSD correlation energy = ",  -0.058503275496498, 1e-8)
    Test("%s" % (sym), "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd]).run()
    execute("rm -rf MRCONEE* MDCINT* DFCOEF")
    execute("rm -rf HINT VINT*")


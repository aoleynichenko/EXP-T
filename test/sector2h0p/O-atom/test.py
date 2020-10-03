#!/usr/bin/env python

# Test:
# (1) Low-lying excited states of the oxygen atom
# (2) FSCC scheme: O2- -> O- -> O (sector 2h0p)
# (3) Basis set: cc-pVDZ (it is just a test!)
# (4) DC - relativistic, 4-comp

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute

print('>>> sector(2,0)/O atom/4c')

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv', 'D2', 'Ci', 'D2h', 'C2h', 'Dinfh']

for sym in symmetries:
    dirac_inp = "TRA.inp"
    if sym == 'Ci' or sym == 'D2h' or sym == 'C2h' or sym == 'Dinfh':
        dirac_inp = "TRAi.inp"
    dirac_mol = "O-%s.mol" % (sym)
    execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
    execute("mv MRCONEE MRCONEE-" + sym)
    execute("mv MDCINT MDCINT-" + sym)
    t1_scf  = Filter("SCF reference energy = ",    -74.235483635912, 1e-7)
    t1_mp2c = Filter("MP2 correlation energy = ",   -0.179872821257, 1e-7)
    t1_mp2  = Filter("Total MP2 energy = ",        -74.415356457169, 1e-7)
    t1_ccsd = Filter("CCSD correlation energy = ",  -0.184200399151, 1e-7)
    t1_e1  = Filter("@    1", -0.5348266381, 1e-7)
    t1_e2  = Filter("@    2", -0.5341280558, 1e-7)
    t1_e3  = Filter("@    3", -0.5337869052, 1e-7)
    t1_e4  = Filter("@    4", -0.4527566128, 1e-7)
    t1_e5  = Filter("@    5", -0.3788442571, 1e-7)
    Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,
        t1_e1,t1_e2,t1_e3,t1_e4,t1_e5]).run()
    execute("rm -rf MRCONEE* MDCINT*")
    execute("rm -rf HINT VINT* modelvectors* HEFF")



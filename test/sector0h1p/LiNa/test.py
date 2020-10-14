#!/usr/bin/env python

# Test: diatomic molecules / LiNa / 4c Dirac-Coulomb / sector (0h,1p)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> sector(0,1)/LiNa+/4c')

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "LiNa-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -168.786710818500950, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",      -0.011109273875, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",          -168.797820092375, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.012008755284061, 1e-7)
	t1_e1 = Filter("@    1", -0.4307257633, 1e-7)
	t1_e2 = Filter("@    2", -0.3590926920, 1e-7)
	t1_e3 = Filter("@    3", -0.3590046164, 1e-7)
	t1_e4 = Filter("@    4", -0.3236142316, 1e-7)
	t1_e5 = Filter("@    5", -0.2387977486, 1e-7)
	t1_e6 = Filter("@    6", -0.2213285449, 1e-7)
	t1_e7 = Filter("@    7", -0.2213071158, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


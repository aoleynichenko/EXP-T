#!/usr/bin/env python

# Test:
# (1) Excited states of the LiNa diatomic molecule
# (2) FSCC scheme: LiNa2+ -> LiNa+ -> LiNa (sector 0h2p)
# (3) Basis set: cc-pVDZ (it is just a test!)
# (4) DC - relativistic, 4-comp

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> sector(0,2)/LiNa/4c')

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "LiNa-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -168.786710818500950, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",      -0.011109273875, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",          -168.797820092375, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.012008755284061, 1e-7)
	t1_e1  = Filter("@    1", -0.6378487411, 1e-7)
	t1_e2  = Filter("@    2", -0.5939084547, 1e-7)
	t1_e3  = Filter("@    3", -0.5939084194, 1e-7)
	t1_e4  = Filter("@    4", -0.5938671497, 1e-7)
	t1_e5  = Filter("@    5", -0.5938258273, 1e-7)
	t1_e6  = Filter("@    6", -0.5742603945, 1e-7)
	t1_e7  = Filter("@    7", -0.5742603817, 1e-7)
	t1_e8  = Filter("@    8", -0.5128750942, 1e-7)
	t1_e9  = Filter("@    9", -0.4955825337, 1e-7)
	t1_e10 = Filter("@   10", -0.4955825118, 1e-7)
	t1_e11 = Filter("@   11", -0.4940053034, 1e-7)
	t1_e12 = Filter("@   12", -0.4851803949, 1e-7)
	t1_e13 = Filter("@   13", -0.4851802936, 1e-7)
	t1_e14 = Filter("@   14", -0.4851386576, 1e-7)
	t1_e15 = Filter("@   15", -0.4850969338, 1e-7)
	t1_e16 = Filter("@   16", -0.4700040772, 1e-7)
	t1_e17 = Filter("@   17", -0.4550044987, 1e-7)
	t1_e18 = Filter("@   18", -0.4470212470, 1e-7)
	t1_e19 = Filter("@   19", -0.3868759477, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,
	    t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,
	    t1_e8,t1_e9,t1_e10,t1_e11,t1_e12,t1_e13,t1_e14,
	    t1_e15,t1_e16,t1_e17,t1_e18,t1_e19
	    ]).run()
	execute("rm -rf MRCONEE* MDCINT* MDPROP")
	execute("rm -rf HINT VINT* modelvectors* HEFF")



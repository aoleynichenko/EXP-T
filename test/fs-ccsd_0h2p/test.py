#!/usr/bin/env python

# Test:
# (1) Excited states of the LiNa diatomic molecule
# (2) FSCC scheme: LiNa2+ -> LiNa+ -> LiNa (sector 0h2p)
# (3) Basis set: cc-pVDZ (it is just a test!)
# (4) DC - relativistic, 4-comp

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "LiNa-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

	t1_scf  = Filter("Total SCF energy = ",       -168.786710818500950, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",      -0.011109273875, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",          -168.797820092375, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.012008755284061, 1e-7)
	t1_e1  = Filter("@    1", -0.6451463611, 1e-7)
	t1_e2  = Filter("@    2", -0.6012962699, 1e-7)
	t1_e3  = Filter("@    3", -0.6012962413, 1e-7)
	t1_e4  = Filter("@    4", -0.6012557395, 1e-7)
	t1_e5  = Filter("@    5", -0.6012151825, 1e-7)
	t1_e6  = Filter("@    6", -0.5816227784, 1e-7)
	t1_e7  = Filter("@    7", -0.5816227662, 1e-7)
	t1_e8  = Filter("@    8", -0.5317074950, 1e-7)
	t1_e9  = Filter("@    9", -0.5219558138, 1e-7)
	t1_e10 = Filter("@   10", -0.5054439772, 1e-7)
	t1_e11 = Filter("@   11", -0.5054438967, 1e-7)
	t1_e12 = Filter("@   12", -0.4935312187, 1e-7)
	t1_e13 = Filter("@   13", -0.4935311510, 1e-7)
	t1_e14 = Filter("@   14", -0.4934900533, 1e-7)
	t1_e15 = Filter("@   15", -0.4934492644, 1e-7)
	t1_e16 = Filter("@   16", -0.4925714975, 1e-7)
	t1_e17 = Filter("@   17", -0.4768345220, 1e-7)
	t1_e18 = Filter("@   18", -0.4688831061, 1e-7)
	t1_e19 = Filter("@   19", -0.4229924575, 1e-7)
	
	ret = Test(sym, "ccsd.inp", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,
	    t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,
	    t1_e8,t1_e9,t1_e10,t1_e11,t1_e12,t1_e13,t1_e14,
	    t1_e15,t1_e16,t1_e17,t1_e18,t1_e19
	    ]).run()
	ret_codes.append(ret)
	
	execute("mv ccsd.inp.test.out ccsd_%s.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP")
	execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(1 if any(ret_codes) else 0)


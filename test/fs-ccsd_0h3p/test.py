#!/usr/bin/env python

#
# Test: CCSD in the 0h3p sector
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

#
# Test 1. F ions
# F5+ (0h0p) -> F4+ (0h1p) -> F3+ (0h2p) -> F2+ (0h3p)
#
for sym in ["C1", "Cs", "C2", "Ci", "C2v", "C2h", "D2", "D2h", "Cinfv"]:
	
	dirac_inp = "TRA_F.inp"
	if sym == "Ci" or sym == "C2h" or sym == "D2h":
		dirac_inp = "TRA_F_i.inp"
	dirac_mol = "F-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	
	t2_scf  = Filter("Total SCF energy = ",       -87.987950437876, 1e-7)
	t2_mp2c = Filter("MP2 correlation energy = ",  -0.060847561030, 1e-7)
	t2_mp2  = Filter("Total MP2 energy = ",       -88.048797998906, 1e-7)
	t2_ccsd = Filter("CCSD correlation energy = ", -0.108826583117, 1e-7)
	t2_e1 = Filter("@    1", -9.6105603807, 1e-7)
	t2_e2 = Filter("@    2", -9.4425250936, 1e-7)
	t2_e3 = Filter("@    3", -9.4424623556, 1e-7)
	t2_e4 = Filter("@    4", -9.3726379865, 1e-7)
	t2_e5 = Filter("@    5", -9.3725415883, 1e-7)
	t2_ip1 = Filter("Ionization potential 0h3p -> 0h2p =", 2.248975900033, 1e-7)
	t2_ip2 = Filter("Ionization potential 0h2p -> 0h1p =", 3.181379904745, 1e-7)
	t2_ip3 = Filter("Ionization potential 0h1p -> 0h0p =", 4.180204596133, 1e-7)

	ret = Test('F5+/F4+/F3+/F2+ %s' % (sym), "ccsd.inp", filters=[t2_scf,t2_mp2,t2_mp2c,t2_ccsd,t2_e1,t2_e2,t2_e3,t2_e4,t2_e5,t2_ip1,t2_ip2,t2_ip3]).run()
	ret_codes.append(ret)
	
	execute("mv ccsd.inp.test.out ccsd_F_%s.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* scratch")
	execute("rm -rf HINT VINT*")

sys.exit(1 if any(ret_codes) else 0)



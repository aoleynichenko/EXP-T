#!/usr/bin/env python

# Test: OpenMP parallelization

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> sector 0h3p')

# Test 1. Ne ions. Source: Kaldor & Hughes, CPL 204, 339 (1993)
# Model: CCSD; basis set: ANO-RCC restricted to [7s7p4d3f]
# Ne6+ (0h0p) -> Ne5+ (0h1p) -> Ne4+ (0h2p) -> Ne3+ (0h3p)
dirac_inp = "TRA_Ne.inp"
dirac_mol = "Ne.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
t1_scf  = Filter("Total SCF energy = ",       -110.087095488765, 1e-7)
t1_mp2c = Filter("MP2 correlation energy = ",   -0.093665954790, 1e-7)
t1_mp2  = Filter("Total MP2 energy = ",       -110.180761443555, 1e-7)
t1_ccsd = Filter("CCSD correlation energy = ",  -0.144802800029, 1e-7)
t1_e1   = Filter("@    1",  -13.9987584068, 1e-7)
t1_e2   = Filter("@    2",  -13.8107014573, 1e-7)
t1_e3   = Filter("@    3",  -13.7284078845, 1e-7)
t1_ip3  = Filter("Ionization potential wrt reference state =", 13.998758406754, 1e-7)
Test('Ne6+/Ne5+/Ne4+/Ne3+/NR/C2v', "input-Ne", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_ip3]).run()
execute("rm -rf MRCONEE* MDCINT* scratch")
execute("rm -rf HINT VINT*")

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
	t2_ip3 = Filter("Ionization potential wrt reference state =", 9.610560380710, 1e-7)

	Test('F5+/F4+/F3+/F2+ %s' % (sym), "input-F", filters=[t2_scf,t2_mp2,t2_mp2c,t2_ccsd,t2_e1,t2_e2,t2_e3,t2_e4,t2_e5,t2_ip3]).run()
	
	execute("cp input-F.test.out input-F_%s.test.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* scratch")
	execute("rm -rf HINT VINT*")


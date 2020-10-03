#!/usr/bin/env python

# Test: Ne atom / nonrelativistic, compared with E. Eliav's result / sector (1h,1p)
# FSCC scheme: Ne(0) -> Ne(1,1) (via 1h0p = Ne+ and 0h1p = Ne-)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute

print('>>> sector(1,1)/Ne/nonrel')

# all symmetries to be tested
symmetries = ['C1', 'C2', 'Cs', 'C2v', 'Ci', 'C2h', 'D2', 'D2h', 'Cinfv', 'Dinfh']

for sym in symmetries:
	if sym == 'C1' or sym == 'C2' or sym == 'Cs' or sym == 'C2v' or sym == 'D2' or sym == 'Cinfv':
		dirac_inp = "TRA.inp"
	else:
		dirac_inp = 'TRAi.inp'
	dirac_mol = "Ne-%s.mol" % (sym)
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	t1_scf  = Filter("Total SCF energy = ",        -128.52777688911166, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",      -0.270756916860, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",          -128.798533805971, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",     -0.271995247737, 1e-7)
	t1_e1  = Filter("@    1", 0.0000000000, 1e-7)
	t1_e2  = Filter("@    2", 0.6046595005, 1e-7)
	t1_e3  = Filter("@    3", 0.6107148345, 1e-7)
	t1_e4  = Filter("@    4", 0.6674738228, 1e-7)
	t1_e5  = Filter("@    5", 0.6770335172, 1e-7)
	t1_e6  = Filter("@    6", 0.6801407839, 1e-7)
	t1_e7  = Filter("@    7", 0.6815529646, 1e-7)
	t1_e8  = Filter("@    8", 0.6815831430, 1e-7)
	t1_e9  = Filter("@    9", 0.6938805723, 1e-7)
	Test(sym, "input", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,
	    t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,t1_e8,t1_e9]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")
	execute("cp input.test.out EXPT_Ne_sym-%s.out" % (sym))


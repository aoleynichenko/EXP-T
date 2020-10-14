#!/usr/bin/env python

# Test: atomic calc-s / Pb atom / quasi-relativistic (AREP)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> atoms/Pb/non-relativistic(AREP)')

# all symmetries to be tested
symmetries = ['C1', 'Ci', 'C2', 'Cs', 'C2v', 'C2h', 'D2', 'D2h']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "Pb-%s.mol" % (sym)
	if sym == 'Ci' or sym == 'D2h' or sym == 'C2h':
		dirac_inp = "TRAi.inp"
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -191.534949024665394, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.111795227863725, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -191.646744252529118, 1e-7)
	t1_mp2  = Filter("CCSD correlation energy = ",  -0.121205718458854, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


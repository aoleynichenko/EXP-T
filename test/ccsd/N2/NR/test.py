#!/usr/bin/env python

# Test: diatomic molecules / N2 / non-relativistic

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> diatomic/N2/nonrelativistic')

# all symmetries to be tested
symmetries = ['C1', 'Ci', 'Cs', 'C2', 'C2v', 'D2', 'D2h']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	if sym == 'Ci' or sym == 'D2h':
		dirac_inp = "TRAi.inp"
	dirac_mol = "N2-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -108.954110777414186, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.310604639864779, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -109.264715417278964, 1e-7)
	t1_mp2  = Filter("CCSD correlation energy = ",  -0.313087787355777, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


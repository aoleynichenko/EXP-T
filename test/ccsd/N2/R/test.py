#!/usr/bin/env python

# Test: diatomic molecules / N2 / relativistic (4c)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))
from minitest import Test, Filter, execute

print '>>> diatomic/N2/relativistic(4c)'

# all symmetries to be tested
symmetries = ['C1', 'Ci', 'Cs', 'C2', 'C2v', 'D2', 'D2h', 'Cinfv', 'Dinfh']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	if sym == 'Ci' or sym == 'D2h' or sym == 'Dinfh':
		dirac_inp = "TRAi.inp"
	dirac_mol = "N2-%s.mol" % (sym)
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -109.016414945143367, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.310687969731266, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -109.327102914874629, 1e-7)
	t1_mp2  = Filter("CCSD correlation energy = ",  -0.313141436233609, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


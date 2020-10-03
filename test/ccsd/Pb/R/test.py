#!/usr/bin/env python

# Test: atomic calc-s / Pb atom / relativistic (SO-REP)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))
from minitest import Test, Filter, execute

print('>>> atoms/Pb/relativistic(SO-REP)')

# all symmetries to be tested
symmetries = ['C1', 'Ci', 'C2', 'Cs', 'C2v', 'C2h', 'D2', 'D2h', 'Cinfv', 'Dinfh']
#symmetries = ['Cinfv', 'Dinfh']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "Pb-%s.mol" % (sym)
	if sym == 'Ci' or sym == 'D2h' or sym == 'C2h' or sym == 'Dinfh':
		dirac_inp = "TRAi.inp"
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -191.730398371209787, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.112225500041208, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -191.842623871250993, 1e-7)
	t1_mp2  = Filter("CCSD correlation energy = ",  -0.121966864239563, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


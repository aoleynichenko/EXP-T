#!/usr/bin/env python

# Test: diatomic molecules / CO / relativistic

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))
from minitest import Test, Filter, execute

print('>>> diatomic/CO/relativistic(4c)')

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "CO-%s.mol" % (sym)
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -112.820480227130517, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.290892396007590, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -113.111372623138109, 1e-7)
	t1_mp2  = Filter("CCSD correlation energy = ",  -0.298117912056231, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


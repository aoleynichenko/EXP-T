#!/usr/bin/env python

# Test: diatomic molecules / CO / relativistic
# Excitation energies in the 1h1p sector
# IH-like shifting technique ("intham-1") is used together with imaginary shifts

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> CO molecule exc energies (1h1p) with imaginary shifts')

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "CO-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")

	t1_scf  = Filter("Total SCF energy = ",       -112.820480227130517, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.290892396007590, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -113.111372623138109, 1e-7)
	t1_mp2  = Filter("CCSD correlation energy = ",  -0.298117912056231, 1e-7)

	ee2 = Filter("@    2", 0.2338425936, 1e-7)
	ee3 = Filter("@    3", 0.2338427622, 1e-7)
	ee4 = Filter("@    4", 0.2340738270, 1e-7)
	ee5 = Filter("@    5", 0.2343060898, 1e-7)
	ee6 = Filter("@    6", 0.3074679705, 1e-7)
	ee7 = Filter("@    7", 0.3074718331, 1e-7)
	ee60= Filter("@   60", 1.0037554050, 1e-7)
	ee61= Filter("@   61", 1.0829364322, 1e-7)
	ee62= Filter("@   62", 1.0829364663, 1e-7)
	ee63= Filter("@   63", 1.1070829157, 1e-7)

	Test(sym, "input", filters=[t1_scf,t1_mp2,t1_mp2c,ee2,ee3,ee4,ee5,ee6,ee7,ee60,ee61,ee62,ee63]).run()
	execute("cp input.test.out %s.test.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


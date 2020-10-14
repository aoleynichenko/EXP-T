#!/usr/bin/env python

# Test: shifts of energy denominators

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> denom shifts/HgH2+/2c-ECP/Cinfv')

sym = 'Cinfv'
shifts = ['real', 'realimag', 'imag']

dirac_inp = "TRA.inp"
dirac_mol = "hgh+.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
execute("mv MRCONEE MRCONEE-" + sym)
execute("mv MDCINT MDCINT-" + sym)

exc_energies = {}
exc_energies['real']     = [-1.0814433294,-0.8090360638,-0.7169341316,-0.6878069165]
exc_energies['realimag'] = [-1.0821336765,-0.8087182315,-0.7159802151,-0.6865209284]
exc_energies['imag']     = [-1.0820083422,-0.8077867956,-0.7147328251,-0.6851998672]

for sh in shifts:
	t1_scf  = Filter("Total SCF energy = ",       -151.402195200361, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.170282605295, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -151.572477805656, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.178257741307, 1e-7)
	t1_e1 = Filter("@    1", exc_energies[sh][0], 1e-7)
	t1_e2 = Filter("@    2", exc_energies[sh][1], 1e-7)
	t1_e3 = Filter("@    3", exc_energies[sh][2], 1e-7)
	t1_e4 = Filter("@    4", exc_energies[sh][3], 1e-7)
	Test('shift type = %s' % sh, "input-%s-%s" % (sym,sh), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_e4]).run()
	execute("rm -rf HINT VINT*")
execute("rm -rf MRCONEE* MDCINT*")


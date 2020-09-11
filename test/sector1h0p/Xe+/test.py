#!/usr/bin/env python

# Test: Xe+ ion/ 2-comp Gatchina RECP / sector (1h,0p)
# FSCC scheme: Xe -> Xe+

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute

print '>>> sector(1,0)/Xe->Xe+/SOREP'

# all symmetries to be tested
symmetries = ['C2v', 'D2h', 'Cinfv', 'Dinfh']

for sym in symmetries:
	if sym == 'C1' or sym == 'C2' or sym == 'Cs' or sym == 'C2v' or sym == 'D2':
		dirac_inp = "TRA.inp"
	else:
		dirac_inp = 'TRAi.inp'
	dirac_mol = "Xe-%s.mol" % (sym)
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -15.298106998368, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",  -0.164624100128, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -15.462731098496, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ", -0.178676372580, 1e-7)
	t1_e1  = Filter("@    1", 0.4370793506, 1e-7)
	t1_e2  = Filter("@    2", 0.4851109206, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


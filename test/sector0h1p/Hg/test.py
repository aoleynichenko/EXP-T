#!/usr/bin/env python

# Test: Hg+ ion / 2-comp Gatchina RECP / sector (0h,1p)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute

print '>>> sector(0,1)/Hg+/SOREP'

# all symmetries to be tested
symmetries = ['C1', 'C2', 'Cs', 'C2v', 'Ci', 'C2h', 'D2', 'D2h']

for sym in symmetries:
	if sym == 'C1' or sym == 'C2' or sym == 'Cs' or sym == 'C2v' or sym == 'D2':
		dirac_inp = "TRA.inp"
	else:
		dirac_inp = 'TRAi.inp'
	dirac_mol = "Hg-%s.mol" % (sym)
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -152.230646138773693, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",      -0.148339079147, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",          -152.378985217920, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.153561572006815, 1e-7)
	t1_e1 = Filter("@    1", -0.6518859983, 1e-7)
	t1_e2 = Filter("@    2", -0.4268435202, 1e-7)
	t1_e3 = Filter("@    3", -0.3884627843, 1e-7)
	Test(sym, "input-%s" % (sym), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


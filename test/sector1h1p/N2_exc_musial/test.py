#!/usr/bin/env python

# Test: N2 molecule excitation energies (R = 2.068 a.u.)
# basis set: cc-pVDZ, 1s electrons are not correlated
# FSCC scheme: N2(0h0h) -> N2(1h1p) (via 1h0p = N2+ and 0h1p = N2-)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH


print('>>> sector(1h1p)/N2/nonrel')

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'Cinfv', 'Dinfh']

for sym in symmetries:
	dirac_inp = 'TRA.inp'
	dirac_mol = "N2-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
	t1_e2  = Filter("@    2", 0.2847991452, 1e-7)
	t1_e3  = Filter("@    3", 0.2970409016, 1e-7)
	t1_e4  = Filter("@    4", 0.3352715553, 1e-7)
	t1_e5  = Filter("@    5", 0.3457793799, 1e-7)
	t1_e6  = Filter("@    6", 0.3714861226, 1e-7)
	t1_e7  = Filter("@    7", 0.3790725681, 1e-7)
	t1_e8  = Filter("@    8", 0.3966063834, 1e-7)
	t1_e9  = Filter("@    9", 0.4202318050, 1e-7)
	t1_e10 = Filter("@   10", 0.5148711303, 1e-7)
	t1_e11 = Filter("@   11", 0.6227912851, 1e-7)
	Test(sym, "input", filters=[t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,t1_e8,t1_e9,t1_e10,t1_e11]).run(options="--no-clean")
	execute("mv input-1h1p.test.out expt_N2_%s_1h1p.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP*")
	execute("rm -rf scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")



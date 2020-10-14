#!/usr/bin/env python

# Test: LiF molecule excitation energies (R = 3.0 a.u.)
# FSCC scheme: LiF(0) -> LiF(1,1) (via 1h0p = LiF+ and 0h1p = LiF-)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

#=======================================================================
# non-relativistic
#=======================================================================

print('>>> sector(1,1)/LiF/nonrel')

# all symmetries to be tested
symmetries = ['C1', 'Cinfv']

for sym in symmetries:
	dirac_inp = 'TRA_NR.inp'
	dirac_mol = "LiF-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
	execute("expt.x --no-clean input-0h0p > expt_LiF_%s_NR_0h0p.out" % (sym))
	t1_scf    = Filter("Total SCF energy = ",          -106.964946434870, 1e-7)
	t1_cc_tot = Filter("Total CCSD energy = ",         -107.194656109696, 1e-7)
	t1_ccsd   = Filter("CCSD correlation energy = ",     -0.229709674826, 1e-7)	
	t1_e2  = Filter("@    2", 0.2298230870, 1e-7)
	t1_e3  = Filter("@    3", 0.2299751759, 1e-7)
	t1_e4  = Filter("@    4", 0.2448229001, 1e-7)
	t1_e5  = Filter("@    5", 0.2463911062, 1e-7)
	Test(sym, "input-1h1p", filters=[t1_scf,t1_cc_tot,t1_ccsd,
	    t1_e2,t1_e3,t1_e4,t1_e5]).run(options="--no-clean")
	execute("mv input-1h1p.test.out expt_LiF_%s_NR_1h1p.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP*")
	execute("rm -rf scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")

#=======================================================================
# relativistic
#=======================================================================

print('>>> sector(1,1)/LiF/rel')

for sym in symmetries:
	dirac_inp = 'TRA_R.inp'
	dirac_mol = "LiF-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
	execute("expt.x --no-clean input-0h0p > expt_LiF_%s_R_0h0p.out" % (sym))
	t1_scf    = Filter("Total SCF energy = ",          -107.057122389929, 1e-7)
	t1_cc_tot = Filter("Total CCSD energy = ",         -107.286923405332, 1e-7)
	t1_ccsd   = Filter("CCSD correlation energy = ",     -0.229801015403, 1e-7)	
	t1_e2  = Filter("@    2", 0.2287013711, 1e-7)
	t1_e3  = Filter("@    3", 0.2287732397, 1e-7)
	t1_e4  = Filter("@    4", 0.2299572197, 1e-7)
	t1_e5  = Filter("@    5", 0.2299636063, 1e-7)
	t1_e6  = Filter("@    6", 0.2300355311, 1e-7)
	t1_e7  = Filter("@    7", 0.2444290430, 1e-7)
	t1_e8  = Filter("@    8", 0.2444293843, 1e-7)
	t1_e9  = Filter("@    9", 0.2459982462, 1e-7)
	Test(sym, "input-1h1p", filters=[t1_scf,t1_cc_tot,t1_ccsd,
	    t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,t1_e8,t1_e9]).run(options="--no-clean")
	execute("mv input-1h1p.test.out expt_LiF_%s_R_1h1p.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP*")
	execute("rm -rf scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


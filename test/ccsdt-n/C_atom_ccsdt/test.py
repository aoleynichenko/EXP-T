#!/usr/bin/env python

# Test: C atom / cc-pvtz / Dirac-Coulomb Hamiltonian / CCSDT

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> ccsdt_0h2p/C/DC')

# all symmetries to be tested
#symmetries = ['C1', 'Cs', 'C2v', 'Cinfv']
symmetries = ['Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "C-%s.mol" % (sym)
	
	# external electric field F=0.1 is added at the integral
	# transformation step (see TRA.inp)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
	
	t_scf         = Filter("Total SCF energy = ",              -36.423980694324, 1e-8)
	t_scf_ref     = Filter("SCF reference energy = ",          -36.423980694324, 1e-8)
	t_mp2         = Filter("Total MP2 energy = ",              -36.474422551054, 1e-8)
	t_ccsdt_corr  = Filter("CCSDT correlation energy = ",       -0.079842226731, 1e-8)
	t_ccsdt       = Filter("Total CCSDT energy = ",            -36.503822921054, 1e-8)
	t_e1          = Filter("@    1", -1.3027424379, 1e-8)
	t_e2          = Filter("@    2", -1.3026583880, 1e-8)
	t_e3          = Filter("@    3", -1.3024915838, 1e-8)
	t_e4          = Filter("@    4", -1.2536383035, 1e-8)
	t_e5          = Filter("@    5", -1.2005053045, 1e-8)

	Test(sym, "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt,t_e1,t_e2,t_e3,t_e4,t_e5]).run()
	
	execute("cp input.test.out input-%s.test.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP* DFCOEF* scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


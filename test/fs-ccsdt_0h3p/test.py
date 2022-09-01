#!/usr/bin/env python

# Test: N atom / cc-pvtz / Dirac-Coulomb Hamiltonian / CCSDT(0h3p)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
#symmetries = ['C1', 'Cs', 'C2v', 'Cinfv']
symmetries = ['Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "N-%s.mol" % (sym)
	
	# external electric field F=0.1 is added at the integral
	# transformation step (see TRA.inp)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
	
	t_scf         = Filter("Total SCF energy = ",              -51.108312781201, 1e-8)
	t_scf_ref     = Filter("SCF reference energy = ",          -51.108312781201, 1e-8)
	t_mp2         = Filter("Total MP2 energy = ",              -51.165832656403, 1e-8)
	t_ccsdt_corr  = Filter("CCSDT correlation energy = ",       -0.092403826589, 1e-8)
	t_ccsdt       = Filter("Total CCSDT energy = ",            -51.200716607789, 1e-8)
	t_e1          = Filter("@    1", -3.3526310493, 1e-7)
	t_e2          = Filter("@    2", -3.2606348402, 1e-7)
	t_e3          = Filter("@    3", -3.2606196835, 1e-7)
	t_e4          = Filter("@    4", -3.2137325896, 1e-7)
	t_e5          = Filter("@    5", -3.2137090068, 1e-7)

	ret = Test(sym, "ccsdt.inp", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt,t_e1,t_e2,t_e3,t_e4,t_e5]).run()
	ret_codes.append(ret)
	
	execute("mv ccsdt.inp.test.out ccsdt_%s.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP* DFCOEF* scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(1 if any(ret_codes) else 0)

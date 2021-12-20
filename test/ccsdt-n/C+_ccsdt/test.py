#!/usr/bin/env python

# Test: C+ ion / cc-pvtz / Dirac-Coulomb Hamiltonian / CCSDT / ext elec field = 0.1 au

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> ccsdt_0h1p/C+/DC/F=0.1')

# all symmetries to be tested
#symmetries = ['C1', 'Cs', 'C2v', 'Cinfv']
#symmetries = ['Cs', 'Cinfv']
symmetries = ['Cinfv']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "C-%s.mol" % (sym)
	
	# external electric field F=0.1 is added at the integral
	# transformation step (see TRA.inp)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=SCF --mol=" + dirac_mol + " --outcmo")
	execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=" + dirac_mol + " --incmo --get=\"MRCONEE MDCINT MDPROP\"")
	
	t_scf         = Filter("Total SCF energy = ",              -36.423980694324, 1e-8)
	t_scf_ref     = Filter("SCF reference energy = ",          -36.423980694324, 1e-8)
	t_mp2         = Filter("Total MP2 energy = ",              -36.490822352282, 1e-8)
	t_ccsdt_corr  = Filter("CCSDT correlation energy = ",       -0.099105449665, 1e-8)
	t_ccsdt       = Filter("Total CCSDT energy = ",            -36.523086143989, 1e-8)
	t_e1          = Filter("@    1", -0.9010908606, 1e-8)
	t_e2          = Filter("@    2", -0.8999174953, 1e-8)
	t_e3          = Filter("@    3", -0.8997145388, 1e-8)
	Test(sym, "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt,t_e1,t_e2,t_e3]).run()
	
	execute("cp input.test.out input-%s.test.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT* MDPROP* DFCOEF* scratch")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


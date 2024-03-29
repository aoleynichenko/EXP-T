#!/usr/bin/env python

# Test: H2O / 3-21G / Dirac-Coulomb Hamiltonian / CCSDT-1b

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
from minitest import Test, Filter, execute, DIRAC_PATH

# list of return codes
ret_codes = []

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2v']

for sym in symmetries:
	dirac_inp = "TRA.inp"
	dirac_mol = "H2O-%s.mol" % (sym)
	execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	t1_scf         = Filter("Total SCF energy = ",             -75.633031179577841, 1e-8)
	t1_mp2         = Filter("Total MP2 energy = ",             -75.758608412401628, 1e-8)
	t1_ccsdt1_corr = Filter("CCSDT-1b correlation energy = ",  -0.134859557360,     1e-8)
	t1_ccsdt1      = Filter("Total CCSDT-1b energy = ",        -75.767890736635,    1e-8)

	ret = Test(sym, "ccsdt-1.inp", filters=[t1_scf,t1_mp2,t1_ccsdt1_corr,t1_ccsdt1]).run()
	ret_codes.append(ret)

	execute("mv ccsdt-1.inp.test.out ccsdt-1_%s.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(1 if any(ret_codes) else 0)


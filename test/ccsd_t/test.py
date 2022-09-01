#!/usr/bin/env python

# Test: H2O / 3-21G / Dirac-Coulomb Hamiltonian / CCSD(T)

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
	t1_scf         = Filter("Total SCF energy = ",         -75.633031179577841, 1e-8)
	t1_mp2         = Filter("Total MP2 energy = ",         -75.758608412401628, 1e-8)
	t1_ccsd_corr   = Filter("CCSD correlation energy = ",   -0.133263334542625, 1e-8)
	t1_ccsd        = Filter("Total CCSD energy = ",        -75.766294514120460, 1e-8)
	t1_ccsd_t4     = Filter("Total CCSD+T(CCSD) energy = ",-75.767932148913530, 1e-8)
	t1_ccsd_t5     = Filter("Total CCSD(T) energy = ",     -75.767833418596567, 1e-8)
	
	ret = Test(sym, "ccsd_t.inp", filters=[t1_scf, t1_mp2, t1_ccsd_corr, t1_ccsd, t1_ccsd_t5, t1_ccsd_t4]).run()
	ret_codes.append(ret)
	
	execute("mv ccsd_t.inp.test.out ccsd_t_%s.out" % (sym))
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(1 if any(ret_codes) else 0)


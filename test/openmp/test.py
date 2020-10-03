#!/usr/bin/env python

# Test: OpenMP parallelization

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute

print('>>> openmp/CO-C2v/R')

sym = 'C2v'
nthreads = [1,2,4,8]

for nth in nthreads:
	dirac_inp = "TRA.inp"
	dirac_mol = "CO-C2v.mol"
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	execute("mv MRCONEE MRCONEE-" + sym)
	execute("mv MDCINT MDCINT-" + sym)
	t1_scf  = Filter("Total SCF energy = ",       -112.820480227130517, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",   -0.290892396007590, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",       -113.111372623138109, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.298117912056231, 1e-7)
	Test('%s - %d threads' % (sym,nth), "input-%s-%s" % (sym,str(nth)), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT*")


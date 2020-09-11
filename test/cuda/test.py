#!/usr/bin/env python

# Test: CUDA parallelization
# Hg atom / 2-comp Gatchina RECP / sector (0h,2p)
# FSCC scheme: Hg2+ -> Hg+ -> Hg

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
from minitest import Test, Filter, execute

print '>>> cuda: Hg in C1 symmetry'

# cuda or no cuda
calctype = ['nocuda', 'cuda']

for ct in calctype:
	dirac_inp = 'TRA.inp'
	dirac_mol = "Hg-C1"
	execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")
	t1_scf  = Filter("Total SCF energy = ",       -152.230646138773693, 1e-7)
	t1_mp2c = Filter("MP2 correlation energy = ",      -0.148339079147, 1e-7)
	t1_mp2  = Filter("Total MP2 energy = ",          -152.378985217920, 1e-7)
	t1_ccsd = Filter("CCSD correlation energy = ",  -0.153561572006815, 1e-7)
	t1_e1  = Filter("@    1", -0.9994119373, 1e-7)
	t1_e2  = Filter("@    2", -0.8488692033, 1e-7)
	t1_e3  = Filter("@    3", -0.8395227548, 1e-7)
	t1_e4  = Filter("@    4", -0.8194611615, 1e-7)
	t1_e5  = Filter("@    5", -0.7277092060, 1e-7)
	t1_e6  = Filter("@    6", -0.6025352885, 1e-7)
	t1_e7  = Filter("@    7", -0.5910759351, 1e-7)
	t1_e8  = Filter("@    8", -0.5784400369, 1e-7)
	t1_e9  = Filter("@    9", -0.5451863419, 1e-7)
	t1_e10 = Filter("@   10", -0.4868173364, 1e-7)
	Test(ct, "input-%s" % (ct), filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,
	    t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,t1_e8,t1_e9,t1_e10]).run()
	execute("rm -rf MRCONEE* MDCINT*")
	execute("rm -rf HINT VINT* modelvectors* HEFF")


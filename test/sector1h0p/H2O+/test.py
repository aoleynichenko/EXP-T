#!/usr/bin/env python

# Test: H2O+ ion / cc-pVTZ / non-relativistic Hamiltonian / sector (1h,0p)
# FSCC scheme: H2O -> H2O+

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> sector(1,0)/H2O->H2O+/cc-pVTZ/nonrel')

dirac_mol = "H2O-C2v.mol"
dirac_inp = "TRA.inp"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")

t1_scf  = Filter("Total SCF energy = ",       -76.057114619775, 1e-7)
t1_mp2c = Filter("MP2 correlation energy = ",  -0.275116992297, 1e-7)
t1_mp2  = Filter("Total MP2 energy = ",       -76.332231612072, 1e-7)
t1_ccsd = Filter("CCSD correlation energy = ", -0.280866129199, 1e-7)
t1_e1  = Filter("@    1", 0.4557505780, 1e-7)
t1_e2  = Filter("@    2", 0.5376795920, 1e-7)
t1_e3  = Filter("@    3", 0.6919460871, 1e-7)
t1_e4  = Filter("@    4", 1.2211489399, 1e-7)
Test("H2O -> H2O+", "input", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_e4]).run()
execute("rm -rf MRCONEE MDCINT")
execute("rm -rf HINT VINT* modelvectors* HEFF")


#!/usr/bin/env python
#
# Test: one-electron property operator (magnetic hyperfine) included
# through the 'oneprop' keyword
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute

print('>>> oneprop/K atom HFS constant')

dirac_inp = "TRA.inp"
dirac_mol = "K.mol"
execute("pam --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
t1_scf  = Filter("Total SCF energy = ",         -601.371625198165, 1e-7)
t1_mp2c = Filter("MP2 correlation energy = ",     -0.330143371715, 1e-7)
t1_mp2  = Filter("Total MP2 energy = ",         -601.701768569880, 1e-7)
t1_ccsd = Filter("CCSD correlation energy = ",    -0.340281231740, 1e-7)
t1_e1 = Filter("@    1", -0.1577455494, 1e-7)
t1_e2 = Filter("@    2", -0.1576629222, 1e-7)
t1_e3 = Filter("@    3", -0.0978697539, 1e-7)
t1_e4 = Filter("@    4", -0.0978601598, 1e-7)
t1_e5 = Filter("@    5", -0.0976169635, 1e-7)
t1_e6 = Filter("@    6", -0.0976147799, 1e-7)
t1_e7 = Filter("@    7", -0.0976125973, 1e-7)
t1_e8 = Filter("@    8", -0.0976104157, 1e-7)
Test("Cinfv", "input", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,t1_e8]).run()
execute("rm -rf MRCONEE* MDCINT* MDPROP* scratch")


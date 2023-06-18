#!/usr/bin/env python

#
# Test: C+ ion electronic states in the 0h1p and 1h2p sectors
# basis set: ANO basis set by Roos & Wodmark, 1990
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

dirac_inp = 'TRA.inp'
dirac_mol = "C.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

t1_e1 = Filter("@    1", -0.6962981331, 1e-7)
t1_e2 = Filter("@    2", -0.5491731540, 1e-7)
t1_e3 = Filter("@    3", -0.4397103407, 1e-7)
t1_e4 = Filter("@    4", -0.3893097627, 1e-7)   
ret = Test("", "ccsd.inp", filters=[t1_e1,t1_e2,t1_e3,t1_e4]).run()
    
execute("mv ccsd.inp.test.out ccsd_1h2p.out")
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(ret)


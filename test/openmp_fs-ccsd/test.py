#!/usr/bin/env python

# Test:
# electronic states of the neutral cesium atom with different number of
# threads (1,2,4,8) (symmetries C1, Cs, C2v, Cinfv)
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
symmetries = ['C1', 'Cs', 'C2v', 'Cinfv']

for sym in symmetries:
    dirac_inp = "TRA.inp"
    dirac_mol = "Cs_%s.mol" % (sym)
    execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

    for nth in [1,2,4,8]:
        t1_e1  = Filter("@    1", -0.1406607633, 1e-6)
        t1_e2  = Filter("@    2", -0.0910861368, 1e-6)
        t1_e3  = Filter("@    3", -0.0886973489, 1e-6)
        t1_e4  = Filter("@    4", -0.0668647242, 1e-6)
        t1_e5  = Filter("@    5", -0.0667297263, 1e-6)

        # algorithm 1: "internal"
        ret = Test(sym, "ccsd_internal_nth%d.inp" % (nth), filters=[t1_e1,t1_e2,t1_e3,t1_e4,t1_e5]).run(options="--no-clean")
	execute("mv ccsd_internal_nth%d.inp.test.out ccsd_%s_internal_nth%d.out" % (nth,sym,nth))
        ret_codes.append(ret)

        # algorithm 2: "external"
        ret = Test(sym, "ccsd_external_nth%d.inp" % (nth), filters=[t1_e1,t1_e2,t1_e3,t1_e4,t1_e5]).run(options="--no-clean")
	execute("mv ccsd_external_nth%d.inp.test.out ccsd_%s_external_nth%d.out" % (nth,sym,nth))
        ret_codes.append(ret)
    
    execute("rm -rf MRCONEE* MDCINT* MDPROP")
    execute("rm -rf scratch")

sys.exit(1 if any(ret_codes) else 0)



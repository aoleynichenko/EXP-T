#!/usr/bin/env python

# Test:
# ground state of the Cs+ ion at the CCSD(T) level
# with different number of threads (1,2,4,8)
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

dirac_inp = "TRA.inp"
dirac_mol = "Cs_Dinfh.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

for nth in [1,2,4,8]:
    t1_sd  = Filter("Total CCSD energy",    -19.949909618500, 1e-6)
    t1_sdt = Filter("Total CCSD(T) energy", -19.956020188732, 1e-6)

    ret = Test("", "ccsd_t_nth%d.inp" % (nth), filters=[t1_sd,t1_sdt]).run(options="--no-clean")
    execute("mv ccsd_t_nth%d.inp.test.out ccsd_t_nth%d.out" % (nth,nth))

    ret_codes.append(ret)

execute("rm -rf MRCONEE* MDCINT* MDPROP")
execute("rm -rf scratch")

print(ret_codes)
sys.exit(1 if any(ret_codes) else 0)




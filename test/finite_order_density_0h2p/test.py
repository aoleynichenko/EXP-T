#!/usr/bin/env python

#
# Test: finite-order calculation of a density matrix for excited state.
# Ca atom, state 3Po_1, Fock space sector 0h2p.
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

#
# obtain transformed integrals
#
dirac_inp = 'moltra.inp'
dirac_mol = "Ca.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")


filter_e1 = Filter("@    1", -0.6573222136, 1e-8)
filter_e2 = Filter("@    2", -0.5874245104, 1e-8)
filter_e3 = Filter("@    3", -0.5871817425, 1e-8)
filter_e4 = Filter("@    4", -0.5866910266, 1e-8)

filter_occ1 = Filter("   1 1/2g+", 0.9959833, 1e-7)
filter_occ2 = Filter("   6 1/2u+", 0.9807309, 1e-7)
filter_occ3 = Filter("  11 3/2u+", 0.2445478, 1e-7)
filter_occ4 = Filter("  16 3/2g+", 0.0105082, 1e-7)
filter_occ5 = Filter("  21 1/2g-", 0.0101876, 1e-7)
filter_occ6 = Filter("  26 3/2u-", 0.0039305, 1e-7)

filter_occ_sum = Filter("sum of occupation numbers:", 10.0000000, 1e-7)


ret = Test("", "ccsd.inp", filters=[
  filter_e1, filter_e2, filter_e3, filter_e4,
  filter_occ1, filter_occ2, filter_occ3, filter_occ4, filter_occ5, filter_occ6,
  filter_occ_sum
]).run(options="--no-clean")
execute("mv ccsd.inp.test.out ccsd_density_matrix.out")


#
# clean up
#
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)

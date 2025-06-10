#!/usr/bin/env python

#
# Test: finite-order calculation of a density matrix
# for the ground state of the H2O+ cation (Fock space sector 1h0p).
#
# basis set: cc-pVDZ-DK.
# hamiltonian: x2cmmf
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

#
# obtain transformed integrals
#
dirac_inp = 'moltra.inp'
dirac_mol = "H2O+.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")


filter_e1 = Filter("@    1", 0.4332660980, 1e-8)
filter_e2 = Filter("@    2", 0.5186236150, 1e-8)
filter_e3 = Filter("@    3", 0.6781790165, 1e-8)

filter_occ1 = Filter("   1 1E", 0.9999608, 1e-7)
filter_occ2 = Filter("   6 1E", 0.9822474, 1e-7)
filter_occ3 = Filter("  11 1E", 0.0146818, 1e-7)
filter_occ4 = Filter("  16 1E", 0.0057380, 1e-7)
filter_occ5 = Filter("  21 1E", 0.0026564, 1e-7)
filter_occ6 = Filter("  26 1E", 0.0011475, 1e-7)
filter_occ7 = Filter("  31 2E", 0.0003846, 1e-7)
filter_occ8 = Filter("  36 1E", 0.0002821, 1e-7)
filter_occ9 = Filter("  41 2E", 0.0002330, 1e-7)
filter_occ10= Filter("  46 2E", 0.0000192, 1e-7)


filter_occ_sum = Filter("sum of occupation numbers:", 9.0000000, 1e-7)

ret = Test("", "ccsd.inp", filters=[
  filter_e1, filter_e2, filter_e3,
  filter_occ1, filter_occ2, filter_occ3, filter_occ4, filter_occ5,
  filter_occ6, filter_occ7, filter_occ8, filter_occ9, filter_occ10,
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

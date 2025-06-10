#!/usr/bin/env python

#
# Test: transition dipole moments in the Xe2^+ molecular ion.
# (direct calculation of TDMs in the 1h0p sector).
# vacuum state: closed-shell vdW dimer Xe2
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH, FSCC_PATH

#
# obtain transformed integrals
#
dirac_inp = 'moltra.inp'
dirac_mol = "Xe2.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

#
# CCSD + direct calculation of matrix elements
#
# up -> down
filter_tdm_fi_1 = Filter("1 (1/2g+) ->  1 (1/2u+)", [15187.21,     0.00, None, None, 1.203814], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_fi_2 = Filter("1 (1/2g+) ->  2 (1/2u+)", [15187.21, 24349.49, None, None, 1.971648], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_fi_3 = Filter("2 (1/2g+) ->  1 (1/2u+)", [31146.73,     0.00, None, None, 1.996281], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_fi_4 = Filter("2 (1/2g+) ->  2 (1/2u+)", [31146.73, 24349.49, None, None, 1.192081], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_fi_5 = Filter("1 (3/2g+) ->  1 (3/2u+)", [ 9437.69, 16548.43, None, None, 2.306687], [1e-1, 1e-1, None, None, 1e-5])

# down -> up
filter_tdm_if_1 = Filter("1 (1/2u+) ->  1 (1/2g+)", [    0.00, 15187.21, None, None, 1.179010], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_if_2 = Filter("1 (1/2u+) ->  2 (1/2g+)", [    0.00, 31146.73, None, None, 1.915530], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_if_3 = Filter("2 (1/2u+) ->  1 (1/2g+)", [24349.49, 15187.21, None, None, 1.989862], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_if_4 = Filter("2 (1/2u+) ->  2 (1/2g+)", [24349.49, 31146.73, None, None, 1.179884], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_if_5 = Filter("1 (3/2u+) ->  1 (3/2g+)", [16548.43,  9437.69, None, None, 2.320246], [1e-1, 1e-1, None, None, 1e-5])


#
# 1. calculate amplitudes
# 2. calculate transition moments
#
execute(FSCC_PATH + " --no-clean ccsd.inp > ccsd.out")

ret = Test("", "ccsd_tdm.inp", filters=[
  filter_tdm_fi_1, filter_tdm_fi_2, filter_tdm_fi_3, filter_tdm_fi_4, filter_tdm_fi_5,
  filter_tdm_if_1, filter_tdm_if_2, filter_tdm_if_3, filter_tdm_if_4, filter_tdm_if_5
]).run(options="--no-clean")
execute("mv ccsd_tdm.inp.test.out ccsd_direct_tdm.out")

#
# clean up
#
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)



 

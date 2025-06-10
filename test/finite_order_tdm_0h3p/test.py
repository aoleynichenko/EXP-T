#!/usr/bin/env python

#
# Test: transition dipole moments in the N atom (2p2 3s -> 2p3).
# (finite-order calculation of TDMs)
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH, FSCC_PATH

#
# obtain transformed integrals
#
dirac_inp = 'moltra.inp'
dirac_mol = "N.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

#
# CCSD + direct calculation of matrix elements
#

filter_tdm_1  = Filter("1 (1/2g+) ->  1 (1/2u-)", [82092.08,     0.00, None, None, 0.116990], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_2  = Filter("1 (1/2g+) ->  2 (1/2u-)", [82092.08, 19325.89, None, None, 0.001266], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_3  = Filter("1 (1/2g+) ->  4 (1/2u-)", [82092.08, 29143.99, None, None, 0.000725], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_4  = Filter("1 (1/2g+) ->  5 (1/2u-)", [82092.08, 29148.36, None, None, 0.001003], [1e-1, 1e-1, None, None, 1e-5])
			 
filter_tdm_5  = Filter("1 (1/2g+) ->  1 (3/2u+)", [82092.08,     0.00, None, None, 0.202632], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_6  = Filter("1 (1/2g+) ->  2 (3/2u+)", [82092.08, 19325.89, None, None, 0.002192], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_7  = Filter("1 (1/2g+) ->  4 (3/2u+)", [82092.08, 29143.99, None, None, 0.001256], [1e-1, 1e-1, None, None, 1e-5])
			 
filter_tdm_8  = Filter("1 (1/2g+) ->  1 (1/2u-)", [82092.08,     0.00, None, None, 0.116990], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_9  = Filter("1 (1/2g+) ->  2 (1/2u-)", [82092.08, 19325.89, None, None, 0.001266], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_10 = Filter("1 (1/2g+) ->  4 (1/2u-)", [82092.08, 29143.99, None, None, 0.000725], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_11 = Filter("1 (1/2g+) ->  5 (1/2u-)", [82092.08, 29148.36, None, None, 0.001003], [1e-1, 1e-1, None, None, 1e-5])

filter_tdm_12 = Filter("1 (1/2g+) ->  1 (3/2u+)", [82092.08,     0.00, None, None, 0.202632], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_13 = Filter("1 (1/2g+) ->  2 (3/2u+)", [82092.08, 19325.89, None, None, 0.002192], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_14 = Filter("1 (1/2g+) ->  4 (3/2u+)", [82092.08, 29143.99, None, None, 0.001256], [1e-1, 1e-1, None, None, 1e-5])

filter_tdm_15 = Filter("1 (1/2g+) ->  1 (1/2u+)", [82092.08,     0.00, None, None, 0.233979], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_16 = Filter("1 (1/2g+) ->  2 (1/2u+)", [82092.08, 19325.89, None, None, 0.002531], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_17 = Filter("1 (1/2g+) ->  4 (1/2u+)", [82092.08, 29143.99, None, None, 0.001451], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_18 = Filter("1 (1/2g+) ->  5 (1/2u+)", [82092.08, 29148.36, None, None, 0.001003], [1e-1, 1e-1, None, None, 1e-5])

filter_tdm_19 = Filter("1 (1/2g+) ->  1 (1/2u-)", [82092.08,     0.00, None, None, 0.116990], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_20 = Filter("1 (1/2g+) ->  1 (3/2u+)", [82092.08,     0.00, None, None, 0.202632], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_21 = Filter("1 (1/2g+) ->  1 (1/2u-)", [82092.08,     0.00, None, None, 0.116990], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_22 = Filter("1 (1/2g+) ->  1 (3/2u+)", [82092.08,     0.00, None, None, 0.202632], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_23 = Filter("1 (1/2g+) ->  1 (1/2u+)", [82092.08,     0.00, None, None, 0.233979], [1e-1, 1e-1, None, None, 1e-5])


#
# 1. calculate amplitudes
# 2. calculate transition moments
#
execute(FSCC_PATH + " --no-clean ccsd.inp > ccsd.out")

ret = Test("", "ccsd_tdm.inp", filters=[
	filter_tdm_1,  filter_tdm_2,  filter_tdm_3,  filter_tdm_4,  filter_tdm_5,
	filter_tdm_6,  filter_tdm_7,  filter_tdm_8,  filter_tdm_9,  filter_tdm_10,
	filter_tdm_11, filter_tdm_12, filter_tdm_13, filter_tdm_14, filter_tdm_15,
	filter_tdm_16, filter_tdm_17, filter_tdm_18, filter_tdm_19, filter_tdm_20,
	filter_tdm_21, filter_tdm_22, filter_tdm_23
]).run(options="--no-clean")
execute("mv ccsd_tdm.inp.test.out ccsd_tdm.out")

#
# clean up
#
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)

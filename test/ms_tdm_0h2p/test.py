#!/usr/bin/env python

#
# Test: transition dipole moments in the Hg atom.
# (model space estimates of TDMs)
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

#
# obtain transformed integrals
#
dirac_inp = 'TRA.inp'
dirac_mol = "Hg.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

#
# CCSD + model-space estimates
#
filter_tdm_1 = Filter("1 (  0g) ->  2 (  0u)", [    0.00, 35725.41, None, None, 0.000000, 0.000000, 0.236732], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_2 = Filter("1 (  0g) ->  4 (  0u)", [    0.00, 53579.09, None, None, 0.000000, 0.000000, 2.033768], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_3 = Filter("1 (  0g) ->  1 ( 1u+)", [    0.00, 35725.41, None, None, 0.167395, 0.167395, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_4 = Filter("1 (  0g) ->  3 ( 1u+)", [    0.00, 53579.09, None, None, 1.438091, 1.438091, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_5 = Filter("1 (  0g) ->  1 ( 1u-)", [    0.00, 35725.41, None, None, 0.167395, 0.167395, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_6 = Filter("1 (  0g) ->  3 ( 1u-)", [    0.00, 53579.09, None, None, 1.438091, 1.438091, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_7 = Filter("2 (  0u) ->  1 (  0g)", [35725.41,     0.00, None, None, 0.000000, 0.000000, 0.251651], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_8 = Filter("4 (  0u) ->  1 (  0g)", [53579.09,     0.00, None, None, 0.000000, 0.000000, 2.060586], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_9 = Filter("1 ( 1u+) ->  1 (  0g)", [35725.41,     0.00, None, None, 0.177944, 0.177944, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_10= Filter("3 ( 1u+) ->  1 (  0g)", [53579.09,     0.00, None, None, 1.457055, 1.457055, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_11= Filter("1 ( 1u-) ->  1 (  0g)", [35725.41,     0.00, None, None, 0.177944, 0.177944, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])
filter_tdm_12= Filter("3 ( 1u-) ->  1 (  0g)", [53579.09,     0.00, None, None, 1.457055, 1.457055, 0.000000], [1e-1, 1e-1, None, None, 1e-5, 1e-5, 1e-5])

ret = Test("", "ccsd.inp", filters=[
	filter_tdm_1,filter_tdm_2,filter_tdm_3,filter_tdm_4,filter_tdm_5,filter_tdm_6,
	filter_tdm_7,filter_tdm_8,filter_tdm_9,filter_tdm_10,filter_tdm_11,filter_tdm_12
]).run(options="--no-clean")
execute("mv ccsd.inp.test.out ccsd_model_space_tdm.out")

#
# clean up
#
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)

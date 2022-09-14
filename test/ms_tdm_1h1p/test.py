#!/usr/bin/env python

#
# Test: transition dipole moments in the CO molecule.
# (model space estimates of TDMs)
#
# This test is based on the TDM calculations for CO published in:
# N. S. Mosyagin, A. V. Oleynichenko, A. Zaitsevskii, A. V. Kudrin,
# E. A. Pazyuk, A. V. Stolyarov.
# J. Quant. Spectrosc. Radiat. Transf. 263, 107532 (2021)
# doi: 10.1016/j.jqsrt.2021.107532
#
# Excited states of the CO molecule are obtained in the 1h1p sector.
# The X1Sigma+ -> A1Pi is studied here (the 0h0p -> 1h1p type transition).
# TDM is obtained using model-space estimates.
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

#
# obtain transformed integrals
#
dirac_inp = 'TRA.inp'
dirac_mol = "CO.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

#
# CCSD + model-space estimates
#
filter_tdm_1 = Filter("1 (   a) ->  1 (   a)", [0.00, 50354.75, 0.000971], [1e-1, 1e-1, 1e-5])
filter_tdm_2 = Filter("1 (   a) ->  2 (   a)", [0.00, 50401.90, 0.002215], [1e-1, 1e-1, 1e-5])
filter_tdm_3 = Filter("1 (   a) ->  4 (   a)", [0.00, 67140.35, 0.009298], [1e-1, 1e-1, 1e-5])
filter_tdm_4 = Filter("1 (   a) ->  5 (   a)", [0.00, 69111.97, 0.873807], [1e-1, 1e-1, 1e-5])
filter_tdm_5 = Filter("1 (   a) ->  8 (   a)", [0.00, 74645.38, 0.005408], [1e-1, 1e-1, 1e-5])

ret = Test("", "ccsd_mstdm.inp", filters=[filter_tdm_1,filter_tdm_2,filter_tdm_3,filter_tdm_4,filter_tdm_5]).run(options="--no-clean")
ret_codes.append(ret)
execute("mv ccsd_mstdm.inp.test.out ccsd_model_space_tdm.out")

#
# CCSD + model-space estimates via the 'mdprop' keywords and selected components (X,Y,Z).
# The alternative code is used to calculate TDMs.
#

# X component, 0h0p -> 1h1p

filter_tdm_1 = Filter("1 (a   ) ->  2 (a   )", [0.00, 50401.90, None, None, 0.002168], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_2 = Filter("1 (a   ) ->  4 (a   )", [0.00, 67140.35, None, None, 0.009362], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_3 = Filter("1 (a   ) ->  5 (a   )", [0.00, 69111.97, None, None, 0.882473], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_4 = Filter("1 (a   ) ->  8 (a   )", [0.00, 74645.38, None, None, 0.005511], [1e-1, 1e-1, None, None, 1e-5])
ret = Test("", "ccsd_xdiplen.inp", filters=[filter_tdm_1,filter_tdm_2,filter_tdm_3,filter_tdm_4]).run(options="--no-clean")
ret_codes.append(ret)
execute("mv ccsd_xdiplen.inp.test.out ccsd_model_space_tdm_x.out")

# Y component, 1h1p -> 0h0p

filter_tdm_1 = Filter("2 (b   ) ->  1 (a   )", [50401.90, 0.00, None, None, 0.002189], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_2 = Filter("4 (b   ) ->  1 (a   )", [67140.35, 0.00, None, None, 0.008720], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_3 = Filter("6 (b   ) ->  1 (a   )", [69111.97, 0.00, None, None, 0.884450], [1e-1, 1e-1, None, None, 1e-5])
filter_tdm_4 = Filter("9 (b   ) ->  1 (a   )", [74645.38, 0.00, None, None, 0.005610], [1e-1, 1e-1, None, None, 1e-5])
ret = Test("", "ccsd_ydiplen.inp", filters=[filter_tdm_1,filter_tdm_2,filter_tdm_3,filter_tdm_4]).run(options="--no-clean")
ret_codes.append(ret)
execute("mv ccsd_ydiplen.inp.test.out ccsd_model_space_tdm_y.out")

# Z component, 1h1p -> 1h1p

filter_tdm_1 = Filter("1 (a   ) ->  1 (a   )", [50354.75, 0.00, None, None, 0.000211], [1e-1, 1e-1, None, None, 1e-5])
ret = Test("", "ccsd_zdiplen.inp", filters=[filter_tdm_1]).run(options="--no-clean")
ret_codes.append(ret)
execute("mv ccsd_zdiplen.inp.test.out ccsd_model_space_tdm_z.out")

#
# clean up
#
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(1 if any(ret_codes) else 0)


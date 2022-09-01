#!/usr/bin/env python

#
# Test: transition dipole moments in the CO molecule.
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

#
# obtain transformed integrals
#
dirac_inp = 'TRA.inp'
dirac_mol = "CO.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

#
# CCSD + model-space estimates
#
filter_tdm = Filter("1 (   a) ->  5 (   a)", [0.00, 69111.97, 0.873807], [1e-5, 1e-5, 1e-4])
ret = Test("", "ccsd.inp", filters=[filter_tdm]).run(options=" --no-clean ")

execute("mv ccsd.inp.test.out ccsd_model_space_tdm.out")
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)



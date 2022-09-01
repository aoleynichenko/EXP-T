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
# TDM is obtained using the finite-field approach.
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
# CCSD + finite-field estimates
#

# Field -F
execute("expt.x --no-clean ccsd_F-.inp > ccsd_F-.out")
execute("mv scratch/HEFF HEFFM")

# Field +F
execute("expt.x --no-clean ccsd_F+.inp > ccsd_F+.out")
execute("mv scratch/HEFF HEFFP")

# TDM calculation
filter_tdm = Filter("1  ->  6", [69111.973, None, None, 0.385949], [1e-1, None, None, 1e-4])
ret = Test("", "ff.inp", filters=[filter_tdm], binary="heffman.x < ").run()

execute("mv ff.inp.test.out finite_field_tdm.out")
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)



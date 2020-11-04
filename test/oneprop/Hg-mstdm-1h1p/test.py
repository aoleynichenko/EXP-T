#!/usr/bin/env python
#
# Test: model-space estimations of transition dipole moments
# mercury atom, 2c-ECP + modest basis set, sector 1h1p, 
# active space: 6s,7s,6d,6p,7p spinors
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> oneprop/Hg atom TDMs (ms-tdm)')

dirac_inp = "TRA.inp"
dirac_mol = "Hg.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")
t1_scf  = Filter("Total SCF energy = ",         -153.178023964342, 1e-7)
t1_mp2c = Filter("MP2 correlation energy = ",     -0.545006422714, 1e-7)
t1_mp2  = Filter("Total MP2 energy = ",         -153.723030387056, 1e-7)
t1_ccsd = Filter("CCSD correlation energy = ",    -0.513654935169, 1e-7)
t1_e1 = Filter("@    1", 0.0000000000, 1e-7)
t1_e2 = Filter("@    2", 0.1900034332, 1e-7)
t1_e3 = Filter("@    3", 0.1975938073, 1e-7)
t1_e4 = Filter("@    4", 0.2176614079, 1e-7)
t1_e5 = Filter("@    5", 0.2548074651, 1e-7)
t1_e6 = Filter("@    6", 0.2846595322, 1e-7)

t1_tdm1 = Filter("1 (  0g) ->  2 (  0u)", [    0.00, 43366.83, 0.462529], [1e-2, 1e-2, 1e-4])
t1_tdm2 = Filter("1 (  0g) ->  4 (  0u)", [    0.00, 55923.77, 2.061929], [1e-2, 1e-2, 1e-4])
t1_tdm3 = Filter("2 (  0u) ->  1 (  0g)", [43366.83,     0.00, 0.450713], [1e-2, 1e-2, 1e-4])
t1_tdm4 = Filter("4 (  0u) ->  1 (  0g)", [55923.77,     0.00, 2.062239], [1e-2, 1e-2, 1e-4])
t1_tdm5 = Filter("1 (  0g) ->  1 ( 1u+)", [    0.00, 43366.83, 0.462529], [1e-2, 1e-2, 1e-4])
t1_tdm6 = Filter("1 (  0g) ->  3 ( 1u+)", [    0.00, 55923.77, 2.061929], [1e-2, 1e-2, 1e-4])
t1_tdm7 = Filter("1 ( 1u+) ->  1 (  0g)", [43366.83,     0.00, 0.450713], [1e-2, 1e-2, 1e-4])
t1_tdm8 = Filter("3 ( 1u+) ->  1 (  0g)", [55923.77,     0.00, 2.062239], [1e-2, 1e-2, 1e-4])
t1_tdm9 = Filter("1 (  0g) ->  1 (  0u)", [62475.55, 41700.93, 1.342872], [1e-2, 1e-2, 1e-4])
t1_tdm10= Filter("1 (  0g) ->  3 (  0u)", [62475.55, 47771.16, 2.293086], [1e-2, 1e-2, 1e-4])
t1_tdm11= Filter("1 (  0g) ->  2 ( 1u+)", [62475.55, 47771.16, 1.985871], [1e-2, 1e-2, 1e-4])

Test("D2h", "input", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,
	t1_tdm1,t1_tdm2,t1_tdm3,t1_tdm4,t1_tdm5,t1_tdm6,t1_tdm7,t1_tdm8,t1_tdm9,t1_tdm10,t1_tdm11]).run()
execute("rm -rf MRCONEE* MDCINT* MDPROP* scratch")


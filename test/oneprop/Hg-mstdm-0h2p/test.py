#!/usr/bin/env python
#
# Test: model-space estimations of transition dipole moments
# mercury atom, 2c-ECP + modest basis set, sector 0h2p, 
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
t1_scf  = Filter("Total SCF energy = ",         -152.233275912798, 1e-7)
t1_mp2c = Filter("MP2 correlation energy = ",     -0.411901595843, 1e-7)
t1_mp2  = Filter("Total MP2 energy = ",         -152.645177508641, 1e-7)
t1_ccsd = Filter("CCSD correlation energy = ",    -0.399992531506, 1e-7)
t1_e1 = Filter("@    1", -1.0604225135, 1e-7)
t1_e2 = Filter("@    2", -0.8902262347, 1e-7)
t1_e3 = Filter("@    3", -0.8820431255, 1e-7)
t1_e4 = Filter("@    4", -0.8611526773, 1e-7)
t1_e5 = Filter("@    5", -0.8129763007, 1e-7)

t1_tdm1 = Filter("1 (  0g) ->  2 (  0u)", [0.00, 39149.75, 0.309501], [1e-2, 1e-2, 1e-4])
t1_tdm2 = Filter("1 (  0g) ->  4 (  0u)", [0.00, 54308.17, 1.868138], [1e-2, 1e-2, 1e-4])
t1_tdm3 = Filter("1 (  0g) ->  1 ( 1u+)", [0.00, 39149.75, 0.309501], [1e-2, 1e-2, 1e-4])
t1_tdm4 = Filter("1 (  0g) ->  3 ( 1u+)", [0.00, 54308.17, 1.868138], [1e-2, 1e-2, 1e-4])
t1_tdm5 = Filter("2 (  0u) ->  1 (  0g)", [39149.75, 0.00, 0.315500], [1e-2, 1e-2, 1e-4])
t1_tdm6 = Filter("4 (  0u) ->  1 (  0g)", [54308.17, 0.00, 1.999312], [1e-2, 1e-2, 1e-4])
t1_tdm7 = Filter("1 ( 1u+) ->  1 (  0g)", [39149.75, 0.00, 0.315500], [1e-2, 1e-2, 1e-4])
t1_tdm8 = Filter("3 ( 1u+) ->  1 (  0g)", [54308.17, 0.00, 1.999312], [1e-2, 1e-2, 1e-4])

Test("D2h", "input", filters=[t1_scf,t1_mp2,t1_mp2c,t1_ccsd,t1_e1,t1_e2,t1_e3,t1_e4,t1_e5,
	t1_tdm1,t1_tdm2,t1_tdm3,t1_tdm4,t1_tdm5,t1_tdm6,t1_tdm7,t1_tdm8]).run()
execute("rm -rf MRCONEE* MDCINT* MDPROP* scratch")


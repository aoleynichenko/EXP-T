#!/usr/bin/env python

# Test: N2 molecule / cc-pVDZ / non-relativistic Hamiltonian / FS-CCSDT(1h1p)
# The Cinfv point group is employed.
#
# This test is borrowed from M. Musial, R. J. Bartlett, 
# J. Chem. Phys. 121, 1670 (2004); doi: 10.1063/1.1765096.
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> ccsdt 1h1p / N2 molecule excitations / Non-rel')

execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=N2 --get=\"MRCONEE MDCINT MDPROP\"")

t_scf         = Filter("Total SCF energy = ",              -108.95453941355765  , 1e-8)
t_scf_ref     = Filter("SCF reference energy = ",          -108.954539413558, 1e-8)
t_mp2         = Filter("Total MP2 energy = ",              -109.259826827630, 1e-8)
t_ccsdt_corr  = Filter("CCSDT correlation energy = ",        -0.320348289208, 1e-8)
t_ccsdt       = Filter("Total CCSDT energy = ",            -109.274887702766, 1e-8)
t_e2          = Filter("@    2", 0.2906764563, 1e-8)
t_e3          = Filter("@    3", 0.2995912335, 1e-8)
t_e4          = Filter("@    4", 0.3382594556, 1e-8)
t_e5          = Filter("@    5", 0.3535514212, 1e-8)
t_e6          = Filter("@    6", 0.3682194407, 1e-8)
t_e7          = Filter("@    7", 0.3795176318, 1e-8)
t_e8          = Filter("@    8", 0.3940401922, 1e-8)
t_e9          = Filter("@    9", 0.4197055451, 1e-8)
t_e10         = Filter("@   10", 0.5065655757, 1e-8)
t_e11         = Filter("@   11", 0.6008743760, 1e-8)
Test("Cinfv", "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt,t_e2,t_e3,t_e4,t_e5,t_e6,t_e7,t_e8,t_e9,t_e10,t_e11]).run()

execute("rm -rf DFCOEF MDPROP* MRCONEE* MDCINT* scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")


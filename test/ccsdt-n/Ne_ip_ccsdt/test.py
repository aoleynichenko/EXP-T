#!/usr/bin/env python

# Test: Ne atom / cc-pVTZ / Dirac-Coulomb Hamiltonian / CCSDT(1h,0p)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> ccsdt 1h0p / Ne atom / DC')

execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=Ne-Cinfv --get=\"MRCONEE MDCINT MDPROP\"")

t_scf         = Filter("Total SCF energy = ",              -128.675592966117, 1e-8)
t_scf_ref     = Filter("SCF reference energy = ",          -128.675592966117, 1e-8)
t_mp2         = Filter("Total MP2 energy = ",              -128.953125745628, 1e-8)
t_ccsdt_corr  = Filter("CCSDT correlation energy = ",        -0.283489980347, 1e-8)
t_ccsdt       = Filter("Total CCSDT energy = ",            -128.959082946465, 1e-8)
t_e1          = Filter("@    1", 0.7794751304, 1e-8)
t_e2          = Filter("@    2", 0.7831702396, 1e-8)
t_e3          = Filter("@    3", 1.7792394091, 1e-8)
Test("Cinfv", "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt,t_e1,t_e2,t_e3]).run()

execute("rm -rf DFCOEF MDPROP* MRCONEE* MDCINT* scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")


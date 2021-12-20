#!/usr/bin/env python

# Test: Be atom / cc-pvdz / Dirac-Coulomb Hamiltonian / CCSDT / ext elec field = 0.1 au

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> ccsdt/Be/DC/F=0.1')

# external electric field F=0.1 is added at the integral
# transformation step (see TRA.inp)
execute(DIRAC_PATH + " --nobackup --noarch --inp=SCF --mol=Be-Cinfv --outcmo")
execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=Be-Cinfv --incmo --get=\"MRCONEE MDCINT MDPROP\"")

t_scf         = Filter("Total SCF energy = ",              -14.575193658320, 1e-8)
t_scf_ref     = Filter("SCF reference energy = ",          -14.575193658320, 1e-8)
t_mp2         = Filter("Total MP2 energy = ",              -14.751448609516, 1e-8)
t_ccsdt_corr  = Filter("CCSDT correlation energy = ",       -0.217623782329, 1e-8)
t_ccsdt       = Filter("Total CCSDT energy = ",            -14.792817440649, 1e-8)
Test("Cinfv", "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt]).run()

execute("rm -rf DFCOEF MDPROP* MRCONEE* MDCINT* scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")


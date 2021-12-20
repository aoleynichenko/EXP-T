#!/usr/bin/env python

# Test: CO / cc-pvdz / Dirac-Coulomb Hamiltonian / CCSDT-3 / ext elec field = 0.01 au

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> ccsdt-3/CO/DC/F=0.01')

# external electric field F=0.01 is added at the integral
# transformation step (see TRA.inp)
execute(DIRAC_PATH + " --nobackup --noarch --inp=SCF --mol=CO-Cinfv --outcmo")
execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA --mol=CO-Cinfv --incmo --get=\"MRCONEE MDCINT\"")

t_scf         = Filter("Total SCF energy = ",             -112.820480227131, 1e-8)
t_scf_ref     = Filter("SCF reference energy = ",         -112.819539718659, 1e-8)
t_mp2         = Filter("Total MP2 energy = ",             -113.112174768284, 1e-8)
t_ccsdt3_corr = Filter("CCSDT-3 correlation energy = ",     -0.311325279213, 1e-8)
t_ccsdt3      = Filter("Total CCSDT-3 energy = ",         -113.130864997872, 1e-8)
Test("Cinfv", "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt3_corr,t_ccsdt3]).run()

execute("rm -rf DFCOEF MRCONEE* MDCINT* scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")


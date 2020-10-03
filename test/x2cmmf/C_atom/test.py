#!/usr/bin/env python

# Test: C atom / cc-pvqz / X2Cmmf Hamiltonian + Gaunt

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute

print('>>> C atom with X2cmmf+Gaunt')

dirac_inp = "TRA.inp"
dirac_mol = "C.mol"

execute("pam --nobackup --noarch --inp=TRA --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

t_scf         = Filter("Total SCF energy = ",              -36.422201451037, 1e-8)
t_scf_ref     = Filter("SCF reference energy = ",          -36.422201451037, 1e-8)
t_mp2         = Filter("Total MP2 energy = ",              -36.490103660871, 1e-8)
t_ccsdt_corr  = Filter("CCSD correlation energy = ",        -0.096747321824, 1e-8)
t_ccsdt       = Filter("Total CCSD energy = ",             -36.518948772861, 1e-8)
t_e1          = Filter("@    1", -1.3048364584, 1e-8)
t_e2          = Filter("@    2", -1.3047707272, 1e-8)
t_e3          = Filter("@    3", -1.3046401154, 1e-8)
t_e4          = Filter("@    4", -1.2600070675, 1e-8)
t_e5          = Filter("@    5", -1.2046505808, 1e-8)

Test("Dinfh", "input", filters=[t_scf,t_scf_ref,t_mp2,t_ccsdt_corr,t_ccsdt,t_e1,t_e2,t_e3,t_e4,t_e5]).run()
	
execute("rm -rf MRCONEE* MDCINT* MDPROP* DFCOEF* scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")


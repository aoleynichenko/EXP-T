#!/usr/bin/env python

# Test: LiF molecule excitation energies (R = 3.0 a.u.)
# FSCC scheme: LiF(0) -> LiF(1,1) (via 1h0p = LiF+ and 0h1p = LiF-)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# obtain transformed integrals
dirac_inp = 'TRA.inp'
dirac_mol = "LiF.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

# solve CCSDT equations in the 0h0p sector
execute("expt.x --no-clean input-0h0p > ccsdt_0h0p.out")

# solve CCSDT equations in the 0h1p, 1h0p, 1h1p sectors
t1_scf    = Filter("Total SCF energy = ",          -107.057122389928, 1e-7)
t1_cc_tot = Filter("Total CCSDT energy = ",        -107.276335630558, 1e-7)
t1_ccsdt  = Filter("CCSDT correlation energy = ",    -0.219213240630, 1e-7)	
t1_e2  = Filter("@    2", 0.2196420989, 1e-7)
t1_e3  = Filter("@    3", 0.2200476676, 1e-7)
t1_e4  = Filter("@    4", 0.2208433071, 1e-7)
t1_e5  = Filter("@    5", 0.2208524433, 1e-7)
t1_e6  = Filter("@    6", 0.2220633019, 1e-7)
t1_e7  = Filter("@    7", 0.2343038379, 1e-7)
t1_e8  = Filter("@    8", 0.2343076953, 1e-7)
t1_e9  = Filter("@    9", 0.2367422911, 1e-7)
t1_ip  = Filter("Ionization potential 0h0p -> 1h0p =", 0.400389335127, 1e-7)
t1_ea  = Filter("Electron affinity    0h0p -> 0h1p =", 0.009458937022, 1e-7)

ret = Test("", "input-1h1p", filters=[t1_scf,t1_cc_tot,t1_ccsdt,
    t1_e2,t1_e3,t1_e4,t1_e5,t1_e6,t1_e7,t1_e8,t1_e9,t1_ip,t1_ea]).run(options="--no-clean")

execute("mv input-1h1p.test.out ccsdt_1h1p.out")
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(ret)

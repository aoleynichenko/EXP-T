#!/usr/bin/env python

# Test:
# (1) Ionization potentials of the magnesium atom
# (2) FSCC scheme: Mg -> Mg+ -> Mg2+ (sector 2h0p)
# (3) Basis set: cc-pVTZ
# (4) DC - relativistic, 4-comp

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
#symmetries = ['C1', 'Cs', 'C2', 'C2v', 'Cinfv', 'D2', 'Ci', 'D2h', 'C2h', 'Dinfh']
symmetries = ['Cs', 'Cinfv', 'Ci', 'D2h']

for sym in symmetries:
    dirac_inp = "TRA.inp"
    dirac_mol = "Mg-%s.mol" % (sym)
    execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")

    filter_scf  = Filter("SCF reference energy = ",          -199.902375277106, 1e-7)
    filter_mp2c = Filter("MP2 correlation energy = ",          -0.024911846482, 1e-7)
    filter_mp2  = Filter("Total MP2 energy = ",              -199.927287123588, 1e-7)
    filter_ccsdt= Filter("CCSDT correlation energy = ",        -0.036362573981, 1e-7)
    filter_ip1  = Filter("Ionization potential 0h0p -> 1h0p =", 0.276972272025, 1e-7)
    filter_ip2  = Filter("Ionization potential 1h0p -> 2h0p =", 0.542081735038, 1e-7)

    ret = Test(sym, "ccsdt.inp", filters=[filter_scf, filter_mp2c, filter_mp2, filter_ccsdt, filter_ip1, filter_ip2]).run()
    ret_codes.append(ret)

    execute("mv ccsdt.inp.test.out ccsdt_%s.out" % (sym))
    execute("rm -rf MRCONEE* MDCINT*")
    execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(1 if any(ret_codes) else 0)

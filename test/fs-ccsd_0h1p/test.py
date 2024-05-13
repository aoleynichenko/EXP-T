#!/usr/bin/env python

# Test:
# (1) HgH molecule, FSCC scheme: HgH+ -> HgH
# (2) hamiltionian: 2-comp gatchina ECP
# (3) symmetries: Cs, Cinfv
# (4) sector (0h,1p)

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

# all symmetries to be tested
symmetries = ['Cs', 'Cinfv']

for sym in symmetries:
    dirac_inp = "moltra.inp"
    dirac_mol = "HgH-%s.mol" % (sym)
    execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

    t1_scf  = Filter("Total SCF energy = ",        -153.399707040104, 1e-7)
    t1_mp2c = Filter("MP2 correlation energy = ",    -0.206580837880, 1e-7)
    t1_mp2  = Filter("Total MP2 energy = ",        -153.606287877984, 1e-7)
    t1_ccsd = Filter("CCSD correlation energy = ",   -0.219213152203, 1e-7)
    t1_e1  = Filter("@    1", -0.2848133503, 1e-7)
    t1_e2  = Filter("@    2", -0.1742116124, 1e-7)
    t1_e3  = Filter("@    3", -0.1588344865, 1e-7)
    
    ret = Test(sym, "ccsd.inp", filters=[
        t1_scf,t1_mp2,t1_mp2c,t1_ccsd,
        t1_e1,t1_e2,t1_e3
        ]).run()
    ret_codes.append(ret)
    
    execute("mv ccsd.inp.test.out ccsd_%s.out" % (sym))
    execute("rm -rf MRCONEE* MDCINT* MDPROP")
    execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(1 if any(ret_codes) else 0)




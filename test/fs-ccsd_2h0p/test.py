#!/usr/bin/env python

#
# Double ionization potential of the ethylene molecule.
# Experimental geometry, cc-pVDZ basis set is used.
# Electronic states of C2H4^{2+} are obtained in the 2h0p FS sector.
#
# See also for ADC calculations:
# E. Ohrendorf, H. Koeppel, L. S. Cederbaum, F. Tarantelli, A. Sgamellotti.
# Doubly ionized states of ethylene: Auger spectrum, potential energy
# surfaces and nuclear dynamics.
# J. Chem. Phys. 91, 1734 (1989).
# DOI: 10.1063/1.457080
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

execute(DIRAC_PATH + " --nobackup --noarch --inp=TRA.inp --mol=C2H4.mol --get=\"MRCONEE MDCINT\"")

filter_e1 = Filter("@    1", 1.1232538467, 1e-7)
filter_e2 = Filter("@    2", 1.1726385103, 1e-7)
filter_e3 = Filter("@    3", 1.1913832631, 1e-7)
filter_e4 = Filter("@    4", 1.2541192076, 1e-7)
filter_e5 = Filter("@    5", 1.2737227646, 1e-7)
filter_e6 = Filter("@    6", 1.2925700154, 1e-7)
filter_e7 = Filter("@    7", 1.3039773489, 1e-7)
filter_e8 = Filter("@    8", 1.3197457532, 1e-7)
filter_e9 = Filter("@    9", 1.3363235318, 1e-7)
filter_e10= Filter("@   10", 1.3452351132, 1e-7)

filter_ip1  = Filter("Ionization potential 0h0p -> 1h0p =", 0.382813727287, 1e-7)
filter_ip2  = Filter("Ionization potential 1h0p -> 2h0p =", 0.740440119418, 1e-7)

ret = Test("dip-ccsd", "dip-ccsd.inp", filters=[
	filter_ip1, filter_ip2,
	filter_e1, filter_e2, filter_e3, filter_e4, filter_e5,
	filter_e6, filter_e7, filter_e8, filter_e9, filter_e10
]).run()

execute("mv dip-ccsd.inp.test.out dip-ccsd.out")

execute("rm -rf MRCONEE* MDCINT*")
execute("rm -rf HINT VINT* modelvectors* HEFF")

sys.exit(ret)

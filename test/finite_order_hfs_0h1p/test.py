#!/usr/bin/env python

#
# Test:
# hyperfine structure constants in the 4s1/2, 4p1/2, 4p3/2 states
# of the potassium atom
# (finite-order calculation using the connected formula)
#
#
# 39K     exptl, MHz    ccsd, MHz     Omega
# S1/2      230.9         226.3        1/2
# P1/2       27.8          26.8        1/2
# P3/2        6.1           6.0        3/2
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

#
# obtain transformed integrals
#
dirac_inp = 'moltra.inp'
dirac_mol = "K.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

#
# CCSD + direct calculation of matrix elements
#
filter_hfs_1  = Filter("1 (1/2g+) ->  1 (1/2g+)", [    0.00,     0.00, None, None, 433.584560], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_2  = Filter("1 (1/2g-) ->  1 (1/2g-)", [    0.00,     0.00, None, None, 433.584560], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_3  = Filter("1 (1/2u+) ->  1 (1/2u+)", [12966.01, 12966.01, None, None,  51.393895], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_4  = Filter("1 (1/2u+) ->  2 (1/2u+)", [12966.01, 13024.86, None, None,   5.008749], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_5  = Filter("2 (1/2u+) ->  1 (1/2u+)", [13024.86, 12966.01, None, None,   5.007758], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_6  = Filter("2 (1/2u+) ->  2 (1/2u+)", [13024.86, 13024.86, None, None,  11.451785], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_7  = Filter("1 (1/2u-) ->  1 (1/2u-)", [12966.01, 12966.01, None, None,  51.393895], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_8  = Filter("1 (1/2u-) ->  2 (1/2u-)", [12966.01, 13024.86, None, None,   5.008749], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_9  = Filter("2 (1/2u-) ->  1 (1/2u-)", [13024.86, 12966.01, None, None,   5.007758], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_10 = Filter("2 (1/2u-) ->  2 (1/2u-)", [13024.86, 13024.86, None, None,  11.451785], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_11 = Filter("1 (3/2u+) ->  1 (3/2u+)", [13024.86, 13024.86, None, None,  34.355354], [1e-1, 1e-1, None, None, 1e-3])
filter_hfs_12 = Filter("1 (3/2u-) ->  1 (3/2u-)", [13024.86, 13024.86, None, None,  34.355354], [1e-1, 1e-1, None, None, 1e-3])

ret = Test("", "ccsd.inp", filters=[
	filter_hfs_1,filter_hfs_2,filter_hfs_3,filter_hfs_4,
	filter_hfs_5,filter_hfs_6,filter_hfs_7,filter_hfs_8,
	filter_hfs_9,filter_hfs_10,filter_hfs_11,filter_hfs_12
]).run(options="--no-clean")
execute("mv ccsd.inp.test.out ccsd_hfs.out")

#
# clean up
#
execute("rm -rf MRCONEE* MDCINT* MDPROP*")
execute("rm -rf scratch")
execute("rm -rf HINT VINT* modelvectors* HEFF*")

sys.exit(ret)


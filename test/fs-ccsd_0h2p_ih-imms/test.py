#!/usr/bin/env python

# Test: simple intermediate Hamiltonian

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

# return codes
ret_codes = []

dirac_inp = "TRA_C.inp"
dirac_mol = "C.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")

#
# manual setting of IH parameters
#

filter_e1  = Filter("@    1", -1.2956321346, 1e-7)   #       0.000000    96.2   1   0g
filter_e2  = Filter("@    2", -1.2955703221, 1e-7)   #      13.566267    96.2   3   0g 1g+ 1g-
filter_e3  = Filter("@    3", -1.2954474169, 1e-7)   #      40.540838    96.2   5   0g 1g+ 1g- 2g+ 2g-
filter_e4  = Filter("@    4", -1.2503517302, 1e-7)   #    9937.900071    94.0   5   0g 1g+ 1g- 2g+ 2g-
filter_e5  = Filter("@    5", -1.1925335317, 1e-7)   #   22627.527853    89.7   1   0g
filter_e6  = Filter("@    6", -1.0155964076, 1e-7)   #   61460.737951    98.6   1   0u
filter_e7  = Filter("@    7", -1.0155144649, 1e-7)   #   61478.722284    98.6   3   0u 1u+ 1u-
filter_e8  = Filter("@    8", -1.0153476400, 1e-7)   #   61515.336123    98.6   5   0u 1u+ 1u- 2u+ 2u-
filter_e9  = Filter("@    9", -1.0015985740, 1e-7)   #   64532.907329    98.6   3   0u 1u+ 1u-
filter_e10 = Filter("@   10", -0.9884834216, 1e-7)   #   67411.350562     0.0   3   0g 1g+ 1g-
filter_e11 = Filter("@   11", -0.9802138583, 1e-7)   #   69226.309921     0.0   3   0g 1g+ 1g-
filter_e12 = Filter("@   12", -0.9801242657, 1e-7)   #   69245.973210     0.0   5   0g 1g+ 1g- 2g+ 2g-
filter_e13 = Filter("@   13", -0.9799871954, 1e-7)   #   69276.056670     0.0   7   0g 1g+ 1g- 2g+ 2g- 3g+ 3g-
filter_e14 = Filter("@   14", -0.9721225068, 1e-7)   #   71002.156296     0.0   3   0g 1g+ 1g-
filter_e15 = Filter("@   15", -0.9685688606, 1e-7)   #   71782.091497     0.9   1   0g
filter_e16 = Filter("@   16", -0.9685194984, 1e-7)   #   71792.925242     0.9   3   0g 1g+ 1g-
filter_e17 = Filter("@   17", -0.9684293190, 1e-7)   #   71812.717325     0.9   5   0g 1g+ 1g- 2g+ 2g-
filter_e18 = Filter("@   18", -0.9560018628, 1e-7)   #   74540.228710     1.8   5   0g 1g+ 1g- 2g+ 2g-
filter_e19 = Filter("@   19", -0.9411826093, 1e-7)   #   77792.678904     3.0   1   0g

filter_list = [
  filter_e1, filter_e2, filter_e3, filter_e4, filter_e5,
  filter_e6, filter_e7, filter_e8, filter_e9, filter_e10,
  filter_e11, filter_e12, filter_e13, filter_e14, filter_e15,
  filter_e16, filter_e17, filter_e18, filter_e19
]

ret = Test('main: 2p^2, 3s2p (manual IH)', "ccsd_ih_manual.inp", filters=filter_list).run()
ret_codes.append(ret)
execute("rm -rf scratch")

#
# automatic determination of IH parameters
#

filter_e1  = Filter("@    1", -1.2956343845, 1e-7)   #        0.000000    96.2   1   0g
filter_e2  = Filter("@    2", -1.2955725720, 1e-7)   #       13.566285    96.2   3   0g 1g+ 1g-
filter_e3  = Filter("@    3", -1.2954496668, 1e-7)   #       40.540866    96.2   5   0g 1g+ 1g- 2g+ 2g-
filter_e4  = Filter("@    4", -1.2503552768, 1e-7)   #     9937.615469    94.0   5   0g 1g+ 1g- 2g+ 2g-
filter_e5  = Filter("@    5", -1.1925510313, 1e-7)   #    22624.180957    89.7   1   0g
filter_e6  = Filter("@    6", -1.0155946503, 1e-7)   #    61461.617433    98.6   1   0u
filter_e7  = Filter("@    7", -1.0155127066, 1e-7)   #    61479.601998    98.6   3   0u 1u+ 1u-
filter_e8  = Filter("@    8", -1.0153458800, 1e-7)   #    61516.216213    98.6   5   0u 1u+ 1u- 2u+ 2u-
filter_e9  = Filter("@    9", -1.0015950164, 1e-7)   #    64534.181929    98.6   3   0u 1u+ 1u-
filter_e10 = Filter("@   10", -0.9884834128, 1e-7)   #    67411.846295     0.0   3   0g 1g+ 1g-
filter_e11 = Filter("@   11", -0.9802125073, 1e-7)   #    69227.100232     0.0   3   0g 1g+ 1g-
filter_e12 = Filter("@   12", -0.9801229175, 1e-7)   #    69246.762912     0.0   5   0g 1g+ 1g- 2g+ 2g-
filter_e13 = Filter("@   13", -0.9799858524, 1e-7)   #    69276.845226     0.0   7   0g 1g+ 1g- 2g+ 2g- 3g+ 3g-
filter_e14 = Filter("@   14", -0.9721224889, 1e-7)   #    71002.654041     0.0   3   0g 1g+ 1g-
filter_e15 = Filter("@   15", -0.9685585861, 1e-7)   #    71784.840290     0.9   1   0g
filter_e16 = Filter("@   16", -0.9685092327, 1e-7)   #    71795.672113     0.9   3   0g 1g+ 1g-
filter_e17 = Filter("@   17", -0.9684190513, 1e-7)   #    71815.464632     0.9   5   0g 1g+ 1g- 2g+ 2g-
filter_e18 = Filter("@   18", -0.9559857281, 1e-7)   #    74544.263674     1.8   5   0g 1g+ 1g- 2g+ 2g-
filter_e19 = Filter("@   19", -0.9411474486, 1e-7)   #    77800.889587     3.0   1   0g

filter_list = [
  filter_e1, filter_e2, filter_e3, filter_e4, filter_e5,
  filter_e6, filter_e7, filter_e8, filter_e9, filter_e10,
  filter_e11, filter_e12, filter_e13, filter_e14, filter_e15,
  filter_e16, filter_e17, filter_e18, filter_e19
]

ret = Test('main: 2p^2, 3s2p (auto IH)', "ccsd_ih_auto.inp", filters=filter_list).run()
ret_codes.append(ret)
execute("rm -rf MRCONEE* MDCINT* scratch")

sys.exit(1 if any(ret_codes) else 0)


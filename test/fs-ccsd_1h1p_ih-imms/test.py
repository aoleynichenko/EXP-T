#!/usr/bin/env python

# Test: simple intermediate Hamiltonian in the 1h1p sector

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

dirac_inp = "TRA.inp"
dirac_mol = "CO.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")


filter_e1  = Filter("@    1", 0.0000000000, 1e-7)    #      0.000000   100.0   1   0
filter_e2  = Filter("@    2", 0.2294092055, 1e-7)    #  50349.500799    99.9   1   0
filter_e3  = Filter("@    3", 0.2294093543, 1e-7)    #  50349.533462    99.9   1   0
filter_e4  = Filter("@    4", 0.2296243638, 1e-7)    #  50396.722596    99.9   2   1+ 1-
filter_e5  = Filter("@    5", 0.2298404121, 1e-7)    #  50444.139714    99.9   2   2+ 2-
filter_e6  = Filter("@    6", 0.3058523610, 1e-7)    #  67126.834177    90.9   2   1+ 1-
filter_e7  = Filter("@    7", 0.3058564311, 1e-7)    #  67127.727469    90.9   1   0
filter_e8  = Filter("@    8", 0.3147979195, 1e-7)    #  69090.157340    97.0   2   1+ 1-
filter_e9  = Filter("@    9", 0.3398211562, 1e-7)    #  74582.122981    92.3   2   3+ 3-
filter_e10 = Filter("@   10", 0.3399374351, 1e-7)    #  74607.643264    92.4   2   2+ 2-
filter_e11 = Filter("@   11", 0.3400658589, 1e-7)    #  74635.829025    92.4   2   1+ 1-
filter_e12 = Filter("@   12", 0.3640361305, 1e-7)    #  79896.695553    93.5   1   0
filter_e13 = Filter("@   13", 0.3640395407, 1e-7)    #  79897.444003    93.5   2   1+ 1-
filter_e14 = Filter("@   14", 0.3680468002, 1e-7)    #  80776.935796    94.0   1   0
filter_e15 = Filter("@   15", 0.3742394991, 1e-7)    #  82136.076107    95.0   2   2+ 2-
filter_e16 = Filter("@   16", 0.3906018700, 1e-7)    #  85727.201419    99.7   1   0
filter_e17 = Filter("@   17", 0.3906019012, 1e-7)    #  85727.208266    99.7   2   1+ 1-
filter_e18 = Filter("@   18", 0.4099573687, 1e-7)    #  89975.242373    99.9   1   0
filter_e19 = Filter("@   19", 0.4225774171, 1e-7)    #  92745.022843    99.6   3   0 1+ 1-
filter_e20 = Filter("@   20", 0.4289705881, 1e-7)    #  94148.161683    99.6   1   0
filter_e21 = Filter("@   21", 0.4317241627, 1e-7)    #  94752.501467    97.4   1   0
filter_e22 = Filter("@   22", 0.4317242134, 1e-7)    #  94752.512590    97.4   1   0
filter_e23 = Filter("@   23", 0.4317288369, 1e-7)    #  94753.527327    97.4   2   1+ 1-
filter_e24 = Filter("@   24", 0.4317334831, 1e-7)    #  94754.547054    97.5   2   2+ 2-
filter_e25 = Filter("@   25", 0.4383585793, 1e-7)    #  96208.587604    98.3   2   1+ 1-
filter_e26 = Filter("@   26", 0.4546895245, 1e-7)    #  99792.815779     5.0   1   0
filter_e27 = Filter("@   27", 0.4546900310, 1e-7)    #  99792.926944     5.0   1   0
filter_e28 = Filter("@   28", 0.4548374342, 1e-7)    #  99825.278209     5.0   2   1+ 1-
filter_e29 = Filter("@   29", 0.4549856663, 1e-7)    #  99857.811381     5.0   2   2+ 2-

filter_list = [
  filter_e1, filter_e2, filter_e3, filter_e4, filter_e5,
  filter_e6, filter_e7, filter_e8, filter_e9, filter_e10,
  filter_e11, filter_e12, filter_e13, filter_e14, filter_e15,
  filter_e16, filter_e17, filter_e18, filter_e19, filter_e20,
  filter_e21, filter_e22, filter_e23, filter_e24, filter_e25,
  filter_e26, filter_e27, filter_e28, filter_e29,
]

ret = Test('auto IH', "ccsd.inp", filters=filter_list).run()
execute("mv ccsd.inp.test.out ccsd.out")
execute("rm -rf MRCONEE* MDCINT* scratch")

sys.exit(ret)


#!/usr/bin/env python
#
# Test: one-electron property operator (magnetic hyperfine) included
# through the 'oneprop' keyword
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> oneprop/Li atom model-space properties')

dirac_inp = "TRA.inp"
dirac_mol = "Li.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT MDPROP\"")

# Transition dipole moments (the 'mstdm' option)
Test("mstdm", "cc_mstdm.inp", filters=[
    Filter("Total SCF energy = ",         -7.23720113394432, 1e-7),
    Filter("MP2 correlation energy = ",     -0.037127930756, 1e-7),
    Filter("Total MP2 energy = ",         -7.27432906470000, 1e-7),
    Filter("CCSD correlation energy = ",    -0.041240544312, 1e-7),
    Filter("@    1", -0.1979693711, 1e-7),
    Filter("@    2", -0.1296541879, 1e-7),
    Filter("@    3", -0.1296513218, 1e-7),
    Filter("1 (  1E) ->  2 (  1E)", [    0.00, 14993.45, 1.942321, 0.171818, 0.000000, 1.373428, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  1E) ->  3 (  1E)", [    0.00, 14994.08, 2.320942, 0.245342, 0.000000, 1.772007, 1.498920], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  1E) ->  4 (  1E)", [    0.00, 14994.08, 1.464459, 0.097678, 0.000000, 1.199532, 0.840098], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("2 (  1E) ->  1 (  1E)", [14993.45,     0.00, 1.942321, 0.171818, 0.000000, 1.373428, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("3 (  1E) ->  1 (  1E)", [14994.08,     0.00, 2.399520, 0.262236, 0.000000, 1.577359, 1.808213], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("4 (  1E) ->  1 (  1E)", [14994.08,     0.00, 1.516990, 0.104812, 0.000000, 0.821270, 1.275451], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  1E) ->  2 (  2E)", [    0.00, 14993.45, 1.373428, 0.085909, 1.373428, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  1E) ->  3 (  2E)", [    0.00, 14994.08, 0.560899, 0.014329, 0.560899, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  1E) ->  4 (  2E)", [    0.00, 14994.08, 1.725165, 0.135552, 1.725165, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("2 (  1E) ->  1 (  2E)", [14993.45,     0.00, 1.373428, 0.085909, 1.373428, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("3 (  1E) ->  1 (  2E)", [14994.08,     0.00, 0.524766, 0.012542, 0.524766, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("4 (  1E) ->  1 (  2E)", [14994.08,     0.00, 1.931791, 0.169967, 1.931791, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  2E) ->  2 (  1E)", [    0.00, 14993.45, 1.373428, 0.085909, 1.373428, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  2E) ->  3 (  1E)", [    0.00, 14994.08, 0.521765, 0.012399, 0.521765, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  2E) ->  4 (  1E)", [    0.00, 14994.08, 1.874665, 0.160063, 1.874665, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("2 (  2E) ->  1 (  1E)", [14993.45,     0.00, 1.373428, 0.085909, 1.373428, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("3 (  2E) ->  1 (  1E)", [14994.08,     0.00, 1.089236, 0.054037, 1.089236, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("4 (  2E) ->  1 (  1E)", [14994.08,     0.00, 2.269609, 0.234610, 2.269609, 0.000000, 0.000000], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  2E) ->  2 (  2E)", [    0.00, 14993.45, 1.942321, 0.171818, 0.000000, 1.373428, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  2E) ->  3 (  2E)", [    0.00, 14994.08, 2.311796, 0.243413, 0.000000, 1.846029, 1.391609], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("1 (  2E) ->  4 (  2E)", [    0.00, 14994.08, 1.637930, 0.122190, 0.000000, 1.552927, 0.520802], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("2 (  2E) ->  1 (  2E)", [14993.45,     0.00, 1.942321, 0.171818, 0.000000, 1.373428, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("3 (  2E) ->  1 (  2E)", [14994.08,     0.00, 2.691322, 0.329895, 0.000000, 1.423896, 2.283798], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5]),
    Filter("4 (  2E) ->  1 (  2E)", [14994.08,     0.00, 1.810669, 0.149321, 0.000000, 0.737208, 1.653797], [1e-2,1e-2,1e-5,1e-5,1e-5,1e-5,1e-5])
]).run(options="--no-clean")

# Transition dipole moments -- X component
# (the 'msprop' option)
Test("dipole moment x", "cc_xdiplen.inp", filters=[
    Filter("1 (1E  ) ->  2 (2E  )", [    0.00, 14993.45,  1.298249, -0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (1E  ) ->  3 (2E  )", [    0.00, 14994.08,  0.090774, -0.553505, 0.560899], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (1E  ) ->  4 (2E  )", [    0.00, 14994.08,  1.667769, -0.441293, 1.725165], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  1 (2E  )", [14993.45,     0.00, -1.298249,  0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  1 (2E  )", [14994.08,     0.00,  0.488066,  0.192797, 0.524766], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  1 (2E  )", [14994.08,     0.00, -1.903168,  0.331310, 1.931791], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  2 (1E  )", [    0.00, 14993.45, -1.298249, -0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  3 (1E  )", [    0.00, 14994.08,  0.123739, -0.506880, 0.521765], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  4 (1E  )", [    0.00, 14994.08, -1.823522, -0.434899, 1.874665], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  1 (1E  )", [14993.45,     0.00,  1.298249,  0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  1 (1E  )", [14994.08,     0.00, -1.089223,  0.005174, 1.089236], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  1 (1E  )", [14994.08,     0.00,  2.257356,  0.235521, 2.269609], [1e-2,1e-2,1e-5,1e-5,1e-5])
]).run("--no-clean")

# Transition dipole moments -- Y component
Test("dipole moment y", "cc_ydiplen.inp", filters=[
    Filter("1 (1E  ) ->  2 (1E  )", [    0.00, 14993.45, -1.298249,  0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (1E  ) ->  3 (1E  )", [    0.00, 14994.08,  1.714009,  0.449645, 1.772007], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (1E  ) ->  4 (1E  )", [    0.00, 14994.08,  1.194041, -0.114637, 1.199532], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  1 (1E  )", [14993.45,     0.00, -1.298249, -0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  1 (1E  )", [14994.08,     0.00,  1.539669, -0.342756, 1.577359], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  1 (1E  )", [14994.08,     0.00,  0.821221, -0.008940, 0.821270], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  2 (2E  )", [    0.00, 14993.45, -1.298249, -0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  3 (2E  )", [    0.00, 14994.08,  1.804475, -0.389478, 1.846029], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  4 (2E  )", [    0.00, 14994.08,  1.544422,  0.162303, 1.552927], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  1 (2E  )", [14993.45,     0.00, -1.298249,  0.448167, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  1 (2E  )", [14994.08,     0.00,  1.395579,  0.282557, 1.423896], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  1 (2E  )", [14994.08,     0.00,  0.735120, -0.055446, 0.737208], [1e-2,1e-2,1e-5,1e-5,1e-5])
]).run("--no-clean")

# Transition dipole moments -- Z component
Test("dipole moment z", "cc_zdiplen.inp", filters=[
    Filter("1 (1E  ) ->  2 (1E  )", [    0.00, 14993.45, -0.448167, -1.298249, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (1E  ) ->  3 (1E  )", [    0.00, 14994.08, -0.025625, -1.498701, 1.498920], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (1E  ) ->  4 (1E  )", [    0.00, 14994.08,  0.668014,  0.509433, 0.840098], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  1 (1E  )", [14993.45,     0.00, -0.448167,  1.298249, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  1 (1E  )", [14994.08,     0.00, -0.110073,  1.804859, 1.808213], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  1 (1E  )", [14994.08,     0.00,  0.922254, -0.881034, 1.275451], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  2 (2E  )", [    0.00, 14993.45, -0.448167,  1.298249, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  3 (2E  )", [    0.00, 14994.08,  0.009847,  1.391575, 1.391609], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  4 (2E  )", [    0.00, 14994.08,  0.519232, -0.040414, 0.520802], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  1 (2E  )", [14993.45,     0.00, -0.448167, -1.298249, 1.373428], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  1 (2E  )", [14994.08,     0.00, -0.393459, -2.249650, 2.283798], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  1 (2E  )", [14994.08,     0.00,  1.151757,  1.186802, 1.653797], [1e-2,1e-2,1e-5,1e-5,1e-5]),
]).run("--no-clean")

# Hyperfine structure matrix elements -- X component
Test("hfs constant x", "cc_xhfs.inp", filters=[
    Filter("1 (1E  ) ->  1 (1E  )", [    0.00,     0.00,-64.874846,  0.000000,64.874846], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  2 (1E  )", [14993.45, 14993.45,  7.426794,  0.000000, 7.426794], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  3 (1E  )", [14993.45, 14994.08, -0.032734, -0.351088, 0.352611], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  4 (1E  )", [14993.45, 14994.08, -1.260803,  0.124312, 1.266917], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  2 (1E  )", [14994.08, 14993.45,  0.269268,  0.230790, 0.354640], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  3 (1E  )", [14994.08, 14994.08, -4.206096, -0.351954, 4.220796], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  4 (1E  )", [14994.08, 14994.08, -1.269284, -0.887748, 1.548928], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  2 (1E  )", [14994.08, 14993.45, -1.288840, -0.208050, 1.305524], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  3 (1E  )", [14994.08, 14994.08, -0.106375,  1.583431, 1.587000], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  4 (1E  )", [14994.08, 14994.08,  1.236061,  0.351954, 1.285192], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  1 (2E  )", [    0.00,     0.00, 64.874846,  0.000000,64.874846], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  2 (2E  )", [14993.45, 14993.45, -7.426794,  0.000000, 7.426794], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  3 (2E  )", [14993.45, 14994.08,  0.180050, -0.333567, 0.379058], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  4 (2E  )", [14993.45, 14994.08,  1.162716,  0.085881, 1.165883], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  2 (2E  )", [14994.08, 14993.45, -0.694674,  0.243503, 0.736115], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  3 (2E  )", [14994.08, 14994.08,  4.606220, -0.949989, 4.703162], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  4 (2E  )", [14994.08, 14994.08,  2.856598, -0.770393, 2.958658], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  2 (2E  )", [14994.08, 14993.45,  1.493977, -0.347348, 1.533825], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  3 (2E  )", [14994.08, 14994.08, -0.527885,  1.933606, 2.004369], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  4 (2E  )", [14994.08, 14994.08, -1.636185,  0.949989, 1.891978], [1e-2,1e-2,1e-5,1e-5,1e-5])
]).run("--no-clean")

# Hyperfine structure matrix elements -- Y component
Test("hfs constant y", "cc_yhfs.inp", filters=[
    Filter("1 (1E  ) ->  1 (2E  )", [    0.00,     0.00, 51.059103, 40.021416,64.874846], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  2 (2E  )", [14993.45, 14993.45,  7.426794, -0.000001, 7.426794], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  3 (2E  )", [14993.45, 14994.08,  1.164670,  0.005429, 1.164683], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  4 (2E  )", [14993.45, 14994.08,  1.103027,  0.069146, 1.105192], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  2 (2E  )", [14994.08, 14993.45,  0.972450,  0.124903, 0.980439], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  3 (2E  )", [14994.08, 14994.08,  0.267749, -1.152614, 1.183304], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  4 (2E  )", [14994.08, 14994.08,  2.457022,  0.928209, 2.626505], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  2 (2E  )", [14994.08, 14993.45,  0.712938, -0.203744, 0.741479], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  3 (2E  )", [14994.08, 14994.08,  1.650339,  1.309549, 2.106784], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  4 (2E  )", [14994.08, 14994.08, -1.219854, -0.969141, 1.557972], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  1 (1E  )", [    0.00,     0.00, 51.059103,-40.021416,64.874846], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  2 (1E  )", [14993.45, 14993.45,  7.426794,  0.000001, 7.426794], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  3 (1E  )", [14993.45, 14994.08,  1.098749,  0.009770, 1.098793], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  4 (1E  )", [14993.45, 14994.08,  0.904769,  0.052743, 0.906305], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  2 (1E  )", [14994.08, 14993.45,  0.824866, -0.258414, 0.864397], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  3 (1E  )", [14994.08, 14994.08, -0.743710,  1.747420, 1.899100], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  4 (1E  )", [14994.08, 14994.08,  3.546802, -1.777853, 3.967439], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  2 (1E  )", [14994.08, 14993.45,  0.703957,  0.224666, 0.738938], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  3 (1E  )", [14994.08, 14994.08,  2.150076, -1.963048, 2.911423], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  4 (1E  )", [14994.08, 14994.08, -2.405185,  2.081532, 3.180832], [1e-2,1e-2,1e-5,1e-5,1e-5])
]).run("--no-clean")

# Hyperfine structure matrix elements -- Z component
Test("hfs constant z", "cc_zhfs.inp", filters=[
    Filter("1 (1E  ) ->  1 (2E  )", [    0.00,     0.00, 40.021416,-51.059103,64.874846], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  2 (2E  )", [14993.45, 14993.45, -0.000000,  7.426794, 7.426794], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  3 (2E  )", [14993.45, 14994.08,  0.339003, -0.984621, 1.041346], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (1E  ) ->  4 (2E  )", [14993.45, 14994.08, -0.016733,  0.059688, 0.061990], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  2 (2E  )", [14994.08, 14993.45,  0.355700, -1.241718, 1.291661], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  3 (2E  )", [14994.08, 14994.08, -0.986222,  0.164788, 0.999895], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (1E  ) ->  4 (2E  )", [14994.08, 14994.08, -0.076613, -1.446027, 1.448055], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  2 (2E  )", [14994.08, 14993.45, -0.411794,  0.575900, 0.707980], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  3 (2E  )", [14994.08, 14994.08, -0.043418, -2.689578, 2.689928], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (1E  ) ->  4 (2E  )", [14994.08, 14994.08,  0.246379, -3.884684, 3.892489], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("1 (2E  ) ->  1 (1E  )", [    0.00,     0.00, 40.021416, 51.059103,64.874846], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  2 (1E  )", [14993.45, 14993.45, -0.000000, -7.426794, 7.426794], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  3 (1E  )", [14993.45, 14994.08,  0.341324,  1.066016, 1.119327], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("2 (2E  ) ->  4 (1E  )", [14993.45, 14994.08, -0.177055, -0.356033, 0.397628], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  2 (1E  )", [14994.08, 14993.45,  0.501924,  1.519540, 1.600290], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  3 (1E  )", [14994.08, 14994.08, -0.934097, -1.116830, 1.455969], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("3 (2E  ) ->  4 (1E  )", [14994.08, 14994.08,  0.421746,  0.631752, 0.759592], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  2 (1E  )", [14994.08, 14993.45, -0.572016, -0.790018, 0.975362], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  3 (1E  )", [14994.08, 14994.08,  0.194002,  2.739081, 2.745943], [1e-2,1e-2,1e-5,1e-5,1e-5]),
    Filter("4 (2E  ) ->  4 (1E  )", [14994.08, 14994.08,  0.068475,  3.939540, 3.940135], [1e-2,1e-2,1e-5,1e-5,1e-5])
]).run("--no-clean")

execute("rm -rf MRCONEE* MDCINT* MDPROP* scratch")

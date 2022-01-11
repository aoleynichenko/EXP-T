#!/usr/bin/env python

# Test: program for rovibrational levels of diatomic molecules

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter, execute, DIRAC_PATH

DIATOMIC_PATH = "expt_diatomic.x"

print('>>> diatomic vib-rot levels')

# AcOH+ molecule: Frank-Condon factors

acoh_filters = [
    Filter("$    0   0   0   0", [ 327.5900, 7607.4633, 7279.8734, None, 0.889149], [1e-4,1e-4,1e-4,None,1e-6]),
    Filter("$    0   0   1   0", [ 979.7793, 7607.4633, 6627.6841, None, 0.105888], [1e-4,1e-4,1e-4,None,1e-6]),
    Filter("$    0   0   2   0", [1628.1902, 7607.4633, 5979.2731, None, 0.004855], [1e-4,1e-4,1e-4,None,1e-6]),
    Filter("$    0   0   3   0", [2273.1511, 7607.4633, 5334.3122, None, 0.000107], [1e-4,1e-4,1e-4,None,1e-6])
]
Test('AcOH+ linear Frank-Condon factors', "AcOH+_FCF.inp", acoh_filters, binary=DIATOMIC_PATH).run()

# KCs molecule: Numerov solver

kcs_num_filters = [
    Filter("@    0   0", [  34.1348,  4.289820], [1e-4,1e-5]),
    Filter("@    0  96", [4058.5663, 10.541161], [1e-4,1e-5]),
    Filter("@    0  97", [4061.5729, 11.060702], [1e-4,1e-5]),
    Filter("@    0  98", [4063.9765, 11.645033], [1e-4,1e-5])
]
Test('KCs ground state, Numerov solver', "KCs_numerov.inp", kcs_num_filters, binary=DIATOMIC_PATH).run()

# KCs molecule: finite-difference (FD2) solver

kcs_fd2_filters = [
    Filter("@    0   0", [  34.1186,  4.289815], [1e-4,1e-5]),
    Filter("@    0  96", [4049.8205,  9.589263], [1e-4,1e-5]),
    Filter("@    0  97", [4054.3493,  9.991194], [1e-4,1e-5]),
    Filter("@    0  98", [4058.1019, 10.450892], [1e-4,1e-5]),
    Filter("@    0  99", [4061.1821, 10.958611], [1e-4,1e-5]),  # note: this is a wrong level!
    Filter("@    0 100", [4063.6584, 11.553705], [1e-4,1e-5])   # note: this is a wrong level!
]
Test('KCs ground state, FD2 solver', "KCs_fd2.inp", kcs_fd2_filters, binary=DIATOMIC_PATH).run()


# KCs molecule: finite-difference (FD2) solver + mapping

kcs_fd2_map_filters = [
    Filter("@    0   0", [  34.1323,  4.289820], [1e-4,1e-5]),
    Filter("@    0  96", [4056.5147, 10.254932], [1e-4,1e-5]),
    Filter("@    0  97", [4059.8778, 10.734295], [1e-4,1e-5]),
    Filter("@    0  98", [4062.6189, 11.295417], [1e-4,1e-5])
]
Test('KCs ground state, FD2 solver + mapping', "KCs_fd2_mapping.inp", kcs_fd2_map_filters, binary=DIATOMIC_PATH).run()

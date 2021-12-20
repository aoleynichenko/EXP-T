#!/usr/bin/env python
#
# Test: CCSDT-1 in 0h0p, 0h1p, 0h2p sectors
#
# Following the paper of S. R. Hughes and U. Kaldor:
# S. R. Hughes, U. Kaldor, Chem. Phys. Lett. V. 204, P. 1993.
#
# FSCC scheme: Ne6+ (0h0p) -> Ne5+ (0h1p) -> Ne4+ (0h2p)
# Basis set: cont-d 7s7p4d3f adapted from ANO-RCC of Widmark, Malmqvist and Roos
# It seems that the basis set employed by Hughes and Kaldor can be slighlty 
# different, however, results are almost identical (the difference is < 0.01 eV)
# Hamiltonian: Schroedinger
# Symmetry: C2v
# Model space: 6 virtual spinors
# All 4 electrons are correlated
#
# Features to be tested:
#  - CCSDT-1 in the 0h1p and 0h2p FS sectors
#  - reuse of cluster amplitudes
#  - storage of amplitudes on disk and in RAM
#
# Stages of the calculation:
#  1. Sector 0h1p: test Ecorr and 1st electron affinity
#  2.
#

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
from minitest import Test, Filter, execute, DIRAC_PATH

print('>>> fs-ccsdt-1\' in 0h1p 0h2p sectors (Ne6+,Ne5+,Ne4+)')

# Hartree-Fock and integral transformation
dirac_inp = "TRA.inp"
dirac_mol = "Ne.mol"
execute(DIRAC_PATH + " --nobackup --noarch --inp=" + dirac_inp + " --mol=" + dirac_mol + " --get=\"MRCONEE MDCINT\"")

t1_scf         = Filter("Total SCF energy = ",             -110.087095488765, 1e-8)
t1_ccsdt1_corr = Filter("CCSDT-1b correlation energy = ",  -0.145396194221,   1e-8)
t1_ccsdt1      = Filter("Total CCSDT-1b energy = ",        -110.232491682986, 1e-8)
t1_ccsdt1_ea   = Filter("@    1 ", -5.7977858157, 1e-7)  # for 0h1p sector
t1_ccsdt1_e1   = Filter("@    1", -10.4350666261, 1e-7)  # for 0h2p sector
t1_ccsdt1_e2   = Filter("@    2", -10.2980224182, 1e-7)
t1_ccsdt1_e3   = Filter("@    3", -10.1513903113, 1e-7)

# I. All diagrams are stored in RAM
# I.1. 0h0p and 0h1p sectors
Test("0h0p+0h1p", "input_0h0p_0h1p", filters=[t1_scf,t1_ccsdt1_corr,t1_ccsdt1,t1_ccsdt1_ea]).run("--no-clean")
# I.2. 0h2p sector; reuse amplitudes from 0h0p and 0h1p
Test("0h2p, reuse 0h0p,0h1p", "input_0h2p", filters=[t1_scf,t1_ccsdt1_corr,t1_ccsdt1,t1_ccsdt1_e1,t1_ccsdt1_e2,t1_ccsdt1_e3]).run("--no-clean")
# I.2. 0h0p, 0h1p, 0h2p sector; reuse all amplitudes
Test("0h2p, reuse all", "input_0h2p_reuse_all", filters=[t1_scf,t1_ccsdt1_corr,t1_ccsdt1,t1_ccsdt1_e1,t1_ccsdt1_e2,t1_ccsdt1_e3]).run("--no-clean")

# II. 6-rank and PPPP diagrams are stored on disk
# II.1. 0h0p and 0h1p sectors
#Test("0h0p+0h1p", "input_0h0p_0h1p_disk", filters=[t1_scf,t1_ccsdt1_corr,t1_ccsdt1,t1_ccsdt1_ea]).run("--no-clean")
# II.2. 0h2p sector; reuse amplitudes from 0h0p and 0h1p
#Test("0h2p, reuse 0h0p,0h1p", "input_0h2p_disk", filters=[t1_scf,t1_ccsdt1_corr,t1_ccsdt1,t1_ccsdt1_e1,t1_ccsdt1_e2,t1_ccsdt1_e3]).run("--no-clean")
# II.2. 0h0p, 0h1p, 0h2p sector; reuse all amplitudes
#Test("0h2p, reuse all", "input_0h2p_reuse_all_disk", filters=[t1_scf,t1_ccsdt1_corr,t1_ccsdt1,t1_ccsdt1_e1,t1_ccsdt1_e2,t1_ccsdt1_e3]).run()


execute("rm -rf MRCONEE* MDCINT*")
execute("rm -rf HINT VINT* modelvectors* HEFF scratch")

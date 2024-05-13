#!/bin/bash

pam --inp=moltra --mol=Rb --noarch --get="MRCONEE MDCINT MDPROP" --gb=2 --ag=5

# field minus
expt.x --no-clean ccsd_F-.inp > ccsd_F-.out
mv scratch/HEFF HEFF1

# field plus
expt.x --no-clean ccsd_F+.inp > ccsd_F+.out
mv scratch/HEFF HEFF2

# calculate TDMs using the finite-difference off-diagonal formula
heffman.x < ff_tdm.inp > ff_tdm.out


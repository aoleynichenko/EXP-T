#!/bin/bash

pam --inp=moltra --mol=Xe2 --noarch --get="MRCONEE MDCINT MDPROP"
expt.x --no-clean ccsd.inp > ccsd_2.out
expt.x --no-clean ccsd_tdm.inp > ccsd_tdm_2.out



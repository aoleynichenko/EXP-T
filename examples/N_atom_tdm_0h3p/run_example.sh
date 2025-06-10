#!/bin/bash

pam --inp=moltra --mol=N --noarch --get="MRCONEE MDCINT MDPROP"
expt.x --no-clean ccsd.inp > ccsd.out
expt.x --no-clean ccsd_tdm.inp > ccsd_tdm.out



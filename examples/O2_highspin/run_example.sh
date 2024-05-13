#!/bin/bash

pam --inp=moltra --mol=O2 --noarch --get="MRCONEE MDCINT MDPROP"
expt.x --no-clean ccsd.inp > ccsd.out




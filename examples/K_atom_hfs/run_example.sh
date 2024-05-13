#!/bin/bash

pam --gb=1 --ag=2 --inp=moltra --mol=K --noarch --get="MRCONEE MDCINT MDPROP"
expt.x --no-clean ccsd.inp > ccsd.out



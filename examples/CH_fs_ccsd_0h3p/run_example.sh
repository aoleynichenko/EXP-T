#!/bin/bash

pam --inp=moltra --mol=CH --noarch --get="MRCONEE MDCINT MDPROP"
expt.x --no-clean ccsd.inp > ccsd.out




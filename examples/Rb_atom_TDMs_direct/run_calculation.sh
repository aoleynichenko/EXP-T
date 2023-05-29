#!/bin/bash

pam --inp=TRA --mol=Rb --noarch --get="MRCONEE MDCINT MDPROP"

expt.x --no-clean ccsd.inp > ccsd.out
expt.x --no-clean tdm.inp  > tdm.out




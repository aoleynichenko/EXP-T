#!/bin/bash

expt2pam.x HgH.exp > HgH.mol
pam --inp=moltra --mol=HgH --noarch --get="MRCONEE MDCINT"
expt.x ccsd.inp > ccsd.out



#!/bin/bash

pam --inp=moltra --mol=CO --noarch --get="MRCONEE MDCINT"

mv MRCONEE MRCONEE-C1
mv MDCINT MDCINT-C1

expt.x ccsd.inp > ccsd.out



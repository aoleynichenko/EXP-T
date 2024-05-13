#!/bin/bash

pam --inp=moltra --mol=Pb --noarch --get="MRCONEE MDCINT" --gb=2 --ag=5
expt.x ccsd.inp > ccsd.out


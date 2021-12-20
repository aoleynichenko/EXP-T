#!/usr/bin/env python

import sys
import os

test_dirs = ['H2O_ccsdt-1', 'CO_ccsdt-3', 'Be_ccsdt', 'C+_ccsdt', 'C_atom_ccsdt', 'N_atom_ccsdt', 'Ne_ip_ccsdt', 'N2_exc_ccsdt_musial']

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


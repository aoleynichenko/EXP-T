#!/usr/bin/env python

import sys
import os

test_dirs = ['N2_exc_musial', 'Ne_uncontracted', 'LiF-FSCC', 'CO_vdz']

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


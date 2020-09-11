#!/usr/bin/env python

import sys
import os

test_dirs = ['Ne_uncontracted', 'LiF-FSCC']

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


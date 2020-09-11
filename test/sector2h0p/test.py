#!/usr/bin/env python

import sys
import os

test_dirs = ['O-atom']

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


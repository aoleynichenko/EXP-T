#!/usr/bin/env python

import sys
import os

test_dirs = ['LiNa', 'Hg']

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


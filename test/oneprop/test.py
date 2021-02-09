#!/usr/bin/env python

import sys
import os

#test_dirs = ['K-hfs', 'Hg-mstdm-0h2p', 'Hg-mstdm-1h1p']
test_dirs = ['K-hfs', 'Li-msprop-x', 'Li-msprop-z']

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


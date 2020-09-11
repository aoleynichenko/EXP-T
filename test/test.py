#!/usr/bin/env python

import sys
import os


# different sets of tests
# can be enabled by the command line parameter: python test.py [<test-set-1> ...]
do_test_quick = False
do_test_omp  = False
do_test_cuda = False
do_test_full = False


if len(sys.argv) > 1:
	for arg in sys.argv[1:]:
		if arg == 'quick':
			do_test_quick = True
		elif arg == 'omp':
			do_test_omp = True
		elif arg == 'cuda':
			do_test_cuda = True
		elif arg == 'full':
			do_test_full= True
		else:
			print 'wrong argument ', arg, ' (allowed arguments: quick, omp, cuda, full)'
else:
	# perform all test: full testing
	do_test_full = True
	
test_dirs = []
if do_test_full:
	test_dirs = [
	'highspin', 'sector1h1p', 'sector2h0p', 'sector0h2p',
    'sector0h1p', 'sector1h0p', 'openmp', 'x2cmmf',
	'efield', 'shifts', 'ccsd', 'cuda', 'oneprop'
	]
if do_test_quick:
	test_dirs += ['oneprop', 'highspin', 'sector0h2p', 'sector1h1p']
if do_test_cuda:
	test_dirs += ['cuda']
if do_test_omp:
	test_dirs += ['openmp']

test_dirs = list(set(test_dirs))  # only unique items

print 'directories with tests to be executed:', test_dirs

for tdir in test_dirs:
	os.chdir(tdir)
	os.system("python test.py")
	os.chdir("..")


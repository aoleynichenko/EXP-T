#!/usr/bin/env python

import sys
import os

# list: directory - enabled/disabled
test_dirs = [
    'highspin',
    'ccsd(t)',
    'sector1h1p',
    'sector2h0p',
    'sector0h2p',
    'diatomic',
    'x2cmmf',
    'efield',
    'ccsdt-n',
    'hughes_kaldor',
    'sector0h3p',
    'sector0h1p',
    'sector1h0p',
    'openmp',
    'shifts',
    'ccsd',
    'cuda',
]

################################################################################

# different sets of tests
# can be enabled by the command line parameter: python test.py [<test-set-1> ...]
do_test_quick = False
do_test_omp  = False
do_test_cuda = False
do_test_full = False
do_test_0h3p = False
do_test_triples = False

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
		elif arg == '0h3p':
			do_test_0h3p = True
		elif arg == 'triples':
			do_test_triples = True
		else:
			print('wrong argument ', arg, ' (allowed arguments: quick, omp, cuda, 0h3p, triples, full)')
else:
	# perform all test: full testing
	do_test_full = True
	
enabled_tests = []
if do_test_full:
	enabled_tests = [
        'highspin',
        'ccsd(t)',
        'sector1h1p',
        'sector2h0p',
        'sector0h2p',
        'diatomic',
        'ccsdt-n',
        'hughes_kaldor',
        'sector0h3p',
        'sector0h1p',
        'sector1h0p',
        'openmp',
        'x2cmmf',
        'efield',
        'shifts',
        'ccsd',
        'cuda',
	]
if do_test_quick:
	enabled_tests += [
        'oneprop',
        'highspin',
        'ccsd(t)',
        'sector0h2p',
        'sector1h1p',
        'diatomic'
    ]
if do_test_cuda:
    enabled_tests += ['cuda']
if do_test_omp:
    enabled_tests += ['openmp']
if do_test_0h3p:
    enabled_tests += ['sector0h3p']
if do_test_triples:
    enabled_tests += ['ccsd(t)', 'ccsdt-n', 'hughes_kaldor']

enabled_tests = list(set(enabled_tests))  # only unique items

################################################################################

print('directories with tests to be executed:')
for tdir in test_dirs:
    if tdir in enabled_tests:
        print("  " + tdir)
print("\n")

################################################################################

for tdir in test_dirs:
    if tdir in enabled_tests:
        os.chdir(tdir)
        os.system("python test.py")
        os.chdir("..")


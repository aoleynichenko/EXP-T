#!/usr/bin/env python

from numpy import *
from scipy.linalg import *
from math import *

N = 2
C = matrix([
 [-3.86843e-02, 9.98522e-01],  
 [ 9.99251e-01, 5.43445e-02]
])

# let vectors be stored by rows
#C = C.T

S = dot(C.H,C)

print 'OVERLAP'
print S

lam, UL, U = eig(S,left=True,right=True)

print 'U = '
print U

UH = matrix(U).H

print 'UH = '
print UH

S12_inv_diag = zeros((N,N))
for i in xrange(0,N):
	S12_inv_diag[i,i] = 1.0/sqrt(lam[i])
print 'S12_inv_diag = '
print S12_inv_diag

# U S12_inv_diag UH

S_12 = dot(U, dot(S12_inv_diag, UH))

print 'S_12'
print S_12

print '=================================='

C_orth = dot(C,S_12)

print dot(C_orth.H, C_orth)

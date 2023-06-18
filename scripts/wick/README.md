
# wick
Automated generation of expressions for matrix elements of arbitrary n-particle operators in the basis of Slater determinants using the Wick theorem.

# Example

Let us evaluate expressions for matrix elements of one- and two-particle effective interaction operators in the (0h2p) Fock space sector. These expressions will be used to construct effective Hamiltonian operator in the Fock-space multireference coupled cluster calculation. We have (in terms of second-quantized operators):

$$ V^{\text{eff}}_{IJ} = \langle \Phi_I |  \hat{V}^{\text{eff}} | \Phi_J \rangle, $$

$$ \hat{V}^{\text{eff}} = \hat{V}^{0h1p} + \hat{V}^{0h2p}, $$

$$  \hat{V}^{0h1p} = \sum_{pq} V^{0h1p}_{pq} \{ a_p^\dagger a_q \}, $$

$$ \hat{V}^{0h2p} = \sum_{pqrs} V^{0h2p}_{pqrs} \{ a_p^\dagger a_q^\dagger a_s a_r \}, $$

$$\langle \Phi_I| = \langle \Phi^{ab} | =  \langle \Phi_0 | a_b a_a, $$

$$|\Phi_J\rangle = |\Phi^{cd} \rangle = a_c^\dagger a_d^\dagger |\Phi_0\rangle, $$

where $V^{\text{eff}}$ stands for the effective interaction operator, $|\Phi_{I}\rangle$ and $|\Phi_{J}\rangle$ stand for model space Slater determinants and $|\Phi_0\rangle$ denotes the vacuum determinant (typically the closed-shell Hartree-Fock one). All indices <i>a</i>, <i>b</i>, <i>p</i>, <i>q</i>, ...  run over active space spin-orbitals (or spinors in the relativistic case).

The input file for <tt>wick</tt> contains definitions of both one- and two-particle interaction operators and bra/ket Slater determinants:
```
#
# comments start with '#'
# Heff operator matrix elements for the Fock-space sector 0h2p
#

# define indices used in determinants and operators
holes i j k
particles a b c d
any p q r s t u

# define determinants
operator bra 1.0 { b a }
operator ket 1.0 { c+ d+ }

# define operators
operator heff1 1.0 { p+ q }
operator heff2 0.5 { p+ q+ s r }

# task
? bra heff1 ket
? bra heff2 ket
```
Note that one can define several tasks in the single input file (we calculate first matrix elements for the one-particle operator and then for the two-particle one). Run:
```
python main.py 02-heff.inp | tee 02-heff.out
```
Output will contain optimized expression for matrix elements:
```
(0) + 1.0 heff1 [ b d ] d_ac
(1) - 1.0 heff1 [ b c ] d_ad
(2) - 1.0 heff1 [ a d ] d_bc
(3) + 1.0 heff1 [ a c ] d_bd
```
and
```
(0) + 1.0 heff2 [ b a d c ]
(1) - 1.0 heff2 [ b a c d ]
```
Note that these expressions do not imply the antisymmetry of the operators. If one uses the Goldstone-Brandow antisymmetrized formalism to deal with diagrams, one have to drop off redundant terms by hands. In this case the last expression should be:
```
(0) + 1.0 heff2 [ b a d c ]
```
Expressions in higher Fock space sectors (0h3p, 1h2p, etc) are evaluated in a completely same manner (see the <tt>examples</tt> directory). All expressions for $\hat{H}^{\text{eff}}$ matrix elements in the EXP-T program package (https://github.com/aoleynichenko/EXP-T) are coded using the <tt>wick</tt> expression generator.

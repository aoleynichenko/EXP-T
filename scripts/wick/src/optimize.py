import copy
import itertools

#
# code from:
# https://stackoverflow.com/questions/1503072/how-to-check-if-permutations-have-equal-parity
#
def parity(permutation):
    permutation = list(permutation)
    length = len(permutation)
    elements_seen = [False] * length
    cycles = 0
    for index, already_seen in enumerate(elements_seen):
        if already_seen:
            continue
        cycles += 1
        current = index
        while not elements_seen[current]:
            elements_seen[current] = True
            current = permutation[current]
    return (length-cycles) % 2 == 0


#
# compares Kroenecker deltas
#
def deltas_are_equivalent(d1, d2):
    d1_sorted = copy.deepcopy(d1)
    d2_sorted = copy.deepcopy(d2)
    return d1_sorted == d2_sorted


#
# checks if terms are equivalent with respect to permutations of electrons
#
def terms_are_equivalent(term1, term2):

    rk1 = len(term1.matr_elem) - 1
    rk2 = len(term2.matr_elem) - 1
    if rk1 != rk2:
        return False

    for d1 in term1.deltas:
        for d2 in term2.deltas:
            if not deltas_are_equivalent(d1, d2):
                return False

    n = int(rk1 / 2)

    bra2 = term2.matr_elem[1:n+1]
    ket2 = term2.matr_elem[n+1:]

    for perm in itertools.permutations([i for i in range(0,n)]):
        permuted = [bra2[perm[i]] for i in range(0,n)] + [ket2[perm[i]] for i in range(0,n)]
        if term1.matr_elem[1:] == permuted:
            return True

    return False


def optimize_expression(expr):

    print('Optimization of analytic expression')
    print('-----------------------------------')
    print('Terms:')
    for i, t in enumerate(expr):
        print('(%d)' % i, t)

    #
    # we will search for all possible permutations of electrons in order to find
    # redundant terms and account for them by increasing the pre-factor
    # of the "main" term
    #

    perm_uniq_expr = []
    factors = []

    for i, t1 in enumerate(expr):
        unique = True
        for j, t2 in enumerate(perm_uniq_expr):
            if terms_are_equivalent(t1, t2):
                factors[j] += 1
                unique = False
                break
        if unique:
            t = copy.deepcopy(t1)
            perm_uniq_expr.append(t)
            factors.append(1)

    for i, expr in enumerate(perm_uniq_expr):
        expr.factor *= factors[i]

    print('\nTerms unique wrt permutations of electrons:')
    for i, t in enumerate(perm_uniq_expr):
        print('(%d)' % i, t)
    print()

    return perm_uniq_expr




import copy

# represents one term in complex matrix element expression
# terms and following simplifications are now implemented for the case
# of only one normal-ordered operator between bra and ket excitations.
class Term:

    def __init__(self, sign, factor, matr_elem, deltas):
        self.sign = sign
        self.factor = factor
        self.matr_elem = copy.deepcopy(matr_elem)
        self.deltas = copy.deepcopy(deltas)

    def __repr__(self):
        s = ''
        s += '+ ' if self.sign == 1 else '- '
        s += str(self.factor) + ' '
        s += self.matr_elem[0]
        s += ' [ '
        for idx in self.matr_elem[1:]:
            s += idx + ' '
        s += ']'
        for d in self.deltas:
            s += ' d_%s%s' % (d)
        return s

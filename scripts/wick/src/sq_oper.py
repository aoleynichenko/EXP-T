#
# Elementary annihilation/creation operator
#
class ElemOperator:

    # name  index name: 'i', 'a', etc
    # typ   '+' for creation, '-' for annihilation
    # hp    'h' for hole index, 'p' for particle index, 'a' - any index
    # op    (type string) name of the normal-ordered operator to which this
    #       elemetary operator belongs
    def __init__(self, name, typ, hp, op):
        self.name = name
        self.typ = typ
        self.hp = hp
        self.op = op

    def __repr__(self):
        s = self.name
        if self.typ == '+':
            s += '+'
        return s

#
# represents second-quantized normal-ordered operator
# for determinants, op == 'bra'/'ket' and factor == 1
#
class Operator:

    # constructor
    #   elem_ops   list with elementary creation/annihilation operators
    #   name       symbolic name of the operator; exceptions only for
    #              names of excited determinants: 'bra' and 'ket'
    #   factor     factor before the sum in this second-quantized normal-
    #              ordered operator
    def __init__(self, elem_ops, name, factor):
        self.elem_ops = elem_ops
        self.name = name
        self.factor = factor

    def __repr__(self):
        if self.name != 'bra' and self.name != 'ket':
            s = str(self.factor) + ' ' + self.name + ' {'
        else:
            s = '{'
        for op in self.elem_ops:
            s += str(op) + ' '
        s += '}'
        return s

import copy

def all_pairs(lst):
    if len(lst) < 2:
        yield copy.deepcopy(lst)
        return
    a = lst[0]
    for i in range(1, len(lst)):
        pair = (a, lst[i])
        for rest in all_pairs(lst[1:i] + lst[i + 1:]):
            yield copy.deepcopy([pair] + rest)


#
# only the following contactions are not equal to zero
# (Bartlett, Shavitt, MBPT in Chem and Phys)
#   i,j - hole indices
#   a,b - particle indices
#   p,q - any indices
# non-zero contractions:
#   i+ j = d_ij (Kronecker delta)
#   a b+ = d_ab
#   i+ q = d_iq
#   p+ j = d_pj
#   a p+ = d_ap
#   q b+ = d_qb
#
def construct_all_possible_contractions(task, operators_dict):

    cntrs = []

    op_str = []
    for op in [operators_dict[op_name] for op_name in task]:
        op_str += op.elem_ops
    n_oper = len(op_str)

    for cntr in all_pairs([i for i in range(0,n_oper)]):
        # check for a-priori equiality to zero
        # loop over elementary contractions
        zero = False
        for p in cntr:
            op1 = op_str[p[0]]
            op2 = op_str[p[1]]
            # eliminate contractions inside single normal-ordered operator
            if op1.op == op2.op:
                zero = True
                break
            # eliminate contractions annihil-annihil and creat-creat
            if op1.typ == op2.typ:
                zero = True
                break
            # eliminate contractions ij+ and a+b
            if op1.typ == '-' and op2.typ == '+' and op1.hp == 'h' and op2.hp == 'h':
                zero = True
                break
            if op1.typ == '+' and op2.typ == '-' and op1.hp == 'p' and op2.hp == 'p':
                zero = True
                break
        if zero:
            continue
        cntrs.append(cntr)

    if len(cntrs) == 0:
        print('----------------------------------------------------------------')
        print('No fully contracted terms possible, matrix element will be zero!')
        exit()

    #
    # beautiful output
    #
    print('\nPossible non-zero contractions:')
    print('-------------------------------\n')
    for nc, cntr in enumerate(cntrs):
        labels = [0 for j in range(0, len(op_str))]
        for i, p in enumerate(cntr):
            for j in p:
                labels[j] = i + 1
        s1, s2 = '', ''
        for i, op in enumerate(op_str):
            s1 += '%-2s ' % (op)
        for lab in labels:
            s2 += '%-2s ' % (str(lab))
        print('(%2d)  ' % nc, s1)
        print('      ', s2, '\n')
    print('-------------------------------')
    print('Total number of possible contractions: ', len(cntrs), '\n')

    return cntrs



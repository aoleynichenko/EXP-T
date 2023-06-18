import copy
import term

def generate_expressions(task, operators_dict, cntrs):

    expr = []

    if len(cntrs) == 0:
        return

    print()
    print('Expression for matrix element:')
    print('------------------------------')

    operators = [operators_dict[op_name] for op_name in task]
    op_str = []
    for op in operators:
        op_str += op.elem_ops

    for contraction in cntrs:
        ops = copy.deepcopy(op_str)
        nop = len(ops)
        sign = +1
        mask = [1 for i in range(0,nop)]
        deltas = []
        for p in contraction:
            i1 = p[0]
            i2 = p[1]
            if i1 > i2:
                i1, i2 = i2, i1
            mask[i1] = 0
            mask[i2] = 0
            nperm = sum(mask[i1 + 1:i2])
            if nperm % 2 != 0:
                sign = -sign
            # put indices to Kronecked delta in lexicographical order
            idx_1 = ops[i1].name
            idx_2 = ops[i2].name
            if idx_2 < idx_1:
                idx_1, idx_2 = idx_2, idx_1
            deltas.append((idx_1, idx_2))

        # construct string of indices for matrix elements (except bra/ket excitation op-s)
        # Note: p+ q+ s r -> Op[p,q,r,s]
        op_ind_strings = []
        op_names = []
        for op in operators:
            if op.name == 'bra' or op.name == 'ket':
                continue

            op_names.append(op.name)
            ind_creat = []
            ind_annih = []
            for elem_op in op.elem_ops:
                if elem_op.typ == '+':
                    ind_creat.append(elem_op.name)
                else:
                    ind_annih.append(elem_op.name)

            op_ind_strings.append([op.factor, op.name] + ind_creat + ind_annih[::-1])

        # extract indices which belong to <bra| and |ket> excitation op-s (if presented)
        bra_ket_indices = []
        if 'bra' in operators_dict:
            bra_ket_indices += [elem_op.name for elem_op in operators_dict['bra'].elem_ops]
        if 'ket' in operators_dict:
            bra_ket_indices += [elem_op.name for elem_op in operators_dict['ket'].elem_ops]
        bra_ket_indices = set(bra_ket_indices)

        # print 'raw expressions' for simple validation
        if sign == 1:
            print('+', end='')
        else:
            print('-', end='')
        for op_ind_s in op_ind_strings:
            print(op_ind_s[0], op_ind_s[1], '[ ', end='')
            for idx in op_ind_s[2:]:
                print(idx, end=' ')
            print(']', end='')
        for d in deltas:
            print(' d_%s%s' % (d), end='')
        print(' => ', end='')

        # now we must remove all redundant deltas
        # example: - f[p,q] d_iq d_ab d_jp => - f[j,i] d_ab
        # remove d_ii type deltas
        new_deltas = []
        for d in deltas:
            if d[0] != d[1]:
                new_deltas.append(d)
        deltas = copy.deepcopy(new_deltas)

        # remove non-bra-ket indices from matrix elements expressions (if possible)
        for op_i, op_ind_s in enumerate(op_ind_strings):
            fact = op_ind_s[0]
            name = op_ind_s[1]
            indices = op_ind_s[2:]
            for i, idx in enumerate(indices):
                deltas = copy.deepcopy(new_deltas)
                new_deltas = []
                if idx in bra_ket_indices:
                    continue
                for d in deltas:
                    if idx == d[0]:
                        indices[i] = d[1]
                    elif idx == d[1]:
                        indices[i] = d[0]
                    else:
                        new_deltas.append(d)
            op_ind_strings[op_i][2:] = indices

        deltas = new_deltas

        # save this term if only one op-r between bra and ket
        if len(op_ind_strings) == 1:
            factor = op_ind_strings[0][0]
            matr_elem = op_ind_strings[0][1:]
            expr.append(term.Term(sign, factor, matr_elem, deltas))

        # print final expression
        if sign == 1:
            print('+', end='')
        else:
            print('-', end='')
        for op_ind_s in op_ind_strings:
            print(op_ind_s[0], op_ind_s[1], '[ ', end='')
            for idx in op_ind_s[2:]:
                print(idx, end=' ')
            print(']', end='')
        for d in deltas:
            print(' d_%s%s' % d, end='')
        print()

    print('------------------------------\n')

    return expr



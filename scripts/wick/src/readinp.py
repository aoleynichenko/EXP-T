import re
import copy
import sq_oper


def echo_input_file(path):

    try:
        print('\necho of input file:')
        print('----------------------------------------------------------')

        with open(path, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                print(line, end='')

        print('----------------------------------------------------------\n')

    except Exception as exc:
        print('Error while reading input file: ', exc)
        exit()


#
# input file parser.
#
def read_input(path):
    names_holes = set()
    names_parts = set()
    names_any = set()
    oper_dict = {}
    tasks = []

    echo_input_file(path)

    try:
        with open(path, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break

                # remove leading and trailing spaces
                line = line.strip()
                if not line:
                    continue

                # skip comments
                if line.startswith('#'):
                    continue

                # letters used to denote indices of holes
                elif line.startswith('holes'):
                    names = line.split()[1:]
                    for nm in names:
                        names_holes.add(nm)
                    print('hole indices are: ', [nm for nm in names_holes])

                # letters used to denote indices of particles
                elif line.startswith('particles'):
                    names = line.split()[1:]
                    for nm in names:
                        names_parts.add(nm)
                    print('particle indices are: ', [nm for nm in names_parts])

                # letters used to denote indices of both holes and particles (any orbitals)
                elif line.startswith('any'):
                    names = line.split()[1:]
                    for nm in names:
                        names_any.add(nm)
                    print('any indices are: ', [nm for nm in names_any])

                # definition of a second-quantized operator
                elif line.startswith('operator'):
                    pattern = re.compile('operator\s+(\w+)\s+(\d+\.\d+)\s+\{\s+((\w+\+?\s+)+)\}')
                    m = pattern.findall(line)[0]
                    name = m[0]
                    factor = float(m[1])
                    elem_op_symbols = m[2].split()
                    elem_op_str = []
                    for op_sym in elem_op_symbols:
                        typ = '-'
                        if op_sym[-1] == '+':
                            typ = '+'
                            op_sym = op_sym[:-1]
                        hp = None
                        if op_sym in names_holes:
                            hp = 'h'
                        elif op_sym in names_parts:
                            hp = 'p'
                        elif op_sym in names_any:
                            hp = 'a'
                        else:
                            raise Exception('unknown type (hole/part/any) for index "' + op_sym + '"')

                        elem_op_str.append(sq_oper.ElemOperator(op_sym, typ, hp, name))

                    oper = sq_oper.Operator(elem_op_str, name, factor)
                    oper_dict[name] = copy.deepcopy(oper)

                elif line.startswith('?'):
                    oper_names = line.split()[1:]
                    tasks.append(oper_names)
                else:
                    raise Exception('wrong line: "' + line + '"')

        # we return dictionary containing:
        # 1. definitions of second-quantized operators
        # 2. list of tasks: combinations of operators to be contracted
        return oper_dict, tasks

    except Exception as exc:
        print('Error while reading input file: ', exc)
        exit()


def print_task(task, operators_dict):

    print('\n\n')
    print('----------------------------------------------------------')
    print('Begin task: ', task)
    print('----------------------------------------------------------\n')

    mat_elem_expr = '<0|'
    oper_list = []
    for op_name in task:
        oper_list.append(operators_dict[op_name])
        mat_elem_expr += str(operators_dict[op_name])
    mat_elem_expr += '|0>'

    print('Matrix element expression: ' + mat_elem_expr)
    op_str = []
    for op in oper_list:
        op_str += op.elem_ops

    print('Elementary operators string: ', op_str)
    print()

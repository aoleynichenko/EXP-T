#
# TOOL FOR AUTOMATIC EVALUATION OF MATRIX ELEMENTS OF ARBITRARY
# NORMAL-ORDERED SECOND-QUANTIZED OPERATORS ON VACUUM DETERMINANTS.
#
# Alexander Oleynichenko, 2017-2022
# NRC "Kurchatov Institute" - PNPI, Gatchina, Russia
# E-mail: alexvoleynichenko@gmail.com
#


import sys
from datetime import datetime

import genexpr
import readinp
import contract
import optimize

#
# parse command-line arguments
#
if len(sys.argv) != 2:
    print('Usage: python %s <input-file>' % (sys.argv[0]))
    exit()

#
# header
#
print("wick v1.1")
print("by a. oleynichenko")
print("run date: ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

#
# parse input file
#
operators, tasks = readinp.read_input(sys.argv[1])

#
# generate analytic expression for all the matrix elements required
#
for task in tasks:

    readinp.print_task(task, operators)

    contractions = contract.construct_all_possible_contractions(task, operators)

    expr = genexpr.generate_expressions(task, operators, contractions)
    optimize.optimize_expression(expr)



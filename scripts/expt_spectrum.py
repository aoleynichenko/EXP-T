#!/usr/bin/env python
#
# EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
# Copyright (C) 2018-2022 The EXP-T developers.
#
# This file is part of EXP-T.
#
# EXP-T is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EXP-T is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
#
# E-mail:        exp-t-program@googlegroups.com
# Google Groups: https://groups.google.com/d/forum/exp-t-program
#

#
# alexander oleynichenko, 21 january 2023
#

import argparse
import sys
import math
import numpy as np

#
# some physical constants
#
AU_TO_CM = 219474.63136320  # atomic units of energy -> cm^-1
AU_TO_EV = 27.211386245988  # atomic units of energy -> eV


class ElectronicState(object):
    """
    information about an electronic state
    """

    def __init__(self, irrep_name, abs_energy):
        self.irrep_name = irrep_name
        self.abs_energy = abs_energy


class ElectronicStatesList(object):
    """
    data structure containing information about all electronic states
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.states_by_irreps = {}

    def add_state(self, irrep_name, abs_energy):
        if not irrep_name in self.states_by_irreps.keys():
            self.states_by_irreps[irrep_name] = []
        self.states_by_irreps[irrep_name].append(ElectronicState(irrep_name, abs_energy))

    def get_states(self, irrep_name):
        return self.states_by_irreps[irrep_name]

    def get_state(self, irrep_name, state_number):
        return self.states_by_irreps[irrep_name][state_number]

    def get_state_energy(self, irrep_name, state_number):
        return self.get_state(irrep_name, state_number).abs_energy

    def min_energy(self):
        emin = 0.0

        for irrep in self.states_by_irreps:
            for state in self.states_by_irreps[irrep]:
                if state.abs_energy < emin:
                    emin = state.abs_energy

        return emin

    def get_irrep_list(self):
        return list(self.states_by_irreps.keys())

    def num_irreps(self):
        return len(self.get_irrep_list())

    def num_states(self, irrep_name=None):
        if irrep_name:
            try:
                return len(self.states_by_irreps[irrep_name])
            except:
                pass
        else:
            count = 0
            for irrep in self.states_by_irreps:
                count += len(self.states_by_irreps[irrep])
            return count


class TransitionMomentsTable(object):
    """
    table of expectation values and transition moments
    """

    def __init__(self, elec_states):
        self.elec_states = elec_states
        irrep_list = elec_states.get_irrep_list()
        self.dm_table = {}
        for irrep_1 in irrep_list:
            self.dm_table[irrep_1] = {}
            for irrep_2 in irrep_list:
                self.dm_table[irrep_1][irrep_2] = [None, None, None]  # x,y,z matrices

    def put_matrix_element(self, bra_irrep, bra_state_no, ket_irrep, ket_state_no, coord, value):
        if coord < 0 or coord > 2:
            raise Exception('wrong coordinate: %d' % coord)

        if self.dm_table[bra_irrep][ket_irrep][coord] is None:
            bra_dim = self.elec_states.num_states(bra_irrep)
            ket_dim = self.elec_states.num_states(ket_irrep)
            self.dm_table[bra_irrep][ket_irrep][coord] = np.zeros((bra_dim, ket_dim), dtype=np.complex128)

        self.dm_table[bra_irrep][ket_irrep][coord][bra_state_no][ket_state_no] = value

    def get_matrix_element(self, bra_irrep, bra_state_no, ket_irrep, ket_state_no, coord):
        try:
            return self.dm_table[bra_irrep][ket_irrep][coord][bra_state_no][ket_state_no]
        except:
            return 0.0 + 0.0j

    def get_dipole_moment(self, bra_irrep, bra_state_no, ket_irrep, ket_state_no):
        dm = np.zeros(3, dtype=np.complex128)
        for coord in [0, 1, 2]:
            dm[coord] = self.get_matrix_element(bra_irrep, bra_state_no, ket_irrep, ket_state_no, coord)
        return dm

    def get_dipole_moment_squared(self, bra_irrep, bra_state_no, ket_irrep, ket_state_no):
        dm = self.get_dipole_moment(bra_irrep, bra_state_no, ket_irrep, ket_state_no)
        return dm.dot(dm.conjugate()).real

    def get_dipole_moment_abs(self, bra_irrep, bra_state_no, ket_irrep, ket_state_no):
        d2 = self.get_dipole_moment_squared(bra_irrep, bra_state_no, ket_irrep, ket_state_no)
        return math.sqrt(d2)

    def get_oscillator_strength(self, bra_irrep, bra_state_no, ket_irrep, ket_state_no):
        energy_1 = self.elec_states.get_state_energy(bra_irrep, bra_state_no)
        energy_2 = self.elec_states.get_state_energy(ket_irrep, ket_state_no)
        d2 = self.get_dipole_moment_squared(bra_irrep, bra_state_no, ket_irrep, ket_state_no)
        return 2.0 / 3.0 * (energy_2 - energy_1) * d2


def print_electronic_states_list(elec_states, upper_thresh=1e10):
    """
    prints a list of electronic states, irrep by irrep
    :param elec_states: object of class ElectronicStatesList
    :param upper_thresh: upper energy limit for the electronic states printed
    """
    emin = elec_states.min_energy()

    for irrep in elec_states.states_by_irreps:
        print('\nirrep %s' % irrep)
        print('state         abs energy [au]     energy [ev]   energy[cm^-1]')

        for i, state in enumerate(elec_states.states_by_irreps[irrep]):
            energy = state.abs_energy
            e_ev = (energy - emin) * AU_TO_EV
            e_cm = (energy - emin) * AU_TO_CM

            if e_cm < upper_thresh:
                print('%5d%24.12f%16.4f%16.0f' % (i + 1, energy, e_ev, e_cm))


def print_transition_moments_table(
        tdm_table,
        irreps_from=None,
        states_from=None,
        irreps_to=None,
        states_to=None,
        energy_thresh=1e10,
        hermitization=False,
        units='cm'
):
    """
    prints a table containing transition dipole moments
    :param tdm_table: object of type TransitionMomentsTable
    :param irreps_from: list of irreps to which initial electronic states should belong
    :param states_from: sequential numbers of initial electronic states
    :param irreps_to: list of irreps to which final electronic states should belong
    :param states_to: sequential numbers of final electronic states
    :param energy_thresh: upper energy limit for electronic states considered (in cm^-1)
    :param hermitization: force hermiticity of the d_if/d_fi moments or not
    :param units: 'cm' for cm^-1, 'ev' for electron-volts
    """
    print('\ntable of transition moments:\n')
    print('initial irreps: ', 'all' if irreps_from is None else irreps_from)
    print('initial states: ', 'all' if states_from is None else states_from)
    print('final irreps  : ', 'all' if irreps_to is None else irreps_to)
    print('final states  : ', 'all' if states_to is None else states_to)

    elec_states = tdm_table.elec_states
    irrep_list_from = irreps_from if irreps_from else elec_states.get_irrep_list()
    irrep_list_to = irreps_to if irreps_to else elec_states.get_irrep_list()
    emin = elec_states.min_energy()

    for irrep_1 in irrep_list_from:
        for irrep_2 in irrep_list_to:

            header_printed = False

            for i, bra_state in enumerate(elec_states.get_states(irrep_1)):
                for j, ket_state in enumerate(elec_states.get_states(irrep_2)):

                    if (not states_from is None) and (not i in states_from):
                        continue
                    if (not states_to is None) and (not j in states_to):
                        continue

                    energy_1 = bra_state.abs_energy
                    energy_2 = ket_state.abs_energy
                    delta_e = energy_2 - energy_1
                    if units == 'cm':
                        energy_1 = AU_TO_CM * (energy_1 - emin)
                        energy_2 = AU_TO_CM * (energy_2 - emin)
                    else:
                        energy_1 = AU_TO_EV * (energy_1 - emin)
                        energy_2 = AU_TO_EV * (energy_2 - emin)

                    if energy_1 > energy_thresh or energy_2 > energy_thresh:
                        continue

                    # i -> f
                    dm_if = tdm_table.get_dipole_moment(irrep_1, i, irrep_2, j)
                    d2_if = tdm_table.get_dipole_moment_squared(irrep_1, i, irrep_2, j)
                    osc_str_if = tdm_table.get_oscillator_strength(irrep_1, i, irrep_2, j)
                    dabs_if = math.sqrt(d2_if)

                    # f -> i
                    dm_fi = tdm_table.get_dipole_moment(irrep_2, j, irrep_1, i)
                    d2_fi = tdm_table.get_dipole_moment_squared(irrep_2, j, irrep_1, i)
                    osc_str_fi = tdm_table.get_oscillator_strength(irrep_2, j, irrep_1, i)
                    dabs_fi = math.sqrt(d2_fi)

                    if dabs_if < 1e-6 and dabs_fi < 1e-6:
                        continue

                    if not header_printed:
                        header_printed = True
                        print('\ntransitions ', irrep_1, ' -> ', irrep_2, '\n')
                        if units == 'cm':
                            print('               e_i,cm^-1   e_f,cm^-1          re(dx)      im(dx)'
                                  '        re(dy)      im(dy)        re(dz)      im(dz)        |d|,a.u.    osc.str.')
                        else:
                            print('                  e_i,ev      e_f,ev          re(dx)      im(dx)'
                                  '        re(dy)      im(dy)        re(dz)      im(dz)        |d|,a.u.    osc.str.')

                    if units == 'cm':
                        format_string = '%4d -> %4d%12.4f%12.4f    %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f    ' \
                                        '%12.6f%12.6f '
                    else:
                        format_string = '%4d -> %4d%12.4f%12.4f    %12.6f%12.6f  %12.6f%12.6f  %12.6f%12.6f    ' \
                                        '%12.6f%12.6f '

                    print(format_string % (i + 1, j + 1, energy_1, energy_2, dm_if[0].real, dm_if[0].imag,
                                           dm_if[1].real, dm_if[1].imag, dm_if[2].real, dm_if[2].imag, dabs_if,
                                           osc_str_if))

                    if hermitization:
                        print(format_string % (j + 1, i + 1, energy_2, energy_1, dm_fi[0].real, dm_fi[0].imag,
                                               dm_fi[1].real, dm_fi[1].imag, dm_fi[2].real, dm_fi[2].imag, dabs_fi,
                                               osc_str_fi))

                        dabs_if_herm = math.sqrt(dabs_if * dabs_fi)
                        osc_str_if_herm = 2.0 / 3.0 * delta_e * dabs_if * dabs_fi
                        print('                                                            '
                              '                                                            '
                              '%12.6f%12.6f' % (dabs_if_herm, osc_str_if_herm))


def parse_expt_output(path):
    """
    parses EXP-T output file.
    :param path: to the exp-t output file
    :return: information about electronic states (ElectronicStatesList)
    and transition dipole moments (TransitionMomentsTable)
    """
    elec_states = ElectronicStatesList()
    tdm_table = None
    coord = 0
    max_irrep_len = 0

    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break

            # begin new sector
            if line.startswith(' solution of amplitude equations'):
                elec_states.reset()
                tdm_table = None

            # information about electronic states energies
            if line.startswith(' Irrep'):

                line_contents = line.replace('(', ' ').replace(')', ' ').split()
                irrep_name = line_contents[2]
                state_number = int(line_contents[4])
                state_energy = float(line_contents[6])

                elec_states.add_state(irrep_name, state_energy)

                if state_number > max_irrep_len:
                    max_irrep_len = state_number

            # information about transition moments
            if line.startswith(' ** MDPROP:'):

                if not tdm_table:
                    tdm_table = TransitionMomentsTable(elec_states)

                if line.startswith(' ** MDPROP: XDIPLEN'):
                    coord = 0
                elif line.startswith(' ** MDPROP: YDIPLEN'):
                    coord = 1
                elif line.startswith(' ** MDPROP: ZDIPLEN'):
                    coord = 2

            if ") -> " in line:
                line_contents = line.replace('(', ' ').replace(')', ' ').split()
                state_1 = int(line_contents[0]) - 1
                state_2 = int(line_contents[3]) - 1
                irrep_1 = line_contents[1]
                irrep_2 = line_contents[4]
                tdm_re = float(line_contents[7])
                tdm_im = float(line_contents[8])

                tdm_table.put_matrix_element(irrep_1, state_1, irrep_2, state_2, coord, tdm_re + 1.0j * tdm_im)

    return elec_states, tdm_table


#
# entry point
#
if __name__ == '__main__':

    # parse command-line arguments

    parser = argparse.ArgumentParser(
        prog='expt_spectrum.py',
        description='parses exp-t output files and collects information '
                    'about electronic transitions and their dipole moments',
        epilog='please report bugs to <exp-t-program@googlegroups.com>.')

    parser.add_argument('filename')  # positional argument
    parser.add_argument('-u', dest='upper', action='store', nargs=1, required=False, type=float,
                        help='upper energy limit for the electronic states considered (cm^-1 or ev)')
    parser.add_argument('-irep', dest='initial_irrep', action='store', nargs=1, required=False,
                        help='irrep of initial electronic states')
    parser.add_argument('-frep', dest='final_irrep', action='store', nargs=1, required=False,
                        help='irrep of final electronic states')
    parser.add_argument('-i', dest='initial_state', action='store', nargs=2, required=False,
                        help='initial electronic state')
    parser.add_argument('-f', dest='final_state', action='store', nargs=2, required=False,
                        help='final electronic state')
    parser.add_argument('-H', dest='do_hermit', action='store_true', required=False,
                        help='force hermiticity of transition moments')
    parser.add_argument('-s', dest='print_spectrum', action='store_true', required=False,
                        help='print detailed information on electronic spectrum')
    parser.add_argument('-ev', dest='units_ev', action='store_true', required=False,
                        help='use electron-volts units for energy')

    args = parser.parse_args()
    energy_units = 'ev' if args.units_ev else 'cm'
    args.upper = float(args.upper[0]) if args.upper else 1e9

    print('a oleynichenko, 2023')
    print()
    print('filename            ', args.filename)
    print('upper limit [%2s]    ' % energy_units, args.upper)
    print('initial irreps      ', args.initial_irrep)
    print('final irreps        ', args.final_irrep)
    print('initial elec state  ', args.initial_state)
    print('final elec state    ', args.final_state)
    print('force hermiticity   ', args.do_hermit)
    print('print elec spectrum ', args.print_spectrum)
    print('energy units        ', energy_units)

    irreps_from = args.initial_irrep
    states_from = None
    if args.initial_state:
        irreps_from = [args.initial_state[0]]
        states_from = [int(args.initial_state[1]) - 1]

    irreps_to = args.final_irrep
    states_to = None
    if args.final_state:
        irreps_to = [args.final_state[0]]
        states_to = [int(args.final_state[1]) - 1]

    # parse EXP-T output file
    elec_states, tdms = parse_expt_output(args.filename)

    # print list of irreducible representations found in the file
    irrep_list = elec_states.get_irrep_list()
    print('irreps: ')
    print(irrep_list)

    # print list of electronic states below the given threshold
    if args.print_spectrum:
        print_electronic_states_list(elec_states, upper_thresh=args.upper)

    # print table of transition moments
    print_transition_moments_table(tdms, irreps_from, states_from, irreps_to, states_to,
                                   energy_thresh=args.upper, hermitization=args.do_hermit,
                                   units=energy_units)

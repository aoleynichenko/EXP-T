!
! EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
! Copyright (C) 2018-2025 The EXP-T developers.
!
! This file is part of EXP-T.
!
! EXP-T is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EXP-T is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
!
! E-mail:        exp-t-program@googlegroups.com
! Google Groups: https://groups.google.com/d/forum/exp-t-program
!

!
! Reads MDPROP unformatted file containing transformed property integrals.
! Prints information about stored property matrices in the human-readable format.
!
! 2023 Alexander Oleynichenko
!


program read_mdprop

    implicit none
    logical :: mdprop_exists
    integer(4) :: mdprop
    character(len = 32) :: achar
    integer(4) :: count_prop

    inquire(file='MDPROP', exist=mdprop_exists)
    if (.not. mdprop_exists) then
        print *, 'Error: no MDPROP file in the working directory'
        stop
    end if

    open(unit=mdprop, file='MDPROP', form='unformatted', status='old')
    rewind (mdprop)

    print *
    print *, "*** MDPROP FILE ***"
    print *

    ! begin loop
    count_prop = 0
    1   continue

    ! read label of the current matrix
    read (mdprop, end = 11, err = 12) achar
    if (achar(25:32) .eq. 'EOFLABEL') goto 11

    count_prop = count_prop + 1
    print '(i4,a,a)', count_prop, '    ', achar(25:)

    ! skip matrix
    read (mdprop, end = 11, err = 12)

    ! go to next record
    goto 1

    12  continue
    print *, 'error while reading MDPROP'
    stop

    11  continue
    close(mdprop)

    print *
    print *, "*** END OF MDPROP FILE ***"
    print *

end program read_mdprop




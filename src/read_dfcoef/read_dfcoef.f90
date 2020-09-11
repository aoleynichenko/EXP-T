!
! EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
! Copyright (C) 2018-2020 The EXP-T developers.
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
! Reads DFCOEF unformatted file with spinors and one-electron energies.
! Prints them in the human-readable format.
!
! See also
! http://diracprogram.org/doc/release-18/programmers/dfcoef.html
! for further details on the DFCOEF files.
!
! 2020 Alexander Oleynichenko
!


program read_dfcoef

    logical :: dfcoef_exists
    integer(8) :: i
    integer(4) :: dfcoef
    integer(8) :: isym, nfsym
    integer(8) :: idim(3,2)
    integer(8) :: npsh(2), nesh(2), nfbas(2)
    character text*74
    real(8) :: toterg
    real(8), dimension(:), allocatable :: eps
    integer(4) :: offset, j

    print *, 'reading DFCOEF'
    inquire(file='DFCOEF', exist=dfcoef_exists)
    if (.not. dfcoef_exists) then
        print *, 'Error: no DFCOEF file in the working directory'
        stop
    end if

    open(unit=dfcoef, file='DFCOEF', form='unformatted', status='old')

    rewind (dfcoef)
    read (dfcoef) text, nfsym, ((idim(i,j),i=1,3),j=1,nfsym),toterg
    nesh = 0
    npsh = 0
    nfbas = 0
    do i = 1, nfsym
        npsh(i) = idim(1,i)
        nesh(i) = idim(2,i)
        nfbas(i) = idim(3,i)
    end do
    print *, 'comment: ', text
    print *, 'nfsym = ', nfsym
    print *, 'number of electronic solutions = ', (nesh(i),i=1,nfsym)
    print *, 'number of positronic solutions = ', (npsh(i),i=1,nfsym)
    print *, 'number of basis functions (L)  = ', (nfbas(i),i=1,nfsym)
    print *, 'total SCF energy = ', toterg

    norbt = nesh(1) + nesh(2) + npsh(1) + npsh(2)
    print *, 'norbt = ', norbt
    allocate(eps(norbt))

    eps = 0
    read (dfcoef)
    read (dfcoef) eps
    read (dfcoef)

    ! print one-electron energies
    ! positronic solutions will be skipped (if they are presented)
    offset = 0
    if (nfsym == 1) then
        print *
        print *, ' one-electron energies '
        print *, '-----------------------'
        print *
    end if
    do isym = 1, nfsym
        if (nfsym == 2) then
            print *
            if (isym == 1) then
                print *, ' one-electron energies (gerade) '
                print *, '--------------------------------'
            end if
            if (isym == 2) then
                print *, ' one-electron energies (ungerade) '
                print *, '----------------------------------'
            end if
            print *
        end if
        offset = offset + npsh(isym)
        j = 0
        do i = offset+1,offset+nesh(isym)
            j = j + 1
            print '(i6,f20.12)', j, eps(i)
        end do
        offset = offset + nesh(isym)
    end do
    print *

    deallocate(eps)
    close(unit=dfcoef)

end program read_dfcoef




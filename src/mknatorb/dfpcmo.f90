!
! EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
! Copyright (C) 2018-2021 The EXP-T developers.
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
! Tools for operating with formatted files DFPCMO.
! DFPCMO files are produced by DIRAC contain coefficients of molecular spinors.
! 2019-2021 Alexander Oleynichenko
!
module dfpcmo

contains

    !
    ! Reads some orbitals expanded in the AO basis from the DFPCMO-type
    ! formatted file named 'filename'
    !
    subroutine read_dfpcmo(filename, nsym, idim, scf_energy, cmo, eps, ibeig)

        implicit none
        character(len = 128), intent(in) :: filename
        integer :: fd
        integer :: i, j, isym, iz
        ! orbital data: output
        character(len = 128) :: comment_line
        integer, intent(out) :: nsym
        integer, intent(out) :: idim(3, 2)
        integer, allocatable, dimension(:, :), intent(inout) :: ibeig
        real(8), allocatable, dimension(:, :), intent(inout) :: eps
        real(8), allocatable, dimension(:, :), intent(inout) :: cmo
        real(8), intent(out) :: scf_energy
        ! local variables
        integer :: npsh(2), nesh(2), nvec(2), nbas(2)
        real(8), allocatable, dimension(:) :: real_buf
        integer, allocatable, dimension(:) :: int_buf
        integer :: nz
        integer :: ncoef(2), offs

        ncoef = 0

        print *, 'READING DFPCMO-LIKE FORMATTED FILE ', trim(filename)

        call estimate_arith(trim(filename), nz)
        if (nz /= 1) then
            print *, 'ERROR! CURRENT VERSION OF THE PROGRAM WORKS ONLY WITH REAL-VALUED SPINORS! (NZ=1)'
            stop
        end if

        fd = 10
        open(unit = fd, file = trim(filename), form = 'formatted', status = 'old', err = 13)
        read (fd, '(a)'), comment_line

        idim = 0
        read (fd, *) nsym, ((idim(i, j), i = 1, 3), j = 1, nsym)
        npsh(1) = idim(1, 1)
        npsh(2) = idim(1, 2)
        nesh(1) = idim(2, 1)
        nesh(2) = idim(2, 2)
        nbas(1) = idim(3, 1)
        nbas(2) = idim(3, 2)
        nvec(1) = npsh(1) + nesh(1)
        nvec(2) = npsh(2) + nesh(2)
        read (fd, *) scf_energy

        write (*, '(a,a)') ' COMMENT LINE = ', trim(comment_line)
        print *, 'NSYM = ', nsym
        write (*, '(a,i0,a,i0)') ' NUM OF POSITRONIC SOLUTIONS = ', idim(1, 1), ' + ', idim(1, 2)
        write (*, '(a,i0,a,i0)') ' NUM OF ELECTRONIC SOLUTIONS = ', idim(2, 1), ' + ', idim(2, 2)
        write (*, '(a,i0,a,i0)') ' NUM OF BASIS FUNCTIONS      = ', idim(3, 1), ' + ', idim(3, 2)
        write (*, '(a,f20.12)')  ' TOTAL SCF ENERGY (A.U.)     = ', scf_energy

        allocate(cmo(nsym, max(nvec(1), nvec(2)) * max(nbas(1), nbas(2)) * nz))
        allocate(eps(nsym, max(nvec(1), nvec(2))))
        allocate(ibeig(nsym, max(nvec(1), nvec(2))))
        allocate(real_buf(nsym * (nvec(1) + nvec(2)) * max(nbas(1), nbas(2)) * nz))
        allocate(int_buf (nsym * (nvec(1) + nvec(2))))

        ! read MO coefficients
        ! total number of coef-s:
        do isym = 1, nsym
            ncoef(isym) = nvec(isym) * nbas(isym)! * nz
        end do
        ! read ALL coef-s to the buffer
        read (fd, *) (real_buf(i), i = 1, nz * (ncoef(1) + ncoef(2)))
        ! redistribution of coef-s by irreps
        offs = 0
        do isym = 1, nsym
            do i = 1, ncoef(isym)
                cmo(isym, i) = real_buf(offs + i)
            end do
            offs = offs + ncoef(isym)     ! to the next irrep if needed
        end do

        ! read orbital energies
        ! read ALL energies to the buffer
        read (fd, *) (real_buf(i), i = 1, nvec(1) + nvec(2))
        ! redistribution of orb energies by irreps
        offs = 0
        do isym = 1, nsym
            do i = 1, nvec(isym)
                eps(isym, i) = real_buf(offs + i)
            end do
            offs = offs + nvec(isym)    ! to the next irrep if needed
        end do

        ! read ibeig
        ! read ALL ibeig to the buffer
        read (fd, *) (int_buf(i), i = 1, nvec(1) + nvec(2))
        offs = 0
        do isym = 1, nsym
            do i = 1, nvec(isym)
                ibeig(isym, i) = int_buf(offs + i)
            end do
            offs = offs + nvec(isym)    ! to the next irrep if needed
        end do

        close(unit = fd)

        ! print orbital energies
        print *, 'ORBITAL ENERGIES'
        do isym = 1, nsym
            print *, 'ISYM = ', isym
            do i = 1, nvec(isym)
                print '(i4,f20.12)', i, eps(isym, i)
            end do
        end do
        print *, 'END OF FILE ', trim(filename)
        print *

        deallocate(real_buf)
        deallocate(int_buf)
        return

        13  continue
        print *, 'ERROR: CANNOT OPEN DFPCMO FILE ', trim(filename)

    end subroutine read_dfpcmo


    !
    ! Write some orbitals expanded in the AO basis to the DFPCMO-type
    ! formatted file named 'filename'
    !
    subroutine write_dfpcmo(filename, nsym, idim, scf_energy, cmo, eps, ibeig)

        implicit none
        integer :: fd
        integer :: i, j, isym, iz
        ! orbital data: input
        character(len = 128), intent(in) :: filename
        integer, intent(in) :: nsym
        integer, intent(in) :: idim(3, 2)
        integer, dimension(:, :), intent(in) :: ibeig
        real(8), dimension(:, :), intent(in) :: eps
        real(8), dimension(:, :, :), intent(in) :: cmo
        real(8), intent(in) :: scf_energy
        ! local variables
        integer :: npsh(2), nesh(2), nvec(2), nbas(2)
        real(8), allocatable, dimension(:) :: real_buf
        integer, allocatable, dimension(:) :: int_buf
        integer :: nz = 2
        integer :: ncoef(2), offs

        print *, 'WRITING DFPCMO-LIKE FORMATTED FILE ', trim(filename)

        fd = 11
        open(unit = fd, file = trim(filename), form = 'formatted', err = 13)
        write (fd, '(a)') 'this DFPCMO file was generated by the MKNATORB program'

        write (fd, '(7(1x,i0))') nsym, ((idim(i, j), i = 1, 3), j = 1, nsym)
        ncoef = 0
        npsh(1) = idim(1, 1)
        npsh(2) = idim(1, 2)
        nesh(1) = idim(2, 1)
        nesh(2) = idim(2, 2)
        nbas(1) = idim(3, 1)
        nbas(2) = idim(3, 2)
        nvec(1) = npsh(1) + nesh(1)
        nvec(2) = npsh(2) + nesh(2)
        write (fd, '(E24.16)') scf_energy

        allocate(real_buf(nsym * (nvec(1) + nvec(2)) * max(nbas(1), nbas(2)) * nz))
        allocate(int_buf (nsym * (nvec(1) + nvec(2))))

        ! copy MO coeff-s to the buffer and then flush to the disk
        offs = 0
        do isym = 1, nsym
            ncoef(isym) = nvec(isym) * nbas(isym) * nz
            do iz = 1, nz
            do i = 1, ncoef(isym)
                real_buf(offs + i) = cmo(isym, i, iz)
            end do
            end do
            offs = offs + ncoef(isym)
        end do
        write (fd, '(6f22.16)') (real_buf(i), i = 1, ncoef(1) + ncoef(2))

        ! copy orbital energies to the buffer and then flush to the disk
        offs = 0
        do isym = 1, nsym
            do i = 1, nvec(isym)
                real_buf(offs + i) = eps(isym, i)
            end do
            offs = offs + nvec(isym)
        end do
        write (fd, '(6E22.12)') (real_buf(i), i = 1, nvec(1) + nvec(2))

        ! the same for the 'ibeig' characteristic
        offs = 0
        do isym = 1, nsym
            do i = 1, nvec(isym)
                int_buf(offs + i) = ibeig(isym, i)
            end do
            offs = offs + nvec(isym)
        end do
        write (fd, '(66(1x,i0))') (int_buf(i), i = 1, nvec(1) + nvec(2))
        ! i don't know why 66

        close(unit = fd)
        print *, 'END OF FILE ', trim(filename)
        print *

        deallocate(real_buf)
        deallocate(int_buf)
        return

        13  continue
        print *, 'ERROR: CANNOT OPEN DFPCMO FILE ', trim(filename)

    end subroutine write_dfpcmo

    !
    ! Write some orbitals expanded in the AO basis to the DFPCMO-type
    ! formatted file named 'filename'
    !
    subroutine privec(nsym, npsh, nesh, nbas, eps, ibeig, cmo)

        implicit none
        integer :: fd
        integer :: i, j, isym
        ! orbital data: input
        integer, intent(in) :: nsym
        integer, intent(in) :: npsh(2), nesh(2), nbas(2)
        integer, dimension(:, :), intent(in) :: ibeig
        real(8), dimension(:, :), intent(in) :: eps
        real(8), dimension(:, :), intent(in) :: cmo
        ! local variables
        integer :: nvec(2)
        real(8), allocatable, dimension(:) :: real_buf
        integer, allocatable, dimension(:) :: int_buf
        integer :: nz = 1
        integer :: ncoef(2), offs, index, ispinor

        nvec(1) = npsh(1) + nesh(1)
        nvec(2) = npsh(2) + nesh(2)

        print *, 'nvec(1) = ', nvec(1)
        print *, 'nvec(2) = ', nvec(2)

        do isym = 1, nsym

            if (isym == 1) then
                print *
                print *, '           GERADE SPINORS '
                print *, '          ================'
                print *
            else
                print *
                print *, '          UNGERADE SPINORS'
                print *, '          ================'
                print *
            end if

            do ispinor = 1, nvec(isym)

                print '(a,i3,a,F22.12)', '[', ispinor, '] Eigenvalue', eps(isym,ispinor)
                print *, '==========================================================='
                do i = 1, nbas(isym)
                    index = nbas(isym) * (ispinor - 1) + i
                    print '(i3,E22.12,E22.12,E22.12,E22.12)', i, cmo(isym,index)
                end do
                print *
            end do

        end do

    end subroutine privec

end module dfpcmo
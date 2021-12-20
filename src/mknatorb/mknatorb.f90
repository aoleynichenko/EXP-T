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

! **********************************************************************
! program mknatorb
!
! Program for transformation of HF spinors to natural spinors.
! Input:
!   (1)  DIRAC formatted file with molecular spinor coefficients
!        (in AO basis) -- usually named DFPCMO
!   (2)  EXPT formatted file with natural orbitals (expanded in the
!        HF spinor basis)
! Output:
!   new DFPCMO-like DIRAC formatted file for subsequent visualization
!
! NOTE: this code is currently limited to real groups only!
!       (C2v, C2h, D2, D2h, Dinfh, Cinfv)
!
! 2019-2021 Alexander Oleynichenko
!
! **********************************************************************

program mknatorb

    use dfpcmo
    implicit none

    character(len = 128) :: dfpcmo_in_file
    character(len = 128) :: dfpcmo_out_file
    character(len = 128) :: natorb_file

    ! input orbitals & symmetry info
    integer :: nsym
    integer :: idim(3, 2)
    real(8) :: scf_energy
    integer :: npsh(2), nesh(2), nvec(2), nbas(2)
    integer, allocatable, dimension(:, :) :: ibeig
    real(8), allocatable, dimension(:, :) :: eps
    real(8), allocatable, dimension(:, :) :: cmo

    ! basis of molecular spinors (active spinors actually)
    integer :: n_act
    real(8), dimension(:), allocatable :: act_eps
    integer, dimension(:), allocatable :: act_nsym

    ! natural orbitals (expanded in the basis of molecular spinors)
    integer :: n_natorb
    real(8), dimension(:), allocatable :: occ_numbers
    complex(8), dimension(:, :), allocatable :: natorb_coef

    ! for NO transformation
    integer, dimension(:, :), allocatable :: map_expt_dirac ! EXP-T numbering -> DIRAC numbering
    real(8), allocatable, dimension(:, :) :: cno  ! NO coef-s in AO basis
    integer :: ino, i, iact, m, isym
    real(8) :: w
    integer :: ino_offs(2)
    integer, dimension(:), allocatable :: no_nsym
    integer :: actorb_beg, actorb_end, natorb_beg, natorb_end

    ! interfaces
    interface

        subroutine read_natorb(filename, n_act, act_eps, act_nsym, n_natorb, occ_numbers, no_nsym, natorb_coef)
            character(len = *), intent(in) :: filename
            integer, intent(out) :: n_act, n_natorb
            real(8), dimension(:), allocatable, intent(out) :: act_eps
            integer, dimension(:), allocatable, intent(out) :: act_nsym
            integer, dimension(:), allocatable, intent(out) :: no_nsym
            real(8), dimension(:), allocatable, intent(out) :: occ_numbers
            complex(8), dimension(:, :), allocatable, intent(out) :: natorb_coef
        end subroutine read_natorb

        subroutine construct_expt_dirac_mapping(map_expt_dirac, nsym, nvec, eps, n_act, act_nsym, act_eps)
            integer, dimension(:, :), allocatable :: map_expt_dirac ! EXP-T numbering -> DIRAC numbering
            integer :: nsym, n_act, nvec(2)
            real(8), dimension(:), allocatable :: act_eps
            integer, dimension(:), allocatable :: act_nsym
            real(8), allocatable, dimension(:, :) :: eps
        end subroutine construct_expt_dirac_mapping

    end interface

    print '(/20x,a/20x,a/)', 'PROGRAM MKNATORB    18 OCT 2021', '==============================='

    ! parse command-line arguments
    call getarg(1, dfpcmo_in_file)
    if (dfpcmo_in_file(1:1) == ' ') then
        print *, 'error: no 1st argument'
        print *, 'required: name of the DIRAC formatted file (input; DFPCMO-like file)'
        stop
    end if
    call getarg(2, natorb_file)
    if (natorb_file(1:1) == ' ') then
        print *, 'error: no 2nd argument'
        print *, 'required: name of the EXPT output -- formatted file with natural spinors'
        stop
    end if
    call getarg(3, dfpcmo_out_file)
    if (dfpcmo_out_file(1:1) == ' ') then
        print *, 'error: no 3rd argument'
        print *, 'required: name of the DIRAC formatted file (output; DFPCMO-like file)'
        stop
    end if
    print *, 'INPUT DFPCMO FILE   : ', trim(dfpcmo_in_file)
    print *, 'OUTPUT DFPCMO FILE  : ', trim(dfpcmo_out_file)
    print *, 'FILE WITH NAT ORB-S : ', trim(natorb_file)
    print *

    ! read input DFPCMO file
    ! = load basis of molecular spinors (expanded in the AO basis)
    call read_dfpcmo(dfpcmo_in_file, nsym, idim, scf_energy, cmo, eps, ibeig)

    npsh(1) = idim(1, 1)
    npsh(2) = idim(1, 2)
    nesh(1) = idim(2, 1)
    nesh(2) = idim(2, 2)
    nbas(1) = idim(3, 1)
    nbas(2) = idim(3, 2)
    nvec(1) = npsh(1) + nesh(1)
    nvec(2) = npsh(2) + nesh(2)

    ! read file with natural orbitals (NO) expansions
    call read_natorb(natorb_file, n_act, act_eps, act_nsym, n_natorb, occ_numbers, no_nsym, natorb_coef)
    call construct_expt_dirac_mapping(map_expt_dirac, nsym, nvec, eps, n_act, act_nsym, act_eps)


    ! мы должны скомбинировать данные спиноры с какими-то коэффициентами
    ! для этого нужно ассоциировать приведенные в этом файле номера спиноров
    ! со спинорами из полного списка со сквозной нумерацией
    ! 1. орбитальные индексы и энергии в DFPCMO определены для крамерс-
    !    ограниченных функций, для одного фактически партнера.
    !    А в файле NATORB орбиталей в два раза больше, т.к. в EXP-T
    !    крамерс-ограниченность снимается

    print *, 'TRANSFORMATION OF NATURAL SPINORS: MO -> AO BASIS...'
    allocate(cno(nsym, max(nvec(1), nvec(2)) * max(nbas(1), nbas(2))))
    cno = 0.0d0

    ino_offs = 1
    do ino = 1, n_natorb
        isym = no_nsym(ino)
        m = nbas(isym)

        natorb_beg = m * (ino_offs(isym) - 1) + 1
        natorb_end = m * ino_offs(isym)

        do iact = 1, n_act
            w = real(natorb_coef(ino, iact))
            i = map_expt_dirac(isym, iact)   ! get DIRAC's MO number 'i'
            if (i == 0) cycle  ! skip spinors with the other g/u symmetry

            actorb_beg = m * (i - 1) + 1
            actorb_end = m * i
            cno(isym, natorb_beg:natorb_end) = cno(isym, natorb_beg:natorb_end) + w * cmo(isym, actorb_beg:actorb_end)
        end do

        ino_offs(isym) = ino_offs(isym) + 1
    end do
    print *, 'DONE'
    print *

    ! write NOs to the output DFPCMO-type file
    ! with their occ numbers instead of orbital energies
    eps = 0.0d0
    do isym = 1, nsym
        i = 1
        do ino = 1, n_natorb
            if (no_nsym(ino) == isym) then
                eps(isym, i) = occ_numbers(ino)
                i = i + 1
            end if
        end do
    end do
    !call write_dfpcmo(dfpcmo_out_file, nsym, idim, scf_energy, cno, eps, ibeig)

    print *, 'MKNATORB TERMINATED NORMALLY'

    ! cleanup
    deallocate(map_expt_dirac)
    deallocate(cno)
    deallocate(cmo)
    deallocate(ibeig)
    deallocate(eps)

end program mknatorb


subroutine construct_expt_dirac_mapping(map_expt_dirac, nsym, nvec, eps, n_act, act_nsym, act_eps)

    implicit none

    ! arguments
    integer, dimension(:, :), allocatable :: map_expt_dirac ! EXP-T numbering -> DIRAC numbering
    integer :: nsym, n_act, nvec(2)
    real(8), dimension(:), allocatable :: act_eps
    integer, dimension(:), allocatable :: act_nsym
    real(8), dimension(:, :), allocatable :: eps

    ! local variables
    integer, dimension(:, :), allocatable :: marks
    integer :: iact, i, isym
    integer :: n_gerade, n_ungerade

    print *, 'CONSTRUCTION OF MAPPING FOR INDICES: EXP-T -> DIRAC...'

    allocate(map_expt_dirac(nsym, n_act))
    allocate(marks(nsym, max(nvec(1), nvec(2))))
    map_expt_dirac = 0
    marks = 0

    ! calculate number of gerade and ungerade spinors (in EXP-T)
    ! for groups with inversion:
    ! in EXP-T and DIRAC: first gerade, then ungerade
    ! we need one or two points where the 'marks' characteristic will be set to zero
    ! if nsym == 1: n_act/2+1
    ! if nsym == 2: n_gerade/2+1; n_gerade+n_ungerade/2+1
    if (nsym == 2) then
        n_gerade = count(act_nsym == 1)
        n_ungerade = count(act_nsym == 2)
        print *, 'N_GERADE   = ', n_gerade
        print *, 'N_UNGERADE = ', n_ungerade
        print *, 'N_ACTIVE   = ', n_act
    else
        n_gerade = n_act
        n_ungerade = 0
        print *, 'N_ACTIVE   = ', n_act
    end if

    ! construct mapping by comparison of one-electron spinor energies
    do iact = 1, n_act
        ! note: orbital energies in DIRAC are for Kramers pairs
        if (nsym == 1 .and. iact == n_act / 2 + 1) then
            marks = 0
        else if (nsym == 2 .and. (iact==n_gerade / 2 + 1 .or. iact==n_gerade + n_ungerade / 2 + 1)) then
            marks = 0
        end if
        isym = act_nsym(iact)
        do i = 1, nvec(isym)
            if (abs(eps(isym, i) - act_eps(iact)) < 1e-6 .and. marks(isym, i) /= 1) then
                map_expt_dirac(isym, iact) = i
                marks(isym, i) = 1
                exit
            end if
        end do
    end do

    do i = 1, n_act
        print '(i3,a,2i3)', i, ' -> ', map_expt_dirac(1:nsym, i)
    end do
    print *

    deallocate(marks)

end subroutine construct_expt_dirac_mapping



!
! subroutine read_natorb
!
! Reads natural orbitals expanded as linear combination of molecular
! spinors.
!
subroutine read_natorb(filename, n_act, act_eps, act_sym, n_natorb, occ_numbers, natorb_sym, natorb_coef)

    implicit none
    ! arguments
    character(len = *), intent(in) :: filename
    integer, intent(out) :: n_act, n_natorb
    real(8), dimension(:), allocatable, intent(out) :: act_eps
    integer, dimension(:), allocatable, intent(out) :: act_sym
    integer, dimension(:), allocatable, intent(out) :: natorb_sym
    real(8), dimension(:), allocatable, intent(out) :: occ_numbers
    complex(8), dimension(:, :), allocatable, intent(out) :: natorb_coef
    ! local variables
    character(len = 128) :: line, buf
    integer :: i, j, isym
    integer :: fd = 15
    real(8) :: real_buf_1, real_buf_2
    real(8) :: w, wmax
    integer :: n_gerade, n_ungerade

    print *, 'READING NATURAL ORBITALS FROM ', trim(filename)
    open(unit = fd, file = filename, status = 'old', err = 13)

    read (fd, *) buf, n_act
    read (fd, *)

    allocate(act_eps(n_act))
    allocate(act_sym(n_act))
    act_sym = 1
    do i = 1, n_act
        read (fd, '(a)'), line
        read (line, *) buf, act_eps(i), buf
        if (index(line(25:), "u") > 0) then
            act_sym(i) = 2
        end if
    end do

    read (fd, *) buf, n_natorb
    allocate(occ_numbers(n_natorb))
    allocate(natorb_coef(n_natorb, n_act))
    occ_numbers = 0
    natorb_coef = (0.0d0, 0.0d0)

    do i = 1, n_natorb
        read (fd, *) buf, occ_numbers(i)
        do j = 1, n_act
            read (fd, *) buf, real_buf_1, real_buf_2
            natorb_coef(i, j) = cmplx(real_buf_1, real_buf_2, kind = 8)
        end do
    end do

    ! try to determine symmetries of natural spinors
    ! symmetry of NO = symmetry of MO with the largest weight
    allocate(natorb_sym(n_natorb))
    natorb_sym = 1
    do i = 1, n_natorb
        wmax = 0.0d0
        do j = 1, n_act
            isym = act_sym(j)
            w = real(natorb_coef(i, j))
            if (w > wmax) then
                wmax = w
                natorb_sym(i) = isym
            end if
        end do
    end do

    n_gerade = count(natorb_sym == 1)
    n_ungerade = count(natorb_sym == 2)

    close(unit = fd)
    print *, 'NUMBER OF NAT SPINORS   = ', n_natorb
    print *, 'BASIS (MOL SPINORS DIM) = ', n_act
    print *, 'ENERGIES AND SYMMETRIES OF BASIS SET FUNCTIONS:'
    do i = 1, n_act
        print '(i3,f16.8,i5)', i, act_eps(i), act_sym(i)
    end do
    print *, 'OCC NUMBERS AND SYMMETRIES OF NATURAL SPINORS:'
    do i = 1, n_natorb
        print '(i3,f16.8,i5)', i, occ_numbers(i), natorb_sym(i)
    end do
    if (n_ungerade > 0 .and. n_gerade > 0) then
        print '(a,i0,5x,a,i0)', ' N(GERADE) = ', n_gerade, 'N(UNGERADE) = ', n_ungerade
    end if
    print *, 'END OF NATORB FILE'

    print *

    return

    13  continue
    print *, 'ERROR: UNABLE TO READ FILE WITH NATURAL ORBITALS ', trim(filename)
    stop

end subroutine read_natorb


! **********************************************************************
! subroutine estimate_arith
!
! Guess arithmetics: NZ = 1 (real), NZ = 2 (complex), NZ = 4 (quatern)
! Uses size of the file in bytes to estimate the number of real coeff-s.
! **********************************************************************
subroutine estimate_arith(filename, nz)

    implicit none
    ! arguments
    character(len = *), intent(in) :: filename
    integer, intent(out) :: nz
    ! local variables
    integer :: file_size
    integer :: fd, i, j
    integer :: nsym, idim(3, 2), nbas(2), nvec(2)
    integer :: ncoef
    real(8) :: approx_nz

    inquire(file = trim(filename), size = file_size)
    print '(3a,i0)', ' ESTIMATE NZ: SIZE OF FILE ', trim(filename), ' = ', file_size

    fd = 20
    open(unit = fd, file = trim(filename), status = 'old', form = 'formatted')

    idim = 0
    read (fd, *)
    read (fd, *) nsym, ((idim(i, j), i = 1, 3), j = 1, nsym)
    nbas(1) = idim(3, 1)
    nbas(2) = idim(3, 2)
    nvec(1) = idim(1, 1) + idim(2, 1)
    nvec(2) = idim(1, 2) + idim(2, 2)
    ncoef = nvec(1) * nbas(1) + nvec(2) * nbas(2)

    ! exclude title (74 chars), line with idim, E(SCF) -- first 3 lines
    ! (just approximately)
    file_size = file_size - 74 - (1 + 2 * nsym) * 3 - 24
    ! exclude lines with orbital energies (1 energy = 22 symbols)
    file_size = file_size - (nvec(1) + nvec(2)) * 22
    ! exclude lines with the 'ibeig' characteristc
    file_size = file_size - (nvec(1) + nvec(2)) * 2
    ! size of each coef-t in bytes = 22
    approx_nz = real(file_size) / 22.0 / ncoef
    ! nearest integer
    nz = nint(approx_nz)
    print '(a,f6.3,a,i3)', ' NZ = ', approx_nz, ' = ', nz

    close(unit = fd)

end subroutine estimate_arith







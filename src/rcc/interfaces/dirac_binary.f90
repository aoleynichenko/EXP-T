!
! EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
! Copyright (C) 2018-2024 The EXP-T developers.
!
! This file is part of EXP-T.
!
! EXP-T is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EXP-T is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
!
! E-mail:        exp-t-program@googlegroups.com
! Google Groups: https://groups.google.com/d/forum/exp-t-program
!

! ******************************************************************************
! Interface to the DIRAC quantum chemistry software.
! Reads DIRAC's unformatted files with transformed MO integrals, expand unique
! integrals to the full list and calculates some quantities from them
! (E(SCF), E(MP2), orbital energies -- all just for check).
!
! Integrals are obtained from the DIRAC MO integrals binary (unformatted) files.
! Required unformatted files (produced by DIRAC):
!   MRCONEE -- spinor info & transformed one-electron integrals
!   MDCINT  -- transformed two-electron integrals
!   (and optional)
!   MDPROP  -- transformed properties integrals
! In fact, they can have other names (but these names are the default ones).
! Binary interface should works much faster than the ASCII one.
! For details on format see source code of the DIRAC program package, files:
! src/moltra/traone.F, src/moltra/traout.F, src/relccsd/ccints.F.
!
! All the symmetries which are used by DIRAC (D2h and subgroups) are supported:
! (relativistic)     C1 Ci C2 Cs C2v D2 C2h D2h C32 D16h
! (non-relativistic) C1 Ci C2 Cs C2v D2 C2h D2h C32 D16h
! Since MRCONEE files don't contain any info about Hamiltonian, the problems
! with the Levy-Leblond Ham-n may occur in quaternion and complex symmetries.
!
! NOTES:
! ======
! (1) all arrays are one-dimensional and we must recompute "compound" indices
!     like (i,j,k,l) to linear indices every time when we want to access array
!     elements. This seems to be very inconvenient and slow, but it has to be
!     done for interoperability with codes written in C (no multidimensional
!     arrays are allowed).
!
! PAPERS:
! ======
! (1) L. Visscher, P.J.C. Aerts, H. Merenga, W.C. Nieuwpoort.
!     Relativistic Quantum Chemistry: the MOLFDIR program package.
!     Computer Physics Communications. 81. 120-144 (1994).
!     DOI: 10.1016/0010-4655(94)90115-5
! (2) Visscher, Lucas. Chapter 6 Post Dirac-Hartree-Fock methods -- electron
!     correlation. Theoretical and Computational Chemistry. 11. 291-331 (2002).
!     DOI: 10.1016/S1380-7323(02)80032-2
! (3) Dyall K., Faegri K. Introduction to Relativistic Quantum Chemistry.
!     Oxford University Press, 2007.
!
! ******************************************************************************


module general
    ! constants
    integer(4), parameter :: CC_MAX_SPINORS = 2048
    integer(4), parameter :: CC_MAX_NUM_IRREPS = 64
    integer(4), parameter :: CC_MAX_NPROP = 256
    integer(4), parameter :: INTCLASS_NO_BARS = 0  ! all indices are unbarred
    integer(4), parameter :: INTCLASS_ONE_BAR = 1  ! one barred index (quaternion groups only)
    integer(4), parameter :: INTCLASS_TWO_BARS = 2  ! two barred indices
    integer(4), parameter :: CC_ARITH_COMPLEX = 1
    integer(4), parameter :: CC_ARITH_REAL = 0
    ! common parameters
    integer(4) :: arith
    logical :: carith
end module general


module c_io

    interface
        integer(4) function io_open(path, mode) bind(C, name = "io_open")
            character, dimension(*) :: path, mode
        end function io_open
        integer(4) function io_close(fd) bind(C, name = "io_close")
            integer(4), value :: fd
        end function io_close
        integer(4) function io_write_compressed(fd, buf, sz) bind(C, name = "io_write_compressed")
            use iso_c_binding
            integer(4), value :: fd
            type(c_ptr), value :: buf
            integer(c_size_t), value :: sz
        end function io_write_compressed
        subroutine c_asctime() bind(C, name = "c_asctime")
        end subroutine c_asctime
    end interface

end module c_io



! ******************************************************************************
! module spinor_blocks
!
! spinors are divided int blocks (1-dimensional arrays of spinors).
! firstly spinors are divided by their symmetry, secondly they divided into
! smaller parts ("tiles") of length no more than 'tile_size'.
!
! for example (here tile_size = 5):
! 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 ...
! |------------ irrep1 ------------| |------- irrep2 ------|
! |---tile1---|  |---tile2---| |t3-| |----tile4---| |--t5--|
! here: n_spinor_blocks == 5
!       nsymrepa == 2
!       spinor_blocks_sizes = [5,5,2,5,3]
!       spinor_blocks_indices =
!         [[1,2,3,4,5],[6,7,8,9,10],[11,12],[13,14,15,16,17],[18,19,20]]
!       spinor_to_spinor_block =
!         [1,1,1,1,1, 2,2,2,2,2, 3,3, 4,4,4,4,4, 5,5,5]
! ******************************************************************************
module spinor_blocks

    use general
    implicit none

    ! total number of spinors
    integer(4) :: n_spinors
    ! max size of the single spinor range ("one tile")
    integer(4) :: tile_size
    ! number of fermion irreps in the Abelian subgroup
    integer(4) :: nsymrepa
    ! number of spinors in each (abelian) irrep
    integer(4), dimension(CC_MAX_NUM_IRREPS) :: repsizes
    ! number of spinor blocks
    integer(4) :: n_spinor_blocks
    ! sizes of spinor blocks
    integer(4), dimension(:), allocatable :: spinor_blocks_sizes
    ! indices of spinors -- for each spinor block
    integer(4), dimension(:, :), allocatable :: spinor_blocks_indices
    ! mapping: number of spinor -> number of spinor block
    integer(4), dimension(:), allocatable :: spinor_to_spinor_block

contains

    ! setup spinor blocks
    subroutine init_spinor_blocks(nspi, nsymrpa, irpamo, repanames, tilesz, multb)
        integer(4) :: nspi, nsymrpa, irpamo(nspi), tilesz
        integer(4) :: i, j, k, isb, sz
        integer(4), dimension(:), allocatable :: rep_indices   ! temporary array
        ! multiplication table for direct products in the Abelian subgroup
        integer(4), dimension(CC_MAX_NUM_IRREPS, CC_MAX_NUM_IRREPS) :: multb
        ! names of these irreps (Abelian subgroup)
        character(len = 4), dimension(CC_MAX_NUM_IRREPS) :: repanames

        ! init variables
        tile_size = tilesz
        n_spinors = nspi
        nsymrepa = nsymrpa
        do i = 1, nsymrepa
            repsizes(i) = 0
        end do
        allocate(rep_indices(n_spinors))
        n_spinor_blocks = 0
        ! calculate sizes of abelian reps
        do i = 1, n_spinors
            repsizes(irpamo(i)) = repsizes(irpamo(i)) + 1
        end do
        ! calculate number of spinor blocks
        do i = 1, nsymrepa
            n_spinor_blocks = n_spinor_blocks + repsizes(i) / tile_size
            if (mod(repsizes(i), tile_size) > 0) then
                n_spinor_blocks = n_spinor_blocks + 1
            end if
        end do
        ! allocate and init arrays
        allocate(spinor_blocks_sizes(n_spinor_blocks))
        allocate(spinor_blocks_indices(n_spinor_blocks, tile_size))
        allocate(spinor_to_spinor_block(n_spinors))
        spinor_blocks_sizes = 0
        spinor_blocks_indices = 0
        spinor_to_spinor_block = 0
        ! calculate sizes of all spinor blocks
        j = 1
        do i = 1, nsymrepa
            do k = 1, repsizes(i) / tile_size
                spinor_blocks_sizes(j) = tile_size
                j = j + 1
            end do
            if (mod(repsizes(i), tile_size) > 0) then
                spinor_blocks_sizes(j) = mod(repsizes(i), tile_size)
                j = j + 1
            end if
        end do
        ! find indices which belong to each spinor block
        isb = 1
        do i = 1, nsymrepa
            ! create list of all indices belonging to this irrep
            k = 1
            do j = 1, n_spinors
                if (irpamo(j) == i) then
                    rep_indices(k) = j
                    k = k + 1
                end if
            end do
            ! cut subarrays of indices (belonging to different spinor blocks)
            do k = 1, repsizes(i) / tile_size
                sz = tile_size
                spinor_blocks_indices(isb, 1:sz) = rep_indices((k - 1) * sz + 1:k * sz)
                isb = isb + 1
            end do
            if (mod(repsizes(i), tile_size) > 0) then
                sz = mod(repsizes(i), tile_size)
                spinor_blocks_indices(isb, 1:sz) = rep_indices((k - 1) * tile_size + 1:(k - 1) * tile_size + sz)
                isb = isb + 1
            end if
        end do
        ! construct mapping "spinor -> spinor block"
        do i = 1, n_spinor_blocks
            do j = 1, spinor_blocks_sizes(i)
                spinor_to_spinor_block(spinor_blocks_indices(i, j)) = i
            end do
        end do
        ! print summary
        !call print_spinor_blocks_info()

        ! cleanup
        deallocate(rep_indices)

    end subroutine init_spinor_blocks

    ! cleanup
    subroutine delete_spinor_blocks()
        deallocate(spinor_blocks_sizes)
        deallocate(spinor_blocks_indices)
        deallocate(spinor_to_spinor_block)
    end subroutine delete_spinor_blocks

    subroutine print_spinor_blocks_info()
        integer(4) :: i
        print *
        print *, '*** spinor blocks ***'
        print *, 'tile size = ', tile_size
        print *, 'number of fermion irreps in the Abelian subgroup = ', nsymrepa
        print '(a,32i3)', ' repsizes = ', (repsizes(i), i = 1, nsymrepa)
        print *, 'number of spinor blocks = ', n_spinor_blocks
        print '(a,32i3)', ' sizes of spinor blocks = ', (spinor_blocks_sizes(i), i = 1, n_spinor_blocks)
        print *, 'indices of spinors for each spinor block:'
        do i = 1, n_spinor_blocks
            print '(a,i3,a,50(20i5/))', '[', i, ']', spinor_blocks_indices(i, 1:spinor_blocks_sizes(i))
        end do
        print *, '*** end spinor blocks ***'
    end subroutine print_spinor_blocks_info

end module spinor_blocks


! ******************************************************************************
! module output_buffers
!
! for all possible non-zero 4d integral blocks we have an individual I/O buffer.
! only really used buffers are stored; if no buffer for integral block IJKL is
! provided, outbuf_idx(I,J,K,L) == 0.
! ******************************************************************************
module output_buffers

    implicit none

    ! integral output:
    ! size of output buffers
    integer(4), parameter :: BUF_SIZE = 16384
    ! number of output buffers
    integer(4) :: n_outbuf
    ! index of output buffers: mapping "sb1 sb2 sb3 sb4 -> seq number of buffer"
    integer(4), dimension(:, :, :, :), allocatable :: outbuf_idx
    ! output buffers:
    ! buffers used to accumulate integrals ...
    complex(8), dimension(:, :), allocatable, target :: buf_vint
    real(8), dimension(:, :), allocatable, target :: buf_vint_re
    ! ... and their indices
    integer(2), dimension(:, :), allocatable, target :: buf_indices
    ! current number of integrals in each buffer
    integer(4), dimension(:), allocatable, target :: buf_nint
    ! == 1 if the buffer was already used; else == 0
    logical, dimension(:), allocatable :: buf_flushed
    ! total number of bytes written
    integer(8) :: nbytes_written = 0
    ! prefix -- used in names of files with integrals
    character(len = 16) :: files_prefix

contains

    ! setup output buffers
    subroutine init_output_buffers(prefix, nspi, nsymrpa, irpamo, repanames, tilesz, multb)
        use general
        use spinor_blocks
        character(len = *), intent(in) :: prefix
        integer(4) :: nspi, nsymrpa, irpamo(nspi), tilesz
        integer(4) :: i1, i2, i3, i4
        integer(4) :: irep1, irep2, irep3, irep4
        ! multiplication table for direct products in the Abelian subgroup
        integer(4), dimension(CC_MAX_NUM_IRREPS, CC_MAX_NUM_IRREPS) :: multb
        ! names of these irreps (Abelian subgroup)
        character(len = 4), dimension(CC_MAX_NUM_IRREPS) :: repanames
        character(len = 4) :: repname1, repname2, repname3, repname4
        logical :: nonrel = .FALSE.

        ! check whether symmetry is relativistic or not
        repname1(1:4) = repanames(1)    ! detect non-rel-c symmetry
        if (repname1(4:4) == 'a') then
            nonrel = .TRUE.
        end if

        files_prefix = prefix

        ! setup output buffers
        n_outbuf = 0
        allocate(outbuf_idx(n_spinor_blocks, n_spinor_blocks, n_spinor_blocks, n_spinor_blocks))
        outbuf_idx = 0
        do i1 = 1, n_spinor_blocks
            do i2 = 1, n_spinor_blocks
                do i3 = 1, n_spinor_blocks
                    do i4 = 1, n_spinor_blocks
                        ! TODO: block selection (space + space + permutational symmetry)
                        ! selection: space symmetry
                        irep1 = irpamo(spinor_blocks_indices(i1, 1))   ! extract rep numbers
                        irep2 = irpamo(spinor_blocks_indices(i2, 1))
                        irep3 = irpamo(spinor_blocks_indices(i3, 1))
                        irep4 = irpamo(spinor_blocks_indices(i4, 1))
                        if (multb(irep1, irep2) /= multb(irep3, irep4)) then
                            cycle
                        end if
                        ! selection: space symmetry (if presented)
                        if (nonrel) then
                            repname1 = repanames(irep1)
                            repname2 = repanames(irep2)
                            repname3 = repanames(irep3)
                            repname4 = repanames(irep4)
                            if (.not. (repname1(4:4) == repname3(4:4) .and. &
                                    repname2(4:4) == repname4(4:4))) then
                                cycle
                            end if
                        end if

                        ! block is suitable, add it
                        n_outbuf = n_outbuf + 1
                        outbuf_idx(i1, i2, i3, i4) = n_outbuf
                    end do
                end do
            end do
        end do  ! end 4-fold loop over spinor blocks

        if (carith) then
            allocate(buf_vint(BUF_SIZE, n_outbuf))      ! I/O buffer for integrals
        else
            allocate(buf_vint_re(BUF_SIZE, n_outbuf))
        end if
        allocate(buf_indices(4 * BUF_SIZE, n_outbuf)) ! I/O buffer for quadruples of indices
        allocate(buf_nint(n_outbuf))               ! number of int-s in each buffer
        allocate(buf_flushed(n_outbuf))
        buf_nint = 0
        buf_flushed = .FALSE.

    end subroutine init_output_buffers


    ! cleanup output buffers
    subroutine delete_output_buffers()
        use general
        use spinor_blocks
        use iso_c_binding
        use c_io
        integer(4) :: fd, status
        integer(4), target :: nint
        integer(4) :: ibuf, sb1, sb2, sb3, sb4
        character(len = 100) :: filename
        integer(8) :: sizeof_int2 = 2
        integer(8) :: sizeof_int4 = 4
        integer(8) :: sizeof_complex8 = 16
        integer(8) :: sizeof_real8 = 8
        integer(4) :: stat
        integer(8) :: n_vint_files = 0

        do sb1 = 1, n_spinor_blocks
            do sb2 = 1, n_spinor_blocks
                do sb3 = 1, n_spinor_blocks
                    do sb4 = 1, n_spinor_blocks
                        ibuf = outbuf_idx(sb1, sb2, sb3, sb4)
                        if (ibuf == 0) then
                            cycle
                        end if

                        write (filename, "(a,a1,i0,a1,i0,a1,i0,a1,i0,a1)") trim(files_prefix), &
                                "-", sb1, "-", sb2, "-", sb3, "-", sb4
                        if (.not. buf_flushed(ibuf)) then
                            open(unit = 1234, iostat = stat, file = filename, status = 'old')
                            if (stat == 0) close(1234, status = 'delete')
                        end if
                        n_vint_files = n_vint_files + 1

                        call flush_buf(sb1, sb2, sb3, sb4)

                        ! terminate file with zero
                        nint = 0
                        fd = io_open(trim(filename) // C_NULL_CHAR, "a" // C_NULL_CHAR)
                        status = io_write_compressed(fd, c_loc(nint), sizeof_int4)
                        status = io_close(fd)
                    end do
                end do
            end do
        end do

        print *, 'number of VINT* files written ', n_vint_files
        print '(a,i0,a,f7.2,a)', ' written to disk: ', nbytes_written, &
                ' bytes = ', nbytes_written / (1024.0**3), ' Gb'

        deallocate(outbuf_idx)
        if (carith) then
            deallocate(buf_vint)
        else
            deallocate(buf_vint_re)
        end if
        deallocate(buf_indices)
        deallocate(buf_nint)
        deallocate(buf_flushed)

    end subroutine delete_output_buffers


    ! puts integral with indices <i1,i2|i3,i4> with value 'val'
    ! to the appropriate file VINT-*
    subroutine put_integral(i1, i2, i3, i4, val)
        use general
        use spinor_blocks
        integer(2) :: i1, i2, i3, i4
        complex(8) :: val
        integer(4) :: sb1, sb2, sb3, sb4
        integer(4) :: ibuf
        integer(4), target :: nint

        sb1 = spinor_to_spinor_block(i1)
        sb2 = spinor_to_spinor_block(i2)
        sb3 = spinor_to_spinor_block(i3)
        sb4 = spinor_to_spinor_block(i4)

        ibuf = outbuf_idx(sb1, sb2, sb3, sb4)
        if (ibuf == 0) then
            return
        end if
        nint = buf_nint(ibuf)

        ! flush buffer if needed
        ! "VINT-*" file must be truncated if opened for the first time
        if (nint == BUF_SIZE) then
            call flush_buf(sb1, sb2, sb3, sb4)
            nint = 0
        end if
        ! write new data to the buffer
        nint = nint + 1
        buf_indices((nint - 1) * 4 + 1, ibuf) = i1
        buf_indices((nint - 1) * 4 + 2, ibuf) = i2
        buf_indices((nint - 1) * 4 + 3, ibuf) = i3
        buf_indices((nint - 1) * 4 + 4, ibuf) = i4
        if (carith) then
            buf_vint(nint, ibuf) = val
        else
            buf_vint_re(nint, ibuf) = real(val)
        end if
        buf_nint(ibuf) = nint

    end subroutine put_integral

    subroutine flush_buf(sb1, sb2, sb3, sb4)
        use general
        use iso_c_binding
        use c_io
        integer(4), intent(in) :: sb1, sb2, sb3, sb4
        integer(4) :: ibuf
        integer(4), target :: nint
        character(len = 100) :: filename
        integer(8), parameter :: sizeof_int2 = 2
        integer(8), parameter :: sizeof_int4 = 4
        integer(8), parameter :: sizeof_complex8 = 16
        integer(8), parameter :: sizeof_real8 = 8
        integer(4) :: stat
        integer(4) :: fd

        ibuf = outbuf_idx(sb1, sb2, sb3, sb4)
        if (ibuf == 0) then
            return
        end if
        nint = buf_nint(ibuf)

        write (filename, "(a,a1,i0,a1,i0,a1,i0,a1,i0,a1)") trim(files_prefix), "-", sb1, "-", sb2, "-", sb3, "-", sb4
        if (.not. buf_flushed(ibuf)) then
            open(unit = 1234, iostat = stat, file = filename, status = 'old')
            if (stat == 0) close(1234, status = 'delete')
        end if
        buf_flushed(ibuf) = .TRUE.
        fd = io_open(trim(filename) // C_NULL_CHAR, "a" // C_NULL_CHAR)
        stat = io_write_compressed(fd, c_loc(nint), sizeof_int4)
        stat = io_write_compressed(fd, c_loc(buf_indices(:, ibuf)), 4 * nint * sizeof_int2)
        nbytes_written = nbytes_written + sizeof_int4 + 4 * nint * sizeof_int2
        if (carith) then
            stat = io_write_compressed(fd, c_loc(buf_vint(:, ibuf)), nint * sizeof_complex8)
            nbytes_written = nbytes_written + nint * sizeof_complex8
        else
            stat = io_write_compressed(fd, c_loc(buf_vint_re(:, ibuf)), nint * sizeof_real8)
            nbytes_written = nbytes_written + nint * sizeof_real8
        end if
        stat = io_close(fd)
        buf_nint(ibuf) = 0

    end subroutine flush_buf

end module output_buffers


module dirac_32_64_compatibility
    use general
    implicit none
    integer(4), parameter :: CC_DIRAC_INT4 = 4
    integer(4), parameter :: CC_DIRAC_INT8 = 8
    integer(4) :: dirac_integer_size

    logical(8) :: breit8
    integer(8) :: nspinors8, invsym8, nz_arith8, is_spinfree8, norb_total8
    integer(8) :: nactive8(8), nstr8(2), nfrozen8(2, 0:2), ndelete8(8)
    integer(8) :: nsymrp8, nsymrpa8, norb8(2), nbsymrp8
    integer(8) :: multb8(CC_MAX_NUM_IRREPS, CC_MAX_NUM_IRREPS)
    integer(8), dimension(CC_MAX_SPINORS) :: irpmo8, irpamo8, ibspi8, iocc8
    integer(8) :: nkr8, ikr8, jkr8, nonzr8
    integer(8), dimension(:), allocatable :: kr8
    integer(8), dimension(:), allocatable :: indk8, indl8

contains

    subroutine test_dirac_integer_size()
        integer(4) :: luone
        character(len = 1024) :: oneel_file   ! passed via the /mrconee_mdcint/ common
        character(len = 1024) :: twoel_file
        character(len = 1024) :: prop_file
        character(len = 1024) :: gaunt_file
        integer(8) :: gaunt_enabled
        integer(8) :: n_twoprop
        character(len = 1024) :: twoprop_files(CC_MAX_NPROP)
        common /mrconee_mdcint/ oneel_file, twoel_file, prop_file, gaunt_file, gaunt_enabled, n_twoprop, twoprop_files
        bind(C) :: /mrconee_mdcint/
        logical(4) :: breit
        integer(4) :: NSPINORS, invsym, nz_arith, is_spinfree, norb_total
        real(8) :: enuc, escf

        luone = 42
        open(unit = luone, file = trim(oneel_file), form = 'unformatted')
        read (luone, err = 2019, end = 2019) nspinors, breit, enuc, invsym, nz_arith, &
                is_spinfree, norb_total, escf
        if (nspinors <= 0) goto 2019
        if (.not. (invsym == 1 .or. invsym == 2)) goto 2019
        if (.not. (nz_arith == 1 .or. nz_arith == 2 .or. nz_arith == 4)) goto 2019
        if (.not. (is_spinfree == -1 .or. is_spinfree == 0 .or. is_spinfree == 1)) goto 2019
        if (norb_total <= 0) goto 2019

        dirac_integer_size = CC_DIRAC_INT4
        close(luone)
        return
        2019 continue
        dirac_integer_size = CC_DIRAC_INT8
        close(luone)
    end subroutine test_dirac_integer_size
end module dirac_32_64_compatibility


! ******************************************************************************
! subroutine dirac_interface_binary
!
! "main" subroutine of this file. subroutine is invoked from the C code (file
! dirac_interface.c).
! reads info about 1-particle functions from the MRCONEE file,
! splits the MDCINT file into integral blocks (files VINT-*)
! ******************************************************************************
subroutine dirac_interface_binary(err) bind(C)

    use iso_c_binding
    use general
    use dirac_32_64_compatibility
    use spinor_blocks
    use output_buffers
    implicit none
    interface
        subroutine flush_fock(nspinors, hint)
            integer(4) :: NSPINORS
            complex(8), dimension(nspinors**2), target :: hint
        end subroutine flush_fock
    end interface

    ! arguments
    integer(4), intent(inout) :: err
    character(len = 1024) :: oneel_file   ! passed via the /mrconee_mdcint/ common
    character(len = 1024) :: twoel_file
    character(len = 1024) :: prop_file
    character(len = 1024) :: gaunt_file
    integer(8) :: gaunt_enabled
    integer(8) :: n_twoprop
    character(len = 1024) :: twoprop_files(CC_MAX_NPROP)
    common /mrconee_mdcint/ oneel_file, twoel_file, prop_file, gaunt_file, gaunt_enabled, n_twoprop, twoprop_files
    bind(C) :: /mrconee_mdcint/

    ! files
    logical :: exists
    integer :: luone   ! 1-el integrals file
    integer :: luint   ! 2-el integrals file

    ! MRCONEE content
    ! basic info
    integer(4) :: NSPINORS
    logical(4) :: breit
    real(8) :: enuc
    integer(4) :: invsym
    integer(4) :: nz_arith
    integer(4) :: is_spinfree
    integer(4) :: norb_total
    real(8) :: escf
    ! NORDER gives the order of the Abelian subgroup used linear molecules.
    ! NORDER=16 (default) uses C_16_h for D_inf_h and C_32 for C_inf_v.
    ! (from DIRAC's moltra/traone.F)
    integer(4), parameter :: NORDER = 16
    ! information about fermion irreps
    integer(4) :: nsymrp      ! number of fermion irreps in parent group
    character(len = 14), dimension(8) :: repnames  ! their names (type: string)
    integer(4), dimension(8) :: nactive          ! number of spinors active in the transf-n (not valence!)
    integer(4), dimension(2) :: nstr             ! total number of orbitals of this ircop
    integer(4), dimension(2, 0:2) :: nfrozen      ! number of occupied frozen (core) spinors
    integer(4), dimension(8) :: ndelete          ! number of deleted spinors
    ! for abelian subgroups
    integer(4) :: nsymrpa                        ! number of fermion irreps in the Abelian subgroup
    character(len = 4), dimension(4 * NORDER) :: repanames ! names of these irreps
    ! Multiplication table for direct products in the Abelian subgroup
    integer(4), dimension(CC_MAX_NUM_IRREPS, CC_MAX_NUM_IRREPS) :: multb
    ! for each spinor:
    integer(4), dimension(CC_MAX_SPINORS) :: irpmo   ! Irrep in parent group (1:gerade, 2:ungerade)
    integer(4), dimension(CC_MAX_SPINORS) :: irpamo  ! Irrep in Abelian subgroup
    real(8), dimension(CC_MAX_SPINORS) :: eorbmo  ! Orbital energy taken from DHF
    integer(4), dimension(CC_MAX_SPINORS) :: ibspi   ! Approximate boson irrep identification (needed in LUCIAREL and LUCITA)
    integer(4) :: nbsymrp                            ! Number of boson symmetry reps (for LUCITA)
    integer(4), dimension(2) :: norb    ! (?) number of Kramers pairs in each u/g irrep
    integer(4), dimension(CC_MAX_SPINORS) :: iocc    ! spinor occupation numbers
    ! core Fock matrix
    complex(8), dimension(:), allocatable :: hint

    ! auxiliary
    integer(4) :: i, j, irp, irpa
    integer(4) :: place4
    integer(8), parameter :: sizeof_int4 = 4
    integer(4) :: tile_sz
    integer(4) :: recommended_carith
    integer(4) :: skip_mdcint
    integer(4), parameter :: STDOUT = 6
    character(len = 256) :: twoprop_file_name

    common /dirac_data/ nspinors, breit, enuc, invsym, nz_arith, is_spinfree, norb_total, escf, &
            nsymrp, repnames, nactive, nstr, nfrozen, ndelete, nsymrpa, repanames, multb, nbsymrp, &
            irpmo, irpamo, ibspi, iocc, place4, eorbmo
    bind(C) :: /dirac_data/
    common /tilesize/ tile_sz
    bind(C) :: /tilesize/
    common /recommendedcarith/ recommended_carith
    bind(C) :: /recommendedcarith/
    common /skip_sorting/ skip_mdcint
    bind(C) :: /skip_sorting/

    print *
    print *, '******************************************************************'
    print *, '**       BINARY INTERFACE TO THE DIRAC PROGRAM PACKAGE          **'
    print *, '**                    version 10 Apr 2020                       **'
    print *, '******************************************************************'
    print *
    print *, 'required unformatted files (produced by DIRAC):'
    print *, '  MRCONEE -- spinor info & transformed one-electron integrals'
    print *, '  MDCINT  -- transformed two-electron integrals'
    print *, 'optional files (DIRAC):'
    print *, '  MDPROP  -- transformed integrals of one-electron properties'
    print *, 'optional files (D. E. Maison, L. V. Skripnikov):'
    print *, '  MGINT   -- transformed Gaunt integrals'
    print *
    print *, 'HINT: in order to obtain these files, add option .4INDEX to the'
    print *, '  **DIRAC section of the DIRAC input file; this will enable the'
    print *, '  code for integrals transformation. Then run:'
    print *, '  $ pam --inp=<inp-file> --mol=<mol-file> --get="MRCONEE MDCINT"'
    print *
    print *

    print *, 'Names of integral files:'
    print *, 'MRCONEE = ', trim(oneel_file)
    print *, 'MDCINT  = ', trim(twoel_file)
    print *, 'MDPROP  = ', trim(prop_file)
    if (gaunt_enabled /= 0) then
        print *, 'MGINT   = ', trim(gaunt_file)
    end if
    if (n_twoprop > 0) then
        do i = 1, n_twoprop
            print *, 'two-electron property = ', trim(twoprop_files(i))
        end do
    end if

    ! check binary files for existance
    inquire(file = trim(oneel_file), exist = exists)
    if (.not. exists) then
        print *, 'MRCONEE unformatted file is not found!'
        err = 1
        return
    end if
    inquire(file = trim(twoel_file), exist = exists)
    if (.not. exists) then
        print *, 'MDCINT unformatted file is not found!'
        err = 1
        return
    end if
    inquire(file = trim(prop_file), exist = exists)
    if (.not. exists) then
        print *, 'MDPROP unformatted file is not found!'
        print *, 'will be continued without properties'
    end if

    ! try to determine size of integer numbers in DIRAC's binary files
    call test_dirac_integer_size()
    if (dirac_integer_size == CC_DIRAC_INT4) then
        print *, 'default integer type in DIRAC              integer(4)'
    else
        print *, 'default integer type in DIRAC              integer(8)'
    end if

    ! read unformatted file with symmetry & spinor information and 1-el integrals.
    print *
    print *, '*** MRCONEE FILE ***'
    luone = 42
    open(unit = luone, file = trim(oneel_file), form = 'unformatted')

    ! <> total number of spinors
    ! <> breit active in DHF
    ! <> core energy (inactive energy + nuclear repulsion)
    ! <> inversion symmetry (yes : 2; no : 1)
    ! <> group type (1 real, 2 complex, 4 quaternion)
    ! <> spinfree formalism
    ! <> total number of orbitals (so including frozen or deleted orbitals)
    ! <> total SCF energy
    if (dirac_integer_size == CC_DIRAC_INT4) then
        read (luone, err = 1299, end = 1299) nspinors, breit, enuc, invsym, nz_arith, &
                is_spinfree, norb_total, escf
    else  ! 8-byte integers in DIRAC
        read (luone, err = 1299, end = 1299) nspinors8, breit8, enuc, invsym8, nz_arith8, &
                is_spinfree8, norb_total8, escf
        nspinors = int(nspinors8, 4)
        breit = logical(breit8, 4)
        invsym = int(invsym8, 4)
        nz_arith = int(nz_arith8, 4)
        is_spinfree = int(is_spinfree8, 4)
        norb_total = int(norb_total8, 4)
    end if
    ! in some cases spinfree maybe == -1 instead of 1 (???)
    if (is_spinfree /= 0) then
        is_spinfree = 1
    end if
    if (is_spinfree == 1 .or. nz_arith == 1) then
        arith = CC_ARITH_REAL
        carith = .FALSE.
    else
        arith = CC_ARITH_COMPLEX
        carith = .TRUE.
    end if
    ! complex arithmetic can be REQUIRED by the user
    if (recommended_carith == 1) then
        arith = CC_ARITH_COMPLEX
        carith = .TRUE.
    end if
    print *, 'NSPINORS                                   ', nspinors
    print *, 'was breit in DHF                           ', breit
    print *, 'nuclear repulsion energy                   ', enuc
    print *, 'inversion symmetry (1-no,2-yes)            ', invsym
    print *, 'group type (1-real,2-cmplx,4-quater)       ', nz_arith
    print *, 'is spinfree                                ', is_spinfree
    if (carith) then
        print *, 'arithmetic                                 complex'
    else
        print *, 'arithmetic                                 real'
    end if
    print *, 'total num of orb-s (+frozen+deleted)       ', norb_total
    print *, 'Total SCF energy =                         ', escf

    ! <> [nsymrp]    number of fermion irreps in parent group
    ! <> [repnames]  names of these irreps (gerade, ungerade)
    ! <> [nactive]   number of spinors active in the transf-n (not valence!)
    ! <> [nstr]      total number of orbitals of this ircop
    ! <> [nfrozen]   number of occupied frozen (core) spinors,
    !                 0: total
    !                 1: positive energy
    !                 2: negative energy
    ! <> [ndelete]   number of deleted spinors.
    if (dirac_integer_size == CC_DIRAC_INT4) then
        read (luone) nsymrp, &
                (repnames(irp), irp = 1, nsymrp), &
                (nactive(irp), irp = 1, nsymrp), &
                (nstr(irp), irp = 1, invsym), &
                (nfrozen(irp, 0), irp = 1, invsym), &
                (nfrozen(irp, 1), irp = 1, invsym), &
                (nfrozen(irp, 2), irp = 1, invsym), &
                (ndelete(irp), irp = 1, invsym)
    else  ! 8-byte integers in DIRAC
        read (luone) nsymrp8, &
                (repnames(irp), irp = 1, nsymrp8), &
                (nactive8(irp), irp = 1, nsymrp8), &
                (nstr8(irp), irp = 1, invsym), &
                (nfrozen8(irp, 0), irp = 1, invsym), &
                (nfrozen8(irp, 1), irp = 1, invsym), &
                (nfrozen8(irp, 2), irp = 1, invsym), &
                (ndelete8(irp), irp = 1, invsym)
        nsymrp = nsymrp8
        nactive = nactive8
        nstr = nstr8
        nfrozen = nfrozen8
        ndelete = ndelete8
    end if
    print *, 'number of fermion irreps in parent group   ', nsymrp
    print *, 'names of these reps (grd, ungrd)           ', (repnames(irp), irp = 1, nsymrp)
    print *, 'number of spinors active in the transf-n   ', (nactive(irp), irp = 1, nsymrp)
    print *, 'total number of orb-s of this ircop        ', (nstr(irp), irp = 1, invsym)
    print *, 'number of occupied frozen (core) spinors   '
    print *, '  - total                                  ', (nfrozen(irp, 0), irp = 1, invsym)
    print *, '  - positive energy                        ', (nfrozen(irp, 1), irp = 1, invsym)
    print *, '  - negative energy                        ', (nfrozen(irp, 2), irp = 1, invsym)
    print *, 'number of deleted spinors                  ', (ndelete(irp), irp = 1, invsym)

    ! number of fermion irreps in the Abelian subgroup,
    ! names of these irreps
    if (dirac_integer_size == CC_DIRAC_INT4) then
        read (luone) nsymrpa, &
                (repanames(irpa), irpa = 1, nsymrpa * 2)
    else  ! 8-byte integers in DIRAC
        read (luone) nsymrpa8, &
                (repanames(irpa), irpa = 1, nsymrpa8 * 2)
        nsymrpa = nsymrpa8
    end if
    print *, 'number of fermion irreps in Abelian subgrp ', nsymrpa
    print *, 'names of these irreps                      ', (repanames(irpa), irpa = 1, nsymrpa * 2)

    ! Multiplication table for direct products in the
    ! Abelian subgroup
    if (dirac_integer_size == CC_DIRAC_INT4) then
        read (luone) ((multb(i, j), i = 1, 2 * nsymrpa), j = 1, 2 * nsymrpa)
    else  ! 8-byte integers in DIRAC
        read (luone) ((multb8(i, j), i = 1, 2 * nsymrpa), j = 1, 2 * nsymrpa)
        multb = multb8
    end if
    ! Information for each spinor :
    !  - Irrep in parent group (1:gerade, 2:ungerade)
    !  - Irrep in Abelian subgroup
    !  - Orbital energy taken from DHF
    !  - Approximate boson irrep identification (needed in LUCIAREL and LUCITA)
    !  - Number of boson symmetry reps (for LUCITA)
    if (dirac_integer_size == CC_DIRAC_INT4) then
        read (luone) (irpmo(i), irpamo(i), eorbmo(i), i = 1, nspinors), &
                (ibspi(i), i = 1, nspinors), &
                (norb(i), i = 1, invsym), &
                nbsymrp
    else  ! 8-byte integers in DIRAC
        read (luone) (irpmo8(i), irpamo8(i), eorbmo(i), i = 1, nspinors), &
                (ibspi8(i), i = 1, nspinors), &
                (norb8(i), i = 1, invsym), &
                nbsymrp8
        irpmo = irpmo8
        irpamo = irpamo8
        ibspi = ibspi8
        norb = norb8
        nbsymrp = nbsymrp8
    end if
    ! set occupation number for every spinor
    iocc = 0
    do i = 1, nspinors
        if (i == 1) then   ! if first rep
            ! make next nactive(irpmo(i)) spinors occupied
            do j = 0, nactive(irpmo(i)) - 1
                iocc(i + j) = 1
            end do
        else if (i > 1 .and. irpmo(i) /= irpmo(i - 1)) then   ! if rep changes
            ! make next nactive(irpmo(i)) spinors occupied
            do j = 0, nactive(irpmo(i)) - 1
                iocc(i + j) = 1
            end do
        end if
    end do

    !print *, 'information for each spinor:'
    !print *, '(i, irrep in parent grp, irrep in abelian subgrp, occ number,'
    !print *,  'orb energy, approx boson irrep identification)'
    !do i = 1, NSPINORS
    !    write (*,'(4i4,f15.8,i4)') i, irpmo(i), irpamo(i), iocc(i), eorbmo(i), ibspi(i)
    !end do
    print *, 'number of g/u Kramers pairs', (norb(i), i = 1, invsym)
    print *, 'number of boson symmetry reps(LUCITA) ', nbsymrp

    ! Matrix elements over the core Fock operator
    allocate(hint(nspinors**2))
    read (luone) ((hint(nspinors * (i - 1) + j), j = 1, nspinors), i = 1, nspinors)
    ! write core Fock operator to disk
    call flush_fock(nspinors, hint)
    deallocate(hint)

    close(unit = luone)
    print *, '*** END OF MRCONEE FILE ***'

    tile_size = tile_sz
    call init_spinor_blocks (nspinors, nsymrpa, irpamo, repanames, tile_size, multb)

    call read_mdprop()
    call expect_dipole(iocc)

    if (skip_mdcint /= 0) then
        print *, 'reading and sorting of MDCINT file will be skipped'
        if (gaunt_enabled /= 0) then
            print *, 'reading and sorting of MGINT file will be skipped'
        end if
    else
        ! Coulomb two-electron integrals
        call init_output_buffers("VINT", nspinors, nsymrpa, irpamo, repanames, tile_size, multb)
        call read_mdcint(nz_arith, is_spinfree, nspinors)
        call delete_output_buffers()

        ! Gaunt two-electron integrals (if required)
        if (gaunt_enabled /= 0) then
            call init_output_buffers("GINT", nspinors, nsymrpa, irpamo, repanames, tile_size, multb)
            call read_two_electron_prop("GAUNT", gaunt_file)
            call delete_output_buffers()
        end if

        do i = 1, n_twoprop
            write (twoprop_file_name, '(a,i0)') 'TWOPROP', i
            call init_output_buffers(trim(twoprop_file_name), nspinors, nsymrpa, irpamo, repanames, tile_size, multb)
            call read_two_electron_prop(trim(twoprop_file_name), trim(twoprop_files(i)))
            call delete_output_buffers()
        end do
    end if

    call delete_spinor_blocks()
    call flush(STDOUT)

    err = 0
    return

    ! error handling
    1299 continue
    print *, 'error reading file MRCONEE'
    err = 1
    return

contains

end subroutine dirac_interface_binary


!*******************************************************************************
! subroutine flush_fock
!
! Flushes core Fock matrix to the binary file (via Unix syscalls)
!*******************************************************************************
subroutine flush_fock(nspinors, hint)

    use iso_c_binding
    use c_io

    integer(4) :: NSPINORS
    complex(8), dimension(nspinors**2), target :: hint

    integer(4) :: fd
    integer(4) :: s   ! return status
    integer(8), parameter :: sizeof_complex8 = 16

    fd = io_open("HINT" // C_NULL_CHAR, "w" // C_NULL_CHAR)
    s = io_write_compressed(fd, c_loc(hint), sizeof_complex8 * nspinors**2)
    s = io_close(fd)

end subroutine flush_fock


!*******************************************************************************
! subroutine read_two_electron_prop
!
! read two-electron property integrals from the unformatted file
! generated by the program of D. E. Maison and L. V. Skripnikov
!
! Maison's order of spinor indices:
!   i1 i4 i3 i2
! then the corresponding integral <12|34>:
!   int dx1 dx2 phi_i1*(1) phi_i2*(2) phi_i3(1) phi_i4(2) / r_12
!*******************************************************************************
subroutine read_two_electron_prop(prop_name, file_name)

    use iso_c_binding
    use dirac_32_64_compatibility
    use spinor_blocks
    use output_buffers
    use c_io
    implicit none

    ! argument
    character(len = *) :: prop_name
    character(len = *) :: file_name
    ! names of integral files
    character(len = 1024) :: oneel_file   ! passed via the /mrconee_mdcint/ common
    character(len = 1024) :: twoel_file
    character(len = 1024) :: prop_file
    character(len = 1024) :: gaunt_file
    integer(8) :: gaunt_enabled
    integer(8) :: n_twoprop
    character(len = 1024) :: twoprop_files(CC_MAX_NPROP)
    common /mrconee_mdcint/ oneel_file, twoel_file, prop_file, gaunt_file, gaunt_enabled, n_twoprop, twoprop_files
    bind(C) :: /mrconee_mdcint/

    ! local varibles
    integer(2) :: i2, j2, k2, l2
    complex(8) :: zgint
    real(8) :: gint
    integer :: luint   ! 2-el property integrals file

    print *
    write (*, '(3a)') ' *** ', trim(prop_name), ' FILE ***'
    call c_asctime()
    luint = 43
    open(unit = luint, file = trim(file_name), form = 'unformatted', access = 'stream')
    rewind(luint)

    10  continue
    read (luint, end = 11) i2, j2, k2, l2, gint
    zgint = cmplx(gint,0.0d0)
    call put_integral(i2, l2, k2, j2, zgint)
    goto 10
    11  continue

    call c_asctime()
    write (*, '(3a)') ' *** END OF ', trim(prop_name), ' FILE ***'
    print *

end subroutine read_two_electron_prop


!*******************************************************************************
! subroutine read_mdcint
!
! read two-electron Coulomb integrals from the MDCINT unformatted file
! (MOLFDIR format by L. Visscher is used).
! For details and more explanations, see DIRAC's source code (relccsd/ccints.F)
!*******************************************************************************
subroutine read_mdcint(nz_arith, is_spinfree, nspinors)

    use iso_c_binding
    use dirac_32_64_compatibility
    use spinor_blocks
    use output_buffers
    use c_io
    implicit none

    ! arguments
    integer(4), intent(in) :: nz_arith, is_spinfree, NSPINORS

    integer(4) :: err
    integer(4) :: ikr, jkr, inz
    integer(8) :: i

    ! names of integral files
    character(len = 1024) :: oneel_file   ! passed via the /mrconee_mdcint/ common
    character(len = 1024) :: twoel_file
    character(len = 1024) :: prop_file
    character(len = 1024) :: gaunt_file
    integer(8) :: gaunt_enabled
    integer(8) :: n_twoprop
    character(len = 1024) :: twoprop_files(CC_MAX_NPROP)
    common /mrconee_mdcint/ oneel_file, twoel_file, prop_file, gaunt_file, gaunt_enabled, n_twoprop, twoprop_files
    bind(C) :: /mrconee_mdcint/
    integer :: luint   ! 2-el integrals file

    ! MDCINT contents
    character(len = 10) :: datex, timex*8
    integer(4) :: nkr
    integer(4), dimension(:), allocatable :: kr
    integer(4) :: nonzr
    integer(4), dimension(:), allocatable :: indk, indl
    real(8), dimension(:, :), allocatable :: cbuf

    print *
    print *, '*** MDCINT FILE ***'
    call c_asctime()
    luint = 43
    open(unit = luint, file = trim(twoel_file), form = 'unformatted')
    rewind(luint)

    allocate(kr(-nspinors:nspinors))
    allocate(indk(nspinors * nspinors))
    allocate(indl(nspinors * nspinors))
    allocate(cbuf(2, nspinors * nspinors))
    ! if 8-byte integers in DIRAC
    if (dirac_integer_size == CC_DIRAC_INT8) then
        allocate(kr8(-nspinors:nspinors))
        allocate(indk8(nspinors * nspinors))
        allocate(indl8(nspinors * nspinors))
    end if


    ! READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,(KR(I),KR(-I),I=1,NKR)
    if (dirac_integer_size == CC_DIRAC_INT4) then
        read (luint, err = 1300, end = 1300) datex, timex, nkr, (kr(i), kr(-i), i = 1, nkr)
    else  ! 8-byte integers in DIRAC
        read (luint, err = 1300, end = 1300) datex, timex, nkr8, (kr8(i), kr8(-i), i = 1, nkr8)
        nkr = nkr8
        kr = kr8
    end if
    if (2 * nkr /= nspinors) then
        print *, 'inconsistent MRCONEE and MDCINT files'
        err = 1
        return
    end if
    print *, 'datex                                 ', datex
    print *, 'timex                                 ', timex
    print *, 'number of Kramers pairs               ', nkr

    !print *, 'Indices of [unbarred] and [barred] spinors'
    !print *, '  kr(-nkr:nkr) = [fmt=i,kr(i),kr(-i)]   '
    !do i = 1, nkr
    !    print *, i, kr(i), kr(-i)
    !end do

    ! read two-electron integrals.
    ! A number of different formats are possible (see DIRAC/src/relccsd/ccints.F).

    ! IF (NZ.EQ.1) THEN
    !  READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
    !       (INDK(INZ),INDL(INZ),INZ=1,NONZR),
    !       (CBUF(1,INZ),INZ=1,NONZR)
    ! real arithmetics
    if (nz_arith == 1 .or. is_spinfree == 1) then
        ! read integral block until (ikr == 0 .and. jkr == 0)
        100     continue
        if (dirac_integer_size == CC_DIRAC_INT4) then
            read (luint, end = 1301, err = 1302) ikr, jkr, nonzr, &
                    (indk(inz), indl(inz), inz = 1, nonzr), &
                    (cbuf(1, inz), inz = 1, nonzr)
        else  ! 8-byte integers in DIRAC
            read (luint, end = 1301, err = 1302) ikr8, jkr8, nonzr8, &
                    (indk8(inz), indl8(inz), inz = 1, nonzr8), &
                    (cbuf(1, inz), inz = 1, nonzr8)
            ikr = ikr8
            jkr = jkr8
            nonzr = nonzr8
            indk = indk8
            indl = indl8
        end if

        ! it was the last entry, break loop
        if (ikr == 0 .and. jkr == 0) goto 110
        ! if no non-zero elements in this batch => goto the next batch
        if (nonzr == 0) goto 100

        ! set imaginary parts strictly to zero
        do inz = 1, nonzr
            cbuf(2, inz) = 0.0d0
        end do

        call expand_ints(ikr, jkr, nonzr, indk, indl, cbuf)

        ! read next batch
        goto 100
        110     continue

        ! ELSE
        !  READ (MDINT,END=10010,ERR=10020) IKR,JKR,NONZR,
        !       (INDK(INZ),INDL(INZ),INZ=1,NONZR),
        !       (CBUF(1,INZ),CBUF(2,INZ),INZ=1,NONZR)
        ! complex arithmetics
    else
        ! read integral block until (ikr == 0 .and. jkr == 0)
        120     continue
        if (dirac_integer_size == CC_DIRAC_INT4) then
            read (luint, end = 1301, err = 1302) ikr, jkr, nonzr, &
                    (indk(inz), indl(inz), inz = 1, nonzr), &
                    (cbuf(1, inz), cbuf(2, inz), inz = 1, nonzr)
        else  ! 8-byte integers in DIRAC
            read (luint, end = 1301, err = 1302) ikr8, jkr8, nonzr8, &
                    (indk8(inz), indl8(inz), inz = 1, nonzr8), &
                    (cbuf(1, inz), cbuf(2, inz), inz = 1, nonzr8)
            ikr = ikr8
            jkr = jkr8
            nonzr = nonzr8
            indk = indk8
            indl = indl8
        end if

        ! it was the last entry, break loop
        if (ikr == 0 .and. jkr == 0) goto 130
        ! if no non-zero elements in this batch => goto the next batch
        if (nonzr == 0) goto 120

        call expand_ints(ikr, jkr, nonzr, indk, indl, cbuf)

        ! read next batch
        goto 120
        130     continue
    end if
    ! NOTE: the case of the Levy-Leblond Hamiltonian is not yet implemented
    !       (NZ != 1 but values are real, to be tested)
    call c_asctime()
    close (unit = luint)
    print *, '*** END OF MDCINT FILE ***'
    print *
    call flush(6)

    ! final cleanup
    !deallocate(kr, indk, indl, cbuf, vint_buf, indices)
    ! if 8-byte integers in DIRAC
    if (dirac_integer_size == CC_DIRAC_INT8) then
        deallocate(kr8, indk8, indl8)
    end if

    err = 0
    return

    ! error handling
    1300 continue
    print *, 'error reading header MDCINT'
    err = 1
    return
    1301 continue
    print *, 'END OF FILE branch taken when reading 2-e integrals (from MDCINT)'
    err = 1
    return
    1302 continue
    print *, 'ERROR branch taken when reading 2-e integrals (from MDCINT)'
    err = 1
    return

contains


    ! expand integrals to the full list (Kramers and permutational symmetries
    ! are taken into account)
    subroutine expand_ints(ikr, jkr, nonzr, indk, indl, cbuf)
        use iso_c_binding
        ! arguments
        integer(4), intent(in) :: ikr, jkr
        integer(4), intent(in) :: nonzr
        integer(4), dimension(nspinors * nspinors), intent(in) :: indk, indl
        real(8), dimension(2, nspinors * nspinors), intent(in) :: cbuf
        ! local variables:
        ! integral class -- in terms of Kramers symmetry, number of barred indices
        integer(4) :: kclass
        ! Kramers pairs indices ("d" is for "Dirac notation")
        integer(4) :: kkr, lkr, ikrd, jkrd, kkrd, lkrd
        ! counter
        integer(4) :: inz
        ! tmp variable for the new integral
        complex(8) :: cint

        if (nonzr == 0) return

        ! determine integral class
        kclass = intclass(ikr, indk(1), jkr, indl(1))

        ! loop over array of non-zero integrals
        do inz = 1, nonzr
            kkr = indk(inz)
            lkr = indl(inz)

            ! pass to Dirac notation: (ij|kl) -> <ik|jl>
            ikrd = ikr
            jkrd = kkr
            kkrd = jkr
            lkrd = lkr

            ! construct integral
            cint = cbuf(1, inz) + cbuf(2, inz) * (0.0d0, 1.0d0)

            if (kclass == INTCLASS_NO_BARS .or. kclass == INTCLASS_TWO_BARS) then
                call perm_symm(cint, kr(ikrd), kr(jkrd), kr(kkrd), kr(lkrd))
                if (is_spinfree == 1) then
                    call perm_symm(cint, kr(ikrd), kr(-lkrd), kr(kkrd), kr(-jkrd))   ! MUST be in NR case
                    call perm_symm(cint, kr(-kkrd), kr(jkrd), kr(-ikrd), kr(lkrd))   ! MUST be in NR case
                end if
                call perm_symm(cint, kr(-kkrd), kr(-lkrd), kr(-ikrd), kr(-jkrd))

            else if (kclass == INTCLASS_ONE_BAR) then
                call perm_symm(cint, kr(ikrd), kr(jkrd), kr(kkrd), kr(lkrd))
                !call perm_symm(cint,kr(ikrd),kr(-lkrd),kr(kkrd),kr(-jkrd))
                !call perm_symm(-cint,kr(-kkrd),kr(jkrd),kr(-ikrd),kr(lkrd))
                call perm_symm(-cint, kr(-kkrd), kr(-lkrd), kr(-ikrd), kr(-jkrd))
            end if
        end do

    end subroutine expand_ints

    ! determine integral class
    integer(4) function intclass(ikr, jkr, kkr, lkr)
        integer(4) :: ikr, jkr, kkr, lkr   ! indices of Kramers pairs, > 0 for unbarred, < 0 for barred
        integer(4) :: sign_product

        sign_product = sign(1,ikr) * sign(1,jkr) * sign(1,kkr) * sign(1,lkr)
        if (sign_product > 0) then
            intclass = INTCLASS_NO_BARS
            ! INTCLASS_TWO_BARS is the same for Coulomb integrals
        else if (sign_product < 0) then
            intclass = INTCLASS_ONE_BAR
        else
            print *, 'unknown Kramers class'
            print *, ikr, jkr, kkr, lkr
            err = 1
            return
        end if

    end function intclass

    subroutine perm_symm(cint, i, j, k, l)
        complex(8), intent(in) :: cint
        integer(4), intent(in) :: i, j, k, l

        call put_integral(int(i, 2), int(j, 2), int(k, 2), int(l, 2), cint)
        call put_integral(int(j, 2), int(i, 2), int(l, 2), int(k, 2), cint)
        call put_integral(int(k, 2), int(l, 2), int(i, 2), int(j, 2), conjg(cint))
        call put_integral(int(l, 2), int(k, 2), int(j, 2), int(i, 2), conjg(cint))

    end subroutine perm_symm

end subroutine read_mdcint


!*******************************************************************************
! subroutine read_mdprop
!
! Reads properties integrals from the MDPROP file (if exists) and writes them
! to formatted files, format for matrix element <i|Op|j>:
! <i>   <j>   Re(<i|Op|j>)   Im(<i|Op|j>)
!
! Note that in DIRAC there are two types of operators:
! (1) "right order": real part is symmetric, imag is antisymmetric
! (2) "wrong order": real part is antisymm, imag is symm, and both are transposed
! the latter case is to be casted to the former one.
!*******************************************************************************
subroutine read_mdprop()

    use spinor_blocks
    implicit none

    character(len = 1024) :: oneel_file   ! passed via the /mrconee_mdcint/ common
    character(len = 1024) :: twoel_file
    character(len = 1024) :: prop_file
    character(len = 1024) :: gaunt_file
    integer(8) :: gaunt_enabled
    integer(8) :: n_twoprop
    character(len = 1024) :: twoprop_files(CC_MAX_NPROP)
    common /mrconee_mdcint/ oneel_file, twoel_file, prop_file, gaunt_file, gaunt_enabled, n_twoprop, twoprop_files
    bind(C) :: /mrconee_mdcint/

    logical :: mdprop_exists
    integer :: mdprop
    integer :: fprop   ! formatted file with integrals
    character(len = 32) :: achar
    real(8), dimension(:, :, :), allocatable :: prop
    integer :: i, j, normal_order

    print *
    print *, "*** MDPROP FILE ***"

    print *, 'path to MDPROP file = ', trim(prop_file)
    inquire(file = prop_file, exist = mdprop_exists)
    if (.not. mdprop_exists) then
        print *, "MDPROP file does not exist, properties integrals cannot be read"
        goto 3001
    end if

    print *, "MDPROP file exists"
    allocate(prop(2, n_spinors, n_spinors))

    mdprop = 10
    open(unit = mdprop, file = prop_file, form = "unformatted")

    1   continue
    read (mdprop, end = 11, err = 12) achar
    if (achar(25:32) .eq. 'EOFLABEL') goto 11
    print *, 'property = ', achar(25:32)
    read (mdprop, end = 12, err = 12) prop

    ! normal order: real part is symmetric, imag is antisymmetric
    ! else: real part is antisymm, imag is symm
    normal_order = 1
    do i = 1, n_spinors
        do j = 1, n_spinors
            ! real parts
            if (abs(prop(1, i, j) - prop(1, j, i)) > 1d-10) then
                normal_order = 0
            end if
            ! imag parts
            if (abs(prop(2, i, j) + prop(2, j, i)) > 1d-10) then
                normal_order = 0
            end if
        end do
    end do
    if (normal_order == 1) then
        print *, 'order of parts: Re + i*Im'
    else
        print *, 'order of parts: Im + i*Re'
    end if

    ! flush integrals to formatted file
    ! together with transposition
    fprop = 12
    open(unit = fprop, file = trim(achar(25:32)))
    do i = 1, n_spinors
        do j = 1, n_spinors
            if (normal_order == 1) then
                write (fprop, *) j, i, prop(1, i, j), prop(2, i, j)
            else
                write (fprop, *) j, i, prop(2, i, j), prop(1, i, j)
            end if
        end do
    end do
    close (unit = fprop)
    goto 1
    12  continue
    print *, 'error while reading MDPROP'
    11  continue
    print *, 'reached end of file MDPROP'

    close(mdprop)
    deallocate(prop)

    3001 continue
    print *, "*** END OF MDPROP FILE ***"

end subroutine read_mdprop


! ******************************************************************************
! subroutine expect_dipole
!
! Calculates DHF expectation value of the dipole moment operator for the given
! list of spinor occupation numbers (iocc parameter)
! NOTE: diagonal elements of the DM operator seems to be incorrect in DIRAC for
!       "infinite" symmetries
!       this fact, however, does not spoil DL-TDM values
! ******************************************************************************
subroutine expect_dipole(iocc)

    use general
    use spinor_blocks
    implicit none
    integer(4), dimension(CC_MAX_SPINORS) :: iocc
    complex(8), dimension(:, :), allocatable :: dx_mat, dy_mat, dz_mat
    complex(8) :: dx = 0.0d0, dy = 0.0d0, dz = 0.0d0
    real(8) :: d_abs
    integer :: i
    integer :: err = 0, no_files = 0

    allocate(dx_mat(n_spinors, n_spinors))
    allocate(dy_mat(n_spinors, n_spinors))
    allocate(dz_mat(n_spinors, n_spinors))

    call read_formatted_matrix('XDIPLEN', dx_mat, n_spinors, err)
    if (err /= 0) no_files = 1
    call read_formatted_matrix('YDIPLEN', dy_mat, n_spinors, err)
    if (err /= 0) no_files = 1
    call read_formatted_matrix('ZDIPLEN', dz_mat, n_spinors, err)
    if (err /= 0) no_files = 1
    if (no_files == 1) then
        print *
        print *, 'expectation value of dipole moment at the SCF level cannot be calculated'
        goto 10
    end if

    print *
    print *, 'ELECTRONIC contribution to dipole moment (DHF level)'

    do i = 1, n_spinors
        if (iocc(i) /= 1) cycle
        dx = dx + dx_mat(i, i)
        dy = dy + dy_mat(i, i)
        dz = dz + dz_mat(i, i)
    end do
    d_abs = sqrt(dx**2 + dy**2 + dz**2)

    print *, 'dx = ', real(dx), ' (re)', aimag(dx), '(im)'
    print *, 'dy = ', real(dy), ' (re)', aimag(dy), '(im)'
    print *, 'dz = ', real(dz), ' (re)', aimag(dz), '(im)'
    print *, '|d| = ', d_abs, ' a.u. = ', d_abs / 0.393430307, ' Debye'
    print *, '1 a.u = 2.54174623 Debye'

    ! cleanup and return
    10  continue
    deallocate(dx_mat)
    deallocate(dy_mat)
    deallocate(dz_mat)
    return

contains

    subroutine read_formatted_matrix(file, mat, n, err)
        implicit none
        character(len = *) :: file
        integer :: n
        complex(8), dimension(n, n) :: mat
        integer :: err
        integer :: fd, i, j
        real(8) :: val_re, val_im

        err = 0
        mat = (0.0d0, 0.0d0)

        fd = 20
        open(unit = fd, file = file, err = 13, status = 'old')
        31  continue
        read (fd, *, err = 30, end = 30) i, j, val_re, val_im
        mat(i, j) = cmplx(val_re, val_im, kind = 8)
        goto 31
        30  continue
        close(unit = fd)
        return

        13  continue
        err = 1
        return

    end subroutine read_formatted_matrix

end subroutine expect_dipole

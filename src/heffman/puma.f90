!***********************************************************************
!                                                                      *
!                             puma                                     *
!                                                                      *
!***********************************************************************
subroutine puma(mlar, msma, zlar, zsma, ii, jj, wtd)
    !
    !     puts square matrix in a cell (ii,jj) of large supermatrix
    !     or gets square matrix from a cell of large supermatrix
    !     all is complex
    !
    use nombres_mod
    implicit none
    integer :: i, j, ii, jj, mlar, msma, ncells, istart, jstart

    !      implicit real*8(a-h,o-x)
    !      implicit complex*16(z)
    !      implicit logical(y)
    character(len = 3) wtd
    complex(kind = pd) :: zlar(mlar, *), zsma(msma, msma)
    ncells = mlar / msma
    if(mlar /= ncells * msma) then
        write(*, '(1x,a)')' subroutine puma:'
        write(*, '(1x,a,2i8)')' wrong supermatrix/matrix sizes', mlar, msma
        stop
    endif
    if(wtd /= 'get' .and. wtd /= 'put') then
        write(*, '(1x,a)')' subroutine puma:'
        write(*, '(1x,a,a)')' do not know what to do:', wtd
        stop
    endif
    istart = (ii - 1) * msma
    jstart = (jj - 1) * msma
    do i = 1, msma
        do j = 1, msma
            if(wtd == 'put') then
                zlar(istart + i, jstart + j) = zsma(i, j)
            else
                zsma(i, j) = zlar(istart + i, jstart + j)
            endif
        enddo
    enddo
    return
end subroutine puma
!***********************************************************************
!                                                                      *
!                             matpade                                  *
!                                                                      *
!***********************************************************************
subroutine matpade(max, ll, mm, zc, zpade, impression)
    !
    !  [ll/mm] (right) rational pade approximant
    !          for a complex (max x max) matrix series zc
    !
    use nombres_mod
    implicit none
    integer :: max, ll, mm, impression, min_mm_i, ii, iii, i, j, k, maxmax, mlar
    complex(kind = pd) :: zc(max * max, *), zpade(max, max), zs
    complex(kind = pd), allocatable :: z(:, :), zb(:, :), za(:, :), zright(:, :), zbt(:, :)
    if(impression > 0) then
        write(6, '(/4x,a,i1,a,i1,a/)') '----- [', ll, '/', mm, '] matrix rational pade approximant -----'
    endif
    !
    !   input data
    !
    if(impression >= 10) then
        do i = 1, ll + mm + 1
            write(6, '(1x,a,i2)') ' series: term', i - 1
            call imprec(zc(1, i), max, max, 5)
        enddo
    endif
    !
    if(mm <= 0) then
        !
        !   trivial [ll/0] pade approximant = taylor partial sum
        !
        ii = 0
        do j = 1, max
            do i = 1, max
                ii = ii + 1
                zpade(i, j) = (0.d0, 0.d0)
                do k = 1, ll + mm + 1
                    zpade(i, j) = zpade(i, j) + zc(ii, k)
                enddo
            enddo
        enddo
        return
        !
    endif
    !
    !
    !   allo
    !
    allocate (z(max * mm, max * mm), zbt(max * mm, max), zright(max * mm, max))
    !
    !   clean z
    !
    z(1:max * mm, 1:max * mm) = (0.d0, 0d0)
    !
    !   fill
    !
    mlar = max * mm
    do i = 1, mm
        call puma(mlar, max, zright, zc(1, ll + i + 1), i, 1, 'put')
        do j = 1, mm
            if(ll + i - j + 1 >= 1 .and. ll + i - j + 1 <= ll + mm + 1) then
                call puma(mlar, max, z, zc(1, ll + i - j + 1), i, j, 'put')
            endif
        enddo
    enddo
    !
    !   rhs sign
    !
    zright(1:mlar, 1:max) = -zright(1:mlar, 1:max)

    if(impression.gt.15)then
        write(*, '(1x,a)')' leq system matrix'
        call imprec(z, mlar, mlar, 5)
        write(*, '(1x,a)')' right hand side'
        call imprec(zright, mlar, max, 5)
    endif
    !
    !   invert z
    !
    call zinver(z, mlar)
    if(impression > 15)then
        write(*, '(1x,a)')' inverted leq system matrix'
        call imprec(z, mlar, mlar, 5)
    endif
    !
    !   compute denominator matrix coefs
    !
    do i = 1, mlar
        do j = 1, max
            zs = (0.d0, 0.d0)
            do k = 1, mlar
                zs = zs + z(i, k) * zright(k, j)
            enddo
            zbt(i, j) = zs
        enddo
    enddo
    !
    !    plus besoin de z,zright
    !
    deallocate (z, zright)
    allocate (za(max * max, ll + 1))
    allocate (zb(max * max, mm))
    !
    !   extract denominator matrix coefs to zb
    !
    do i = 1, mm
        call puma(mlar, max, zbt, zb(1, i), i, 1, 'get')
        if(impression > 10) then
            write(*, '(1x,a,i1)')' b', i
            call imprec(zb(1, i), max, max, 5)
        endif
    enddo
    !
    deallocate (zbt)
    !
    !   compute numerator and put it into za
    !
    maxmax = max * max
    do i = 0, ll
        !
        call zcopy(maxmax, zc(1, i + 1), 1, za(1, i + 1), 1)
        if(i > 0) then
            min_mm_i = i
            if(mm < min_mm_i) min_mm_i = mm
            if(min_mm_i .ge. 1) then
                do k = 1, min_mm_i
                    call summat(za(1, i + 1), zc(1, i + 1 - k), zb(1, k), max)
                enddo
            endif
        endif
        if(impression > 10) then
            write(*, '(1x,a,i1)')' a', i
            call imprec(za(1, i + 1), max, max, 5)
        endif
    enddo
    !
    !     compute total numerator for z=1 in za( ,1)
    !
    if(ll > 0) then
        do i = 2, ll + 1
            do j = 1, maxmax
                za(j, 1) = za(j, 1) + za(j, i)
            enddo
        enddo
    endif
    !
    !     compute total denominator for z=1 in zb( ,1)
    !
    if (mm > 1) then
        do i = 2, mm
            do j = 1, maxmax
                zb(j, 1) = zb(j, 1) + zb(j, i)
            enddo
        enddo
    endif
    !    add unit
    iii = 0
    do i = 1, max
        do j = 1, max
            iii = iii + 1
            if(i == j) zb(iii, 1) = zb(iii, 1) + (1.d0, 0.d0)
        enddo
    enddo
    if(impression > 10) then
        write(*, '(1x,a)')' total numerator'
        call imprec(za(1, 1), max, max, 5)
        write(*, '(1x,a)')' total denominator'
        call imprec(zb(1, 1), max, max, 5)
    endif
    !  invert
    call zinver(zb, max)
    !
    !     compute [ll/mm] value at z=0 and put it to zpade
    !
    zpade = (0.d0, 0.d0)
    call summat(zpade, za, zb, max)
    if(impression > 10)then
        write(*, '(1x,a)')' approximant at z=1'
        call imprec(zpade, max, max, 5)
    endif
    !
    !   deallouer
    !
    deallocate (zb, za)
    !
    return
end subroutine matpade
!***********************************************************************
!                                                                      *
!                             summat                                   *
!                                                                      *
!***********************************************************************
subroutine summat(za, zb, zc, max)
    !
    use nombres_mod
    implicit real*8(a-h, o-x)
    integer, intent(in) :: max
    integer :: i, j, k
    complex(kind = pd), dimension(max, *), intent(in) :: zb, zc
    complex(kind = pd), dimension(max, *), intent(inout) :: za
    !
    !
    !     a=a+b*c
    !
    do i = 1, max
        do j = 1, max
            do k = 1, max
                za(i, j) = za(i, j) + zb(i, k) * zc(k, j)
            enddo
        enddo
    enddo
    return
end subroutine summat
            

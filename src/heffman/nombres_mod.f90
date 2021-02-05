module nombres_mod
    implicit none

    public

    integer, parameter :: ps = selected_real_kind(6, 30)         ! precision simple
    integer, parameter :: pd = selected_real_kind(15, 305)       ! precision double
    real(kind = pd), parameter :: pi = 3.141592653589793238462643d0
    real(kind = pd), parameter :: pi234 = 0.356352735177495080176699228d0
    !     pi234 = (2*(pi**3))**(-1/4)
    real(kind = pd), parameter :: zero = 0.d0, two = 2.d0, half = 0.5d0, twothirds = 0.666666666666667d0
    real(kind = pd), parameter :: a_dau_eau = 21420.007784d0
    ! einsten coeff conversion from au        (luuk visscher's value)
    real(kind = pd), parameter :: au2cm = 219474.631280634d0

end module nombres_mod


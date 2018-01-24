! Code to carry out variational quantum monte carlo calculations
! for simple systems, we work in units of Angstroms/Electronvolts

! The constants used by the program
module constants
implicit none
    ! The floating point precision we use
    integer,    parameter :: prec = selected_real_kind(15,307)
    real(prec), parameter :: bohrRadius = 0.529177
end module

module basisSets
use constants
implicit none

contains

    ! A Hydrogen-like atomic orbital with the given
    ! qunatum numbers
    function atomicOrbital(n,l,m,x) result(ret)
        real(prec) :: x(3), ret
        integer    :: n, l, m

    end function

    ! The spherical harmonics
    function spericalHarmonic(l,m,x) result(ret)
        real(prec) :: x(3), ret
        integer    :: l, m
        
    end function

    ! The associated Legendre polynomials
    function associatedLegendrePolynomial(l,m,x) result(ret)
        real(prec) :: x, ret
        integer    :: l, m

    end function

    ! The legendre polynomials (defined recursively)
    recursive function legendrePolynomial(n,x) result(ret)
        real(prec) :: x, ret
        integer    :: n
        if (n == 0) then
            ret = 1
            return
        else if (n == 1) then
            ret = x
            return
        endif
        ret = legendrePolynomial(n-1,x)*x*(2*n-1)/n - &
              legendrePolynomial(n-1,x)*(n-1)/n
    end function

end module

program main
use basisSets
implicit none
    integer :: i
    do i=0,10
        print *, legendrePolynomial(i,real(0, kind=prec))
    enddo
end program
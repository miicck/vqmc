! Atomic basis functions
module atomicBasis
use constants
use vqmc, only: basisState
implicit none

    type, extends(basisState) :: atomicState
    private
        integer :: n, l, m, z
    contains
        procedure :: value
        procedure :: printDebugInfo
    end type

    interface atomicState
        module procedure constructor
    end interface

contains

    ! Implementation of atomicBasis value
    function value(this, x)
    implicit none
    class(atomicState) :: this
    complex(prec) :: value
    real(prec)    :: x(3)
        value = atomicWavefunction(this%n,this%l,this%m,this%z,x)
    end function

    ! Implementation of the debug print subroutine for atomic states
    subroutine printDebugInfo(this)
    implicit none
    class(atomicState) :: this
        print *, this%n , this%l, this%m, this%z
    end subroutine

    ! Atomic state constructor
    function constructor(n, l, m, z) result(this)
    implicit none
        type(atomicState), pointer :: this
        integer :: n, l, m, z
        allocate(this)

        if (n<1) print *, "ERROR: n < 1 in atomic basis state!"
        if (l>=n) print *, "ERROR: l > n-1 in atomic basis state!"
        if (abs(m)>l) print *, "ERROR: abs(m) > l in atomic basis state!"

        this%n = n
        this%l = l
        this%m = m
        this%z = z
    end function

    ! Returns L_n^a(x) via recursion relation
    recursive function laguerrePolynomial(n, a, x) result(ret)
    implicit none  
    integer    :: n, a, k
    real(prec) :: x, ret
        !print *, n,a,x
        if (n == 0) then
            ret = 1
        else if (n == 1) then
            ret = 1 + a - x
        else
            k = n - 1
            ret = ((2*k+1+a-x) * laguerrePolynomial(k,a,x) - &
                   (k+a) * laguerrePolynomial(k-1,a,x)) &
                   /(k+1)
        endif
    end function

    ! Returns P_l^m(x) via recursion relation
    recursive function associatedLegendrePolynomial(l, m, x) result(ret)
    implicit none
    integer    :: l, m
    real(prec) :: x, ret
        if (l==0 .and. m==0) then
            ret = 1
        else if (l == m) then
            ret = -(2*(l-1)+1)*sqrt(1-x**2)*associatedLegendrePolynomial(l-1,l-1,x)
        else if (l-1 > m) then
            ret = (2*(l-1)+1)*x*associatedLegendrePolynomial(l-1,m,x) - &
                (l-1+m)*associatedLegendrePolynomial(l-2,m,x)
            ret = ret / (l-m)
        else if (m == l-1) then
            ret = x*(2*(l-1)+1)*associatedLegendrePolynomial(l-1,l-1,x)
        else
            print *, "Error, recursion for associatedLegendrePolynomial broke..."
            ret = 0
        endif
    end function

    ! Returns the factorial of n
    recursive function factorial(n) result(fact)
    implicit none        
    integer :: n, fact
        if (n == 0) then
            fact = 1
        else
            fact = n*factorial(n-1)
        endif
    end function

    ! Get the radial part of the n, l atomic orbital
    ! of an atom with charge z
    function radialPart(n, l, r, z)
    implicit none        
    integer    :: n, l, z
    real(prec) :: r, radialPart
        radialPart = laguerrePolynomial(n-1-l,2*l+1,2*z*r/(n*a0))
        radialPart = radialPart * (2*z*r/(n*a0))**l
        radialPart = radialPart * exp(-z*r/(a0*n))
    end function

    ! Y_l^m
    function sphericalHarmonic(l, m, theeta, phi) result(ret)
    implicit none
    integer    :: l, m
    real(prec) :: theeta, phi
    complex(prec)    :: ret
        ret = exp((0,1)*m*phi)
        ret = ret * associatedLegendrePolynomial(l,m,cos(theeta))      
    end function

    ! The n, l, m atomic wavefunction with nuclear charge z
    function atomicWavefunction(n, l, m, z, x) result(ret)
    implicit none
    integer    :: n, l, m, z
    real(prec) :: x(3), theeta, phi
    complex(prec)    :: ret
        phi = atan(x(2)/x(1))
        theeta = acos(x(3)/norm2(x))
        ret = sphericalHarmonic(l, m, theeta, phi)
        ret = ret*radialPart(n, l, norm2(x), z)        
    end function

    ! Log wavefunctions for debugging
    subroutine debugAtomicState(state)
    implicit none
    class(atomicState) :: state
    integer    :: xi, yi
    integer, parameter :: grid = 100
    real(prec) :: r(3), a
        print *, "Debugging atomic state: "
        call state%printDebugInfo()
        open(unit=1,file="wavefunctionDebug")
        do xi=-grid,grid
            do yi=-grid,grid
                r(1) = 5*state%n*angstrom*xi/real(grid)
                r(3) = 5*state%n*angstrom*yi/real(grid)
                a = abs(state%value(r))**2
                if (.not. isnan(a)) then
                    write(1,*) r(1)/angstrom,",",r(3)/angstrom,",",a
                endif
            enddo
        enddo
        close(unit=1)
    end subroutine

    ! Log the radial parts for debugging
    subroutine debugRadialPart(unit)
    implicit none
    integer    :: unit, n, l, i
    real(prec) :: r
        do n=1,3
            do l=0,n-1
                write(2,*) "#"
                do i=1,1000
                    r = real(i)/1000
                    r = r * 15
                    write(2,*) r,",",radialPart(n,l,r*angstrom,1)
                enddo
            enddo
        enddo
    end subroutine

    ! Log the associatedLegendrePolynomials for debugging
    subroutine debugAssociatedLegendrePolynomial(unit)
    implicit none
    integer    :: unit, l, m, i
    real(prec) :: x
        do l=0,4
            do m=0,l
                write(2,*) "#"
                do i=1,1000
                    x = real(i)/1000
                    x = 2*(x-0.5)
                    write(2,*) x,",",associatedLegendrePolynomial(l,m,x)
                enddo
            enddo
        enddo
    end subroutine

end module

! ----- OLD SYSTEMS ----- !

! A system with a nucleus of the given charge
! at the origin
module hydrogenicSystem
use constants
implicit none

    real(prec) :: nuclearCharge = 1

contains

    ! hydrogen ground state
    ! Wavefunction (not normalized)
    function groundState(x)
    implicit none
        real(prec)    :: x(3)
        complex(prec) :: groundState
        groundState = exp(-nuclearCharge*norm2(x)/a0)
    end function

    ! Hydrogen potential
    function potential(x)
    implicit none
        real(prec) :: x(3), potential
        potential = -nuclearCharge*qElectron**2/(4*pi*epsNaught*norm2(x))
    end function
    
end module

! Things to do with the particle in a 3D box
module particleInBox
use constants
implicit none

    real(prec), parameter :: boxSize = angstrom
    real(prec), parameter :: infiniteEnergy = huge(prec)

contains

    ! The particle in a box ground state
    function groundState(x)
    implicit none
        real(prec)    :: x(3)
        complex(prec) :: groundState
        if (outsideBounds(x(1))) then
            groundState = 0
        else if (outsideBounds(x(2))) then
            groundState = 0
        else if (outsideBounds(x(3))) then
            groundState = 0
        else
            groundState = sin(pi*x(1)/boxSize)*sin(pi*x(2)/boxSize)*sin(pi*x(3)/boxSize);
        endif
    end function

    ! The particle in a box potential
    function potential(x)
    implicit none
        real(prec) :: x(3), potential
        if (outsideBounds(x(1))) then
            potential = infiniteEnergy
        else if (outsideBounds(x(2))) then
            potential = infiniteEnergy
        else if (outsideBounds(x(3))) then
            potential = infiniteEnergy
        else
            potential = 0
        endif
    end function

    ! Returns true if the single coordinate x is outside the box
    function outsideBounds(x)
    implicit none
        real(prec) :: x
        logical :: outsideBounds
        outsideBounds = or(x<0,x>boxSize) 
    end function

end module
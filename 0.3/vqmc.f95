! Version 0.2
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
implicit none

    ! The number of monte-carlo itterations
    integer, parameter :: MC_ITTER = 1000000
    integer, parameter :: MC_INIT_STEPS = 1000

    ! The minimum distance from the divergent origin
    ! our trial electrons are allowed
    real(prec), parameter :: minR = angstrom/10000

    ! An object that can be used as a basis function
    type, abstract :: basisState
    contains
        procedure(basisWfn), deferred :: value
    end type

    abstract interface
        ! A wavefunction
        function wvfn(x)
            import
            real(prec) :: x(3)
            complex(prec) :: wvfn
        end function

        ! Interface for basis function
        function basisWfn(this, x)
            import
            class(basisState) :: this
            real(prec)    :: x(3)
            complex(prec) :: basisWfn
        end function
    end interface

contains

    ! Optimize the coefficeints in the given basis to
    ! minimize the energy in the given potential
    subroutine optimizeBasis(basis, potential)
    implicit none
        class(basisState)     :: basis(:)
        procedure(real(prec)) :: potential
    end subroutine

    ! Carry out a monte carlo integration of the local energy of the
    ! wavefunction in the given potential, using importance sampling
    ! of the wavefunction
    function energy(wavefunction, potential) result(ret)
    implicit none
        procedure(real(prec)) :: potential
        procedure(wvfn)       :: wavefunction
        complex(prec)         :: ret
        real(prec) :: x(3)
        integer    :: i

        x = 0
        ret = 0

        ! Get rid of any dependence on initial conditions
        ! by carring out MC_INIT_STEPS trial moves and
        ! ignoring them
        do i=1,MC_INIT_STEPS
            call metro(x,wavefunction)
        enddo

        ! Actually sample
        do i=1, MC_ITTER
            call metro(x, wavefunction)
            ret = ret + localEnergy(wavefunction, potential, x)
        enddo
        ret = ret/MC_ITTER

    end function

    ! Calculate the second derivative of the wavefunction
    ! in the given direction using finite differences
    function secondDerivative(wavefunction, x, dir)
        procedure(wvfn) :: wavefunction
        complex(prec)   :: secondDerivative
        real(prec)      :: x(3), dir(3)
        real(prec), parameter :: eps = 0.00001 * angstrom
        dir = eps * dir/norm2(dir)
        secondDerivative = wavefunction(x+dir) - 2*wavefunction(x) + wavefunction(x-dir)
        secondDerivative = secondDerivative / (eps**2)
    end function

    ! Calculate psi^-1(x) T psi(x)
    function localKineticEnergy(wavefunction, x)
        procedure(wvfn) :: wavefunction
        real(prec)      :: x(3)
        complex(prec)   :: dx, dy, dz, localKineticEnergy
        dx = secondDerivative(wavefunction,x,real((/1,0,0/),kind=prec))
        dy = secondDerivative(wavefunction,x,real((/0,1,0/),kind=prec))
        dz = secondDerivative(wavefunction,x,real((/0,0,1/),kind=prec))
        localKineticEnergy = -(hbar**2/(2*mElectron))*(dx + dy + dz)
        localKineticEnergy = localKineticEnergy/wavefunction(x)
    end function

    ! Calculate the local energy of the given wavefunction
    ! with the given potential at x
    function localEnergy(wavefunction, potential, x)
    implicit none
        procedure(wvfn) :: wavefunction
        procedure(real(prec)) :: potential
        complex(prec)         :: localEnergy        
        real(prec) :: x(3)
        localEnergy = localKineticEnergy(wavefunction,x) + potential(x)
    end function

    ! Make a trial metropolis move on x
    subroutine metro(x, wavefunction)
    implicit none
        real(prec) :: x(3), newX(3) , dx(3), r, theeta, phi
        procedure (wvfn) :: wavefunction

        ! Move some distance r in any direction
        theeta = rand()*pi
        phi    = rand()*2*pi
        r      = rand()*angstrom/10

        dx(1) = r*sin(theeta)*cos(phi)
        dx(2) = r*sin(theeta)*sin(phi)
        dx(3) = r*cos(theeta)

        newX = x + dx
        if (rand() < abs(wavefunction(newX))**2 / & 
                     abs(wavefunction(x))**2) then
            x = newX ! Accept the move via metropolis criteria
        endif

        ! Ensure x doesn't get too close to a divergent origin
        if (norm2(x) < minR) then
            x(1) = minR
        endif
    end subroutine

end module vqmc
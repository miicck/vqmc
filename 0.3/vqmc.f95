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

    ! Optimize the coefficients in the given basis to
    ! minimize the energy in the given potential
    subroutine optimizeBasis(basis, potential)
    implicit none
        class(basisState)          :: basis(:)
        procedure(real(prec))      :: potential
        complex(prec), allocatable :: coefficients(:)
        complex(prec)              :: lowestEnergy, currentEnergy
        real(prec)                 :: inf
        integer                    :: i, c

        ! Initialize our best energy to infinity
        ! and allocate space for our coefficients
        lowestEnergy = huge(inf)
        allocate(coefficients(size(basis)))

        ! Run optimization cycles
        c = 0
        do while(.true.)
            c = c + 1
            
            call randomizeCoefficients(coefficients)

            ! Calculate energy with new coefficients
            print *, ""
            print *, ""
            print *, "Optimization cycle: ", c
            currentEnergy = energy(basis, coefficients, potential)

            if (realpart(currentEnergy) < realpart(lowestEnergy)) then

                ! Improved ground state estimate
                lowestEnergy = currentEnergy
                print *, "Decreased energy to: ", lowestEnergy/electronVolt
                do i=1,size(basis)
                    print *, "    Coefficient: ", i, ":", coefficients(i)
                enddo

            else

                ! Didn't improve ground state estimate
                print *, "Failed to decrease energy."

            endif

        enddo

    end subroutine

    subroutine randomizeCoefficients(coefficients)
    implicit none
        complex(prec) :: coefficients(:)
        integer       :: i
        do i=1,size(coefficients)
            coefficients(i) = real(2*(rand()-0.5),kind=prec)
        enddo
    end subroutine

    ! A wavefucntion in a given basis, with given
    ! coefficients evaluated at x
    function wavefunction(basis, coefficients, x) result(ret)
    implicit none
        class(basisState) :: basis(:)
        complex(prec)     :: ret, coefficients(:)
        real(prec)        :: x(3)
        integer           :: i
        if (size(basis) /= size(coefficients)) then
            print *, "ERROR: Basis/Coefficients mismatch!"
        endif
        ret = 0
        do i=1,size(basis)
            ret = ret + coefficients(i)*basis(i)%value(x)
        enddo
    end function

    ! Calculate the second derivative of the wavefunction
    ! in the given direction using finite differences
    function secondDerivative(basis, coefficients, x, dir)
        class(basisState) :: basis(:)
        complex(prec)     :: secondDerivative, coefficients(:)
        real(prec)        :: x(3), dir(3)
        real(prec), parameter :: eps = 0.00001 * angstrom
        dir = eps * dir/norm2(dir)
        secondDerivative = wavefunction(basis, coefficients, x+dir) - &
                           2*wavefunction(basis, coefficients, x) + &
                           wavefunction(basis, coefficients, x-dir)
        secondDerivative = secondDerivative / (eps**2)
    end function

    ! Calculate psi^-1(x) T psi(x)
    function localKineticEnergy(basis, coefficients, x)
        class(basisState) :: basis(:)
        real(prec)        :: x(3)
        complex(prec)     :: dx, dy, dz, localKineticEnergy, coefficients(:)
        dx = secondDerivative(basis, coefficients, x, real((/1,0,0/), kind=prec))
        dy = secondDerivative(basis, coefficients, x, real((/0,1,0/), kind=prec))
        dz = secondDerivative(basis, coefficients, x, real((/0,0,1/), kind=prec))
        localKineticEnergy = -(hbar**2/(2*mElectron))*(dx + dy + dz)
        localKineticEnergy = localKineticEnergy/wavefunction(basis, coefficients, x)
    end function

    ! Calculate the local energy of the given wavefunction
    ! with the given potential at x
    function localEnergy(basis, coefficients, potential, x)
    implicit none
        class(basisState)     :: basis(:)
        procedure(real(prec)) :: potential
        complex(prec)         :: localEnergy, coefficients(:)      
        real(prec) :: x(3)
        localEnergy = localKineticEnergy(basis, coefficients, x) + potential(x)
    end function

    ! Carry out a monte carlo integration of the local energy of the
    ! wavefunction in the given potential, using importance sampling
    ! of the wavefunction
    function energy(basis, coefficients, potential) result(ret)
    implicit none
        procedure(real(prec)) :: potential
        class(basisState)     :: basis(:)
        complex(prec)         :: ret, coefficients(:)
        real(prec) :: x(3)
        integer    :: i

        x = 0
        ret = 0

        ! Get rid of any dependence on initial conditions
        ! by carring out MC_INIT_STEPS trial moves and
        ! ignoring them
        do i=1,MC_INIT_STEPS
            call metro(x, basis, coefficients)
        enddo

        ! Actually sample
        do i=1, MC_ITTER
            call metro(x, basis, coefficients)
            ret = ret + localEnergy(basis, coefficients, potential, x)
        enddo
        ret = ret/MC_ITTER

    end function

    ! Make a trial metropolis move on x
    subroutine metro(x, basis, coefficients)
    implicit none
        real(prec) :: x(3), newX(3) , dx(3), r, theeta, phi
        class(basisState) :: basis(:)
        complex(prec)     :: coefficients(:)

        ! Move some distance r in any direction
        theeta = rand()*pi
        phi    = rand()*2*pi
        r      = rand()*angstrom/10

        dx(1) = r*sin(theeta)*cos(phi)
        dx(2) = r*sin(theeta)*sin(phi)
        dx(3) = r*cos(theeta)

        newX = x + dx
        if (rand() < abs(wavefunction(basis,coefficients,newX))**2 / & 
                     abs(wavefunction(basis,coefficients,x))**2) then
            x = newX ! Accept the move via metropolis criteria
        endif

        ! Ensure x doesn't get too close to a divergent origin
        if (norm2(x) < minR) then
            x(1) = minR
        endif
    end subroutine

end module vqmc
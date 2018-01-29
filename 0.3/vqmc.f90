! Version 0.2
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
!use IFPORT
implicit none

    ! The minimum distance from the divergent origin
    ! our trial electrons are allowed
    real(prec), parameter :: minR = angstrom/10000

    ! An object that can be used as a basis function
    type, abstract :: basisState
    contains
        procedure(basisWfn), deferred :: value
        procedure(printStateDebug), deferred :: printDebugInfo
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

        ! Interface for print statements
        subroutine printStateDebug(this)
            import
            class(basisState) :: this
        end subroutine

    end interface

    interface
        ! A potential
        function pot(x)
            import
            real(prec) :: x(3), pot
        end function
    end interface

contains

    ! Optimize the coefficients in the given basis to
    ! minimize the energy in the given potential
    subroutine optimizeBasis(basis, potential)
    implicit none
        class(basisState)          :: basis(:)
        procedure(pot)             :: potential
        complex(prec), allocatable :: coefficients(:), newCoefficients(:)
        real(prec),    allocatable :: gradReal(:), gradImag(:)
        real(prec),    parameter   :: eps = 0.0001
        real(prec)                 :: samples(3,10000), shunt, itterVar, newVar, workVar
        complex(prec)              :: itterEnergy, workE
        integer                    :: i, n, m

        ! Allocate space
        allocate(coefficients(size(basis)))
        allocate(newCoefficients(size(basis)))
        allocate(gradReal(size(basis)))
        allocate(gradImag(size(basis)))

        do i=1,size(basis)
            coefficients(i) = (1,1)
        enddo

        shunt = 0.1

        do n=1,100

            ! Generate samples in coordinate space using
            ! the metropolis algorithm for the current
            ! trial wavefunction
            call metro(samples, basis, coefficients)        
                
            ! Calculate the energetics with the current coefficients
            call energetics(basis, coefficients, potential, samples, itterEnergy, itterVar)
            
            ! Calculate the real and imaginary coefficient
            ! gradients of the real part of the energy variance
            do i=1,size(basis)

                ! Real gradients
                coefficients(i) = coefficients(i) + (eps, 0)
                call energetics(basis, coefficients, potential, samples, workE, workVar)
                gradReal(i) = (workVar-itterVar) / eps
                coefficients(i) = coefficients(i) - (eps, 0)

                ! Imaginary gradients
                coefficients(i) = coefficients(i) + (0, eps)
                call energetics(basis, coefficients, potential, samples, workE, workVar)
                gradImag(i) = (workVar-itterVar) / eps
                coefficients(i) = coefficients(i) - (0, eps)

            enddo

            ! Shunt the coefficients closer to the minimum
            ! by about the size of the minimum coefficient
            call makeOrderUnityReal(gradReal)
            call makeOrderUnityReal(gradImag)
            shunt = minimumModulus(coefficients)*0.5
            gradReal = gradReal * shunt
            gradImag = gradImag * shunt

            ! Move along steepest decent
            newCoefficients = coefficients - gradReal - gradImag*(0,1)
            call makeOrderUnityComplex(newCoefficients)
            call energetics(basis, newCoefficients, potential, samples, workE, newVar)

            ! Accept the new coefficients if they reduced the variance
            if (newVar < itterVar) coefficients = newCoefficients

            print *, ""
            print *, ""
            print *, ""
            print *, "Itteration:", n
            print *, "    Energy:  ", itterEnergy/electronVolt, "eV"
            print *, "    Variance:", min(newVar, itterVar)/(electronVolt**2)
            do m=1, size(basis)
                print *, "    Basis", m, "character", abs(coefficients(m))
            enddo

        enddo

    end subroutine

    ! Get the minimum modulus of a set of complex numbers
    function minimumModulus(set) result(min)
    implicit none
        complex(prec) :: set(:)
        real(prec) :: min
        integer    :: i
        min = huge(min)
        do i=1,size(set)
            if (abs(set(i))<min) min = abs(set(i))
        enddo
    end function

    ! Make the values in a set O(1)
    subroutine makeOrderUnityReal(set)
    implicit none
        real(prec) :: set(:), max
        integer    :: i
        max = 0
        do i=1,size(set)
            if (abs(set(i))>max) max = abs(set(i))
        enddo
        set = set / max
    end subroutine

    ! Make the values in a set O(1)
    subroutine makeOrderUnityComplex(set)
    implicit none
        complex(prec) :: set(:)
        real(prec)    :: max
        integer       :: i
        max = 0
        do i=1,size(set)
            if (abs(set(i))>max) max = abs(set(i))
        enddo
        set = set / max
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
        class(basisState) :: basis(:)
        procedure(pot)    :: potential
        complex(prec)     :: localEnergy, coefficients(:)      
        real(prec) :: x(3)
        localEnergy = localKineticEnergy(basis, coefficients, x) + potential(x)
    end function

    ! Carry out a monte carlo integration of the local energy of the
    ! wavefunction in the given potential, using the sample positions
    ! provided
    subroutine energetics(basis, coefficients, potential, samples, energy, variance)
    implicit none
        procedure(pot)             :: potential
        class(basisState)          :: basis(:)
        complex(prec)              :: energy, coefficients(:)
        complex(prec), allocatable :: localEnergies(:)
        real(prec)                 :: samples(:,:), variance
        integer                    :: i

        allocate(localEnergies(size(samples,2)))

        ! Calculate the average local energy of our samples
        energy = 0
        do i=1, size(localEnergies)
            localEnergies(i) = localEnergy(basis, coefficients, potential, samples(:,i))
            energy = energy + localEnergies(i)
        enddo
        energy = energy/size(localEnergies)

        ! Calculate the local energy variance of our samples
        variance = 0
        do i=1, size(localEnergies)
            variance = variance + realpart(localEnergies(i)-energy)**2
        enddo
        variance = variance/size(localEnergies)

    end subroutine

    ! Get a set of samples using the metropolis algorithm
    subroutine metro(samples, basis, coefficients)
    implicit none
        real(prec) :: newX(3), r, theeta, phi
        real(prec), allocatable :: x(:), dx(:)
        class(basisState) :: basis(:)
        complex(prec)     :: coefficients(:)
        real(prec) :: samples(:,:)
        integer    :: i

        ! Number of metropolis steps to take and discard
        ! to remove dependance on initial position
        integer, parameter :: INIT_STEPS = 1000

        ! Initial position is 0
        allocate(x(size(samples,1)))
        allocate(dx(size(samples,1)))
        x = 0

        do i=1,size(samples,2)+INIT_STEPS

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

            if (i>INIT_STEPS) samples(:,i-INIT_STEPS) = x
        enddo
    end subroutine

end module vqmc
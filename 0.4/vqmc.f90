! Version 0.3
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
!use IFPORT
implicit none

    ! Parameters controlling our vqmc simulation
    
    integer,    parameter :: metroSamples = 10000      ! The number of metropolis electron configurations that are used
    real(prec), parameter :: minR = angstrom/10000     ! The minimum distance e's are alowed from nucleii
    real(prec), parameter :: maxMetroJump = 4*angstrom ! The maximum distance an electron can move in a metropolis trial move

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
    class(basisState) :: basis(:)
    procedure(pot)    :: potential
    real(prec)        :: randCoeffChange
    complex(prec), allocatable :: coefficients(:)

        allocate(coefficients(size(basis)))

        call optimizeBasisRandomSearch(basis, coefficients, potential, 50, randCoeffChange)
        !call optimizeBasisSteepestDecent(basis, coefficients, potential, 100, randCoeffChange)
        call sampleWavefunctionToFile(basis, coefficients)
        call writeWavefunctionToFile(basis, coefficients)

    end subroutine

    ! Optimize the coefficients of our basis using a random search
    ! of coefficient space. lastCoeffChange will contain the amount
    ! that |coefficients| changed at the last step in energy (useful
    ! to know for further optimization)
    subroutine optimizeBasisRandomSearch(basis, coefficients, potential, itterations, lastCoeffChange)
    implicit none
    class(basisState)          :: basis(:)
    procedure(pot)             :: potential
    complex(prec), allocatable :: newCoefficients(:)
    real(prec)                 :: samples(3,metroSamples), workVar, lastCoeffChange
    complex(prec)              :: workE, minEnergy, coefficients(:)
    integer                    :: i, n, itterations
    real(prec), parameter      :: relaxation = 0

        ! Allocate space
        allocate(newCoefficients(size(basis)))
        do i=1,size(basis)
            newCoefficients(i) = 1
        enddo

        open(unit=1,file="output/randomOptimizationEnergy")
        open(unit=2,file="output/randomOptimizationBasis1Char")
        open(unit=3,file="output/randomOptimizationBasis2Char")

        ! Initialize coefficients to the lowest energy
        ! of a random set of coefficients
        print *, "Initializing basis..."
        minEnergy = huge(workVar)
        do n=1,itterations
            do i=1,size(basis)
                coefficients(i) = rand()*(1-relaxation) + newCoefficients(i)*relaxation
            enddo
            call normalizeCoefficients(coefficients)
            call metro(samples, basis, coefficients)  
            call energetics(basis, coefficients, potential, samples, workE, workVar)
            if (realpart(workE)<realpart(minEnergy)) then
                minEnergy = realpart(workE)
                lastCoeffChange = norm2(abs(newCoefficients-coefficients))
                newCoefficients = coefficients
            endif
            print *, "    Itteration", n, "of", itterations, "energy:", realpart(minEnergy)/electronVolt, "eV"
            write(1,*) n,",", realpart(minEnergy)/electronVolt
            write(2,*) n,",", abs(newCoefficients(1))**2
            write(3,*) n,",", abs(newCoefficients(2))**2
        enddo
        coefficients = newCoefficients

        close(unit=1)
        close(unit=2)
        close(unit=3)
    end subroutine

    ! Optimize the basis using a steepest decent method, initially modifying
    ! coefficients by +/- startingShunt, then reducing the shunts to narrow
    ! in on the optimal set
    subroutine optimizeBasisSteepestDecent(basis, coefficients, potential, itterations, startingShunt)
    implicit none
    class(basisState)          :: basis(:)
    procedure(pot)             :: potential
    complex(prec), allocatable :: newCoefficients(:)
    real(prec),    allocatable :: gradReal(:), gradImag(:)
    real(prec),    parameter   :: eps = 0.001
    real(prec)                 :: samples(3,metroSamples), shunt, itterVar, newVar, workVar, startingShunt
    complex(prec)              :: itterEnergy, workE, coefficients(:)
    integer                    :: i, n, m, itterations
    logical                    :: energyDecrease

        ! TODO - remove unnecassary calls to energetics(), should use results
        ! of previous itteration to calculate gradients

        ! Allocate space
        allocate(newCoefficients(size(basis)))
        allocate(gradReal(size(basis)))
        allocate(gradImag(size(basis)))

        open(unit=1,file="output/steepestDecentOptimizationEnergy")
        open(unit=2,file="output/steepestDecentOptimizationBasis1Char")
        open(unit=3,file="output/steepestDecentOptimizationBasis2Char")

        shunt = startingShunt
        
        do n=1,itterations

            call normalizeCoefficients(coefficients) ! Normalize our starting coefficients
            newCoefficients = coefficients ! Reset working set used to calculate gradients

            ! Generate samples in coordinate space using
            ! the metropolis algorithm for the current
            ! trial wavefunction
            call metro(samples, basis, coefficients)    
                
            ! Calculate the energetics with the current coefficients
            call energetics(basis, coefficients, potential, samples, itterEnergy, itterVar)           

            ! Calculate the real and imaginary coefficient
            ! gradients of the cost function
            do i=1,size(basis)

                ! Real gradients
                newCoefficients(i) = coefficients(i) + (eps, 0)
                call energetics(basis, newCoefficients, potential, samples, workE, workVar)
                gradReal(i) = realpart(workE-itterEnergy) / eps
                newCoefficients(i) = coefficients(i) ! Reset
                            
                ! Imaginary gradients
                newCoefficients(i) = coefficients(i) + (0, eps)               
                call energetics(basis, newCoefficients, potential, samples, workE, workVar)
                gradImag(i) = realpart(workE-itterEnergy) / eps
                newCoefficients(i) = coefficients(i) ! Reset                

            enddo

            ! Shunt the coefficients closer to the minimum
            ! by about the size of the minimum coefficient
            call makeOrderUnityReal(gradReal)
            call makeOrderUnityReal(gradImag)
            gradReal = gradReal * shunt
            gradImag = gradImag * shunt

            ! Move along steepest decent
            newCoefficients = coefficients - gradReal! - gradImag*(0,1)
            call energetics(basis, newCoefficients, potential, samples, workE, newVar)

            ! Accept the new coefficients if they reduced the cost function
            if (realpart(workE) < realpart(itterEnergy)) then
                energyDecrease = .true.
                coefficients = newCoefficients
            else
                energyDecrease = .false.
                shunt = shunt * 0.5 ! Reduce the coefficient shunt, to search more finely
            endif

            write(1,*) n,",", min(realpart(workE),realpart(itterEnergy))/electronVolt
            write(2,*) n,",", abs(coefficients(1))**2
            write(3,*) n,",", abs(coefficients(2))**2

            print *, ""
            print *, ""
            print *, ""
            print *, "Steepest decent basis opimization"
            print *, "Itteration:", n
            print *, "Energy decrease:", energyDecrease
            print *, "    Energy:  ", realpart(itterEnergy)/electronVolt, "eV"
            print *, "    Variance:", min(newVar, itterVar)/(electronVolt**2), "eV^2"
            print *, ""
            print *, "Basis character:"
            do m=1, size(basis)
                print *, "        Basis", m, "character", abs(coefficients(m))**2, coefficients(m)
            enddo
            print *, ""
            print *, "Basis changes (shunt =",shunt,")"
            do m=1, size(basis)
                print *, "        Real, imag change", m, -gradReal(m), -gradImag(m)
            enddo

        enddo

        close(unit=1)
        close(unit=2)
        close(unit=3)
        
    end subroutine

    ! Get data for plots of quantities vs basis character
    ! for two basis states
    subroutine characterPlots(basis, potential)
    implicit none
        class(basisState) :: basis(2)
        procedure(pot)    :: potential        
        complex(prec)     :: coefficients(2), energy
        real(prec)        :: variance, samples(3,metroSamples)
        integer           :: s1, grid

        print *, "Creating character plot data..."
        open(unit=1,file="output/energyVsCharacter")
        open(unit=2,file="output/varianceVsCharacter")
        grid = 100
        do s1=1,grid
            print *, "Progress: ", s1, "/", grid
            coefficients(1) = s1/real(grid)
            coefficients(2) = 1-coefficients(1)
            call normalizeCoefficients(coefficients)
            call metro(samples, basis, coefficients)
            call energetics(basis, coefficients, potential, samples, energy, variance)
            write(1,*)  abs(coefficients(1))**2,",",realpart(energy)/electronVolt
            write(2,*)  abs(coefficients(1))**2,",",variance/(electronVolt**2)
        enddo
        close(unit=1)
        close(unit=2)

    end subroutine

    ! Output wavefunction values to a file for plotting etc
    subroutine writeWavefunctionToFile(basis, coefficients)
    implicit none
    class(basisState) :: basis(:)
    complex(prec)     :: coefficients(:)
    integer    :: xi, yi
    integer, parameter :: grid = 100
    real(prec) :: r(3), a
        print *, ""
        print *, "Writing wavefunction to file..."
        open(unit=1,file="output/wavefunctionValues")
        do xi=-grid,grid
            do yi=-grid,grid
                r(1) = 5*angstrom*xi/real(grid)
                r(3) = 5*angstrom*yi/real(grid)
                a = abs(wavefunction(basis, coefficients, r))**2
                if (.not. isnan(a)) then
                    write(1,*) r(1)/angstrom,",",r(3)/angstrom,",",a
                endif
            enddo
        enddo
        close(unit=1)
    end subroutine

    ! Output a sampled wavefuntion to a file for plotting etc
    subroutine sampleWavefunctionToFile(basis, coefficients)
    implicit none
        class(basisState) :: basis(:)    
        complex(prec)     :: coefficients(:)
        real(prec) :: samples(3,metroSamples)
        integer    :: i
        open(unit=1,file="output/wavefunctionSamples")
        print *, ""
        print *, "Sampling wavefunction to file..."
        call metro(samples,basis,coefficients)
        do i=1,size(samples,2)
            write(1,*) samples(1,i)/angstrom, ",", samples(3,i)/angstrom
        enddo
        close(unit=1)
    end subroutine

    ! Normalize a coefficient set
    subroutine normalizeCoefficients(set)
    implicit none
        complex(prec) :: set(:)
        real(prec)    :: norm
        integer       :: i
        norm = 0
        do i=1,size(set)
            norm = norm + abs(set(i))**2
        enddo
        set = set/sqrt(norm)
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
            r      = rand()*maxMetroJump

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
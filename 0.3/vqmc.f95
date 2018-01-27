! Version 0.2
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
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

contains

    ! Optimize the coefficients in the given basis to
    ! minimize the energy in the given potential
    subroutine optimizeBasis(basis, potential)
    implicit none
        class(basisState)          :: basis(:)
        procedure(real(prec))      :: potential
        complex(prec), allocatable :: coefficients(:), designMatrix(:,:), dTd(:,:)
        integer                    :: i, n, m
        real(prec)                 :: samples(3,10), varE
        complex(prec), parameter   :: oneCplx = 1 ! FORTRAN

        ! Allocate space for our coefficients
        ! initialize them all to 1
        allocate(coefficients(size(basis)))
        do i=1,size(basis)
            coefficients(i) = 1
        enddo

        ! Allocate our matricies for least-squares
        allocate(designMatrix(size(basis),size(samples,2)))
        allocate(dTd(size(basis),size(basis)))

        ! Get a set of samples using the metropolis algorithm
        call metro(samples, basis, coefficients)

        ! Calculate our variational energy <E_L>
        varE = energy(basis, coefficients, potential, samples)
        
        ! Calculate the entries of our design matrix
        print *, "singleBASIS"
        do n=1,size(basis)
            do m=1,size(samples,2)
                designMatrix(m,n) = &
                    localEnergy(basis((/n/)), (/oneCplx/), &
                    potential, samples(:,m)) / &
                    wavefunction(basis,coefficients,samples(:,m))
            enddo
        enddo

        ! Calculate D^T*D, our least squares projection matrix
        dTd = matmul(transpose(designMatrix), designMatrix)
        !call potrf(dTd)
        !call potri(dTd)

        ! Print the energy with the given coefficients etc..
        print *, energy(basis, coefficients, potential, samples)/electronVolt

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
        call basis(1)%printDebugInfo()
        localEnergy = localKineticEnergy(basis, coefficients, x) + potential(x)
    end function

    ! Carry out a monte carlo integration of the local energy of the
    ! wavefunction in the given potential, using the sample positions
    ! provided
    function energy(basis, coefficients, potential, samples) result(ret)
    implicit none
        procedure(real(prec)) :: potential
        class(basisState)     :: basis(:)
        complex(prec)         :: ret, coefficients(:)
        real(prec) :: samples(:,:)
        integer    :: i

        ret = 0
        do i=1, size(samples,2)
            ret = ret + localEnergy(basis, coefficients, potential, samples(:,i))
        enddo
        ret = ret/size(samples,2)

    end function

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
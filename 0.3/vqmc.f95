! Version 0.2
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
implicit none

    ! Flag to output the sampled points
    logical :: DEBUG_LOG = .false.

    ! The number of monte-carlo itterations
    integer, parameter :: MC_ITTER = 1000000
    integer, parameter :: MC_INIT_STEPS = 1000

    ! The minimum distance from the divergent origin
    ! our trial electrons are allowed
    real(prec), parameter :: minR = angstrom/10000

contains

    ! Carry out a monte carlo integration of the local energy of the
    ! wavefunction in the given potential, using importance sampling
    ! of the wavefunction
    function energy(wavefunction, potential) result(ret)
    implicit none
        procedure(real(prec)) :: wavefunction, potential
        real(prec) :: x(3), ret, localE
        integer    :: i

        if (DEBUG_LOG) open(unit=1,file="log")
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
            localE = localEnergy(wavefunction, potential, x)
            if (DEBUG_LOG) write(1,*) localE,","
            ret = ret + localE
        enddo
        ret = ret/MC_ITTER

    end function

    ! Calculate the second derivative of the wavefunction
    ! in the given direction using finite differences
    function secondDerivative(wavefunction, x, dir)
        procedure(real(prec)) :: wavefunction
        real(prec) :: x(3), dir(3), secondDerivative
        real(prec), parameter :: eps = 0.00001 * angstrom
        dir = eps * dir/norm2(dir)
        secondDerivative = wavefunction(x+dir) - 2*wavefunction(x) + wavefunction(x-dir)
        secondDerivative = secondDerivative / (eps**2)
    end function

    ! Calculate psi^-1(x) T psi(x)
    function localKineticEnergy(wavefunction, x)
        procedure(real(prec)) :: wavefunction
        real(prec) :: x(3), localKineticEnergy, dx, dy, dz
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
        procedure(real(prec)) :: wavefunction, potential
        real(prec) :: x(3), localEnergy
        localEnergy = localKineticEnergy(wavefunction,x) + potential(x)
    end function

    ! Make a trial metropolis move on x
    subroutine metro(x, wavefunction)
    implicit none
        real(prec) :: x(3), newX(3) , dx(3), r, theeta, phi
        procedure (real(prec)) :: wavefunction

        ! Move some distance r in any direction
        theeta = rand()*pi
        phi    = rand()*2*pi
        r      = rand()*angstrom/10

        dx(1) = r*sin(theeta)*cos(phi)
        dx(2) = r*sin(theeta)*sin(phi)
        dx(3) = r*cos(theeta)

        newX = x + dx
        if (rand() < importance(newX, wavefunction)/importance(x, wavefunction)) then
            x = newX ! Accept the move via metropolis criteria
        endif

        ! Ensure x doesn't get too close to a divergent origin
        if (norm2(x) < minR) then
            x(1) = minR
        endif
    end subroutine

    ! The importance function used in importance sampling
    ! in our case this is our electronic wavefunction
    function importance(x, wavefunction)
    implicit none
        real(prec) :: x(3), importance, wavefunction
        importance = wavefunction(x)**2
    end function

end module vqmc

program main
use vqmc
use particleInBox
use atomicBasis
implicit none

    real(prec) :: start, end, x
    integer :: i, n, l

    !open(unit=2,file="toPlot")
    !call debugRadialPart(2)
    !call debugAssociatedLegendrePolynomial(2)

    open(unit=3,file="toPlot2D")
    call debugWavefunctions(3)

    return

    DEBUG_LOG = .false.

    ! Calculate the ground state energy using vqmc
    call cpu_time(start)
    print *, "Calculated system energy (eV):", energy(groundState, potential)/electronVolt
    call cpu_time(end)
    print *, "Monte-carlo itterations:", MC_ITTER  
    print *, "Elapsed time:           ", end-start

end program main
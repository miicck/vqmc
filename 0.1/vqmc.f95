! Version 0.1
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units. Manages 7
! signifigant figures of the rydberg in ~ 0.5s!

! Constants used by the program
module constants
implicit none

    integer, parameter :: prec = selected_real_kind(15,307)
    real(prec), parameter :: pi = 3.14159265358979
    real(prec), parameter :: hbar = 1.054571800E-34
    real(prec), parameter :: mElectron = 9.10938356E-31
    real(prec), parameter :: qElectron = 1.6021766208E-19
    real(prec), parameter :: mProton = 1.672621898E-27
    real(prec), parameter :: epsNaught = 8.854187817620E-12
    real(prec), parameter :: electronVolt = qElectron
    real(prec), parameter :: angstrom = 1E-10
    real(prec), parameter :: a0 = 5.29177E-11 ! Bohr radius

end module

! Things to do with the hydrogent ground state
module hydrogenGroundState
use constants
implicit none

contains

    ! hydrogen ground state
    ! Wavefunction normalized s.t
    ! int 4*pi*r^2*|wavefunction|^2 = 1
    function wavefunction(x)
    implicit none
        real(prec) :: x(3), wavefunction
        wavefunction = exp(-norm2(x)/a0)/sqrt(pi*a0**3)
    end function

    ! Calculate psi^-1(x) laplacian psi(x)
    function localKineticEnergy(x)
    implicit none
        real(prec) :: x(3),r, localLaplacian, localKineticEnergy
        r = norm2(x)
        localLaplacian = (r-2*a0)*exp(-r/a0)/(r*sqrt(pi*(a0**7)))
        localLaplacian = localLaplacian / wavefunction(x)
        localKineticEnergy = - (hbar**2/(2*mElectron)) * localLaplacian
    end function

    ! Calculate V(x)
    function localPotentialEnergy(x)
    implicit none
        real(prec) :: x(3), r, localPotentialEnergy
        r = norm2(x)
        localPotentialEnergy = - (qElectron**2)/(4*pi*epsNaught*r)
    end function    

    ! Calculate psi^-1(x) H psi(x)
    function localEnergy(x)
        implicit none
            real(prec) :: x(3), localEnergy
            localEnergy = localKineticEnergy(x) + localPotentialEnergy(x)
    end function
    
end module

! Module that carries out the monte-carlo integration
module vqmc
use hydrogenGroundState
implicit none

    ! The number of monte-carlo itterations
    integer,    parameter :: MC_ITTER = 1000000

    ! The minimum distance from the divergent origin
    ! our trial electrons are allowed
    real(prec), parameter :: minR = angstrom/10000

contains

    ! Carry out a monte carlo integration of the local
    ! function f, using the given importance sampling
    function integrate(f) result(ret)
    implicit none
        real(prec) :: x(3), ret, f
        integer    :: i
        ret = 0
        x = 0
        do i=1, MC_ITTER
            call metro(x)
            ret = ret + f(x)
        enddo
        ret = ret/MC_ITTER
    end function

    ! Make a trial metropolis move on x
    subroutine metro(x)
    implicit none
        real(prec) :: x(3), newX(3) , dx(3), r, theeta, phi

        ! Move some distance r in any direction
        theeta = rand()*pi
        phi    = rand()*2*pi
        r      = rand()*angstrom/10

        dx(1) = r*sin(theeta)*cos(phi)
        dx(2) = r*sin(theeta)*sin(phi)
        dx(3) = r*cos(theeta)

        newX = x + dx
        if (rand() < importance(newX)/importance(x)) then
            x = newX ! Accept the move via metropolis criteria
        endif

        ! Ensure x doesn't get too close to a divergent origin
        if (norm2(x) < minR) then
            x(1) = minR
        endif
    end subroutine

    ! The importance function used in importance sampling
    ! in our case this is our electronic wavefunction
    function importance(x)
    implicit none
        real(prec) :: x(3), importance
        importance = wavefunction(x)**2
    end function

end module vqmc

program main
use vqmc
implicit none

    real(prec) :: start, end

    ! Calculate the ground state energy of the
    ! hydrogen atom using vqmc
    call cpu_time(start)
    print *, "Calculated rydberg (eV):", integrate(localEnergy)/electronVolt
    call cpu_time(end)
    print *, "Monte-carlo itterations:", MC_ITTER    
    print *, "Elapsed time:           ", end-start

end program main
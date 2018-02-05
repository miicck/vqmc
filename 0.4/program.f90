program main
use vqmc
use atomicBasis
implicit none

    integer, parameter :: ITTER = 10000
    real(prec) :: startT, endT, time

    call cpu_time(startT)

    call hydrogen()
    call helium()
    call beryllium()
    !call neon() ! Takes a while...

    call cpu_time(endT)
    time = endT-startT  

    print *, ""
    print *, "Performance information, total time: ",time, "seconds"
    print *, "    CPU calculating atomic states: ", 100*atomicCPUtime/time, "%"
    print *, "    CPU calculating permutations : ", 100*permutationCPUtime/time, "%"

contains

    ! ----- HYDROGEN ----- !

    function hydrogenPotential(x)
    real(prec) :: hydrogenPotential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        hydrogenPotential = -qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

    subroutine hydrogen()
    implicit none

        type(atomicState), allocatable :: upBasis(:), downBasis(:)
        real(prec)    :: startT, endT
        complex(prec) :: energy

        ! Basis for hydrogen is 1s (it's ground state)
        allocate(downBasis(0))     

        allocate(upBasis(1))
        upBasis(1) = atomicState(1,0,0,1)

        call cpu_time(startT)
        call monteCarloEnergetics(upBasis, downBasis, slaterDeterminant, hydrogenPotential, ITTER, energy)
        call cpu_time(endT)
        print *, "Hydrogen energy:", realpart(energy)/electronVolt, &
                "eV Calculation took:", (endT-startT), "seconds"

    end subroutine

    ! ----- HELIUM ----- !

    function heliumPotential(x)
    real(prec) :: heliumPotential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        heliumPotential = -2*qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

    subroutine helium()
    implicit none

        type(atomicState), allocatable :: upBasis(:), downBasis(:)
        real(prec)    :: startT, endT
        complex(prec) :: energy

        ! Basis for helium 1s^2
        allocate(upBasis(1))
        upBasis(1)   = atomicState(1,0,0,2)

        allocate(downBasis(1))
        downBasis(1) = atomicState(1,0,0,2)

        call cpu_time(startT)
        call monteCarloEnergetics(upBasis, downBasis, slaterDeterminant, heliumPotential, ITTER, energy)
        call cpu_time(endT)
        print *, "Helium energy:", realpart(energy)/electronVolt, &
                "eV Calculation took:", (endT-startT), "seconds"

    end subroutine

    ! ----- BERYLLIUM ----- !

    function berylliumPotential(x)
    real(prec) :: berylliumPotential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        berylliumPotential = -4*qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

    subroutine beryllium()
    implicit none
        type(atomicState), allocatable :: upBasis(:), downBasis(:)
        real(prec)    :: startT, endT
        complex(prec) :: energy

        ! Basis for beryllium 1s^2 2s^2
        allocate(upBasis(2))
        upBasis(1)   = atomicState(1,0,0,4)
        upBasis(2)   = atomicState(2,0,0,4)

        allocate(downBasis(2))
        downBasis(1) = atomicState(1,0,0,4)
        downBasis(2) = atomicState(2,0,0,4)

        call cpu_time(startT)
        call monteCarloEnergetics(upBasis, downBasis, slaterDeterminant, berylliumPotential, ITTER, energy)
        call cpu_time(endT)
        print *, "Beryllium energy:", realpart(energy)/electronVolt, &
                 "eV Calculation took:", (endT-startT), "seconds"

    end subroutine

    ! ----- NEON ----- !

    function neonPotential(x)
    real(prec) :: neonPotential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        neonPotential = -10*qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

    subroutine neon()
    implicit none
        type(atomicState), allocatable :: upBasis(:), downBasis(:)
        real(prec)    :: startT, endT
        complex(prec) :: energy

        ! Basis for neon 1s^2 2s^2 2p^6
        allocate(upBasis(5))
        upBasis(1)   = atomicState(1,0,0,10)
        upBasis(2)   = atomicState(2,0,0,10)
        upBasis(3)   = atomicState(2,1,-1,10)
        upBasis(4)   = atomicState(2,1,0,10)
        upBasis(5)   = atomicState(2,1,1,10)

        allocate(downBasis(5))
        downBasis(1) = atomicState(1,0,0,10)
        downBasis(2) = atomicState(2,0,0,10)
        downBasis(3) = atomicState(2,1,-1,10)
        downBasis(4) = atomicState(2,1,0,10)
        downBasis(5) = atomicState(2,1,1,10)

        call cpu_time(startT)
        call monteCarloEnergetics(upBasis, downBasis, slaterDeterminant, neonPotential, ITTER, energy)
        call cpu_time(endT)
        print *, "Neon energy:", realpart(energy)/electronVolt, &
                "eV Calculation took:", (endT-startT), "seconds"

    end subroutine

end program main

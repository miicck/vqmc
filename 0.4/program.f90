program main
use vqmc
use atomicBasis
implicit none

    integer, parameter :: ITTER = 10000
    real(prec) :: startT, endT, time

    call cpu_time(startT)
  
    call hydrogen()
    call H2plusIon()
    call H2molecule()
    call beryllium()

    call cpu_time(endT)
    time = endT-startT  

    print *, ""
    print *, "Performance information, total time: ",time, "seconds"
    print *, "    CPU calculating atomic states: ", 100*atomicCPUtime/time, "%"
    print *, "    CPU calculating permutations : ", 100*permutationCPUtime/time, "%"

contains

    ! ----- HYDROGEN ----- !

    ! Hydrogen atomic potential
    function hydrogenPotential(x)
    real(prec) :: hydrogenPotential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        hydrogenPotential = -qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

    subroutine hydrogen()
    implicit none
        complex(prec) :: energy
 
        ! Setup our basis and many body wavefunction
        allocate(basis(1))
        basis(1)%state => atomicState(1,0,0,1)

        ! Describe our system
        upElectrons    = 1
        potential      => hydrogenPotential

        call initialize()
        call monteCarloEnergetics(ITTER, energy)
        call cleanUp()
        print *, "Hydrogen energy:", realpart(energy)/electronVolt, "(eV)"

    end subroutine

    ! Hydrogen molecular ion potential
    function H2plusPotential(x)
    real(prec) :: H2plusPotential, x(3)
    real(prec), parameter :: zero = 0, atom2(3) = (/zero,zero,angstrom*1.05687/)
        H2plusPotential = hydrogenPotential(x) + hydrogenPotential(x-atom2)
    end function

    subroutine H2plusIon()
    implicit none
        complex(prec) :: energy

        ! Setup our basis and many body wavefunction
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,1)
        basis(2)%state => atomicState(1,0,0,1)
        basis(2)%centre(3) = angstrom*1.05687

        ! Describe our system
        upElectrons      = 1
        potential        => H2plusPotential
        initialCharacter => explicitCharacter ! Use explicit electronic characters
        
        allocate(upCharacters(2,1))
        upCharacters(1,1) = 1
        upCharacters(2,1) = -1
        
        allocate(downCharacters(2,0))

        call initialize()
        call monteCarloEnergetics(ITTER, energy)
        !call sampleElectronPositionsToFile(ITTER,1,.true.)
        call cleanUp()        
        print *, "Hydrogen molecular ion energy:", realpart(energy)/electronVolt, "(eV)"

    end subroutine

    ! Hydrogen molecular potential
    function H2potential(x)
    real(prec) :: H2potential, x(3)
    real(prec), parameter :: zero = 0, atom2(3) = (/zero,zero,angstrom*0.741/)
        H2potential = hydrogenPotential(x) + hydrogenPotential(x-atom2)
    end function

    subroutine H2molecule()
    implicit none
        complex(prec) :: energy

        ! Setup our basis and many body wavefunction
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,1)
        basis(2)%state => atomicState(1,0,0,1)
        basis(2)%centre(3) = angstrom*0.741

        ! Describe our system
        upElectrons      = 1
        downElectrons    = 1
        potential        => H2potential
        initialCharacter => explicitCharacter ! Use explicit electronic characters
        
        allocate(upCharacters(2,1))
        upCharacters(1,1) = 1
        upCharacters(2,1) = -1
        
        allocate(downCharacters(2,1))
        downCharacters(1,1) = 1
        downCharacters(2,1) = -1

        call initialize()        
        call monteCarloEnergetics(ITTER, energy)
        call cleanUp()        
        print *, "Hydrogen molecule energy:", realpart(energy)/electronVolt, "(eV)"

    end subroutine

    ! ----- BERYLLIUM ----- !

    function berylliumPotential(x)
    real(prec) :: berylliumPotential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        berylliumPotential = -4*qElectron**2/(4*pi*epsNaught*norm2(x))
    end function
    
    subroutine beryllium()
    implicit none
        complex(prec) :: energy

        ! Setup our basis and many body wavefunction
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,4)
        basis(2)%state => atomicState(2,0,0,4)

        ! Describe our system
        upElectrons    = 2
        downElectrons  = 2
        potential      => berylliumPotential

        call initialize()
        call monteCarloEnergetics(ITTER, energy)
        call cleanUp()
        print *, "Beryllium energy:", realpart(energy)/electronVolt, "(eV)"

    end subroutine

end program main

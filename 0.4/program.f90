program main
use vqmc
use atomicBasis
implicit none

    real(prec) :: startT, endT, time
    real(prec) :: h2plusBondLength = angstrom*1.05687

    energyUnit = hartree
    metroSamples = 10000

    call cpu_time(startT)
  
    !call hydrogen()
    !call hydrogen1s2sMixing()
    !call H2plusIon()
    !call H2plusIonBondLength()
    call helium()
    !call beryllium()

    call cpu_time(endT)
    time = endT-startT  

    print *, ""
    print *, "Performance information, total time: ",time, "seconds"
    print *, "    CPU calculating atomic states: ", 100*atomicCPUtime/time, "%"
    print *, "    CPU calculating permutations : ", 100*permutationCPUtime/time, "%"

contains

    ! ----- HYDROGEN ----- !

    ! Hydrogen atomic potential
    subroutine hydrogen()
    implicit none
        complex(prec) :: energy
 
        ! Setup our basis
        allocate(basis(1))
        basis(1)%state => atomicState(1,0,0,1)

        ! Describe our system
        upElectrons    = 1
        allocate(nucleii(1))
        nucleii(1)%centre = 0
        nucleii(1)%charge = 1

        call initialize()
        call monteCarloEnergetics()
        call cleanUp()
        print *, "Hydrogen"
        call printLastEnergetics()

    end subroutine

    ! Get data for hydrogen exited state mixing
    subroutine hydrogen1s2sMixing()
    implicit none

        real(prec) :: c = 0
        integer    :: i, grid, n, nMax, b

        open(unit=2,file="hydrogenMixingEnergies")
        open(unit=3,file="hydrogenMixingVariances")
        open(unit=4,file="hydrogenMixingErrors")

        grid = 20
        nMax = 10
        do n=1,nMax-1
            do i=0,grid

                c = i/real(grid)

                ! Setup our basis
                allocate(basis(nMax))
                do b=1,nMax
                    basis(b)%state => atomicState(b,0,0,1)
                enddo

                ! Describe our system
                upElectrons      = 1
                initialCharacter => explicitCharacter ! Use explicit electronic characters
                allocate(nucleii(1))
                nucleii(1)%centre = 0
                nucleii(1)%charge = 1

                allocate(upCharacters(nMax,1))
                upCharacters = 0
                upCharacters(n,1) = 1-c
                upCharacters(n+1,1) = c

                call initialize()
                call monteCarloEnergetics()
                call cleanUp()
                print *, "Hydrogen mixing progress:",100*((n-1+c)/real(nMax-1)),"%"
                write(2,*) n+c,",",realpart(mcEnergyLast)/energyUnit
                write(3,*) n+c,",",realpart(mcReblockedVarianceLast)/(energyUnit**2)
                write(4,*) n+c,",",realpart(mcEnergyErrorLast)/(energyUnit**2)

            enddo
        enddo
        close(unit=2)
        close(unit=3)
        close(unit=4)
    end subroutine

    subroutine H2plusIon()
    implicit none
        complex(prec) :: energy
        integer :: i, grid

        ! Setup our basis and many body wavefunction
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,1)
        basis(2)%state => atomicState(1,0,0,1)
        basis(2)%centre(3) = h2plusBondLength

        ! Describe our system
        upElectrons      = 1
        allocate(nucleii(2))
        nucleii(1)%centre    = 0
        nucleii(1)%charge    = 1
        nucleii(2)%centre(3) = h2plusBondLength
        nucleii(2)%charge    = 1

        ! Use explicit electronic characters
        initialCharacter => explicitCharacter
        allocate(upCharacters(2,1))
        upCharacters(1,1) = 1 ! 1 for bonding or antibonding
        upCharacters(2,1) = -1 ! 1 for bonding, -1 for antibonding
        allocate(downCharacters(2,0))

        call initialize()
        call monteCarloEnergetics()
        call sampleElectronPositionsToFile(1,.true.)        
        call cleanUp()
        call printLastEnergetics()

    end subroutine

    subroutine H2plusIonBondLength()
    implicit none
        complex(prec) :: energy
        integer :: i, grid

        open(unit=2,file="H2ionEnergyVsBondLength")
        open(unit=3,file="H2ionVarianceVsBondLength")
        grid = 500
        do i=0,grid

            print *, "Ion bond length progress:",100*i/real(grid),"%"

            h2plusBondLength = 10*angstrom * i/real(grid)

            ! Setup our basis and many body wavefunction
            allocate(basis(2))
            basis(1)%state => atomicState(1,0,0,1)
            basis(2)%state => atomicState(1,0,0,1)
            basis(2)%centre(3) = h2plusBondLength

            ! Describe our system
            upElectrons      = 1
            allocate(nucleii(2))
            nucleii(1)%centre    = 0
            nucleii(1)%charge    = 1
            nucleii(2)%centre(3) = h2plusBondLength
            nucleii(2)%charge    = 1

            ! Use explicit electronic characters
            initialCharacter => explicitCharacter
            allocate(upCharacters(2,1))
            upCharacters(1,1) = 1 ! 1 for bonding or antibonding
            upCharacters(2,1) = -1 ! 1 for bonding, -1 for antibonding  
            allocate(downCharacters(2,0))

            call initialize()
            call monteCarloEnergetics()
            call cleanUp()

            write (2,*) h2plusBondLength/angstrom, ",", realpart(mcEnergyLast)/energyUnit
            write (3,*) h2plusBondLength/angstrom, ",", realpart(mcReblockedVarianceLast)/energyUnit
        enddo
        close(unit=2)
        close(unit=3)
    
    end subroutine

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
        initialCharacter => explicitCharacter ! Use explicit electronic characters
        
        allocate(upCharacters(2,1))
        upCharacters(1,1) = 1
        upCharacters(2,1) = -1
        
        allocate(downCharacters(2,1))
        downCharacters(1,1) = 1
        downCharacters(2,1) = -1

        call initialize()        
        call monteCarloEnergetics()
        call cleanUp()        
        print *, "H2 Molecule"
        call printLastEnergetics()

    end subroutine

    ! ----- HELIUM ----- !

    subroutine helium()
    implicit none
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,2)
        basis(2)%state => atomicState(1,0,0,2)

        upElectrons = 1
        downElectrons = 1

        allocate(nucleii(1))
        nucleii(1)%charge = 2

        !manyBodyMethod => slaterDeterminant ! Don't use a Jastow factor

        call initialize()
        !call optimizeJastrow()
        call monteCarloEnergetics()
        print *, "Helium"
        call printLastEnergetics()
        call cleanup()
    end subroutine

    ! ----- BERYLLIUM ----- !
    
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

        call initialize()
        call monteCarloEnergetics()
        call cleanUp()
        print *, "Beryllium"
        call printLastEnergetics()

    end subroutine

end program main

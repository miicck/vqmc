program main
use vqmc
use atomicBasis
use mpi
implicit none

    real(prec) :: startT, endT, time
    integer    :: ierr, i, rank, nproc, clock

    call mpi_init(ierr)

    ! Set a different random seed for each mpi process.
    ! Also randomize the time via the clock so different
    ! runs will give (hopefully slightly) different answers
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, nproc, ierr)
    call system_clock(clock)
    call srand(rank+clock)
    if (rank == 0) call cpu_time(startT)

    ! Carry out our calculation
    energyUnit = hartree
    metroSamples = 100000/nproc
    !printMetroProgress = .true.

    call hydrogen()
    !call helium()
    !call heliumExtendedBasis()
    !call lithium()
    !call beryllium()
    !call boron()

    !call hydrogen1s2sMixing()
    !call H2plusIon()
    !call H2plusIonBondLength()
    !call H2molecule()
    !call H2moleculeBondLength()

    ! Output performance info for rank 0
    if (rank == 0) then
        call cpu_time(endT)
        time = endT-startT    
        print *, ""
        print *, "Performance information for mpi rank 0, total time: ",time, "seconds"
        print *, "    CPU calculating atomic states: ", 100*atomicCPUtime/time, "%"
        print *, "    CPU calculating permutations : ", 100*permutationCPUtime/time, "%"      
    endif

    ! End mpi session
    call mpi_finalize(ierr)
contains

    ! ----- HYDROGEN ----- !

    ! Hydrogen atomic potential
    subroutine hydrogen()
    implicit none
        
 
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
        call printLastEnergetics()

    end subroutine

    ! ----- HELIUM ----- !

    subroutine helium()
    implicit none
        allocate(basis(1))
        basis(1)%state => atomicState(1,0,0,2)

        upElectrons = 1
        downElectrons = 1

        allocate(nucleii(1))
        nucleii(1)%charge = 2

        !manyBodyMethod => slaterDeterminant ! Don't use a Jastow factor

        call initialize()
        !call optimizeJastrow()
        call monteCarloEnergetics()
        call printLastEnergetics()
        call cleanup()

    end subroutine

    subroutine heliumExtendedBasis()
    implicit none
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,2)
        basis(2)%state => atomicState(1,0,0,1)

        upElectrons = 1
        downElectrons = 1

        allocate(nucleii(1))
        nucleii(1)%charge = 2

        call initialize()
        call setAllJastrowParams(real(1.30672,kind=prec))
        call optimizeCharacters()
        call monteCarloEnergetics()
        call printLastEnergetics()
        call cleanup()
    end subroutine

    ! ----- LITHIUM ----- !

    subroutine lithium()
    implicit none

        ! Setup our basis 1s2 2s
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,3)
        basis(2)%state => atomicState(2,0,0,3)

        ! Describe the system
        upElectrons   = 2
        downElectrons = 1

        allocate(nucleii(1))
        nucleii(1)%charge = 3

        call initialize()
        call monteCarloEnergetics()
        call printLastEnergetics()
        call cleanUp()    

    end subroutine

    ! ----- BERYLLIUM ----- !
    
    subroutine beryllium()
    implicit none

        ! Setup our basis 1s2 2s2
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,4)
        basis(2)%state => atomicState(2,0,0,4)

        ! Describe our system
        upElectrons    = 2
        downElectrons  = 2

        allocate(nucleii(1))
        nucleii(1)%charge = 4

        call initialize()
        call monteCarloEnergetics()
        call printLastEnergetics()
        call cleanUp()

    end subroutine

    ! ----- BORON ----- !

    subroutine boron()
    implicit none

        ! Setup our basis 1s2 2s2
        allocate(basis(3))
        basis(1)%state => atomicState(1,0,0,4)
        basis(2)%state => atomicState(2,0,0,4)
        basis(3)%state => atomicState(2,1,-1,4)

        ! Describe our system
        upElectrons    = 3
        downElectrons  = 2

        allocate(nucleii(1))
        nucleii(1)%charge = 5

        call initialize()
        call monteCarloEnergetics()
        call printLastEnergetics()
        call cleanUp()

    end subroutine

    ! ----- MOLECULES ----- !

    ! Get data for hydrogen exited state mixing
    subroutine hydrogen1s2sMixing()
    implicit none

        real(prec) :: c = 0
        integer    :: i, grid, n, nMax, b, rank, ierr

        call mpi_comm_rank(mpi_comm_world, rank, ierr)

        if (rank==0) then
            open(unit=2,file="hydrogenMixingEnergies")
            open(unit=3,file="hydrogenMixingVariances")
            open(unit=4,file="hydrogenMixingErrors")
        endif

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

                if (rank == 0) then
                    print *, "Hydrogen mixing progress:",100*((n-1+c)/real(nMax-1)),"%"
                    write(2,*) n+c,",",mcEnergyLast/energyUnit
                    write(3,*) n+c,",",mcReblockedVarianceLast/(energyUnit**2)
                    write(4,*) n+c,",",mcEnergyErrorLast/(energyUnit**2)
                endif

            enddo
        enddo
        
        if (rank == 0) then
            close(unit=2)
            close(unit=3)
            close(unit=4)
        endif
    end subroutine

    subroutine H2plusIon()
    implicit none
        
        real(prec), parameter :: h2plusBondLength = angstrom*1.240
        integer :: i, grid

        ! Setup our basis and many body wavefunction
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,1)
        basis(2)%state => atomicState(1,0,0,1)
        basis(2)%centre = (/0.0D0,0.0D0,h2plusBondLength/)

        ! Describe our system
        upElectrons          = 1
        allocate(nucleii(2))
        nucleii(1)%centre    = 0
        nucleii(1)%charge    = 1
        nucleii(2)%centre    = (/0.0D0,0.0D0,h2plusBondLength/)

        print *, nucleii(2)%centre

        ! Use explicit electronic characters
        initialCharacter => explicitCharacter
        allocate(upCharacters(2,1))
        upCharacters(1,1) = 1 ! 1 for bonding or antibonding
        upCharacters(2,1) = 1 ! 1 for bonding, -1 for antibonding
        allocate(downCharacters(2,0))

        call initialize()
        call monteCarloEnergetics()
        call sampleElectronPositionsToFile(1,.true.)
        call printLastEnergetics()     
        call cleanUp()     

    end subroutine

    subroutine H2plusIonBondLength()
    implicit none
        
        integer :: i, grid, rank, ierr
        real(prec) :: h2plusBondLength = 1
        call mpi_comm_rank(mpi_comm_world, rank, ierr)

        if (rank == 0) then
            open(unit=2,file="H2ionEnergyVsBondLength")
            open(unit=3,file="H2ionVarianceVsBondLength")
        endif

        grid = 500
        do i=0,grid

            print *, "Ion bond length progress:",100*i/real(grid),"%"

            h2plusBondLength = 10*angstrom * i/real(grid)

            ! Setup our basis and many body wavefunction
            allocate(basis(2))
            basis(1)%state => atomicState(1,0,0,1)
            basis(2)%state => atomicState(1,0,0,1)
            basis(2)%centre = (/0.0D0,0.0D0,h2plusBondLength/)

            ! Describe our system
            upElectrons      = 1
            allocate(nucleii(2))
            nucleii(1)%centre    = 0
            nucleii(1)%charge    = 1
            nucleii(2)%centre    = (/0.0D0,0.0D0,h2plusBondLength/)
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

            if (rank == 0) then
                write (2,*) h2plusBondLength/angstrom, ",", mcEnergyLast/energyUnit
                write (3,*) h2plusBondLength/angstrom, ",", mcReblockedVarianceLast/energyUnit
            endif
        enddo

        if (rank == 0) then
            close(unit=2)
            close(unit=3)
        endif
    
    end subroutine

    subroutine H2molecule()
    implicit none

        ! Setup our basis
        allocate(basis(2))
        basis(1)%state => atomicState(1,0,0,1)
        basis(2)%state => atomicState(1,0,0,1)
        basis(2)%centre(3) = angstrom*0.741

        ! Describe our system
        upElectrons      = 1
        downElectrons    = 1

        allocate(nucleii(2))
        nucleii(1)%centre = 0
        nucleii(2)%centre = (/0.0D0,0.0D0,angstrom*0.741/)

        ! Use explicit electronic characters
        initialCharacter => explicitCharacter 
        
        allocate(upCharacters(2,1))
        upCharacters(1,1) = 1
        upCharacters(2,1) = 1
        
        allocate(downCharacters(2,1))
        downCharacters(1,1) = 1
        downCharacters(2,1) = 1

        !manyBodyMethod => slaterDeterminant ! Don't use a Jastow factor

        call initialize()        
        call monteCarloEnergetics()
        call printLastEnergetics()
        call cleanUp()

    end subroutine

    subroutine H2moleculeBondLength()
    implicit none
        integer    :: i, grid, rank, ierr
        real(prec) :: length

        call mpi_comm_rank(mpi_comm_world, rank, ierr)
        if (rank == 0) then
            open(unit=2,file="output/H2EnergyVsBondLength")
            open(unit=3,file="output/H2VarianceVsBondLength")
        endif

        grid = 100
        do i=1,grid

            length = 10*angstrom*i/real(grid)

            ! Setup our basis
            allocate(basis(2))
            basis(1)%state  => atomicState(1,0,0,1)
            basis(2)%state  => atomicState(1,0,0,1)
            basis(2)%centre = (/0.0D0,0.0D0,length/)

            ! Describe our system
            upElectrons      = 1
            downElectrons    = 1

            allocate(nucleii(2))
            nucleii(1)%centre = 0
            nucleii(2)%centre = (/0.0D0,0.0D0,length/)

            ! Use explicit electronic characters
            initialCharacter => explicitCharacter 
            
            allocate(upCharacters(2,1))
            upCharacters(1,1) = 1
            upCharacters(2,1) = 1
            
            allocate(downCharacters(2,1))
            downCharacters(1,1) = 1
            downCharacters(2,1) = 1

            manyBodyMethod => slaterDeterminant ! Don't use a Jastow factor

            call initialize()        
            call monteCarloEnergetics()
            if (rank == 0) then
                write(2,*) length/angstrom, ",", mcEnergyLast/energyUnit
                write(3,*) length/angstrom, ",", mcReblockedVarianceLast/(energyUnit**2)
            endif
            call cleanUp()

        enddo

        if (rank==0) then
            close(unit=2)
            close(unit=3)
        endif

    end subroutine

end program main

! Version 0.3
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
!use IFPORT
implicit none

    real(prec) :: permutationCPUtime = 0

    ! An object that can be used as a basis function
    type, abstract :: basisState
    contains
        procedure(basisWfn), deferred :: value
        procedure(printStateDebug), deferred :: printDebugInfo
    end type

    type :: nucleus
        real(prec) :: centre(3) = 0
        real(prec) :: charge = 1
    contains
        procedure :: potential => nuclearPotential
    end type

    interface

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

        ! A potential
        function pot(x)
            import
            real(prec) :: x(3), pot
        end function

        ! Takes our electron positions (up, down) and
        ! returns a complex amplitude
        function manyBodyWfn(upConfig, downConfig)
            import
            real(prec) :: upConfig(:,:), downConfig(:,:)
            complex(prec) :: manyBodyWfn
        end function

        ! How we assign electron character
        subroutine charAssign()
            import
        end subroutine
    end interface

    ! Container for basisState pointers !FORTRAN
    type basisListElm
        class(basisState), pointer :: state
        real(prec) :: centre(3) = 0
    end type

    ! Elements definining the calculation to carry out
    type(basisListElm), allocatable :: basis(:)
    class(nucleus), allocatable     :: nucleii(:)
    integer                         :: upElectrons = 0, downElectrons = 0
    real(prec), allocatable         :: upCharacters(:,:), downCharacters(:,:)
    procedure(manyBodyWfn), pointer :: manyBodyMethod   => slaterJastrow
    procedure(charAssign), pointer  :: initialCharacter => sequentialCharacter
    integer                         :: metroSamples = 10000
    integer                         :: initMetroSamples = 1000
    integer                         :: reblockingLength = 1000
    real(prec)                      :: maxMetroJump = 4*angstrom ! The maximum distance an electron can move in a metropolis trial move
    complex(prec)                   :: mcEnergyLast = 0
    complex(prec)                   :: mcKineticEnergyLast = 0
    complex(prec)                   :: mcPotentialEnergyLast = 0
    complex(prec)                   :: mcReblockedVarianceLast = 0
    complex(prec)                   :: mcEnergyErrorLast = 0
    complex(prec)                   :: mcNuclearEnergyLast = 0
    real(prec)                      :: energyUnit = electronVolt

contains

    ! A nuclear potential
    function nuclearPotential(this, x) result(ret)
    implicit none
        class(nucleus) :: this
        real(prec)     :: x(3), ret
        ret = -this%charge*(qElectron**2)/(4*pi*epsNaught*norm2(x-this%centre))
    end function

    ! The potential of our system
    function potential(x) result(ret)
    implicit none
        real(prec) :: ret, x(3)
        integer    :: i
        ret = 0
        do i=1,size(nucleii)
            ret = ret + nucleii(i)%potential(x)
        enddo
    end function

    ! Initialize the program
    subroutine initialize()
    implicit none
        ! Carry out initial checks
        if (upElectrons < 0) print *, "Error: upElectrons < 0!"
        if (downElectrons < 0) print *, "Error: downElectrons < 0!"
        if (upElectrons + downElectrons == 0) print *, "Error: no electrons!"

        ! Carry out initialization
        call initialCharacter()
    end subroutine

    ! Clean up memeory after run
    subroutine cleanUp()
    implicit none
        ! Clean up memory from previous run
        if (allocated(upCharacters)) deallocate(upCharacters)
        if (allocated(downCharacters)) deallocate(downCharacters)
        if (allocated(basis)) deallocate(basis)
        if (allocated(nucleii)) deallocate(nucleii)

        ! Reset default parameters
        initialCharacter => sequentialCharacter
        manyBodyMethod   => slaterJastrow
        upElectrons = 0
        downElectrons = 0
    end subroutine

    ! Print the results of the last monteCarloEnergetics call
    subroutine printLastEnergetics()
    implicit none
        print *, ""
        print *, "Energy:", realpart(mcEnergyLast)/energyUnit
        print *, "     of which electron kinetic:  ", realpart(mcKineticEnergyLast)/energyUnit
        print *, "     of which electron Potential:", realpart(mcPotentialEnergyLast)/energyUnit
        print *, "     of which nuclear:           ", realpart(mcNuclearEnergyLast)/energyUnit
        print *, "Standard error in energy:", realpart(mcEnergyErrorLast)/(energyUnit)
        print *, "Reblocked variance:", realpart(mcReblockedVarianceLast)/(energyUnit**2) 
    end subroutine

    ! Use monte carlo integration to get the energetics of the given wavefunction in the given potential
    subroutine monteCarloEnergetics()
    implicit none
        real(prec), allocatable    :: upPos(:,:,:), downPos(:,:,:), allPos(:,:,:)
        integer                    :: i, j, i2, j2, N_u, N_d
        complex(prec)              :: energy, kineticEnergy, potentialEnergy, nuclearEnergy
        complex(prec), allocatable :: localEnergies(:), localKineticEnergies(:), localPotentialEnergies(:)
        real(prec)                 :: dr(3)      

        ! Allocate/fill our arrays
        call metropolis(upPos, downPos)     ! Fill array of up/down electron positions

        N_u = size(upPos, 2)                            ! # up electrons
        N_d = size(downPos, 2)                          ! # down electrons
        allocate(allPos(3, N_u + N_d, metroSamples))    ! Array of all electron positions
        if (N_d > 0) allPos(:,1:N_d,:) = downPos        ! Fill the down part of allPos
        if (N_u > 0) allPos(:,N_d+1:N_d+N_u,:) = upPos  ! Fill the up part of allPos
        allocate(localEnergies(metroSamples))
        allocate(localKineticEnergies(metroSamples))
        allocate(localPotentialEnergies(metroSamples))

        ! Perform monte carlo integration
        energy = 0
        potentialEnergy = 0
        kineticEnergy = 0
        do i=1,metroSamples

            localEnergies(i) = 0
            localKineticEnergies(i) = 0
            localPotentialEnergies(i) = 0

            ! Sum potentials for all electrons
            do j=1,N_u + N_d
                localPotentialEnergies(i) = localPotentialEnergies(i) + potential(allPos(:,j,i))
            enddo

            ! Sum electron - electron coulomb interactions
            do i2=1, N_u + N_d - 1
                do j2=i2+1, N_u + N_d
                    dr = allPos(:,i2,i) - allPos(:,j2,i)
                    localPotentialEnergies(i) = localPotentialEnergies(i) + qElectron**2/(4*pi*epsNaught*norm2(dr))
                enddo
            enddo
            
            ! Add kinetic energy of up and down electrons
            localKineticEnergies(i) = localKinetic(upPos(:,:,i), downPos(:,:,i))
            localEnergies(i) = localKineticEnergies(i) + localPotentialEnergies(i)

            potentialEnergy = potentialEnergy + localPotentialEnergies(i)
            kineticEnergy = kineticEnergy + localKineticEnergies(i)
            energy = energy + localEnergies(i)

        enddo

        mcKineticEnergyLast = kineticEnergy/real(metroSamples)
        mcPotentialEnergyLast = potentialEnergy/real(metroSamples)
        mcEnergyLast = energy/real(metroSamples)
        mcReblockedVarianceLast = reblockedVariance(localEnergies, mcEnergyLast)
        mcEnergyErrorLast = sqrt(mcReblockedVarianceLast/real(metroSamples/reblockingLength))

        ! Add the nuclear energy, it's the same regardless of electron config
        do i=1,size(nucleii)-1
            do j=i+1,size(nucleii)
                ! N.B nuclear potental assumes it's given an electron, modify it for nuclear-nuclear interactions
                nuclearEnergy = -nucleii(i)%potential(nucleii(j)%centre)*nucleii(j)%charge
            enddo
        enddo

        mcEnergyLast = mcEnergyLast + nuclearEnergy
        mcNuclearEnergyLast = nuclearEnergy

    end subroutine

    ! Get the reblocked variance of a set of local energies
    function reblockedVariance(localEnergies, energy)
    implicit none
        complex(prec) :: localEnergies(:), reblockedVariance, energy
        complex(prec), allocatable :: blockAverages(:)
        integer :: i, j, numBlocks, k

        numBlocks = metroSamples/reblockingLength
        allocate(blockAverages(numBlocks))

        do i=1,numBlocks
            blockAverages(i) = 0
            do j=1,reblockingLength
                k = (i-1)*reblockingLength + j
                blockAverages(i) = blockAverages(i) + localEnergies(k)
            enddo
            blockAverages(i) = blockAverages(i)/reblockingLength
        enddo

        reblockedVariance = 0
        do i=1,numBlocks
            reblockedVariance = reblockedVariance + &
                (energy - blockAverages(i))**2
        enddo
        reblockedVariance = reblockedVariance/(numBlocks-1)
        
    end function

    ! Calculate the local kinetic energy of wavefunction
    ! made from the given single particle basis states
    ! evaluated at the given electron configuration
    function localKinetic(upConfig, downConfig) result(ret)
    implicit none
        real(prec)              :: upConfig(:,:), downConfig(:,:)
        real(prec), parameter   :: eps = angstrom / 10E4
        real(prec), allocatable :: deltaUp(:,:), deltaDown(:,:)
        complex(prec)           :: laplacian, valueAtConfig, ret
        integer                 :: i, j

        allocate(deltaUp(size(upConfig,1),size(upConfig,2)))
        allocate(deltaDown(size(downConfig,1),size(downConfig,2)))

        valueAtConfig = manyBodyMethod(upConfig, downConfig)

        ! Calculate laplacian = sum_i(laplacian_i)
        ! where laplacian_i is the laplacian of
        ! the wavefunction w.r.t the ith
        ! electronic position
        laplacian = 0
        deltaUp = 0
        deltaDown = 0

        do i=1,size(upConfig,2) ! i <=> up electrons                 
            do j=1,3 ! j <=> directions
                ! Evaluate d^2/dx_j^2 using finite differences and add it to the laplacian
                deltaUp(j,i) = eps ! Create a displacement of the ith electron in the jth direction
                laplacian = laplacian + manyBodyMethod(upConfig + deltaUp, downConfig) - &
                                        2*valueAtConfig + &
                                        manyBodyMethod(upConfig - deltaUp, downConfig)         
                deltaUp(j,i) = 0   ! Reset
            enddo
        enddo

        do i=1,size(downConfig,2) ! i <=> down electrons                 
            do j=1,3 ! j <=> directions
                ! Evaluate d^2/dx_j^2 using finite differences and add it to the laplacian
                deltaDown(j,i) = eps ! Create a displacement of the ith electron in the jth direction
                laplacian = laplacian + manyBodyMethod(upConfig, downConfig + deltaDown) - &
                                        2*valueAtConfig + &
                                        manyBodyMethod(upConfig, downConfig - deltaDown)         
                deltaDown(j,i) = 0   ! Reset
            enddo
        enddo

        laplacian = laplacian / (eps**2)
        ret = -laplacian*(hbar**2)/(2*mElectron)
        ret = ret / valueAtConfig ! Calculating *local* kinetic energy

    end function

    ! Sample metropolis x, z coordinates to a file for plotting
    subroutine sampleElectronPositionsToFile(electron, isUp)
    implicit none
        integer :: i, electron
        logical :: isUp
        real(prec), allocatable :: upConfigs(:,:,:), downConfigs(:,:,:)

        call metropolis(upConfigs, downConfigs)
        open(unit=2,file="sampledElectronPositions")

        if (isUp) then
            do i=1,size(upConfigs,3)
                write(2,*) upConfigs(1,electron,i),",",upConfigs(3,electron,i)
            enddo
        else        
            do i=1,size(downConfigs,3)
                write(2,*) downConfigs(1,electron,i),",",downConfigs(3,electron,i)
            enddo
        endif

    end subroutine

    ! Sample n configurations from the given wavefunction
    ! built out of the given basis
    subroutine metropolis(upConfigs, downConfigs)
    implicit none
        real(prec), allocatable :: upConfigs(:,:,:), downConfigs(:,:,:) ! # dim, # electrons, # samples
        real(prec), allocatable :: upConfig(:,:), newUpConfig(:,:), downConfig(:,:), newDownConfig(:,:)
        integer                 :: i, N_d, N_u
        real(prec)              :: oldProb, newProb

        N_u = upElectrons
        N_d = downElectrons

        allocate(upConfigs(3,N_u,metroSamples))
        allocate(upConfig(3,N_u))
        allocate(newUpConfig(3,N_u))

        allocate(downConfigs(3,N_d,metroSamples))
        allocate(downConfig(3,N_d))
        allocate(newDownConfig(3,N_d))

        ! Initial configuration is all electrons at the origin
        upConfig = 0
        downConfig = 0
        
        do i=1,metroSamples + initMetroSamples

            ! Make a metropolis trial move (move one electron
            ! by up to maxMetroJump in a random direction)
            newUpConfig = upConfig
            newDownConfig = downConfig

            ! Move up or down electron with probability 0.5
            ! Also deal with the N_u = 0, N_d = 0 cases
            if ((N_u > 0) .and. ((rand()<0.5) .or. (N_d==0))) then
                call makeMetropolisMove(newUpConfig)      
            else if (N_d > 0) then
                call makeMetropolisMove(newDownConfig)
            endif
            
            ! Apply the metropolis rejection
            oldProb = abs(manyBodyMethod(upConfig, downConfig))**2
            newProb = abs(manyBodyMethod(newUpConfig, newDownConfig))**2
            if ((oldProb == 0) .or. (newProb/oldProb > rand())) then
                upConfig = newUpConfig
                downConfig = newDownConfig
            endif

            ! Sample the current config if it's not an intitalization step
            if (i>initMetroSamples) then
                upConfigs(:,:,i-initMetroSamples) = upConfig
                downConfigs(:,:,i-initMetroSamples) = downConfig
            endif
        enddo
    end subroutine

    ! Make a metropolis move on the given electron configuration
    subroutine makeMetropolisMove(config)
    implicit none
        real(prec) :: config(:,:), r, theeta, phi, dx(3)
        integer    :: n

        r      = rand()*maxMetroJump
        theeta = rand()*pi
        phi    = rand()*2*pi        

        dx(1) = r*sin(theeta)*cos(phi)
        dx(2) = r*sin(theeta)*sin(phi)
        dx(3) = r*cos(theeta)

        n = int(rand()*size(config,2)+1)
        config(:,n) = config(:,n) + dx
    end subroutine

    ! Assign each electron to a basis sequentially
    ! kind of like the aufbau principle, but in the
    ! order which the basis states were given to the
    ! program
    subroutine sequentialCharacter()
    implicit none
    integer :: j

        if (size(basis)<upElectrons) print *, "Error: could not construct sequential character, too many up electrons!"
        if (size(basis)<downElectrons) print *, "Error: could not construct sequential character, too many down electrons!"

        allocate(upCharacters(size(basis),upElectrons))     !M_ij = jth electrons ith basis character
        allocate(downCharacters(size(basis),downElectrons))
        
        upCharacters = 0
        downCharacters = 0
        
        ! Put the jth up electron in the jth basis state (Aufbau-like)
        do j=1, upElectrons
            upCharacters(j,j) = 1
        enddo

        ! Same for down electrons
        do j=1, downElectrons
            downCharacters(j,j) = 1
        enddo        
    end subroutine

    ! Our electronic characters have been explicitly constructed
    subroutine explicitCharacter()
    implicit none
        ! Don't need to do anything!
    end subroutine

    ! The single particle electron state of the nth electron
    ! as described by the electron characters given
    function singleParticleState(n, characters, x) result(ret)
    implicit none
        integer    :: n, i
        logical    :: isUp
        real(prec) :: x(3), characters(:,:)
        complex(prec) :: ret

        ret = 0
        do i=1,size(basis)
            ret = ret + characters(i,n)*basis(i)%state%value(x-basis(i)%centre)
        enddo

    end function

    ! A simple exponential jastrow factor
    function jastrowFactor(upConfig, downConfig)
    implicit none
        real(prec) :: upConfig(:,:), downConfig(:,:), jastrowFactor, r
        integer    :: i, j, N_u, N_d

        if (size(upConfig,1) /= 3) print *, "Jastrow factor error: wrong spatial dimension!"
        if (size(downConfig,1) /= 3) print *, "Jastrow factor error: wrong spatial dimension!"
    
        N_u = size(upConfig,2)
        N_d = size(downConfig,2)

        jastrowFactor = 0
        
        ! Calculate same-spin jastrow contributions
        do i=1,N_u-1
            do j=i+1,N_u ! i,j <=> up spin pair
                r = norm2(upConfig(:,i) - upConfig(:,j))/angstrom
                jastrowFactor = jastrowFactor + r/(1+r)
            enddo 
        enddo
        do i=1,N_d-1
            do j=i+1,N_d ! i,j <=> down spin pair
                r = norm2(downConfig(:,i) - downConfig(:,j))/angstrom
                jastrowFactor = jastrowFactor + r/(1+r)
            enddo 
        enddo

        ! Calculate different-spin jastrow contributions
        do i=1,N_u
            do j=1,N_d
                r = norm2(upConfig(:,i) - downConfig(:,j))/angstrom
                jastrowFactor = jastrowFactor + 0.5*r/(1+r)
            enddo
        enddo

        jastrowFactor = exp(jastrowFactor)

    end function

    ! A slater determinant combined with a jastrow factor
    function slaterJastrow(upConfig, downConfig) result(ret)
    implicit none
        real(prec)        :: upConfig(:,:), downConfig(:,:)
        complex(prec)     :: ret
            ret = slaterDeterminant(upConfig, downConfig) * &
                  jastrowFactor(upConfig, downConfig)
    end function

    ! A slater determinant of two spin species seperates into the product
    ! of slater determinants of each species (no pauli exclusion between spins)
    function slaterDeterminant(upConfig, downConfig) result(ret)
    implicit none
        real(prec)        :: upConfig(:,:), downConfig(:,:)
        complex(prec)     :: ret
        ret = 1
        if (size(upConfig,2)>0) ret = ret * slaterDeterminantSpinless(upConfig, upCharacters)
        if (size(downConfig,2)>0) ret = ret * slaterDeterminantSpinless(downConfig, downCharacters)
    end function

    ! A slater determinant of the given single particle states, evaluated
    ! at the given electron positions
    function slaterDeterminantSpinless(config, characters) result(ret)
    implicit none
        real(prec)           :: config(:,:), characters(:,:)
        complex(prec)        :: ret, prod, elm
        integer              :: M, i, j
        integer, allocatable :: perm(:,:), sgn(:)

        M = size(config,2)

        if (size(config,1) /= 3) then
            print *, "Error in calculating a slater determinant: Spatial dimension incorrect!"
            ret = 0
            return
        endif

        if (M == 0) then
            print *, "Error in slaterDeterminantSpinless: no electrons!"
        endif

        ! Deal with the single electron case (i.e a 1x1 matrix)
        if (M == 1) then
            ret = singleParticleState(1, characters, config(:,1))
            return
        endif

        ! Fill our permutation matrix and parity array
        call permutations(M,perm,sgn)

        ! Calculate ret = sum_{permutations} sign(permutation) product_i s_{i,permutation(i)}
        ! where s is our slater matrix
        ret = 0
        do i=1,factorial(M) ! Sum over permutations
            prod = 1
            do j=1,M ! Product over entries of our slater matrix
                elm = singleParticleState(j, characters, config(:,perm(i,j)))
                prod = prod*elm
            enddo
            ret = ret + sgn(i)*prod
        enddo
        ret = ret/sqrt(real(M)) ! Normalization (unneccasary as so far everything else isn't normalized anyway but whatever)

    end function
    
    ! A subroutine that generates the permutaions of a and
    ! stores the result in the factorial(size(a)) x size(a)
    ! matrix p. Uses Heap's algorithm. As an extension it
    ! also stores the signs of the permutations in s
    subroutine permutations(n, p, s)
    implicit none
        integer :: n, i, temp, currentRow, sgn
        integer, allocatable :: a(:), c(:), p(:,:), s(:)
        real(prec) :: startT, endT

        call cpu_time(startT)
        
        allocate(a(n))
        allocate(c(n))
        allocate(p(factorial(n),n))
        allocate(s(factorial(n)))

        do i=1,size(a)
            c(i) = 1
            a(i) = i
        enddo

        p(1,:)     = a(:)
        s(1)       = 1
        currentRow = 2

        i = 1
        sgn = 1
        do while(i < n+1)
            if (c(i) < i) then
                if (modulo(i,2)/=0) then
                    ! Swap a(i), a(1)
                    temp = a(i)
                    a(i) = a(1)
                    a(1) = temp
                    sgn  = -sgn
                else
                    ! Swap a(c(i)), a(i)
                    temp = a(i)
                    a(i) = a(c(i))
                    a(c(i)) = temp
                    sgn  = -sgn
                endif
                
                ! Generated a permutation
                p(currentRow,:) = a(:)
                s(currentRow) = sgn
                currentRow = currentRow + 1        

                c(i) = c(i) + 1
                i = 1
            else
                c(i) = 1
                i = i + 1
            endif
        end do

        call cpu_time(endT)
        permutationCPUtime = permutationCPUtime + (endT-startT)
    
    end subroutine
    
    ! Return the factorial of n
    recursive function factorial(n) result(ret)
    implicit none
        integer :: n, ret

        if (n < 0) then
            print *, "Error: tried to calculate the factorial of a negative integer!"
            ret = 0
            return
        endif

        if (n < 2) then
            ret = 1
            return
        endif

        ret = n * factorial(n-1)

    end function

end module vqmc

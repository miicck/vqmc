! Version 0.3
! A program that carries out variational qunatum monte carlo
! for the hydrogen atom, we work in SI units.

! Module that carries out the monte-carlo integration
module vqmc
use constants
!use IFPORT
implicit none

    ! Parameters controlling our vqmc simulation
    integer,    parameter :: metroSamples = 100000     ! The number of metropolis electron configurations that are used
    real(prec), parameter :: minR = angstrom/10000     ! The minimum distance e's are alowed from nucleii
    real(prec), parameter :: maxMetroJump = 4*angstrom ! The maximum distance an electron can move in a metropolis trial move
    integer, parameter    :: metroInit = 1000          ! Number of metropolis steps to take and discard to remove dependance on initial position

    real(prec) :: permutationCPUtime = 0

    ! An object that can be used as a basis function
    type, abstract :: basisState
    contains
        procedure(basisWfn), deferred :: value
        procedure(printStateDebug), deferred :: printDebugInfo
    end type

    abstract interface

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

        ! A method of taking a single particle basis and
        ! electronic positions and returning
        ! a complex amplitude (slater, slater-jastrow etc...)
        function manyBodyWfn(upBasis, downBasis, upConfig, downConfig)
            import
            real(prec) :: upConfig(:,:), downConfig(:,:)
            complex(prec) :: manyBodyWfn
            class(basisState) :: upBasis(:), downBasis(:)
        end function
    end interface

contains

    ! Use monte carlo integration to get the energetics of the given wavefunction in the given potential
    subroutine monteCarloEnergetics(upBasis, downBasis, wavefunction, potential, numConfigs, energy)
    implicit none
        class(basisState)       :: upBasis(:), downBasis(:)
        procedure(manyBodyWfn)  :: wavefunction
        procedure(pot)          :: potential
        real(prec), allocatable :: upPos(:,:,:), downPos(:,:,:), allPos(:,:,:)
        integer                 :: i, j, i2, j2, numConfigs, N_u, N_d
        complex(prec)           :: energy
        real(prec)              :: dr(3)

        N_u = size(upBasis)   ! # up electrons   = # up states
        N_d = size(downBasis) ! # down electrons = # down states
        energy = 0

        if ((N_u == 0) .and. (N_d == 0)) then
            print *, "Error, no basis states given for monte carlo energetics!"
        endif

        ! Allocate/fill our arrays
        call metropolis(upBasis, downBasis, wavefunction, numConfigs, upPos, downPos) ! Fill array of up/down electron positions
        allocate(allPos(3, N_u + N_d, numConfigs))      ! Array of all electron positions
        if (N_d > 0) allPos(:,1:N_d,:) = downPos        ! Fill the down part of allPos
        if (N_u > 0) allPos(:,N_d+1:N_d+N_u,:) = upPos  ! Fill the up part of allPos

        ! Perform monte carlo integration
        do i=1,numConfigs

            ! Sum potentials for all electrons
            do j=1,N_u + N_d
                energy = energy + potential(allPos(:,j,i))
            enddo
            
            ! Add kinetic energy of up and down electrons
            energy = energy + localKinetic(upBasis, downBasis, wavefunction, upPos(:,:,i), downPos(:,:,i))

            ! Sum electron - electron coulomb interactions
            do i2=1, N_u + N_d - 1
                do j2=i2+1, N_u + N_d
                    dr = allPos(:,i2,i) - allPos(:,j2,i)
                    energy = energy + qElectron**2/(4*pi*epsNaught*norm2(dr))
                enddo
            enddo

        enddo

        energy = energy/real(numConfigs)

    end subroutine

    ! Calculate the local kinetic energy of wavefunction
    ! made from the given single particle basis states
    ! evaluated at the given electron configuration
    function localKinetic(upBasis, downBasis, wavefunction, upConfig, downConfig) result(ret)
    implicit none
        class(basisState)       :: upBasis(:), downBasis(:)
        procedure(manyBodyWfn)  :: wavefunction
        real(prec)              :: upConfig(:,:), downConfig(:,:)
        real(prec), parameter   :: eps = angstrom / 10E4
        real(prec), allocatable :: deltaUp(:,:), deltaDown(:,:)
        complex(prec)           :: laplacian, valueAtConfig, ret
        integer                 :: i, j

        if (size(upBasis)+size(downBasis)==0) Print *, "Error in localKinetic, no basis!"
        if (size(upConfig,1) /= 3) print *, "Error in localKinetic, wrong spatial dimension for up electrons!"
        if (size(downConfig,1) /= 3) print *, "Error in localKinetic, wrong spatial dimension for down electrons!"

        allocate(deltaUp(size(upConfig,1),size(upConfig,2)))
        allocate(deltaDown(size(downConfig,1),size(downConfig,2)))

        valueAtConfig = wavefunction(upBasis, downBasis, upConfig, downConfig)

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
                laplacian = laplacian + wavefunction(upBasis, downBasis, upConfig + deltaUp, downConfig) - &
                                        2*valueAtConfig + &
                                        wavefunction(upBasis, downBasis, upConfig - deltaUp, downConfig)         
                deltaUp(j,i) = 0   ! Reset
            enddo
        enddo

        do i=1,size(downConfig,2) ! i <=> down electrons                 
            do j=1,3 ! j <=> directions
                ! Evaluate d^2/dx_j^2 using finite differences and add it to the laplacian
                deltaDown(j,i) = eps ! Create a displacement of the ith electron in the jth direction
                laplacian = laplacian + wavefunction(upBasis, downBasis, upConfig, downConfig + deltaDown) - &
                                        2*valueAtConfig + &
                                        wavefunction(upBasis, downBasis, upConfig, downConfig - deltaDown)         
                deltaDown(j,i) = 0   ! Reset
            enddo
        enddo

        laplacian = laplacian / (eps**2)
        ret = -laplacian*(hbar**2)/(2*mElectron)
        ret = ret / valueAtConfig ! Calculating *local* kinetic energy

    end function

    ! Sample n configurations from the given wavefunction
    ! built out of the given basis
    subroutine metropolis(upBasis, downBasis, wavefunction, n, upConfigs, downConfigs)
    implicit none
        class(basisState)       :: upBasis(:), downBasis(:)
        procedure(manyBodyWfn)  :: wavefunction
        real(prec), allocatable :: upConfigs(:,:,:), downConfigs(:,:,:) ! # dim, # electrons, # samples
        real(prec), allocatable :: upConfig(:,:), newUpConfig(:,:), downConfig(:,:), newDownConfig(:,:)
        integer                 :: n, i, N_d, N_u
        real(prec)              :: oldProb, newProb

        ! Get electron counts
        N_u = size(upBasis)
        N_d = size(downBasis)
        if (N_u + N_d == 0) print *, "Error in metropolis algorithm, no basis!"

        allocate(upConfigs(3,N_u,n))
        allocate(upConfig(3,N_u))
        allocate(newUpConfig(3,N_u))

        allocate(downConfigs(3,N_d,n))
        allocate(downConfig(3,N_d))
        allocate(newDownConfig(3,N_d))

        ! Initial configuration is all electrons at the origin
        upConfig = 0
        downConfig = 0
        
        do i=1,n + metroInit

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
            oldProb = abs(wavefunction(upBasis, downBasis, upConfig, downConfig))**2
            newProb = abs(wavefunction(upBasis, downBasis, newUpConfig, newDownConfig))**2
            if ((oldProb == 0) .or. (newProb/oldProb > rand())) then
                upConfig = newUpConfig
                downConfig = newDownConfig
            endif

            ! Sample the current config if it's not an intitalization step
            if (i>metroInit) then
                upConfigs(:,:,i-metroInit) = upConfig
                downConfigs(:,:,i-metroInit) = downConfig
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
    function slaterJastrow(upBasis, downBasis, upConfig, downConfig) result(ret)
    implicit none
        class(basisState) :: upBasis(:), downBasis(:)
        real(prec)        :: upConfig(:,:), downConfig(:,:)
        complex(prec)     :: ret
            ret = slaterDeterminant(upBasis, downBasis, upConfig, downConfig) * &
                  jastrowFactor(upConfig, downConfig)
    end function

    ! A slater determinant of two spin species seperates into the product
    ! of slater determinants of each species (no pauli exclusion between spins)
    function slaterDeterminant(upBasis, downBasis, upConfig, downConfig) result(ret)
    implicit none
        class(basisState) :: upBasis(:), downBasis(:)
        real(prec)        :: upConfig(:,:), downConfig(:,:)
        complex(prec)     :: ret
        if (size(upBasis) + size(downBasis)==0) print *, "Error in slaterDeterminant, no basis!"
        ret = 1
        if (size(upBasis)>0) ret = ret * slaterDeterminantSpinless(upBasis, upConfig)
        if (size(downBasis)>0) ret = ret * slaterDeterminantSpinless(downBasis, downConfig)
    end function

    ! A slater determinant of the given single particle states, evaluated
    ! at the given electron positions
    function slaterDeterminantSpinless(singleParticleStates, electronPositions) result(ret)
    implicit none
        class(basisState)    :: singleParticleStates(:)
        real(prec)           :: electronPositions(:,:)
        complex(prec)        :: ret, prod, elm
        integer              :: M, i, j
        integer, allocatable :: perm(:,:), sgn(:)

        M = size(singleParticleStates)

        if (size(electronPositions,1) /= 3) then
            print *, "Error in calculating a slater determinant: Spatial dimension incorrect!"
            ret = 0
            return
        endif

        if (size(electronPositions,2) /= M) then
            print *, "Error in calculating a slater determinant: # Electrons /= # Basis States!"
        endif

        if (M == 0) then
            print *, "Error in slaterDeterminantSpinless: no basis states!"
        endif

        ! Deal with the single basis case (i.e a 1x1 matrix)
        if (M == 1) then
            ret = singleParticleStates(1)%value(electronPositions(:,1))
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
                elm = singleParticleStates(j)%value(electronPositions(:,perm(i,j)))
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

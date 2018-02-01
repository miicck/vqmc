program main
use vqmc
use atomicBasis
implicit none
    type(atomicState), allocatable :: basis(:), basis2(:)
    complex(prec) :: coefficients2(1)

    allocate(basis(3))
    basis(1) = atomicState(1,0,0,1) ! A hydrogenic 1s state (n=1,l=0,m=0,z=1)
    basis(2) = atomicState(2,0,0,1) ! A hydrogenic 2s state (n=2,l=0,m=0,z=1)
    basis(3) = atomicState(2,1,0,1) ! A hydrogenic 2p state (n=2,l=1,m=0,z=1)
    !basis(3) = atomicState(3,0,0,1) ! A hydrogenic 3s state (n=2,l=0,m=0,z=1)

    !call characterPlots(basis, potential) 
    !call optimizeBasis(basis, potential)

    allocate(basis2(1))
    basis2(1) = atomicState(5,3,0,1)
    coefficients2(1) = 1
    
    call sampleWavefunctionToFile(basis2, coefficients2)
    call debugAtomicState(basis2(1))

contains

    function potential(x)
    real(prec) :: potential, x(3)
        potential = -qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

end program main
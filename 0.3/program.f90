program main
use vqmc
use atomicBasis
implicit none
    type(atomicState), allocatable :: basis(:)
    complex(prec) :: coeff(2)

    allocate(basis(2))
    basis(1) = atomicState(1,0,0,1) ! A hydrogenic 1s state (n=1,l=0,m=0,z=1)
    basis(2) = atomicState(2,0,0,1) ! A hydrogenic 2s state (n=2,l=0,m=0,z=1)
    !basis(3) = atomicState(3,0,0,1) ! A hydrogenic 3s state (n=2,l=0,m=0,z=1)

    coeff(1) = 1
    coeff(2) = 0
    call sampleWavefunctionToFile(basis,coeff)

    call characterPlots(basis, potential) 
    !call optimizeBasis(basis, potential)
contains

    function potential(x)
    real(prec) :: potential, x(3)
        potential = -qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

end program main
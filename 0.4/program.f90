program main
use vqmc
use atomicBasis
implicit none
    type(atomicState), allocatable :: basis(:)
    complex(prec) :: energy

    allocate(basis(2))
    basis(1) = atomicState(1,0,0,1) ! A hydrogenic 1s state (n=1,l=0,m=0,z=1)
    basis(2) = atomicState(5,2,0,1) ! A hydrogenic 2s state (n=2,l=0,m=0,z=1)
    !basis(3) = atomicState(2,1,0,1) ! A hydrogenic 2p state (n=2,l=1,m=0,z=1)
    !basis(3) = atomicState(3,0,0,1) ! A hydrogenic 3s state (n=2,l=0,m=0,z=1)

    !call characterPlots(basis, potential)
    !call optimizeBasis(basis, potential)

    call energeticsSlater(basis, potential, 1000, energy)
    print *, energy/electronVolt
    !call testMetroSlater(basis,1000000)

contains

    function potential(x)
    real(prec) :: potential, x(3)
        if (norm2(x)==0) x = (/1,0,0/)*angstrom/10000
        potential = -qElectron**2/(4*pi*epsNaught*norm2(x))
    end function

end program main

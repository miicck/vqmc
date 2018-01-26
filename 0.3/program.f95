program main
use vqmc
!use particleInBox, only: wavefunction=>groundState, potential
use hydrogenicSystem, only: wavefunction=>groundState, potential 
use atomicBasis
implicit none
    real(prec) :: start, end, x(3)

    open(unit=2,file="toPlot")
    call debugRadialPart(2)

    open(unit=3,file="toPlot2D")
    call debugWavefunctions(3)

    ! Calculate the ground state energy using vqmc
    call cpu_time(start)
    print *, "Calculated system energy (eV):", energy(wavefunction, potential)/electronVolt
    call cpu_time(end)
    print *, "Monte-carlo itterations:", MC_ITTER  
    print *, "Elapsed time:           ", end-start

end program main
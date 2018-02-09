! Constants used by the program
module constants
    implicit none
    
        integer, parameter :: prec = selected_real_kind(15,307)
        real(prec), parameter :: pi = 3.14159265358979
        real(prec), parameter :: hbar = 1.054571800E-34
        real(prec), parameter :: mElectron = 9.10938356E-31
        real(prec), parameter :: qElectron = 1.6021766208E-19
        real(prec), parameter :: mProton = 1.672621898E-27
        real(prec), parameter :: epsNaught = 8.854187817620E-12
        real(prec), parameter :: electronVolt = qElectron
        real(prec), parameter :: hartree = 4.35974465054E-18
        real(prec), parameter :: angstrom = 1E-10
        real(prec), parameter :: a0 = 5.29177E-11 ! Bohr radius
    
end module
! A system with a nucleus of the given charge
! at the origin
module hydrogenicSystem
    use constants
    implicit none
    
        real(prec) :: nuclearCharge = 1
    
    contains
    
        ! hydrogen ground state
        ! Wavefunction (not normalized)
        function groundState(x)
        implicit none
            real(prec) :: x(3), groundState
            groundState = exp(-nuclearCharge*norm2(x)/a0)
        end function
    
        ! Hydrogen potential
        function potential(x)
        implicit none
            real(prec) :: x(3), potential
            potential = -nuclearCharge*qElectron**2/(4*pi*epsNaught*norm2(x))
        end function
        
    end module
    
    ! Things to do with the particle in a 3D box
    module particleInBox
    use constants
    implicit none
    
        real(prec), parameter :: boxSize = angstrom
        real(prec), parameter :: infiniteEnergy = huge(prec)
    
    contains
    
        ! The particle in a box ground state
        function groundState(x)
        implicit none
            real(prec) :: x(3), groundState
            if (outsideBounds(x(1))) then
                groundState = 0
            else if (outsideBounds(x(2))) then
                groundState = 0
            else if (outsideBounds(x(3))) then
                groundState = 0
            else
                groundState = sin(pi*x(1)/boxSize)*sin(pi*x(2)/boxSize)*sin(pi*x(3)/boxSize);
            endif
        end function
    
        ! The particle in a box potential
        function potential(x)
        implicit none
            real(prec) :: x(3), potential
            if (outsideBounds(x(1))) then
                potential = infiniteEnergy
            else if (outsideBounds(x(2))) then
                potential = infiniteEnergy
            else if (outsideBounds(x(3))) then
                potential = infiniteEnergy
            else
                potential = 0
            endif
        end function
    
        ! Returns true if the single coordinate x is outside the box
        function outsideBounds(x)
        implicit none
            real(prec) :: x
            logical :: outsideBounds
            outsideBounds = or(x<0,x>boxSize) 
        end function
    
    end module
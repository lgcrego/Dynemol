module MM_input

use constants_m
use MM_parms_module 

contains

!================================
 subroutine Define_MM_Environment
!================================
implicit none

!------------------------------------------------------------------------------
! NOT USED ANY MORE ...
!------------------------------------------------------------------------------

end subroutine Define_MM_Environment

!
!
!
!================================
 subroutine allocate_species( N )
!================================
implicit none
integer , intent(in)    :: N

! local variables ...
integer :: i

allocate( species ( N ) )

do i = 1 , N
    species(i) % my_species     = 0
    species(i) % N_of_atoms     = 0
    species(i) % N_of_molecules = 0
    species(i) % cm(3)          = 0.0d0
    species(i) % mass           = 0.0d0
    species(i) % flex           = .false.
    species(i) % residue        = "XXX"
    species(i) % nr             = 0
    species(i) % Nbonds         = 0
    species(i) % Nangs          = 0
    species(i) % Ndiheds        = 0
    species(i) % Nharm          = 0
    species(i) % Nbonds14       = 0
    species(i) % NintraIJ       = 0
end do

end subroutine allocate_species
!
!
!
end module MM_input


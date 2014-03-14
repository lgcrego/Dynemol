module MM_input

use MM_types , only : MM_molecular , MM_system

type(MM_system)                  :: MM
type(MM_molecular) , allocatable :: species(:) 

real*8  :: temperature, pressure, cutoff_radius, temperature_coupling, pressure_coupling, damping_Wolf
integer :: read_velocities, gmx_input_format, MM_log_step, MM_frame_step

logical , parameter :: T_ = .true. , F_ = .false. 

contains

!================================
 subroutine Define_MM_Environment
!================================
implicit none

!------------------------------------------------------------------------------
! SYSTEM  INFO
!
  MM % N_of_molecules = 1                   ! <== total number of molecules
  MM % N_of_species   = 1                   ! <== total number of species

  CALL allocate_species( MM % N_of_species )
!------------------------------------------------------------------------------
! repeat the following information filling for all the different species ...
!

  species(1) % residue         = "Ru2"      ! <== Residue label for species i ; character(len3)
  species(1) % N_of_molecules  = 1          ! <== Number of molecules of species i
  species(1) % N_of_atoms      = 61         ! <== Number of atoms comprosing a single molecule of species i
  species(1) % flex            = T_         ! <== Flexible : T_ , F_

!------------------------------------------------------------------------------
! ENVIRONMENT parameters ...
!
  temperature            = 300.d0           ! <== Bath Temperature (K)
  pressure               = 1.d0             ! <== Pressure
  temperature_coupling   = 1.d-1            ! <== Temperature coupling term with the bath
  pressure_coupling      = 2.d-1            ! <== Pressure coupling term 

  cutoff_radius          = 16.d0            ! <== Cut off radius (Angs) for electrostatic interactions
  damping_Wolf           = 6.d0             ! <== Wolf's method damping parameter (length^{-1}) ; (J. Chem. Phys. 1999; 110(17):8254)

!------------------------------------------------------------------------------
! GENERAL INFO ...
!
  read_velocities        = F_               ! <== reads the initial velocities : T_ , F_
  gmx_input_format       = T_               ! <== reads FF parameters from gmx input files : T_ , F_  

  MM_log_step            = 50               ! <== step for saving MM results & parameters
  MM_frame_step          = 01               ! <== step for saving MM results & parameters

! =====================================================================================

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
    species(i) % N_of_molecules = "XXX"
    species(i) % cm(:)          = 0.0d0
    species(i) % mass           = 0.0d0
    species(i) % flex           = .false.
    species(i) % residue        = "XXX"
    species(i) % nr             = 0
    species(i) % Nbonds         = 0
    species(i) % bonds          = 0.0d0
    species(i) % Nangs          = 0
    species(i) % Ndiheds        = 0
    species(i) % Nharm          = 0
    species(i) % Nbonds14       = 0
    species(i) % fact14         = 0.0d0
end do

end subroutine allocate_species
!
!
!
end module MM_input


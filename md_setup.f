module VV_Parent

    type, public  :: VV
        real*8    :: Kinetic
        real*8    :: Temperature
        real*8    :: Pressure
        real*8    :: Density
        character(len=:), allocatable :: thermostat_type
    contains
        procedure :: VV1
        procedure :: VV2
    end type  

contains

    subroutine VV1( me , dt )
        class(VV) , intent(inout) :: me
        real*8    , intent(in)    :: dt
    end subroutine

    subroutine VV2( me , dt )
        class(VV) , intent(inout) :: me
        real*8    , intent(in)    :: dt
    end subroutine

end module VV_Parent
!
!
!
!
module setup_m

    use parameters_m , only : PBC
    use constants_m  , only : imol, half
    use MD_read_m    , only : MM , atom , molecule

    implicit none

    public :: Molecular_CM , move_to_box_CM

contains
!
!
!
!========================
 subroutine Molecular_CM
!========================
    implicit none
   
   ! local variables ...
    integer              :: i, xyz, k, l 
    real*8               :: atom_mass 
    real*8, dimension(3) :: ref_xyz
    real*8, dimension(3) :: intra_xyz
    real*8, dimension(3) :: mass_weighted_xyz
    real*8, dimension(3) :: displacement
    logical              :: SmallMolecule

    if ( any(molecule% mass <= tiny(1.d0)) ) then
        error stop "Molecular_CM: zero total mass for at least one molecule"
    end if
    
   ! calculates the center of mass of molecule i ...

    l = 1
    do i = 1 , MM % N_of_molecules

         ref_xyz = atom(l)% xyz 
         ! atom(l) % mass = massa do primeiro atomo da molecula i
         mass_weighted_xyz = atom(l)% xyz * atom(l) % mass

         SmallMolecule = NINT(real(molecule(i)%N_of_atoms) / real(MM%N_of_atoms)) == 0

         do k = 1 , molecule(i) % N_of_atoms - 1
              intra_xyz = atom(l+k)% xyz
              ! wrapped displacement
              displacement = intra_xyz - ref_xyz
              do xyz = 1 , 3
                     if ( abs( displacement(xyz) ) > MM % box(xyz) * HALF .AND. SmallMolecule ) then
                        ! unwrap atom k of molecule i ...
                        intra_xyz(xyz) = intra_xyz(xyz) - sign( MM % box(xyz) , displacement(xyz) ) * PBC(xyz)
                     endif
              end do
              ! atom(l+k) % mass = massa do atomo k da molecula i
              atom_mass = atom(l+k) % mass
              mass_weighted_xyz = mass_weighted_xyz + atom_mass*intra_xyz
         end do  

         molecule(i)% cm = imol * mass_weighted_xyz / molecule(i)% mass
         l = l + molecule(i) % N_of_atoms

    end do

end subroutine Molecular_CM
!
!
!
!================================
 subroutine move_to_box_CM(frame)
!================================
    implicit none
    integer, intent(in) :: frame
    
    ! local varibales ...
    integer :: j
    real*8  :: atom_mass, masstot
    real*8, dimension(3) :: mass_times_position
    real*8, dimension(3) :: mass_times_velocity
    real*8, dimension(3) :: rcm, vcm
    
    ! determine the center of mass and center-of-mass velocity of the entire simulation box

    If( mod(frame,10) == 0 ) then
    
        mass_times_position = 0.d0
        mass_times_velocity = 0.d0
        masstot             = 0.d0
        
        do j = 1 , MM%N_of_atoms
        
            atom_mass = atom(j)%mass
        
            mass_times_position = mass_times_position + atom_mass * atom(j)%xyz
            mass_times_velocity = mass_times_velocity + atom_mass * atom(j)%vel
        
            masstot = masstot + atom_mass
        
        end do
        
        if (masstot <= tiny(1.d0)) then
            error stop "move_to_box_CM: zero total mass"
        end if
        
        rcm = mass_times_position / masstot
        vcm = mass_times_velocity / masstot
        
        ! transform atomic coordinates to CM frame ...
        ! subtract VCM from atomic velocities...
        do j = 1 , MM%N_of_atoms
            atom(j)%xyz = atom(j)%xyz - rcm
            atom(j)%vel = atom(j)%vel - vcm
        end do

    end if

end subroutine move_to_box_CM
!
!
!
end module setup_m

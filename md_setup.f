module VV_Parent

    type, public  :: VV
        integer   :: thermostat_type
        real*8    :: Kinetic
        real*8    :: Temperature
        real*8    :: Pressure
        real*8    :: Density
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

    use constants_m
    use MD_read_m   , only : MM , atom , molecule

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
    real*8               :: massa 
    real*8, dimension(3) :: t0, t, t1, dr
    logical              :: SmallMolecule

   ! calculates the center of mass of molecule i ...

    l = 1
    do i = 1 , MM % N_of_molecules

         t0(1:3) = atom(l) % xyz(1:3) 
         t(1:3)  = atom(l) % xyz(1:3) * atom(l) % mass ! massa do primeiro átomo da molécula i

         SmallMolecule = NINT(float(molecule(i)%N_of_atoms) / float(MM%N_of_atoms)) == 0

         do k = 1 , molecule(i) % N_of_atoms - 1
              t1(1:3) = atom(l+k) % xyz(1:3)
              dr(1:3) = t1(1:3) - t0(1:3)
              do xyz = 1 , 3
                     if ( abs( dr(xyz) ) > MM % box(xyz) * HALF .AND. SmallMolecule ) then
                        ! fix devided molecule i ...
                        t1(xyz) = t1(xyz) - sign( MM % box(xyz) , dr(xyz) ) 
                     endif
              end do
              massa = atom(l+k) % mass ! massa do átomo k da molécula i
              t(1:3) = t(1:3) + massa * t1(1:3)
         end do  

         molecule(i) % cm(1:3) = imol * t(1:3) / molecule(i) % mass
         l = l + molecule(i) % N_of_atoms

    end do

end subroutine Molecular_CM
!
!
!
!==========================
 subroutine move_to_box_CM
!==========================
    implicit none
    
    ! local varibales ...
    integer :: i, j, l
    real*8  :: massa, masstot
    real*8, dimension(3) :: t, p, rcm, vcm
    
    ! determines the atomic Center of Mass and its velocity for the box of atoms ... 
    p(:) = 0.d0
    t(:) = 0.d0
    masstot = 0.d0
    
    l = 1
    do i = 1 , MM % N_of_molecules 
        do j = l , l + molecule(i) % N_of_atoms - 1
           massa = atom(j) % mass
           p(:) = p(:) + massa * atom(j) % xyz(:)
           t(:) = t(:) + massa * atom(j) % vel(:)
           masstot = masstot + massa
        end do
        l = l + molecule(i) % N_of_atoms
    end do
    rcm(:) = p(:) / masstot  ! <== center of mass of the box of atoms
    vcm(:) = t(:) / masstot  ! <== velocity of the center of mass
    
    ! transform atomic coordinates and velocities to CM frame ...
    l = 1
    do i = 1 , MM % N_of_molecules
        do j = l , l + molecule(i) % N_of_atoms -1
           atom(j) % xyz(:) = atom(j) % xyz(:) - rcm(:)  ! <== atomic coordinates with the origin at CM
           atom(j) % vel(:) = atom(j) % vel(:) - vcm(:)  ! <== atomic velocities measured with respect with vcm
        end do
        l = l + molecule(i) % N_of_atoms
    end do

end subroutine move_to_box_CM
!
!
!
end module setup_m

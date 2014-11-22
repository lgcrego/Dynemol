module MM_ERG_class_m

    use type_m
    use constants_m
    use MD_read_m       , only : atom , MM 
    use MM_types        , only : MM_system , MM_atomic , debug_MM
    use MM_input        , only : driver_MM
    use MD_dump_m       , only : saving_MM_frame
    use F_intra_m       , only : ForceIntra, Pot_Intra                                     

    public :: CG_OPT

    private

    type :: CG_OPT
        integer                 :: N_of_Freedom
        real*8  , allocatable   :: p(:)
        character (len=11)      :: driver
        logical                 :: profiling = .false.
    contains
        procedure :: cost => energy
        procedure :: cost_variation => forces
        procedure :: output => dump_geometry
    end type CG_OPT


    interface CG_OPT
        module procedure  constructor
    end interface

    ! module variables ...

contains
!
!
!
!====================================
 function constructor( ) result( me )
!====================================
implicit none

type(CG_OPT) :: me 

!local variable ...
integer :: i 

me % driver = driver_MM

If( driver_MM == "MM_Optimize" .OR. driver_MM == "NormalModes" ) me % profiling = .true.

! number of degrees of freedom ...
me % N_of_Freedom = 3 * MM % N_of_atoms

allocate( me % p( me % N_of_Freedom ) , source = [(atom(:)%xyz(i) , i=1,3)] )

end function constructor
!
!
!
!
!=====================
 function Energy( me )
!=====================
implicit none
class(CG_OPT) , intent(inout)  :: me
real*8                         :: Energy

!local variables ...
integer :: i

do i = 1 , 3

    atom(:) % ftotal(i) = D_zero
    atom(:) % xyz(i)    = me % p( (i-1)*MM%N_of_atoms+1 : i*MM%N_of_atoms ) 

end do

CALL ForceIntra

Energy = Pot_Intra * mol * micro / MM % N_of_molecules

end function Energy
!
!
!
!=================================
 subroutine Forces( me , f_intra )
!=================================
implicit none
class(CG_OPT)   , intent(in)    :: me
real*8          , intent(out)   :: f_intra(:)

! local variables ...
integer :: i , GeneSize

f_intra = - [(atom(:)%ftotal(i) , i=1,3)] * 1.d10

end subroutine Forces
!
!
!
!====================================
 subroutine dump_geometry( me , iter)
!====================================
implicit none
class(CG_OPT)   , intent(in) :: me
integer         , optional , intent(in) :: iter

! local variables ...

if( iter == 0 ) call system( "rm frames-MM.pdb" )

call saving_MM_frame( iter , D_zero )

end subroutine dump_geometry
!
!
!
end module MM_ERG_class_m

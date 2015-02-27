module MM_ERG_class_m

    use type_m
    use constants_m
    use MD_read_m               , only : atom , MM 
    use MM_types                , only : MM_system , MM_atomic , debug_MM
    use MM_input                , only : driver_MM
    use MD_dump_m               , only : saving_MM_frame
    use F_intra_m               , only : ForceIntra, Pot_Intra                                     
    use F_inter_m               , only : ForceInter
    use for_force               , only : Pot_Total 
    use OPT_Parent_class_m      , only : OPT_Parent

    implicit none

    private 

    public :: MM_OPT

    type, extends(OPT_Parent)    :: MM_OPT
        integer                  :: ITMAX_MM = 4000             ! <== 100-300 is a good compromise of accuracy and safety
        real*8                   :: BracketSize_MM = 1.d-2      ! <== this value may vary between 1.0d-2 and 1.0d-3
        logical                  :: profiling_MM = .true.
    contains
        procedure :: cost => energy
        procedure :: cost_variation => forces
        procedure :: output => dump_geometry
    end type MM_OPT


    interface MM_OPT
        module procedure  constructor
    end interface

    ! module variables ...
    integer :: N_of_free

contains
!
!
!
!====================================
 function constructor( ) result( me )
!====================================
implicit none

type(MM_OPT) :: me 

!local variable ...
integer :: i , j

! Cause you are my kind / You're all that I want ...
me % ITMAX       = me % ITMAX_MM
me % BracketSize = me % BracketSize_MM
me % profiling   = me % profiling_MM

me % driver = driver_MM

If( driver_MM == "MM_Optimize" .OR. driver_MM == "NormalModes" ) me % profiling = .true.

! number of degrees of freedom allowed to relax ...
N_of_free = count( atom % flex )
me % N_of_Freedom = 3 * N_of_free

allocate( me % p( me % N_of_Freedom ) )
forall(i=1:3) me % p( (i-1)*N_of_free+1 : i*N_of_free ) = pack( atom(:)%xyz(i) , atom(:)%flex , me%p)

end function constructor
!
!
!
!
!=====================
 function Energy( me )
!=====================
implicit none
class(MM_OPT) , intent(inout)  :: me
real*8                         :: Energy

!local variables ...
integer :: i , j

do i = 1 , 3
    atom(:) % ftotal(i) = D_zero
    where( atom % flex ) atom(:) % xyz(i) = me % p( (i-1)*N_of_free+1 : i*N_of_free ) 
end do

If( MM % N_of_molecules == 1 ) then

    CALL ForceIntra
    Energy = Pot_Intra * mol * micro / MM % N_of_molecules

else

    CALL ForceInter
    CALL ForceIntra
    Energy = Pot_Total

end if

end function Energy
!
!
!
!================================
 subroutine Forces( me , vector )
!================================
implicit none
class(MM_OPT)   , intent(in)     :: me
real*8          , intent(inout)  :: vector(:)

! local variables ...
integer :: i , GeneSize

do i = 1 , 3 
    where( atom % flex ) vector( (i-1)*N_of_free+1 : i*N_of_free ) = - atom(:)%ftotal(i) * mts_2_Angs
end do

end subroutine Forces
!
!
!
!====================================
 subroutine dump_geometry( me , iter)
!====================================
implicit none
class(MM_OPT)   , intent(in) :: me
integer         , optional , intent(in) :: iter

! local variables ...

if( iter == 0 ) call system( "rm frames-MM.pdb" )

call saving_MM_frame( iter , D_zero )

end subroutine dump_geometry
!
!
!
end module MM_ERG_class_m

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
        character(len=72)        :: my_message
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
!===================================
 function constructor() result( me )
!===================================
use setup_m, only : setup
implicit none
type(MM_OPT) :: me 

!local variable ...
integer :: i 

! Cause you are my kind / You're all that I want ...
me % ITMAX       = me % ITMAX_MM
me % BracketSize = me % BracketSize_MM
me % profiling   = me % profiling_MM
me % accuracy    = mid_prec
me % driver      = driver_MM
me % message     = me % my_message

If( driver_MM == "Parametrize" ) me % profiling = .false.

! setup cutoff parameters for LJ and Coulomb interactions ...
CALL setup

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
integer :: i , j , k

k = 1
do i = 1 , 3
    do j = 1 , size(atom)
         atom(j) % ftotal(i) = D_zero
         If( atom(j)%flex ) then
              atom(j) % xyz(i) = me % p( k ) 
              k = k + 1
          end If     
    end do
end do

! this may not work for very large systems ...
!do i = 1 , 3
!    atom(:) % ftotal(i) = D_zero
!    where( atom % flex ) atom(:) % xyz(i) = me % p( (i-1)*N_of_free+1 : i*N_of_free ) 
!end do

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
integer :: i 
real*8  :: dumb


do i = 1 , 3 
    where( atom % flex ) vector( (i-1)*N_of_free+1 : i*N_of_free ) = - atom(:)%ftotal(i) * mts_2_Angs
end do

! to avoid compiler warnings ...
dumb = me%p(1)

end subroutine Forces
!
!
!
!====================================
 subroutine dump_geometry( me , iter)
!====================================
implicit none
class(MM_OPT)              , intent(in) :: me
integer         , optional , intent(in) :: iter

! local variables ...
real*8 :: dumb

if( iter == 0 ) call system( "rm frames-MM.pdb" )

call saving_MM_frame( iter , D_zero )

! to avoid compiler warnings ...
dumb = me%p(1)

end subroutine dump_geometry
!
!
!
end module MM_ERG_class_m

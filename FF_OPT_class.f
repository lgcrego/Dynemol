module FF_OPT_class_m

    use type_m
    use constants_m
    use MD_read_m               , only : atom , molecule , MM 
    use MM_types                , only : MM_system , MM_atomic , MM_molecular , debug_MM
    use MM_input                , only : driver_MM
!    use MD_dump_m               , only : saving_MM_frame
    use F_intra_m               , only : ForceIntra, Pot_Intra                                     
    use OPT_Parent_class_m      , only : OPT_Parent

    public :: FF_OPT

    private

    type, extends(OPT_Parent)   :: FF_OPT
        integer                 :: ITMAX_FF = 20            ! <== 100-300 is a good compromise of accuracy and safety
        real*8                  :: BracketSize_FF = 5.d-3   ! <== this value may vary between 1.0d-2 and 1.0d-3
        logical                 :: profiling_FF = .true.
    contains
        procedure :: cost => energy
        procedure :: cost_variation => Generalized_Forces
        procedure :: output => dump_geometry
    end type FF_OPT


    interface FF_OPT
        module procedure  constructor
    end interface

    ! module variables ...
    integer               :: bonds , angs , diheds
    logical , allocatable :: bonds_mask(:,:) , angs_mask(:,:) , diheds_mask(:,:)

contains
!
!
!
!====================================
 function constructor( ) result( me )
!====================================
implicit none

type(FF_OPT) :: me 

!local variable ...
integer :: i , j 

! Cause you are my kind / You're all that I want ...
me % ITMAX       = me % ITMAX_FF
me % BracketSize = me % BracketSize_FF
me % profiling   = me % profiling_FF

me % driver = driver_MM

If( driver_MM == "Parametrize" ) me % profiling = .true.

allocate( bonds_mask  ( size(molecule(1)%kbond0 (:,1)) , size(molecule(1)%kbond0 (1,:)) ) , source = (abs(molecule(1)%kbond0)  > low_prec) )
allocate( angs_mask   ( size(molecule(1)%kang0  (:,1)) , size(molecule(1)%kang0  (1,:)) ) , source = (abs(molecule(1)%kang0)   > low_prec) )
allocate( diheds_mask ( size(molecule(1)%kdihed0(:,1)) , size(molecule(1)%kdihed0(1,:)) ) , source = (abs(molecule(1)%kdihed0) > low_prec) )

! number of degrees of freedom in optimization space ...
bonds  = count( bonds_mask  ) 
angs   = count( angs_mask   ) 
diheds = count( diheds_mask )

me % N_of_Freedom = bonds + angs + diheds

allocate( me % p( me % N_of_Freedom ) , source = D_zero )

me % p( 1 :            ) = pack( molecule(1)%kbond0  , abs(molecule(1) % kbond0  ) > low_prec , me%p )

me % p( bonds+1 :      ) = pack( molecule(1)%kang0   , abs(molecule(1) % kang0   ) > low_prec , me%p ) 

me % p( bonds+angs+1 : ) = pack( molecule(1)%kdihed0 , abs(molecule(1) % kdihed0 ) > low_prec , me%p ) 

end function constructor
!
!
!
!=====================
 function Energy( me )
!=====================
implicit none
class(FF_OPT) , intent(inout)  :: me
real*8                         :: Energy

!local variables ...
integer :: i , j , i1, i2 , indx


!reset forces ...
forall( i=1:3 ) atom(:) % ftotal(i) = D_zero

molecule(1) % kbond0  = unpack( me%p(:bonds)             , bonds_mask  , molecule(1)%kbond0  )

molecule(1) % kang0   = unpack( me%p(bonds+1:bonds+angs) , angs_mask   , molecule(1)%kang0   )

molecule(1) % kdihed0 = unpack( me%p(bonds+angs+1:)      , diheds_mask , molecule(1)%kdihed0 )

CALL ForceIntra

Energy = Pot_Intra * mol * micro / MM % N_of_molecules

end function Energy
!
!
!
!
!============================================
 subroutine Generalized_Forces( me , vector )
!============================================
implicit none
class(FF_OPT)   , intent(in)    :: me
real*8          , intent(inout) :: vector(:)

! local parameters ...
real*8  , parameter :: small = 1.d-4

! local variables ...
integer         :: i 
type(FF_OPT)    :: before , after

do i = 1 , me % N_of_Freedom

    after  = me
    before = me

    after  % p(i) = me % p(i) * (D_one + small)
    before % p(i) = me % p(i) * (D_one - small)

    vector(i) = ( after%cost() - before%cost() ) / (two*small*me%p(i))


print*, i , vector(i)

end do
 
print*, " "

end subroutine Generalized_Forces
!
!
!
!====================================
 subroutine dump_geometry( me , iter)
!====================================
 implicit none
 class(FF_OPT)   , intent(in) :: me
 integer         , optional , intent(in) :: iter

 ! local variables ...

print*, "here"
 end subroutine dump_geometry
 !
 !
 !
end module FF_OPT_class_m

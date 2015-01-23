module FF_OPT_class_m

    use type_m
    use constants_m
    use MD_read_m               , only : atom , molecule , MM 
    use MM_types                , only : MM_system , MM_atomic , MM_molecular , debug_MM
    use MM_input                , only : driver_MM
    use F_intra_m               , only : ForceIntra, Pot_Intra                                     
    use OPT_Parent_class_m      , only : OPT_Parent

    public :: FF_OPT

    private

    type, extends(OPT_Parent)   :: FF_OPT
        integer                 :: ITMAX_FF = 100           ! <== 100-300 is a good compromise of accuracy and safety
        real*8                  :: BracketSize_FF = 1.d-4   ! <== this value may vary between 1.0d-3 and 1.0d-4
        logical                 :: profiling_FF = .FALSE.
    contains
        procedure :: cost => energy
        procedure :: cost_variation => Generalized_Forces
        procedure :: output => dump_parameters
    end type FF_OPT


    interface FF_OPT
        module procedure  constructor
    end interface

    ! module variables ...
    integer                                     :: bonds , angs , diheds
    integer             , allocatable           :: bonds_indx(:) , angs_indx(:) , diheds_indx(:)
    logical             , allocatable           :: bonds_mask(:,:) , angs_mask(:,:) , diheds_mask(:,:)
    real*8              , allocatable , target  :: bond_target(:,:) , ang_target(:,:) , dihed_target(:,:)
    type(real_pointer)  , allocatable           :: bond(:,:) , ang(:,:) , dihed(:,:)

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

! catch the distinct parameters and make (bond, ang, dihed) point to them ...
allocate(bond_target , source=Select_Different( molecule(1)%kbond0  , instance="bonds"  ) )
allocate(bond        , source=Define_Pointer  ( molecule(1)%kbond0  , bond_target       ) )

allocate(ang_target  , source=Select_Different( molecule(1)%kang0   , instance="angs"   ) )
allocate(ang         , source=Define_Pointer  ( molecule(1)%kang0   , ang_target        ) )

allocate(dihed_target, source=Select_Different( molecule(1)%kdihed0 , instance="diheds" ) )
allocate(dihed       , source=Define_Pointer  ( molecule(1)%kdihed0 , dihed_target      ) )

! the parameters = zero are not optimized ...
allocate( bonds_mask  ( size(bond_target (:,1)) , size(bond_target (1,:)) ) , source = (abs(bond_target)  > low_prec) )
allocate( angs_mask   ( size(ang_target  (:,1)) , size(ang_target  (1,:)) ) , source = (abs(ang_target)   > low_prec) )
allocate( diheds_mask ( size(dihed_target(:,1)) , size(dihed_target(1,:)) ) , source = (abs(dihed_target) > low_prec) )


bonds_mask(:,1) = .false.
bonds_mask(:,3) = .false.

angs_mask(:,1) = .false.

diheds_mask(:,:) = .false.


! number of degrees of freedom in optimization space ...
bonds  = count( bonds_mask  ) 
angs   = count( angs_mask   ) 
diheds = count( diheds_mask )

me % N_of_Freedom = bonds + angs + diheds

allocate( me % p( me % N_of_Freedom ) , source = D_zero )

me % p( 1 :            ) = pack( bond_target  , bonds_mask  , me%p )
me % p( bonds+1 :      ) = pack( ang_target   , angs_mask   , me%p ) 
me % p( bonds+angs+1 : ) = pack( dihed_target , diheds_mask , me%p ) 

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

! local variables ...
integer :: i , j  

! reset forces ...
forall( i=1:3 ) atom(:) % ftotal(i) = D_zero

bond_target  = unpack( me%p(:bonds)             , bonds_mask  , bond_target  )
ang_target   = unpack( me%p(bonds+1:bonds+angs) , angs_mask   , ang_target   )
dihed_target = unpack( me%p(bonds+angs+1:)      , diheds_mask , dihed_target )

forall( i=1:molecule(1)%nbonds  , j=1:size(molecule(1)%kbond0 (1,:)) ) molecule(1)%kbond0(i,j)  =  bond(i,j)%PTR

forall( i=1:molecule(1)%nangs   , j=1:size(molecule(1)%kang0  (1,:)) ) molecule(1)%kang0(i,j)   =  ang (i,j)%PTR

forall( i=1:molecule(1)%ndiheds , j=1:size(molecule(1)%kdihed0(1,:)) ) molecule(1)%kdihed0(i,j) =  dihed(i,j)%PTR

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

end do
 
end subroutine Generalized_Forces
!
!
!
!======================================
 subroutine dump_parameters( me , iter)
!======================================
implicit none
class(FF_OPT)   , intent(in) :: me
integer         , optional , intent(in) :: iter

! local variables ...
integer :: i , at1 , at2

do i = 1 , size(bonds_indx)
    at1 = molecule(1)%bonds(bonds_indx(i),1)
    at2 = molecule(1)%bonds(bonds_indx(i),2)

    write(*,'(2A4,3F15.4)') atom(at1)%MMSymbol                    , &
                            atom(at2)%MMSymbol                    , &
                            bond_target(i,2) / nano_2_angs        , &
                            bond_target(i,1) / ( factor2 *  imol) , &
                            bond_target(i,3)

end do




print*, "there"
end subroutine dump_parameters
!
!
!
!
!==========================================
 function Select_Different( a , instance ) result( b )
!==========================================
implicit none
real*8        , intent(in)  :: a(:,:)
character(*)  , intent(in)  :: instance
real*8        , allocatable :: b(:,:)

! local variables ...
integer                 :: j , k , diff , diff_size
integer , allocatable   ::  indx(:)
real*8  , allocatable   ::  tmp(:,:)

allocate( tmp  (size(a(:,1)),size(a(1,:))), source=D_zero )
allocate( indx (size(a(:,1))             ), source=I_zero )

! select different elements of a ...
diff_size = 1
indx(1)   = 1
tmp(1,:)  = a(1,:)

do j = 2 , size(a(:,1))
    
    do k = 1 , diff_size
        if( all(tmp(k,:) == a(j,:)) )  exit
    end do

    if( k > diff_size ) then
        tmp(k,:)  = a(j,:)
        indx(k)   = j
        diff_size = k
    end if

end do

allocate( b , source=tmp(1:diff_size,:) )

select case (instance)
    case( "bonds" )
        allocate( bonds_indx , source=indx(1:diff_size) )
    case( "angs" )
        allocate( angs_indx , source=indx(1:diff_size) )
    case( "diheds" )
        allocate( diheds_indx , source=indx(1:diff_size) )
end select

deallocate( tmp , indx )

end function Select_Different
!
!
!
!============================================
 function Define_Pointer( a , b ) result( c )
!============================================
implicit none
real*8 , intent(in)                 :: a(:,:)
real*8 , intent(in) , target        :: b(:,:)
type(real_pointer)  , allocatable   :: c(:,:) 

! local variables ...
integer :: i , j , k 

allocate( c(size(a(:,1)),size(a(1,:))) )

! creates a copy of a that points to b, i.e., c=a and c => b ...
do i = 1 , size(b(:,1))
    do j = 1  , size( a(:,1) )
        if( all(a(j,:) == b(i,:)) ) forall( k=1:size(a(1,:)) ) c(j,k)%PTR => b(i,k)
    end do
end do

end function Define_Pointer
!
!
!
end module FF_OPT_class_m

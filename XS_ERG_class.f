module XS_ERG_class_m

    use type_m             
    use constants_m
    use MD_read_m           , only : atom , MM 
    use MM_types            , only : MM_system , MM_atomic , debug_MM
    use MM_input            , only : driver_MM
    use MD_dump_m           , only : saving_MM_frame
    use F_intra_m           , only : ForceIntra, Pot_Intra                                     
    use F_inter_m           , only : ForceInter
    use for_force           , only : Pot_Total 
    use OPT_Parent_class_m  , only : OPT_Parent
    use parameters_m        , only : electron_state , hole_state
    use EH_parms_module     , only : HFP_Forces , verbose
    use QCModel_Huckel      , only : EigenSystem
    use Structure_Builder   , only : Unit_Cell ,               &
                                     system => Extended_Cell , &
                                     Generate_Structure ,      &
                                     Basis_Builder 
    use HuckelForces_m      , only : HuckelForces ,            &
                                     Force_QM => Force

    implicit none

    private 

    public :: XS_OPT

    type, extends(OPT_Parent)    :: XS_OPT
        integer                  :: ITMAX_XS = 40000             ! <== 100-300 is a good compromise of accuracy and safety
        real*8                   :: BracketSize_XS = 1.d-5       ! <== this value may vary between 1.0d-2 and 1.0d-3
        logical                  :: profiling_XS = .true.
        character(len=120)       :: my_message
    contains
        procedure :: cost => energy
        procedure :: cost_variation => forces
        procedure :: output => dump_geometry
    end type XS_OPT


    interface XS_OPT
        module procedure  constructor
    end interface

    ! module variables ...
    type(STO_basis) , allocatable  :: basis(:)
    type(R_eigen)                  :: UNI
    real*8          , allocatable  :: f_XS(:,:)
    integer                        :: N_of_free , eh_PES(2)

contains
!
!
!
!===================================
 function constructor() result( me )
!===================================
use setup_m, only : setup
implicit none
type(XS_OPT) :: me 

!local variable ...
integer :: i 

! Cause you are my kind / You're all that I want ...
me % ITMAX       = me % ITMAX_XS
me % BracketSize = me % BracketSize_XS
me % profiling   = me % profiling_XS
me % accuracy    = high_prec!mid_prec
me % driver      = driver_MM
me % message     = me % my_message

If( driver_MM == "Parametrize" ) me % profiling = .false.

open( unit=32, file='opt.trunk/XS_ERG.dat', status='unknown' )

! setup cutoff parameters for LJ and Coulomb interactions ...
CALL setup

!-----------------------------------------------------
! setup QM environment for XS calculations ...
 verbose    = .false.
 HFP_Forces = .true.
 CALL Generate_Structure(1)
 CALL Basis_Builder( system, basis )
 eh_PES = [ electron_state , hole_state ]
 allocate( f_XS( system% atoms , 3 ) , source = d_zero )
!-----------------------------------------------------

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
class(XS_OPT) , intent(inout)  :: me
real*8                         :: Energy

!local variables ...
integer :: i , j , k

k = 1
do i = 1 , 3
    do j = 1 , size(atom)
         atom(j) % ftotal(i) = D_zero
         If( atom(j)%flex ) &
         then
             atom(j) % xyz(i) = me % p( k ) 
             k = k + 1
         end If     
    end do
end do

If( MM % N_of_molecules == 1 ) &
then
    CALL ForceIntra
    Energy = Pot_Intra * mol * micro * kJmol_2_eV  !<== eV units
else
    CALL ForceInter
    CALL ForceIntra
    Energy = Pot_Total * kJmol_2_eV * MM% N_of_Molecules !<== eV units
end if

CALL EigenSystem( system , basis , UNI )
!----------------------------------
! approach, type character :
! HFP = Hellmann-Feynman-Pulay
! FDM = Finite Difference Method
!----------------------------------
CALL HuckelForces( system , basis , UNI , approach="HFP" , eh_PES=eh_PES )

Energy = Energy + ( UNI%erg(electron_state) - UNI%erg(hole_state) )
 
end function Energy
!
!
!
!================================
 subroutine Forces( me , vector )
!================================
implicit none
class(XS_OPT)   , intent(in)     :: me
real*8          , intent(inout)  :: vector(:)

! local parameters ...
real*8 , parameter :: Newton_2_eVAngs = 6.2415093433d8

! local variables ...
integer :: i , xyz , indx
real*8  :: dumb

do concurrent ( i=1:system% atoms , ( system% QMMM(i) == "QM" .AND. system% flex(i) == .true. ) )
     do xyz = 1 , 3
            indx = (i-1)*3 + xyz
            ! f_XS units = eV/Angs
            f_XS(i,xyz) = atom(i)% ftotal(xyz)*Newton_2_eVAngs + (Force_QM(indx,1) - Force_QM(indx,2)) 
     end do 
end do

do i = 1 , 3 
    where( atom % flex ) vector( (i-1)*N_of_free+1 : i*N_of_free ) = - f_XS(:,i)
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
class(XS_OPT)              , intent(in) :: me
integer         , optional , intent(in) :: iter

! local variables ...
real*8  :: dumb
logical :: exist

inquire(file="velocity_MM.pdb", EXIST=exist)

if( exist .and. iter == 0 ) call systemQQ( "rm frames.pdb" )

call saving_MM_frame( iter , D_zero )

! to avoid compiler warnings ...
dumb = me%p(1)

end subroutine dump_geometry
!
!
!
end module XS_ERG_class_m

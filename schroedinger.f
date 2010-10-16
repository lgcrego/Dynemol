 module Schroedinger_m

 use type_m
 use constants_m
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube , DP_Moment , initial_state
 use mkl95_precision
 use mkl95_blas
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use FMO_m                      , only : orbital
 use Multipole_Core             , only : Dipole_Moment
 use Data_Output                , only : Populations
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format

 public :: Huckel_dynamics , DeAllocate_QDyn

 private

 contains
!
! 
!===========================================================
 subroutine Huckel_dynamics(system, basis, UNI, FMO , QDyn )
!===========================================================
 implicit none
 type(structure) , intent(inout)    :: system
 type(STO_basis) , intent(in)       :: basis(:)
 type(C_eigen)   , intent(in)       :: UNI 
 type(C_eigen)   , intent(in)       :: FMO 
 type(f_time)    , intent(inout)    :: QDyn

! local variables ...
integer                             :: it , j
real*8                              :: t , t_rate
real*8                              :: Total_DP(3)
real*8          , ALLOCATABLE       :: Pops(:,:)
complex*16      , ALLOCATABLE       :: MO_bra(:,:)   , MO_ket(:,:)
complex*16      , ALLOCATABLE       :: AO_bra(:,:)   , AO_ket(:,:) 
complex*16      , ALLOCATABLE       :: DUAL_ket(:,:) , DUAL_bra(:,:) 
complex*16      , ALLOCATABLE       :: phase(:)      , bra(:)        , ket(:)

! preprocessing stuff ..........................................................

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 ) ) 

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(UNI%L(1,:))     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

MO_bra = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
MO_ket = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
!...............................................................................

!=========================================================
!                       Q-DYNAMICS  
!=========================================================
t = t_i              

t_rate = (t_f - t_i) / float(n_t)

DO it = 1 , n_t    

   phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

   If( t == t_i ) phase = C_one

   forall(j=1:n_part)   
      MO_bra(:,j) = conjg(phase(:)) * MO_bra(:,j) 
      MO_ket(:,j) =       phase(:)  * MO_ket(:,j) 
   end forall

!--------------------------------------------------------------------------
! . LOCAL representation for film STO production ...

! coefs of <k(t)| in AO basis 
   CALL gemm(UNI%L,MO_bra,AO_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in AO basis 
   CALL gemm(UNI%L,MO_ket,AO_ket,'T','N',C_one,C_zero)

   bra(:) = AO_bra(:,1)
   ket(:) = AO_ket(:,1)
 
   if ( GaussianCube ) CALL Gaussian_Cube_Format(bra,ket,it,t)

!--------------------------------------------------------------------------
! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
   DUAL_bra = AO_bra

! coefs of |k(t)> in DUAL basis ...
   CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

   Pops(it,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

   if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , bra , ket , Dual_ket(:,1) , Total_DP )

   t = t + t_rate

END DO

! sum population dynamics over frames ...
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops )

include 'formats.h'

end subroutine Huckel_dynamics
!
!
!
!=========================================
 subroutine DeAllocate_QDyn( QDyn , flag )
!=========================================
implicit none
type(f_time)  , intent(inout) :: QDyn
character(*)  , intent(in)    :: flag

! local variable ...
integer      :: i , N_of_fragments
character(1) :: first_in_line

select case( flag )

    case( "alloc" )

        if( allocated(trj) ) then

            CALL Coords_from_Universe( Unit_Cell, trj(2) )          ! <== use number 2 to avoid verbose
            CALL Generate_Structure( 2 )
            N_of_fragments = size( Extended_Cell%list_of_fragments ) 

        else

            CALL Generate_Structure( 2 )                            ! <== use number 2 to avoid verbose
            N_of_fragments = size( Extended_Cell%list_of_fragments )

        end if

        ! for the sake of having the donor survival probability in the first column at output ...
        first_in_line = Extended_Cell%list_of_fragments(1)
        where( Extended_Cell%list_of_fragments == "D" ) Extended_Cell%list_of_fragments = first_in_line
        Extended_Cell%list_of_fragments(1) = "D"

        ! QDyn%dyn = ( time ; fragments ; all fragments ) ...
        allocate( QDyn%fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )
        allocate( QDyn%dyn      ( n_t , 0:N_of_fragments+1 )                , source = 0.d0                              )

        ! cleaning the mess ...
        CALL DeAllocate_Structures( Extended_Cell )

    case( "dealloc" )

        deallocate( QDyn%dyn , QDyn%fragments )

end select

end subroutine DeAllocate_QDyn
!
!
end module Schroedinger_m

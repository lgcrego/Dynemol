 module Schroedinger_m

 use type_m
 use constants_m
 use mkl95_precision
 use mkl95_blas
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use FMO_m                      , only : orbital
 use Data_Output                , only : get_Populations
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format

 public :: Huckel_dynamics

 contains
!
! 
!===========================================================
 subroutine Huckel_dynamics(system, basis, UNI, FMO , QDyn )
!===========================================================
 implicit real*8      (a-h,o-y)
 implicit complex*16  (z)

 type(structure) , intent(in)               :: system
 type(STO_basis) , intent(in)               :: basis(:)
 type(C_eigen)   , intent(in)               :: UNI 
 type(C_eigen)   , intent(in)               :: FMO 
 real*8          , intent(inout)            :: QDyn(:,:)

! local variables ...
real*8     , ALLOCATABLE :: Pops(:,:)
complex*16 , ALLOCATABLE :: zG_L(:,:)     , MO_bra(:,:) 
complex*16 , ALLOCATABLE :: zG_R(:,:)     , MO_ket(:,:)
complex*16 , ALLOCATABLE :: AO_bra(:,:)   , AO_ket(:,:) 
complex*16 , ALLOCATABLE :: DUAL_ket(:,:) , DUAL_bra(:,:) 
complex*16 , ALLOCATABLE :: phase(:)      , bra(:)        , ket(:)

complex*16 , parameter   :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)


! preprocessing stuff ..........................................................

allocate( Pops( n_t , size(system%list_of_fragments)+1 ) ) 

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(UNI%L(1,:))     ,      &
                         zG_L     , zG_R     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

zG_L = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
zG_R = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
!...............................................................................

!=========================================================
!                       Q-DYNAMICS  
!=========================================================
t = t_i              

t_rate = (t_f - t_i) / float(n_t)

DO it = 1 , n_t    

   phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

   If( t == t_i ) phase = one

   forall(j=1:n_part)   
      MO_bra(:,j) = conjg(phase(:)) * zG_L(:,j) 
      MO_ket(:,j) =       phase(:)  * zG_R(:,j) 
   end forall

!--------------------------------------------------------------------------
! . LOCAL representation for film STO production ...

! coefs of <k(t)| in AO basis 
   CALL gemm(UNI%L,MO_bra,AO_bra,'T','N',one,zero)

! coefs of |k(t)> in AO basis 
   CALL gemm(UNI%L,MO_ket,AO_ket,'T','N',one,zero)

   bra(:) = AO_bra(:,1)
   ket(:) = AO_ket(:,1)
   
   if( GaussianCube ) CALL Gaussian_Cube_Format(bra,ket,it,t)

!--------------------------------------------------------------------------
! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
   CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',one,zero)

! coefs of |k(t)> in DUAL basis ...
   CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',one,zero)

    
   Pops(it,:) = get_Populations(system,basis,DUAL_bra,DUAL_ket)

   t = t + t_rate

   zG_L = MO_bra       ! <== updating expansion coefficients at t 
   zG_R = MO_ket       ! <== updating expansion coefficients at t

END DO

! sum population dynamics over frames ...
QDyn = QDyn + Pops

deallocate( Pops )

include 'formats.h'

end subroutine Huckel_dynamics
!
!
!
!==========================================================
 subroutine DeAllocate_QDyn( QDyn , QDyn_fragments , flag )
!==========================================================
implicit none
real*8          , optional      , allocatable   , intent(inout) :: QDyn(:,:)
character(1)    , optional      , allocatable   , intent(inout) :: QDyn_fragments(:)
character(*)                                    , intent(in)    :: flag

! local variable ...
integer :: i , N_of_fragments

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

        If( (survival) .AND. (Extended_Cell%list_of_fragments(1) /= "D") ) pause ">>> list_of_fragments(1) /= 'D' <<<"

        allocate( QDyn_fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )
        allocate( QDyn          (n_t,N_of_fragments+1)                      , source = 0.d0                              )

        ! cleaning the mess ...
        CALL DeAllocate_Structures( Extended_Cell )

    case( "dealloc" )

        deallocate(QDyn , QDyn_fragments)

end select

end subroutine DeAllocate_QDyn
!
!
end module Schroedinger_m

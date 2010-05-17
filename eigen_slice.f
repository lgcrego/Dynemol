module Eigen_Slice_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Allocation_m                , only : Allocate_Brackets
    use QCModel_Huckel              , only : EigenSystem
    use FMO_m                       , only : FMO_analysis ,                 &
                                             orbital
    use Data_Output                 , only : Populations
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format


    public :: Huckel_Slice_Dynamics , Preprocess_Huckel_Slice

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE        , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Complex*16 , ALLOCATABLE , save , dimension(:)   :: phase 
    Complex*16 , ALLOCATABLE        , dimension(:)   :: bra , ket

    integer     , save  :: it = 1
    real*8      , save  :: t , t_rate =  MD_dt * frame_step

contains
!
!
!
!===========================================================
 subroutine preprocess_Huckel_Slice( system , basis , Qdyn , zG_L , zG_R )
!===========================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
type(f_time)    , intent(inout) :: QDyn
complex*16      , ALLOCATABLE   , intent(inout) :: zG_L(:,:) , zG_R(:,:)

!local variables ...
type(C_eigen) :: UNI , FMO

! prepare  DONOR  state ...
CALL EigenSystem  ( system , basis , UNI         )
CALL FMO_analysis ( system , basis , UNI%R , FMO )

! allocate work vectors ...
CALL Allocate_Brackets( size(basis) , zG_L , zG_R , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

zG_L = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
zG_R = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 

Print 56 , initial_state     ! <== initial state of the isolated molecule 

! DUAL representation for efficient calculation of survival probabilities ...
! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,ZG_L,DUAL_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,ZG_R,DUAL_ket,'N','N',C_one,C_zero)

QDyn%dyn(it,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

it = it + 1

deallocate( UNI%L , UNI%R , UNI%erg )

include 'formats.h'

end subroutine preprocess_Huckel_Slice
!
!
!
!
!=============================================================
 subroutine Huckel_Slice_Dynamics( system , basis , QDyn , zG_L , zG_R , t )
!=============================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
type(f_time)    , intent(inout) :: QDyn
complex*16      , ALLOCATABLE   , intent(inout) :: zG_L(:,:) , zG_R(:,:)
real*8          , intent(out)   :: t

! local variables ...
integer       :: j
type(C_eigen) :: UNI 


CALL EigenSystem( system , basis, UNI )

phase(:) = exp(- zi * UNI%erg(:) * t_rate / h_bar)

! peak from previous phase ...
forall( j=1:n_part )   
    MO_bra(:,j) = conjg(phase(:)) * ZG_L(:,j) 
    MO_ket(:,j) =       phase(:)  * ZG_R(:,j) 
end forall

! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

t = (it-1)*t_rate

QDyn%dyn(it,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )
print*, t , Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )
stop


ZG_L(:,n_part) = MO_bra(:,n_part)       ! <== updating expansion coefficients at t 
ZG_R(:,n_part) = MO_ket(:,n_part)       ! <== updating expansion coefficients at t


deallocate( UNI%L , UNI%R , UNI%erg )

include 'formats.h'

end subroutine Huckel_Slice_Dynamics
!
!
!
end module Eigen_Slice_m

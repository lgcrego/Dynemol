 module Schroedinger_m

 use constants_m
 use mkl95_precision
 use mkl95_blas
 use Allocation_m
 use FMO_m
 use Data_Output
 use Psi_Squared_Cube_Format

 public :: Huckel_dynamics

 contains
!
! 
!==================================================
 subroutine Huckel_dynamics(system, basis, UNI, FMO)
!==================================================
 implicit real*8      (a-h,o-y)
 implicit complex*16  (z)

 type(structure) , intent(in)               :: system
 type(STO_basis) , intent(in)               :: basis(:)
 type(eigen)     , intent(in)               :: UNI 
 type(eigen)     , intent(in)               :: FMO 

! . local variables
 complex*16 , ALLOCATABLE :: zG_L(:,:)     , MO_bra(:,:) 
 complex*16 , ALLOCATABLE :: zG_R(:,:)     , MO_ket(:,:)
 complex*16 , ALLOCATABLE :: AO_bra(:,:)   , AO_ket(:,:) 
 complex*16 , ALLOCATABLE :: DUAL_ket(:,:) , DUAL_bra(:,:) 
 complex*16 , ALLOCATABLE :: phase(:)      , bra(:)        , ket(:)

 complex*16 , parameter   :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)


 OPEN(unit=26,file='survival.dat',status='unknown')   
     
!========================================================
!     define the initial states
!========================================================

 Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
 CALL Allocate_Brackets( size(UNI%L(1,:))    ,      &
                         zG_L     , zG_R     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

 zG_L = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
 zG_R = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 

!=========================================================
!             DYNAMIC  QUANTITIES
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
! . DUAL representation for efficient calculation of survival probabilities ...

! . coefs of <k(t)| in DUAL basis
   CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',one,zero)

! . coefs of |k(t)> in DUAL basis 
   CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',one,zero)

   if( Survival ) CALL Dump_Populations(system,basis,DUAL_bra,DUAL_ket,t)

   t = t + t_rate

   zG_L = MO_bra       ! <== updating expansion coefficients at t 
   zG_R = MO_ket       ! <== updating expansion coefficients at t

 END DO

!=============================================== 

 CLOSE(26)

 include 'formats.h'

 end subroutine Huckel_dynamics

 end module Schroedinger_m

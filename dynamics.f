 module dynamics_m

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
!==============================================================================
 subroutine Huckel_dynamics(system, basis, zL, zR, FMO_L, FMO_R, erg)
!==============================================================================
 implicit real*8      (a-h,o-y)
 implicit complex*16  (z)

 type(structure) , intent(in)               :: system
 type(STO_basis) , intent(in)               :: basis(:)
 complex*16      , intent(in) , allocatable :: zL(:,:)
 complex*16      , intent(in) , allocatable :: zR(:,:)
 complex*16      , intent(in) , allocatable :: FMO_L(:,:) 
 complex*16      , intent(in) , allocatable :: FMO_R(:,:) 
 real*8          , intent(in) , allocatable :: erg(:) 

! . local variables
 complex*16 , ALLOCATABLE :: zG_L(:,:)     , zGtL(:,:) 
 complex*16 , ALLOCATABLE :: zG_R(:,:)     , zGtR(:,:)
 complex*16 , ALLOCATABLE :: AO_bra(:,:)   , AO_ket(:,:) 
 complex*16 , ALLOCATABLE :: DUAL_ket(:,:) , DUAL_bra(:,:) 
 complex*16 , ALLOCATABLE :: phase(:)      , bra(:)        , ket(:)

 complex*16 , parameter   :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)


 OPEN(unit=26,file='survival.dat',status='unknown')   
     
!========================================================
!     define the initial states
!========================================================

 Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
 CALL Allocate_Brackets( size(zL(1,:))       ,      &
                         zG_L     , zG_R     ,      &
                         zGtL     , zGtR     ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

 zG_L = FMO_L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
 zG_R = FMO_R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 

!=========================================================
!             DYNAMIC  QUANTITIES
!=========================================================
 t = t_i              

 t_rate = (t_f - t_i) / float(n_t)
  
 DO it = 1 , n_t    

   phase(:) = cdexp(- zi * erg(:) * t_rate / h_bar)

   If( t == t_i ) phase = one

   forall(j=1:n_part)   
      zGtL(:,j) = conjg(phase(:)) * zG_L(:,j) 
      zGtR(:,j) =       phase(:)  * zG_R(:,j) 
   end forall

! -----------  LOCAL representation  ----------------- c

! coefs of <k(t)| in AO basis 
   CALL gemm(zL,zGtL,AO_bra,'T','N',one,zero)

! coefs of |k(t)> in AO basis 
   CALL gemm(zL,zGtR,AO_ket,'T','N',one,zero)

   bra(:) = AO_bra(:,1)
   ket(:) = AO_ket(:,1)
   
   if( GaussianCube ) CALL Gaussian_Cube_Format(bra,ket,it,t)

! -----------  DUAL representation  ----------------- c

! coefs of <k(t)| in DUAL basis
   CALL gemm(zL,zGtL,DUAL_bra,'T','N',one,zero)

! coefs of |k(t)> in DUAL basis 
   CALL gemm(zR,zGtR,DUAL_ket,'N','N',one,zero)

   if( Survival ) CALL Dump_Populations(system,basis,DUAL_bra,DUAL_ket,t)

   t = t + t_rate

   zG_L = zGtL       ! <== updating expansion coefficients at t 
   zG_R = zGtR       ! <== updating expansion coefficients at t

 END DO

!=============================================== 

 CLOSE(26)

 include 'formats.h'

 end subroutine Huckel_dynamics

 end module dynamics_m

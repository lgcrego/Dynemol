 Program MACH95

 use type_m
 use constants_m
 use mkl95_precision
 use mkl95_blas
 use Allocation_m
 use FMO_m
 use QCModel_Huckel
 use EHT_parameters
 use projectors
 use DOS_m
 use Structure_Builder
 use Data_Output
 use Psi_Squared_Cube_Format
 use Multipole_Core
 use Oscillator_m


 implicit real*8      (a-h,o-y)
 implicit complex*16  (z)

 complex*16 , ALLOCATABLE :: zG_L(:,:)     , FMO_L(:,:)    , zGtL(:,:) 
 complex*16 , ALLOCATABLE :: zG_R(:,:)     , FMO_R(:,:)    , zGtR(:,:)
 complex*16 , ALLOCATABLE :: AO_bra(:,:)   , AO_ket(:,:) 
 complex*16 , ALLOCATABLE :: DUAL_ket(:,:) , DUAL_bra(:,:) 
 complex*16 , ALLOCATABLE :: zL(:,:)       , zR(:,:)
 complex*16 , ALLOCATABLE :: phase(:)      , bra(:)        , ket(:)
 real*8     , ALLOCATABLE :: erg(:)        , erg_FMO(:) 
 
 complex*16 , parameter   :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

 OPEN(unit=26,file='survival.dat',status='unknown')   
     
!========================================================
!     starting up 
!========================================================
 
 CALL read_EHT_parameters

 CALL Read_Structure

 CALL Generate_Structure(t_i)

 CALL Basis_Builder(Extended_Cell,ExCell_basis)

 CALL eigen(Extended_Cell,ExCell_basis,zL,zR,erg)

 CALL TDOS(erg)

 CALL PDOS(extended_cell,zL,zR,erg)

 CALL FMO_analysis(Extended_Cell,zR,FMO_L,FMO_R,erg_FMO)

! CALL Dipole_Matrix(Extended_Cell,ExCell_basis,zL,zR)  

! CALL Optical_Transitions(Extended_Cell, ExCell_basis, M_matrix, zL, zR, erg)

!========================================================
!     define the initial states
!========================================================

 Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
 CALL Allocate_Brackets( size(ExCell_basis)  ,      &
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

   if( Survival ) CALL Dump_Populations(extended_cell,DUAL_bra,DUAL_ket,t)

   t = t + t_rate

   zG_L = zGtL       ! <== updating expansion coefficients at t 
   zG_R = zGtR       ! <== updating expansion coefficients at t

 END DO

!=============================================== 

 CLOSE(26)

 include 'formats.h'

 END

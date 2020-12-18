module ElHl_Chebyshev_m

    use type_m              , g_time => f_time  
    use blas95
    use lapack95
    use constants_m
    use ifport
    use parameters_m        , only : t_i , frame_step , Coulomb_ , EnvField_ , n_part, driver , QMMM , CT_dump_step , HFP_Forces
    use Structure_Builder   , only : Unit_Cell 
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis , eh_tag 
    use Data_Output         , only : Populations
    use Hamiltonians        , only : X_ij , even_more_extended_Huckel
    use Taylor_m            , only : Propagation, dump_Qdyn
    use Ehrenfest_Builder   , only : store_Hprime
    use Matrix_Math

    public  :: ElHl_Chebyshev , preprocess_ElHl_Chebyshev 

    private

! module parameters ...
    integer       , parameter   :: order       = 25
    real*8        , parameter   :: error       = 1.0d-12
    real*8        , parameter   :: norm_error  = 1.0d-12

! module variables ...
    logical       , save        :: necessary_  = .true.
    logical       , save        :: first_call_ = .true.
    real*8        , save        :: save_tau(2)
    real*8        , pointer     :: h(:,:)
    real*8, target, allocatable :: h0(:,:)
    real*8        , allocatable :: S_matrix(:,:) , H_prime(:,:)
    complex*16    , allocatable :: Psi_t_bra(:,:) , Psi_t_ket(:,:)

    interface preprocess_ElHl_Chebyshev
        module procedure preprocess_ElHl_Chebyshev
        module procedure preprocess_from_restart
    end interface
    
#ifdef USE_GPU
#  define _electron_  -1
#  define _hole_      +1
#endif

contains
!
!
!==========================================================================================================
 subroutine preprocess_ElHl_Chebyshev( system , basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , it )
! used for normal start (see interface) ...
!==========================================================================================================
implicit none
type(structure) , intent(inout) :: system
type(STO_basis) , intent(inout) :: basis(:)
complex*16      , intent(out)   :: AO_bra(:,:)
complex*16      , intent(out)   :: AO_ket(:,:)
complex*16      , intent(out)   :: DUAL_bra(:,:) 
complex*16      , intent(out)   :: DUAL_ket(:,:) 
type(g_time)    , intent(inout) :: QDyn
integer         , intent(in)    :: it

!local variables ...
integer                         :: N
type(R_eigen)   :: AO

N = size(basis)

#ifdef USE_GPU
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    allocate( S_matrix(N,N) )
    call GPU_Pin(S_matrix, N*N*8)
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#endif

! MUST compute S_matrix before FMO analysis ...
CALL Overlap_Matrix( system , basis , S_matrix )

allocate( h0(N,N) , source = D_zero )

If( EnvField_ ) then
    h0(:,:) = even_more_extended_Huckel( system , basis , S_matrix )
else
    h0(:,:) = Build_Huckel( basis , S_matrix )
end If

! for a rigid structure once is enough ...
If( driver == 'q_dynamics' ) necessary_ = .false.

!========================================================================
! prepare electron state ...
  CALL FMO_analysis( system , basis, AO=AO , instance="E" )

  AO_bra(:,1) = dcmplx(AO%L(:,1))
  AO_ket(:,1) = dcmplx(AO%R(:,1))
  deallocate( AO%L , AO%R )

! prepare hole state ...
  CALL FMO_analysis( system , basis, AO=AO , instance="H" )
  
  AO_bra(:,2) = dcmplx(AO%L(:,2))
  AO_ket(:,2) = dcmplx(AO%R(:,2))
  deallocate( AO%L , AO%R  )
!========================================================================

!==============================================
! prepare DUAL basis for local properties ...
! DUAL_bra = (C*)^T    ;    DUAL_ket = S*C ...
  DUAL_bra = AO_bra
  call op_x_ket( DUAL_ket, S_matrix , AO_ket )
!==============================================

!==============================================
! vector states to be propagated ...
! Psi_bra = C^T*S       ;      Psi_ket = C ...
  allocate( Psi_t_bra(N,n_part) )
  allocate( Psi_t_ket(N,n_part) )
  call bra_x_op( Psi_t_bra, AO_bra , S_matrix ) 
  Psi_t_ket = AO_ket
!==============================================

!==============================================
! preprocess stuff for EhrenfestForce ...
  CALL QuasiParticleEnergies(AO_bra, AO_ket, H0)

  CALL syInvert( S_matrix, return_full ) ! <== S_matrix content is destroyed and S_inv is returned
#define S_inv S_matrix

  allocate( H(N,N) )
  h => h0

! allocate and compute H' = S_inv * H ...
  allocate( H_prime(N,N) , source = D_zero )
  CALL syMultiply( S_inv , h , H_prime )
 
  CALL store_Hprime( N , H_prime )
!==============================================

! save populations(time=t_i) ...
QDyn%dyn(it,:,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )
CALL dump_Qdyn( Qdyn , it )

!==============================================
! clean and exit ...
If( driver /= 'q_dynamics' ) then
    deallocate( h0 )
    nullify( h )
end if
#ifndef USE_GPU
deallocate( H_prime )
#endif
#undef S_inv
!==============================================

end subroutine preprocess_ElHl_Chebyshev
!
!
!
!=============================================================================================================
 subroutine ElHl_Chebyshev( system , basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , t , delta_t , it )
!=============================================================================================================
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)
complex*16       , intent(inout) :: AO_bra(:,:)
complex*16       , intent(inout) :: AO_ket(:,:)
complex*16       , intent(inout) :: Dual_bra(:,:)
complex*16       , intent(inout) :: Dual_ket(:,:)
type(g_time)     , intent(inout) :: QDyn
real*8           , intent(inout) :: t
real*8           , intent(in)    :: delta_t
integer          , intent(in)    :: it

! local variables...
integer :: j , N , it_sync
real*8  :: t_max , tau_max , tau(2) , t_init

t_init = t

! max time inside slice ...
t_max = delta_t*frame_step*(it-1)  
! constants of evolution ...
tau_max = delta_t / h_bar

! trying to adapt time step for efficient propagation ...
tau(:) = merge( tau_max , save_tau(:) * 1.15d0 , first_call_ )
! but tau should be never bigger than tau_max ...
tau(:) = merge( tau_max , tau(:) , tau(:) > tau_max )

N = size(basis)

if(first_call_) then           ! allocate matrices
#ifdef USE_GPU
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    allocate( H(N,N) )         ! no need of H_prime in the cpu: will be calculated in the gpu
    call GPU_Pin( H, N*N*8 )
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#else
    allocate( H(N,N) , H_prime(N,N) )
#endif
end if

If ( necessary_ ) then ! <== not necessary for a rigid structures ...

    ! compute S and S_inverse ...
    CALL Overlap_Matrix( system , basis , S_matrix )

    allocate( h0(N,N) , source = D_zero )
    
    If( EnvField_ ) then
        it_sync = it-1  ! <== for synchronizing EnvSetUp call therein ...
        h0(:,:) = even_more_extended_Huckel( system , basis , S_matrix , it_sync)
    else
        h0(:,:) = Build_Huckel( basis , S_matrix )
    end If

    CALL syInvert( S_matrix, return_full )   ! S_matrix content is destroyed and S_inv is returned
#define S_inv S_matrix

end If

!=======================================================================
!           Electron Hamiltonian : upper triangle ...
!=======================================================================

h => h0

#ifndef USE_GPU
! allocate and compute H' = S_inv * H ...
If (.not. allocated(H_prime) ) allocate( H_prime(N,N) )
call syMultiply( S_inv , h , H_prime )
#endif

If( eh_tag(1) == "el" ) then  ! <== proceed evolution of ELECTRON wapacket with best tau ...
#ifdef USE_GPU
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    call PropagationElHl_gpucaller(_electron_, Coulomb_, N, S_inv(1,1), h(1,1), Psi_t_bra(1,1), Psi_t_ket(1,1), t_init, t_max, tau(1), save_tau(1))
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#else
    CALL Propagation( N , H_prime , Psi_t_bra(:,1) , Psi_t_ket(:,1) , t_init , t_max, tau(1) , save_tau(1) )
#endif
End If

!=======================================================================
!            Hole Hamiltonian : lower triangle ...
!=======================================================================

If( eh_tag(2) == "hl" ) then  ! <==  proceed evolution of HOLE wapacket with best tau ...
#ifdef USE_GPU
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    call PropagationElHl_gpucaller(_hole_, Coulomb_, N, S_inv, h, Psi_t_bra(1,2), Psi_t_ket(1,2), t_init, t_max, tau(2), save_tau(2))
    !GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#else
    CALL Propagation( N , H_prime , Psi_t_bra(:,2) , Psi_t_ket(:,2) , t_init , t_max , tau(2) , save_tau(2) )
#endif
End If

!=======================================================================

t = t_init + (delta_t*frame_step)

! prepare DUAL basis for local properties ...
DUAL_bra = dconjg(Psi_t_ket)
DUAL_ket = Psi_t_bra

! prepare Slater basis for FORCE properties ...
call op_x_ket( AO_bra, S_inv, Psi_t_bra )
AO_bra = dconjg(AO_bra)
AO_ket = Psi_t_ket

CALL QuasiParticleEnergies( AO_bra , AO_ket , H )
CALL store_Hprime( N , H_prime ) 

! save populations(time) ...
QDyn%dyn(it,:,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

if( mod(it,CT_dump_step) == 0 ) CALL dump_Qdyn( Qdyn , it )

! clean and exit ...
If( driver /= 'q_dynamics' ) then
    deallocate( h0 )
    nullify( h )
end if
#ifndef USE_GPU
deallocate( H_prime )
#endif
#undef S_inv

first_call_ = .false.

include 'formats.h'

end subroutine ElHl_Chebyshev
!
!
!
!===================================================
 function Build_Huckel( basis , S_matrix ) result(h)
!===================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

! local variables ... 
integer               :: i , j , N
real*8  , allocatable :: h(:,:)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
!----------------------------------------------------------

N = size(basis)
ALLOCATE( h(N,N) , source = D_zero )

do j = 1 , N
  do i = 1 , j

        h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j)

        h(j,i) = h(i,j)

    end do
end do

end function Build_Huckel
!
!
!
!
!=======================================================
 subroutine QuasiParticleEnergies( AO_bra , AO_ket , H )
!=======================================================
implicit none
complex*16 , intent(in) :: AO_bra(:,:)
complex*16 , intent(in) :: AO_ket(:,:)
real*8     , intent(in) :: H(:,:)

!local variables ...
integer :: i , j , mm
complex*16 :: erg_el , erg_hl

mm = size(AO_bra(:,1))

erg_el = (0.d0,0.d0)
erg_hl = (0.d0,0.d0)

If( eh_tag(1) == "el" ) then  
    !$OMP parallel do private(i,j) default(shared) reduction(+ : erg_el )
    do j = 1 , mm
        do i = 1 , mm
            erg_el = erg_el + AO_bra(i,1)*H(i,j)*AO_ket(j,1)
        end do
    end do
    !$OMP end parallel do  
End If


If( eh_tag(2) == "hl" ) then  
    !$OMP parallel do private(i,j) default(shared) reduction(+ : erg_hl)
    do j = 1 , mm
        do i = 1 , mm
            erg_hl = erg_hl + AO_bra(i,2)*H(i,j)*AO_ket(j,2)
        end do
    end do
    !$OMP end parallel do  
End If

Unit_Cell% QM_wp_erg(1) = erg_el
Unit_Cell% QM_wp_erg(2) = erg_hl

! QM_erg = E_occ - E_empty ; to be used in MM_dynamics energy balance ...
Unit_Cell% QM_erg = erg_el - erg_hl

end subroutine QuasiParticleEnergies
!
!
!
!
!=================================================================================
 subroutine preprocess_from_restart( system , basis , DUAL_ket , AO_bra , AO_ket )
!=================================================================================
implicit none
type(structure) , intent(inout) :: system
type(STO_basis) , intent(inout) :: basis(:)
complex*16      , intent(in)    :: DUAL_ket (:,:)
complex*16      , intent(in)    :: AO_bra   (:,:)
complex*16      , intent(in)    :: AO_ket   (:,:)

!local variables ...
integer :: N

N = size(basis)

#ifdef USE_GPU
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
    allocate( S_matrix(N,N) )
    call GPU_Pin(S_matrix, N*N*8)
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#endif

CALL Overlap_Matrix( system , basis , S_matrix )

allocate( h0(N,N) , source = D_zero )

If( EnvField_ ) then
    h0(:,:) = even_more_extended_Huckel( system , basis , S_matrix )
else
    h0(:,:) = Build_Huckel( basis , S_matrix )
end If

!vector states to be propagated ...
allocate( Psi_t_bra(N,n_part) )
allocate( Psi_t_ket(N,n_part) )

Psi_t_bra = DUAL_ket
Psi_t_ket = AO_ket

CALL QuasiParticleEnergies(AO_bra, AO_ket, h0)

!==============================================
! preprocess stuff for EhrenfestForce ...

  CALL syInvert( S_matrix, return_full ) ! <== S_matrix content is destroyed and S_inv is returned
#define S_inv S_matrix

  allocate( H(N,N) )
  h => h0

! allocate and compute H' = S_inv * H ...
  allocate( H_prime(N,N) , source = D_zero )
  CALL syMultiply( S_inv , h , H_prime )

  CALL store_Hprime( N , H_prime )
!==============================================

!==============================================
! clean and exit ...
deallocate( h0 )
nullify( h )
#ifndef USE_GPU
deallocate( H_prime )
#endif
#undef S_inv
!==============================================

! IF QM_erg < 0 => turn off QMMM ; IF QM_erg > 0 => turn on QMMM ...
QMMM = (.NOT. (Unit_Cell% QM_erg < D_zero)) .AND. (HFP_Forces == .true.)

end subroutine preprocess_from_restart
!
!
!
end module ElHl_Chebyshev_m


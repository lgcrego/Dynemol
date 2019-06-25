module ElHl_Chebyshev_m

    use type_m              , g_time => f_time  
    use blas95
    use lapack95
    use constants_m
    use ifport
    use parameters_m        , only : t_i , frame_step , Coulomb_ , n_part, driver , QMMM , CT_dump_step , HFP_Forces
    use Structure_Builder   , only : Unit_Cell 
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis , eh_tag    
    use Data_Output         , only : Populations 
    use Coulomb_SMILES_m    , only : Build_Coulomb_potential
!    use Chebyshev_m         , only : Propagation, dump_Qdyn
    use Taylor_m            , only : Propagation, dump_Qdyn
    use Ehrenfest_Builder   , only : store_Hprime                               
    use Matrix_Math

    public  :: ElHl_Chebyshev , preprocess_ElHl_Chebyshev 

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-12
    real*8      , parameter :: norm_error   = 1.0d-12

! module variables ...
    real*8      ,   save          :: save_tau(2)
    logical     ,   save          :: ready = .false.
    logical     ,   save          :: necessary_  = .true.
    logical     ,   save          :: first_call_ = .true.
    real*8, target, allocatable   :: h0(:,:)
    real*8      ,   allocatable   :: S_matrix(:,:)
    real*8      ,   allocatable   :: V_Coul_El(:) , V_Coul_Hl(:)
    complex*16  ,   allocatable   :: V_Coul(:,:)
    complex*16  ,   allocatable   :: Psi_t_bra(:,:) , Psi_t_ket(:,:)

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
!
!==========================================================================================================
 subroutine preprocess_ElHl_Chebyshev( system , basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , it )
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
integer                         :: li , N 
real*8          , allocatable   :: wv_FMO(:)
complex*16      , allocatable   :: ElHl_Psi(:,:)
type(R_eigen)                   :: FMO

!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#ifdef USE_GPU
allocate( S_matrix(size(basis),size(basis)) )
call GPU_Pin(S_matrix, size(basis)*size(basis)*8)
#endif
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

! MUST compute S_matrix before FMO analysis ...
CALL Overlap_Matrix( system , basis , S_matrix )
CALL Huckel( basis , S_matrix , h0 )
! for a rigid structure once is enough ...
If( driver == 'q_dynamics' ) necessary_ = .false.

!========================================================================
! prepare electron state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="E" )

! place the electron state in Structure's Hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%El )
N  = size(wv_FMO)
allocate( ElHl_Psi( size(basis) , n_part ) , source=C_zero )
ElHl_Psi(li:li+N-1,1) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )

!========================================================================
! prepare hole state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="H" )

! place the hole state in Structure's Hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%Hl )
N  = size(wv_FMO)
ElHl_Psi(li:li+N-1,2) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )
!========================================================================

!==============================================
! prepare DUAL basis for local properties ...
! DUAL_bra = (C*)^T    ;    DUAL_ket = S*C ...
DUAL_bra = dconjg( ElHl_Psi )
call op_x_ket( DUAL_ket, S_matrix , ElHl_Psi )
!==============================================

!==============================================
!vector states to be propagated ...
!Psi_bra = C^T*S       ;      Psi_ket = C ...
allocate( Psi_t_bra(size(basis),n_part) )
allocate( Psi_t_ket(size(basis),n_part) )
call bra_x_op( Psi_t_bra, ElHl_Psi , S_matrix ) 
Psi_t_ket = ElHl_Psi
!==============================================

!==============================================
AO_bra = ElHl_Psi 
AO_ket = ElHl_Psi 
CALL QuasiParticleEnergies(AO_bra, AO_ket, H0)
!==============================================

! save populations(time=t_i) ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )
CALL dump_Qdyn( Qdyn , it )

! leaving S_matrix allocated

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
integer                     :: j , N
real*8                      :: t_max , tau_max , tau(2) , t_init
real*8      , pointer       :: h(:,:)
real*8      , allocatable   :: V_coul_El(:) , V_coul_Hl(:) , H_prime(:,:)

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

    CALL Huckel( basis , S_matrix , h0 )

    CALL syInvert( S_matrix, return_full )   ! S_matrix content is destroyed and S_inv is returned
#define S_inv S_matrix

end If

If( Coulomb_ ) then

    If( ready ) then
        CALL Build_Coulomb_Potential( system , basis , AO_bra , AO_ket , V_coul , V_coul_El , V_coul_Hl )
        deallocate( V_Coul )
    else
        allocate( V_coul_El(N) , source = D_zero )
        allocate( V_coul_Hl(N) , source = D_zero )
        ready = .true.
    end If

end If    

!=======================================================================
!           Electron Hamiltonian : upper triangle of V_coul ...
!=======================================================================

if( Coulomb_ ) then
    allocate( h(N,N), source = h0 )   ! h = h0
    forall( j=1:N ) h(j,j) = h(j,j) + V_coul_El(j)
else
    h => h0
end if

#ifndef USE_GPU
! allocate and compute H' = S_inv * H ...
If (.not. allocated(H_prime) ) allocate( H_prime(N,N) )
call syMultiply( S_inv , h , H_prime )
#endif

! proceed evolution of ELECTRON wapacket with best tau ...
#ifdef USE_GPU
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
call PropagationElHl_gpucaller(_electron_, Coulomb_, N, S_inv(1,1), h(1,1), Psi_t_bra(1,1), Psi_t_ket(1,1), t_init, t_max, tau(1), save_tau(1))
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#else
CALL Propagation( N , H_prime , Psi_t_bra(:,1) , Psi_t_ket(:,1) , t_init , t_max, tau(1) , save_tau(1) )
#endif

!=======================================================================
!            Hole Hamiltonian : lower triangle of V_coul ...
!=======================================================================

if( Coulomb_ ) forall(j=1:N) h(j,j) = h0(j,j) + V_coul_Hl(j)

#ifdef USE_GPU
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
call PropagationElHl_gpucaller(_hole_, Coulomb_, N, S_inv, h, Psi_t_bra(1,2), Psi_t_ket(1,2), t_init, t_max, tau(2), save_tau(2))
!GGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#endif

if( Coulomb_ ) then
#ifndef USE_GPU
    call syMultiply( S_inv , h , H_prime )
#endif
    deallocate( h, V_coul_El, V_coul_Hl)
end if

#ifndef USE_GPU
! proceed evolution of HOLE wapacket with best tau ...
CALL Propagation( N , H_prime , Psi_t_bra(:,2) , Psi_t_ket(:,2) , t_init , t_max , tau(2) , save_tau(2) )
#endif

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
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

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
!=========================================
subroutine Huckel( basis , S_matrix , h0 )
!=========================================
implicit none
type(STO_basis)               , intent(in)  :: basis(:)
real*8                        , intent(in)  :: S_matrix(:,:)
real*8          , allocatable , intent(out) :: h0(:,:)

! local variables ... 
real*8  :: k_eff , k_WH , c1 , c2 , c3 , c4
integer :: i , j

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

ALLOCATE( h0(size(basis),size(basis)) )

do j = 1 , size(basis)

    do i = 1 , j - 1

       c1 = basis(i)%IP - basis(j)%IP
       c2 = basis(i)%IP + basis(j)%IP
    
       c3 = (c1/c2)*(c1/c2)
    
       c4 = (basis(i)%V_shift + basis(j)%V_shift) * HALF
    
       k_WH = (basis(i)%k_WH + basis(j)%k_WH) * HALF
    
       k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)
    
       h0(i,j) = (k_eff*c2*HALF + c4) * S_matrix(i,j)

    end do

    h0(j,j) = basis(j)%IP + basis(j)%V_shift

end do

call Matrix_Symmetrize( h0, 'U' )

end subroutine Huckel
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
!$OMP parallel do private(i,j) default(shared) reduction(+ : erg_el , erg_hl)
do j = 1 , mm
    do i = 1 , mm

        erg_el = erg_el + AO_bra(i,1)*H(i,j)*AO_ket(j,1)
        erg_hl = erg_hl + AO_bra(i,2)*H(i,j)*AO_ket(j,2)

    end do
end do
!$OMP end parallel do  

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

!vector states to be propagated ...
allocate( Psi_t_bra(size(basis),n_part) )
allocate( Psi_t_ket(size(basis),n_part) )

Psi_t_bra = DUAL_ket
Psi_t_ket = AO_ket

CALL Overlap_Matrix( system , basis , S_matrix )
CALL Huckel( basis , S_matrix , h0 )

CALL QuasiParticleEnergies(AO_bra, AO_ket, h0)

! IF QM_erg < 0 => turn off QMMM ; IF QM_erg > 0 => turn on QMMM ...
QMMM = (.NOT. (Unit_Cell% QM_erg < D_zero)) .AND. (HFP_Forces == .true.)

end subroutine preprocess_from_restart
!
!
!
end module ElHl_Chebyshev_m

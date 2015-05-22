module ElHl_Chebyshev_m

    use type_m              , g_time => f_time  
    use blas95
    use lapack95
    use constants_m
    use ifport
    use parameters_m        , only : t_i , frame_step ,         &
                                     Coulomb_ , n_part,         &
                                     driver , restart
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis , eh_tag    
    use Data_Output         , only : Populations 
    use Coulomb_SMILES_m    , only : Build_Coulomb_potential
    use Chebyshev_m         , only : Propagation, dump_Qdyn
    use Matrix_Math

    public  :: ElHl_Chebyshev , preprocess_ElHl_Chebyshev

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-12
    real*8      , parameter :: norm_error   = 1.0d-12

! module variables ...
    real*8      ,   save          :: save_tau(2)
    logical     ,   save          :: done = .false.
    logical     ,   save          :: necessary_  = .true.
    logical     ,   save          :: first_call_ = .true.
    real*8, target, allocatable   :: h0(:,:)
    real*8      ,   allocatable   :: S_inv(:,:) 
    real*8      ,   allocatable   :: V_Coul_El(:) , V_Coul_Hl(:)
    complex*16  ,   allocatable   :: V_Coul(:,:)
    
#ifdef USE_GPU
#  define _electron_  -1
#  define _hole_      +1
#endif

contains
!
!
!
!============================================================================================================
 subroutine preprocess_ElHl_Chebyshev( system , basis , Psi_bra , Psi_ket , Dual_bra , Dual_ket , QDyn , it )
!============================================================================================================
implicit none
type(structure) , intent(inout) :: system
type(STO_basis) , intent(inout) :: basis(:)
complex*16      , intent(out)   :: Psi_bra(:,:)
complex*16      , intent(out)   :: Psi_ket(:,:)
complex*16      , intent(out)   :: DUAL_bra(:,:) 
complex*16      , intent(out)   :: DUAL_ket(:,:) 
type(g_time)    , intent(inout) :: QDyn
integer         , intent(in)    :: it

!local variables ...
integer                         :: li , N 
real*8          , allocatable   :: wv_FMO(:) , S_matrix(:,:)
complex*16      , allocatable   :: ElHl_Psi(:,:)
type(R_eigen)                   :: FMO


!========================================================================
! prepare electron state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="E" )

! place the electron state in Structure's hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%El )
N  = size(wv_FMO)

allocate( ElHl_Psi( size(basis) , n_part ) , source=C_zero )
ELHl_Psi(li:li+N-1,1) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )
!========================================================================
! prepare hole state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="H" )

! place the hole state in Structure's hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%Hl )
N  = size(wv_FMO)

ElHl_Psi(li:li+N-1,2) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )
!========================================================================

! prepare DUAL basis for local properties ...
CALL Overlap_Matrix( system , basis , S_matrix )
DUAL_bra = dconjg( ElHl_Psi )
call op_x_ket( DUAL_ket, S_matrix , ElHl_Psi )

call bra_x_op( Psi_bra, ElHl_Psi , S_matrix )
Psi_ket = ElHl_Psi

If( .not. restart ) then
    ! save populations(time=t_i) ...
    QDyn%dyn(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

    CALL dump_Qdyn( Qdyn , it )
end If    

! clean and exit ...
deallocate( S_matrix )

end subroutine preprocess_ElHl_Chebyshev
!
!
!
!===================================================================================================================
 subroutine ElHl_Chebyshev( system , basis , Psi_t_bra , Psi_t_ket , Dual_bra , Dual_ket , QDyn , t , delta_t , it )
!===================================================================================================================
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)
complex*16       , intent(inout) :: Psi_t_bra(:,:)
complex*16       , intent(inout) :: Psi_t_ket(:,:)
complex*16       , intent(inout) :: Dual_bra(:,:)
complex*16       , intent(inout) :: Dual_ket(:,:)
type(g_time)     , intent(inout) :: QDyn
real*8           , intent(inout) :: t
real*8           , intent(in)    :: delta_t
integer          , intent(in)    :: it

! local variables...
integer                     :: i, j , N
real*8                      :: t_max , tau_max , tau(2) , t_init
real*8      , pointer       :: h(:,:)
real*8      , allocatable   :: V_coul_El(:) , V_coul_Hl(:) , H_prime(:,:) , S_matrix(:,:) 
complex*16  , allocatable   :: AO_bra(:,:) , AO_ket(:,:) , V_coul(:,:) 

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

If ( necessary_ ) then

    ! compute S and S_inverse ...
#define S_matrix S_inv
    CALL Overlap_Matrix( system , basis , S_matrix )

    CALL Huckel( basis , S_matrix , h0 )

    call GPU_Pin( S_matrix, N*N*8 )
    call syInvert( S_matrix, return_full )   ! S_matrix content is destroyed and S_inv is returned
#undef S_matrix

    ! for a rigid structure once is enough ...
    If( driver == 'q_dynamics' ) necessary_ = .false.

end If


If( Coulomb_ ) then

    If( done ) then

        allocate( AO_bra(N,n_part) )
        allocate( AO_ket(N,n_part) )
        call op_x_ket( AO_bra, S_inv, Psi_t_bra )
        forall( i=1:N, j=1:n_part ) AO_bra(i,j) = dconjg(AO_bra(i,j))
        AO_ket = Psi_t_ket

        CALL Build_Coulomb_Potential( system , basis , AO_bra , AO_ket , V_coul , V_coul_El , V_coul_Hl )
        deallocate( V_Coul , AO_bra , AO_ket )

    else

        allocate( V_coul_El(N) , source = D_zero )
        allocate( V_coul_Hl(N) , source = D_zero )

        done = .true.

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
allocate( H_prime(N,N) )
call syMultiply( S_inv , h , H_prime )
#endif

! proceed evolution of ELECTRON wapacket with best tau ...
#ifdef USE_GPU
call PropagationElHl_gpucaller(_electron_, Coulomb_, N, S_inv, h, Psi_t_bra(1,1), Psi_t_ket(1,1), t_init, t_max, tau(1), save_tau(1))
#else
CALL Propagation( N , H_prime , Psi_t_bra(:,1) , Psi_t_ket(:,1) , t_init , t_max, tau(1) , save_tau(1) )
#endif

!=======================================================================
!            Hole Hamiltonian : lower triangle of V_coul ...
!=======================================================================

if( Coulomb_ ) forall(j=1:N) h(j,j) = h0(j,j) + V_coul_Hl(j)

#ifdef USE_GPU
call PropagationElHl_gpucaller(_hole_, Coulomb_, N, S_inv, h, Psi_t_bra(1,2), Psi_t_ket(1,2), t_init, t_max, tau(2), save_tau(2))
#endif

if( Coulomb_ ) then
#ifndef USE_GPU
    call syMultiply( S_inv , h , H_prime )
#endif
    deallocate( h, V_coul_El, V_coul_Hl)
end if

If( driver /= 'q_dynamics' ) then
    call GPU_Unpin( S_inv )
    deallocate( h0 , S_inv )
    nullify( h )
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

! save populations(time) ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

CALL dump_Qdyn( Qdyn , it )

! clean and exit ...
#ifndef USE_GPU
deallocate( H_prime )
#endif

Print 186, t

include 'formats.h'

first_call_ = .false.

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
real*8  :: k_eff , k_WH , c1 , c2 , c3
integer :: i , j

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

ALLOCATE( h0(size(basis),size(basis)) )

do j = 1 , size(basis)

    do i = 1 , j - 1

        c1 = basis(i)%IP - basis(j)%IP
        c2 = basis(i)%IP + basis(j)%IP

        c3 = (c1/c2)*(c1/c2)

        k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

        k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

        h0(i,j) = k_eff * S_matrix(i,j) * c2 / two

    end do

    h0(j,j) = basis(j)%IP

end do

call Matrix_Symmetrize( h0, 'U' )

end subroutine Huckel
!
!
!
!
end module ElHl_Chebyshev_m

module Taylor_m

    use type_m              , g_time => f_time  
    use constants_m
    use blas95
    use lapack95
    use ifport
    use parameters_m        , only : t_i, frame_step,               &
                                     DP_Field_, driver,             &
                                     n_part, restart, CT_dump_step                  
    use Structure_Builder   , only : Unit_Cell                                      
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis, eh_tag                  
    use QCmodel_Huckel      , only : Huckel,                        &
                                     even_more_extended_Huckel
    use Data_Output         , only : Populations
    use Matrix_Math

    public  :: Taylor, preprocess_Taylor, Propagation, dump_Qdyn

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-8
    real*8      , parameter :: norm_error   = 1.0d-8

! module variables ...
    real*8  , save :: save_tau 
    logical , save :: necessary_  = .true.
    logical , save :: first_call_ = .true.
    real*8  , allocatable , save :: H(:,:), H_prime(:,:), S_matrix(:,:)

contains
!
!
!
!=====================================================================================================
 subroutine preprocess_Taylor( system , basis , Psi_bra , Psi_ket , Dual_bra , Dual_ket , QDyn , it )
!! Excatly equal to preprocess_Chebyshev
!=====================================================================================================
implicit none
type(structure) , intent(inout) :: system
type(STO_basis) , intent(inout) :: basis(:)
complex*16      , intent(out)   :: Psi_bra(:)
complex*16      , intent(out)   :: Psi_ket(:)
complex*16      , intent(out)   :: Dual_bra(:)
complex*16      , intent(out)   :: Dual_ket(:)
type(g_time)    , intent(inout) :: QDyn
integer         , intent(in)    :: it

!local variables ...
integer                         :: li , N 
real*8          , allocatable   :: wv_FMO(:)
complex*16      , allocatable   :: Psi(:)
type(R_eigen)                   :: FMO

! prepare  DONOR  state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="E" )

! place the  DONOR  state in Structure's hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%El )
N  = size(wv_FMO)

allocate( Psi(size(basis)) , source=C_zero )
Psi(li:li+N-1) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )

#ifdef USE_GPU
allocate( S_matrix(size(basis),size(basis)) )
call GPU_Pin(S_matrix, size(basis)*size(basis)*8)
#endif

! prepare DUAL basis for local properties ...
CALL Overlap_Matrix( system , basis , S_matrix )
DUAL_bra = dconjg( Psi )
call op_x_ket( DUAL_ket, S_matrix , Psi )

call bra_x_op( Psi_bra, Psi , S_matrix )
Psi_ket = Psi

If( .not. restart ) then
    ! save populations(time=t_i) ...
    QDyn%dyn(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

    CALL dump_Qdyn( Qdyn , it )
end If

! leaving S_matrix allocated

end subroutine preprocess_Taylor
!
!
!
!==============================================================================================================
 subroutine Taylor( system , basis , Psi_t_bra , Psi_t_ket , Dual_bra , Dual_ket , QDyn , t , delta_t , it )
!! Equal to Chebyshev
!==============================================================================================================
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)
complex*16       , intent(inout) :: Psi_t_bra(:)
complex*16       , intent(inout) :: Psi_t_ket(:)
complex*16       , intent(inout) :: Dual_bra(:)
complex*16       , intent(inout) :: Dual_ket(:)
type(g_time)     , intent(inout) :: QDyn
real*8           , intent(inout) :: t
real*8           , intent(in)    :: delta_t
integer          , intent(in)    :: it

! local variables... 
integer                          :: N
real*8                           :: tau , tau_max , t_init, t_max 

t_init = t

! max time inside slice ...
t_max = delta_t*frame_step*(it-1)  
! constants of evolution ...
tau_max = delta_t / h_bar

! trying to adapt time step for efficient propagation ...
tau = merge( tau_max , save_tau * 1.15d0 , first_call_ )
! but tau should be never bigger than tau_max ...
tau = merge( tau_max , tau , tau > tau_max )

N = size(basis)

if(first_call_) then           ! allocate matrices
    allocate( H(N,N) , H_prime(N,N) )
end if

If ( necessary_ ) then
    
    call Overlap_Matrix( system , basis , S_matrix )
    call Huckelx( basis , S_matrix , H )
    call syInvert( S_matrix )   ! S_matrix content is destroyed and S_inv is returned
    call syMultiply( S_matrix , H , H_prime )

    ! for a rigid structure once is enough ...
    If( driver ==  "q_dynamics" ) necessary_ = .false.

end If

! Propagation(..., type=Taylor)? 
call Propagation( N, H_prime, Psi_t_bra, Psi_t_ket, t_init, t_max, tau, save_tau )

t = t_init + (delta_t*frame_step)

! prepare DUAL basis for local properties ...
DUAL_bra = dconjg(Psi_t_ket)
DUAL_ket = Psi_t_bra

! save populations(time) ...
QDyn%dyn(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

if( mod(it,CT_dump_step) == 0 ) CALL dump_Qdyn( Qdyn , it )

Print 186, t

include 'formats.h'

first_call_ = .false.

end subroutine Taylor
!
!
!
!=============================================================================================
 subroutine Propagation( N , H_prime , Psi_t_bra , Psi_t_ket , t_init , t_max , tau , save_tau )
!=============================================================================================
implicit none
integer         , intent(in)    :: N
real*8          , intent(in)    :: H_prime(:,:)
complex*16      , intent(inout) :: Psi_t_bra(:)
complex*16      , intent(inout) :: Psi_t_ket(:)
real*8          , intent(in)    :: t_init
real*8          , intent(in)    :: t_max
real*8          , intent(inout) :: tau    
real*8          , intent(out)   :: save_tau

! local variables...
complex*16                  :: r
complex*16  , allocatable   :: bra(:,:), ket(:,:), tmp_bra(:), tmp_ket(:), C(:)
real*8                      :: norm_ref, norm_test, t
integer                     :: k, k_ref
logical                     :: OK


! complex variables ...
allocate( bra     ( N , order ) )   ! redundant, only need n=2 vectors to generate the series, not n=order vectors
allocate( ket     ( N , order ) )
allocate( tmp_bra ( N         ) )
allocate( tmp_ket ( N         ) )
allocate( C       (     order ) )

norm_ref = abs(dotc(Psi_t_bra, Psi_t_ket))

! first convergence: best tau-parameter for k_ref ...
do
    call Convergence( Psi_t_bra, Psi_t_ket, C, k_ref, tau, H_prime, norm_ref, OK )

    if( OK ) exit
    tau = tau * 0.9d0
end do
save_tau = tau

t = t_init + tau*h_bar

if( t_max-t < tau*h_bar ) then
    tau = ( t_max - t ) / h_bar
    C = coefficient(tau,order)
end if

! proceed evolution with best tau ...
do while( t < t_max )
    
    ! Ѱ₁ = c₁|Ѱ⟩ = |Ѱ⟩
    bra(:,1) = Psi_t_bra
    ket(:,1) = Psi_t_ket

    tmp_bra = bra(:,1)
    tmp_ket = ket(:,1)

    do k = 2, k_ref

        ! Ѱₙ = (cₙ/cₙ₋₁) H'Ѱₙ₋₁
        r = c(k)/c(k-1)
        call bra_x_op( bra(:,k), bra(:,k-1), H_prime, r )
        call op_x_ket( ket(:,k), H_prime, ket(:,k-1), r )

        ! add term Ѱₙ to series expansion
        tmp_bra = tmp_bra + bra(:,k)
        tmp_ket = tmp_ket + ket(:,k)
    end do

    ! convergence criteria
    norm_test = abs(dotc( tmp_bra, tmp_ket ))
    if (abs( norm_test - norm_ref ) < norm_error) then
        Psi_t_bra = tmp_bra
        Psi_t_ket = tmp_ket
    else
        OK = .false.
        do while( .not. OK )
            tau = tau * 0.975d0
            print*, "rescaling tau", tau
            call Convergence( Psi_t_bra, Psi_t_ket, C, k_ref, tau, H_prime, norm_ref, OK )
        end do
    end if

    t = t + (tau * h_bar)

    if( t_max-t < tau*h_bar ) then
        tau = ( t_max-t ) / h_bar
        C = coefficient(tau, order)
    end if

end do

deallocate( bra, ket, tmp_bra, tmp_ket, C )

end subroutine Propagation
!
!
!
!===============================================================================
subroutine Convergence( Psi_bra, Psi_ket, C, k_ref, tau, H_prime, norm_ref, OK )
!===============================================================================
implicit none
complex*16  , intent(inout) :: Psi_bra(:)
complex*16  , intent(inout) :: Psi_ket(:)
complex*16  , intent(out)   :: C(:)
integer     , intent(inout) :: k_ref
real*8      , intent(in)    :: tau
real*8      , intent(in)    :: H_prime(:,:)
real*8      , intent(in)    :: norm_ref
logical     , intent(out)   :: OK

! local variables...
integer                     :: k, k_max, N  
real*8                      :: norm_tmp
complex*16                  :: r
complex*16  , allocatable   :: bra(:,:), ket(:,:), tmp_bra(:,:), tmp_ket(:,:)

#define old 1
#define new 2

N = size(Psi_bra)

allocate( bra      ( N , order ) , source=C_zero )
allocate( ket      ( N , order ) , source=C_zero )
allocate( tmp_bra  ( N , 2     ) )
allocate( tmp_ket  ( N , 2     ) )

OK = .false.

! get C_k coefficients ...
C = coefficient(tau, order)

k_max = order
do k = 2, order
    if( abs( c(k) ) < 1.0d-16 ) then
        k_max = k
        exit
    end if
end do

k_ref = k_max

! Ѱ₁ = c₁|Ѱ⟩ = |Ѱ⟩
bra(:,1) = Psi_bra(:)
ket(:,1) = Psi_ket(:)

tmp_bra(:,old) = bra(:,1)
tmp_ket(:,old) = ket(:,1)

do k = 2, k_max

    ! Ѱₙ = (cₙ/cₙ₋₁) H'Ѱₙ₋₁
    r = c(k)/c(k-1)
    call bra_x_op( bra(:,k), bra(:,k-1), H_prime, r )
    call op_x_ket( ket(:,k), H_prime, ket(:,k-1), r )

    ! add term Ѱₙ to the expansion
    tmp_bra(:,new) = tmp_bra(:,old) + bra(:,k)
    tmp_ket(:,new) = tmp_ket(:,old) + ket(:,k)

!   convergence criteria...
    if( isConverged( tmp_bra(:,new), tmp_bra(:,old), error ) ) then
    if( isConverged( tmp_ket(:,new), tmp_ket(:,old), error ) ) then

        norm_tmp = abs(dotc( tmp_bra(:,new), tmp_ket(:,new) ))

        if( abs( norm_tmp - norm_ref ) < norm_error ) then
            Psi_bra = tmp_bra(:,new)
            Psi_ket = tmp_ket(:,new)
            ok = .true.
            exit
        end if

    end if ! ket conv.
    end if ! bra conv.

    tmp_bra(:,old) = tmp_bra(:,new)
    tmp_ket(:,old) = tmp_ket(:,new)

end do !k

deallocate( bra, ket, tmp_bra, tmp_ket )

#undef old
#undef new

end subroutine Convergence
!
!
!
!==================================
 function coefficient(tau , k_max ) 
!==================================
implicit none
complex*16  , dimension(k_max)  :: coefficient
real*8      , intent(in)        :: tau
integer     , intent(in)        :: k_max

!local variables ...
integer :: k

coefficient(1) = C_one
do k = 2 , k_max
   coefficient(k) = -zi * coefficient(k-1) * (tau/(k-1))
end do

end function coefficient
!
!
!
!=================================
 subroutine dump_Qdyn( Qdyn , it )
!=================================
implicit none
type(g_time)    , intent(in) :: QDyn
integer         , intent(in) :: it 

! local variables ...
integer :: nf , n

do n = 1 , n_part

    If( it == 1 ) then
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
        write(52,15) "#" ,( nf+1 , nf=0,size(QDyn%fragments)+1 )  ! <== numbered columns for your eyes only ...
        write(52,12) "#" , QDyn%fragments , "total"

        open( unit = 53 , file = "tmp_data/"//eh_tag(n)//"_wp_energy.dat" , status = "replace" , action = "write" , position = "append" )

    else
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat"  , status = "unknown", action = "write" , position = "append" )
        open( unit = 53 , file = "tmp_data/"//eh_tag(n)//"_wp_energy.dat" , status = "unknown", action = "write" , position = "append" )        
    end If

    ! dumps el-&-hl populations ...    
    write(52,13) ( QDyn%dyn(it,nf,n) , nf=0,size(QDyn%fragments)+1 )

    ! dumps el-&-hl wavepachet energies ...
    write(53,14) QDyn%dyn(it,0,n) , real( Unit_Cell% QM_wp_erg(n) ) , dimag( Unit_Cell% QM_wp_erg(n) )

    close(52)
    close(53)

end do

12 FORMAT(/15A10)
13 FORMAT(F11.6,14F10.5)
14 FORMAT(3F12.6)
15 FORMAT(A,I9,14I10)

end subroutine dump_Qdyn
!
!
!
!==============================
function isConverged( a, b, tol )
! returns true if abs(a-b)<tol
!==============================
    logical                :: isConverged
    complex*16, intent(in) :: a(:), b(:)
    real*8,     intent(in) :: tol
    integer :: i
    
    isConverged = .false.
    do i = 1, size(a)
        if( abs(a(i)-b(i)) > tol ) return  ! allow earlier return if not converged
    end do
    isConverged = .true.
end function isConverged
!
!
!
!=========================================
subroutine Huckelx( basis , S_matrix , H )
!=========================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)
real*8          , intent(out)   :: H(:,:)

! local variables ... 
real*8  :: k_eff , k_WH , c1 , c2 , c3
real*8  :: basis_j_IP, basis_j_k_WH
integer :: i, j, n

n = size(basis)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

!$omp parallel private(i,j,basis_j_IP,basis_j_k_WH,c1,c2,c3,k_WH,k_eff) default(shared)
!$omp do schedule(dynamic,1)
do j = 1, n

    basis_j_IP   = basis(j)%IP
    basis_j_k_WH = basis(j)%k_WH

    do i = 1, j - 1

        c1 = basis(i)%IP - basis_j_IP
        c2 = basis(i)%IP + basis_j_IP
        
        c3 = (c1/c2)**2

        k_WH = (basis(i)%k_WH + basis_j_k_WH) * half

        k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

        H(i,j) = k_eff * S_matrix(i,j) * c2 * half

    end do

    H(j,j) = basis_j_IP

end do
!$omp end do

! call Matrix_Symmetrize( H, 'U' )
!$omp do
do i = 1, n
    H( i+1:n, i ) = H( i, i+1:n )
end do
!$omp end do nowait
!$omp end parallel

end subroutine Huckelx


end module Taylor_m




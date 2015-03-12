module Chebyshev_m

    use type_m              , g_time => f_time  
    use constants_m
    use blas95
    use lapack95
    use ifport
    use parameters_m        , only : t_i , t_f , n_t ,              &
                                     frame_step , DP_Field_ ,       &
                                     file_type , n_part                   
    use Babel_m             , only : MD_dt
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis , eh_tag                  
    use QCmodel_Huckel      , only : Huckel ,                       &
                                     Huckel_with_FIELDS
    use Data_Output         , only : Populations
    use Matrix_Math

    public  :: Chebyshev , preprocess_Chebyshev

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-12
    real*8      , parameter :: norm_error   = 1.0d-12

! module variables ...
    real*8  , save :: save_tau 
    logical , save :: necessary = .true.
    real*8  , allocatable , save :: H_prime(:,:)

contains
!
!
!
!=================================================================================
 subroutine preprocess_Chebyshev( system , basis , Psi_bra , Psi_ket , QDyn , it )
!=================================================================================
implicit none
type(structure)                 , intent(inout) :: system
type(STO_basis)                 , intent(inout) :: basis(:)
complex*16      , allocatable   , intent(out)   :: Psi_bra(:)
complex*16      , allocatable   , intent(out)   :: Psi_ket(:)
type(g_time)                    , intent(inout) :: QDyn
integer                         , intent(in)    :: it

!local variables ...
integer                         :: li , N 
real*8          , allocatable   :: wv_FMO(:) 
real*8          , allocatable   :: S_matrix(:,:)
complex*16      , allocatable   :: DUAL_bra(:) , DUAL_ket(:) , Psi(:)
type(R_eigen)                   :: FMO
integer :: i

! prepare  DONOR  state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="E" )

! place the  DONOR  state in Structure's hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%El )
N  = size(wv_FMO)

allocate( Psi(size(basis)) , source=C_zero )
Psi(li:li+N-1) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )

! prepare DUAL basis for local properties ...
allocate( Dual_bra (size(basis)) )
allocate( Dual_ket (size(basis)) )
CALL Overlap_Matrix( system , basis , S_matrix )
DUAL_bra(:) = dconjg( Psi )
call op_x_ket( DUAL_ket, S_matrix , Psi )

allocate( Psi_bra (size(basis)) )
allocate( Psi_ket (size(basis)) )
call bra_x_op( Psi_bra, Psi , S_matrix )
Psi_ket = Psi

! save populations(time=t_i) ...
QDyn%dyn(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

CALL dump_Qdyn( Qdyn , it )

! clean and exit ...
deallocate( DUAL_ket , S_matrix , DUAL_bra )

end subroutine preprocess_Chebyshev
!
!
!
!==============================================================================
 subroutine Chebyshev( system , basis , Psi_t_bra , Psi_t_ket , QDyn , t , it )
!==============================================================================
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)
complex*16       , intent(inout) :: Psi_t_bra(:)
complex*16       , intent(inout) :: Psi_t_ket(:)
type(g_time)     , intent(inout) :: QDyn
real*8           , intent(out)   :: t
integer          , intent(in)    :: it

! local variables...
complex*16  , allocatable   :: C_Psi_bra(:,:) , C_Psi_ket(:,:)
complex*16  , allocatable   :: Psi_tmp_bra(:) , Psi_tmp_ket(:) , C_k(:) , DUAL_bra(:) , DUAL_ket(:)
real*8                      :: delta_t , tau , tau_max , norm_ref , norm_test
real*8                      :: t_max 
integer                     :: j , k , k_ref , N
logical                     :: OK

! max time inside slice ...
t_max = MD_dt*frame_step*(it-1)  

call start_clock

If ( necessary ) then
    
    ! building  S_matrix  and  H'= S_inv * H ...
    CALL Build_Hprime( system , basis )

    ! for a rigid structure once is enough ...
    If( file_type == 'structure' ) necessary = .false.

end If

N = size(basis)

allocate( C_Psi_bra   (N , order ) , source=C_zero )
allocate( C_Psi_ket   (N , order ) , source=C_zero )
allocate( C_k         (order     ) , source=C_zero )
allocate( Psi_tmp_bra (N         ) )
allocate( Psi_tmp_ket (N         ) )

delta_t = merge( (t_f) / dfloat(n_t) , MD_dt * frame_step , MD_dt == epsilon(1.0) )
tau_max = delta_t / h_bar

! constants of evolution ...

! trying to adapt time step for efficient propagation ...
tau = merge( delta_t / h_bar , save_tau * 1.15d0 , it == 2 )
! but tau should be never bigger than tau_max ...
tau = merge( tau_max , tau , tau > tau_max )


#ifdef USE_GPU
call chebyshev_gpucaller( N, tau, save_tau, t_max), t, Psi_t_bra, Psi_t_ket, H_prime )
#else

norm_ref = abs(dotc( Psi_t_bra , Psi_t_ket ))
k_ref = 0

! first convergence: best tau-parameter for k_ref ...
do
    CALL Convergence( Psi_t_bra , Psi_t_ket , C_k , k_ref , tau , H_prime , norm_ref , OK )

    if( OK ) exit       
    tau = tau * 0.9d0
end do

t = t + tau * h_bar 
save_tau = tau

if( t_max-t < tau*h_bar ) then
    tau = ( t_max-t ) / h_bar
    C_k = coefficient(tau,order)
end if

! proceed evolution with best tau ...
do while( t < t_max )

    C_Psi_bra(:,1) = Psi_t_bra(:)
    C_Psi_ket(:,1) = Psi_t_ket(:)

    call bra_x_op( C_Psi_bra(:,2), C_Psi_bra(:,1), H_prime )
    call op_x_ket( C_Psi_ket(:,2), H_prime, C_Psi_ket(:,1) )

    do j = 3 , k_ref
        call bra_x_op( C_Psi_bra(:,j), C_Psi_bra(:,j-1), H_prime, C_two )
        C_Psi_bra(:,j) = C_Psi_bra(:,j) - C_Psi_bra(:,j-2)
        call op_x_ket( C_Psi_ket(:,j), H_prime, C_Psi_ket(:,j-1), C_two )
        C_Psi_ket(:,j) = C_Psi_ket(:,j) - C_Psi_ket(:,j-2)
    end do

    do j = 1, N
        Psi_tmp_bra(j) = sum( C_k(1:k_ref) * C_Psi_bra(j,1:k_ref) )
        Psi_tmp_ket(j) = sum( C_k(1:k_ref) * C_Psi_ket(j,1:k_ref) )
    end do

!   convergence criteria ...
    norm_test = abs(dotc( Psi_tmp_bra , Psi_tmp_ket ))
    if( abs( norm_test - norm_ref ) < norm_error ) then
        Psi_t_bra(:) = Psi_tmp_bra(:)
        Psi_t_ket(:) = Psi_tmp_ket(:)
    else
        OK = .false.
        do while( .not. OK )
            tau = tau * 0.975d0
            print*, "rescaling tau" , tau
            CALL Convergence( Psi_t_bra , Psi_t_ket , C_k , k_ref , tau , H_prime , norm_ref , OK )
        end do
    end if

    t = t + (tau * h_bar)

    if( t_max-t < tau*h_bar ) then
        tau = (t_max - t) / h_bar
        C_k = coefficient(tau,order)
    end if

end do
#endif

! prepare DUAL basis for local properties ...
allocate( Dual_bra (N) )
allocate( Dual_ket (N) )
DUAL_bra = dconjg(Psi_t_ket)
DUAL_ket = Psi_t_bra

! save populations(time) ...
QDyn%dyn(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

CALL dump_Qdyn( Qdyn , it )

! clean and exit ...
deallocate( C_k , C_Psi_bra , C_Psi_ket , Psi_tmp_bra , Psi_tmp_ket , DUAL_bra , DUAL_ket )
If( file_type == 'trajectory' ) deallocate( H_prime )

Print 186, t

include 'formats.h'

call stop_clock('Chebyshev')

end subroutine Chebyshev
!
!
!
!========================================================================================
subroutine Convergence( Psi_bra , Psi_ket , C_k , k_ref , tau , H_prime , norm_ref , OK )
!========================================================================================
implicit none
complex*16  , intent(inout) :: Psi_bra(:)
complex*16  , intent(inout) :: Psi_ket(:)
complex*16  , intent(out)   :: C_k(:)
integer     , intent(inout) :: k_ref
real*8      , intent(in)    :: tau
real*8      , intent(in)    :: H_prime(:,:)
real*8      , intent(in)    :: norm_ref
logical     , intent(out)   :: OK

! local variables...
integer                     :: j , l , k , N  
real*8                      :: norm_tmp
complex*16  , allocatable   :: C_Psi_bra(:,:) , C_Psi_ket(:,:) , Psi_tmp_bra(:,:) , Psi_tmp_ket(:,:)

N = size(Psi_bra)

allocate( C_Psi_bra    ( N , order ) , source=C_zero )
allocate( C_Psi_ket    ( N , order ) , source=C_zero )
allocate( Psi_tmp_bra  ( N , 2     ) )
allocate( Psi_tmp_ket  ( N , 2     ) )

OK = .false.

! get C_k coefficients ...
 C_k = coefficient(tau, order)

 C_Psi_bra(:,1) = Psi_bra(:)
 C_Psi_ket(:,1) = Psi_ket(:)
call bra_x_op( C_Psi_bra(:,2), C_Psi_bra(:,1), H_prime )
call op_x_ket( C_Psi_ket(:,2), H_prime, C_Psi_ket(:,1) )

Psi_tmp_bra(:,1) = C_k(1)*C_Psi_bra(:,1) + C_k(2)*C_Psi_bra(:,2)
Psi_tmp_ket(:,1) = C_k(1)*C_Psi_ket(:,1) + C_k(2)*C_Psi_ket(:,2)

do k = 3, order

    k_ref = k

    call bra_x_op( C_Psi_bra(:,k), C_Psi_bra(:,k-1), H_prime, C_two )
    C_Psi_bra(:,k) = C_Psi_bra(:,k) - C_Psi_bra(:,k-2)
    call op_x_ket( C_Psi_ket(:,k), H_prime, C_Psi_ket(:,k-1), C_two )
    C_Psi_ket(:,k) = C_Psi_ket(:,k) - C_Psi_ket(:,k-2)

    Psi_tmp_bra(:,2) = Psi_tmp_bra(:,1) + C_k(k)*C_Psi_bra(:,k)
    Psi_tmp_ket(:,2) = Psi_tmp_ket(:,1) + C_k(k)*C_Psi_ket(:,k)

!   convergence criteria...
    if( isConverged( Psi_tmp_bra(:,2), Psi_tmp_bra(:,1), error ) ) then
    if( isConverged( Psi_tmp_ket(:,2), Psi_tmp_ket(:,1), error ) ) then
        norm_tmp = abs(dotc( Psi_tmp_bra(:,2) , Psi_tmp_ket(:,2) ))
        if( abs( norm_tmp - norm_ref ) < norm_error ) then
            Psi_bra = Psi_tmp_bra(:,2)
            Psi_ket = Psi_tmp_ket(:,2)
            OK = .true.
            exit
        end if
    end if ! ket conv.
    end if ! bra conv.

    Psi_tmp_bra(:,1) = Psi_tmp_bra(:,2)
    Psi_tmp_ket(:,1) = Psi_tmp_ket(:,2)

end do !k = 3, order

deallocate( C_Psi_bra , C_Psi_ket , Psi_tmp_bra , Psi_tmp_ket )

end subroutine Convergence

!
!
!=========================================
 subroutine Build_Hprime( system , basis )
!=========================================
implicit none
type(structure)                 , intent(in)  :: system
type(STO_basis)                 , intent(in)  :: basis(:)

! local variables...
real*8  , allocatable   :: Hamiltonian(:,:)
real*8  , allocatable   :: S(:,:)
integer                 :: i, j, n

n = size(basis)
CALL Overlap_Matrix( system , basis , S )

allocate( Hamiltonian(n,n) )
call GPU_Pin( Hamiltonian, n*n*8)

If( DP_field_ ) then
 
    do j = 1 , n
        do i = 1 , j
     
            Hamiltonian(i,j) = huckel_with_FIELDS(i,j,S(i,j),basis)
            Hamiltonian(j,i) = Hamiltonian(i,j)

        end do
    end do  

else

    do j = 1 , n
        do i = 1 , j
     
            Hamiltonian(i,j) = huckel(i,j,S(i,j),basis)
            Hamiltonian(j,i) = Hamiltonian(i,j)

        end do
    end do  

end If

! compute S_inverse...
call GPU_Pin( S, n*n*8)
call Matrix_syInvert( S, 'U' )

! allocate and compute H' = S_inv * H ...
allocate( H_prime(n,n) )

call syMultiply( 'L', 'U', S, Hamiltonian, H_prime )

call GPU_Unpin(S)
call GPU_Unpin(Hamiltonian)
deallocate( S , Hamiltonian )

end subroutine Build_Hprime
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

!local parameters ...
complex*16  , parameter :: zi_k(0:3) = [ zi , C_one , -zi , -C_one ]

 coefficient(1) = dcmplx( dbesjn(0,tau), 0.d0 )
do k = 2 , k_max
   coefficient(k) = (two * dbesjn(k-1,tau)) * zi_k(mod(k,4))
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
        write(52,12) "#" , QDyn%fragments , "total"
    else
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "unknown", action = "write" , position = "append" )
    end If

    write(52,13) ( QDyn%dyn(it,nf,n) , nf=0,size(QDyn%fragments)+1 ) 

    close(52)

end do

12 FORMAT(10A10)
13 FORMAT(10F10.5)

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
        if( abs(a(i) - b(i) ) >= tol ) return  ! allow earlier return if not convverged
    end do
    isConverged = .true.
end function isConverged

end module Chebyshev_m

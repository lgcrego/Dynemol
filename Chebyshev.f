module Chebyshev_m

    use type_m              , g_time => f_time  
    use constants_m
    use mkl95_blas
    use mkl95_lapack
    use ifport
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis                   
    use QCmodel_Huckel      , only : Huckel 
    use Data_Output         , only : Populations

    public  :: Chebyshev , preprocess_Chebyshev

    private

! module parameters ...
    integer     , parameter :: order     = 15
    real*8      , parameter :: delta_t   = 2.5d0
    real*8      , parameter :: E_range   = 2.0d0
    real*8      , parameter :: error     = 1.0d-6
    real*8      , parameter :: norm_error= 1.0d-6

! module variables ...
    integer , save :: it = 1
    real*8  , save :: save_tau 

contains
!
!
!
!=============================================================
 subroutine preprocess_Chebyshev( system , basis , MO , QDyn )
!=============================================================
implicit none
type(structure)                 , intent(in)    :: system
type(STO_basis)                 , intent(in)    :: basis(:)
complex*16      , allocatable   , intent(out)   :: MO(:)
type(g_time)                    , intent(inout) :: QDyn

!local variables ...
integer                         :: li , N 
real*8          , allocatable   :: wv_FMO(:) 
real*8          , allocatable   :: S_matrix(:,:)
complex*16      , allocatable   :: DUAL_bra(:) , DUAL_ket(:)
type(C_eigen)                   :: FMO

! prepare  DONOR  state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO )

! place the  DONOR  state in structure Hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%fragment == "D" )
N  = size(wv_FMO)
allocate( MO(size(basis)) , source=C_zero )
MO(li:li+N-1) = cmplx( wv_FMO(:) )
deallocate( wv_FMO )

! prepare DUAL basis for local properties ...
allocate( Dual_bra (size(basis)), source=C_zero )
allocate( Dual_ket (size(basis)), source=C_zero )
CALL Overlap_Matrix( system , basis , S_matrix )
DUAL_bra(:) = conjg( MO )
DUAL_ket(:) = matmul( S_matrix , MO )

! save populations(time=t_i) ...
QDyn%dyn(it,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

! clean and exit ...
it = it  + 1

deallocate( DUAL_ket , S_matrix , DUAL_bra )

end subroutine preprocess_Chebyshev
!
!
!
!=========================================================
 subroutine Chebyshev( system , basis , Psi_t , QDyn , t )
!=========================================================
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)
complex*16       , intent(inout) :: Psi_t(:)
type(g_time)     , intent(inout) :: QDyn
real*8           , intent(out)   :: t

! local variables...
complex*16  , allocatable   :: C_Psi(:,:) , Psi_temp(:) , C_k(:) , DUAL_bra(:) , DUAL_ket(:) 
real*8      , allocatable   :: H_prime(:,:) , S_matrix(:,:) 
real*8                      :: tau , inv , norm_ref , norm_test
integer                     :: i , j , k_ref , N
logical                     :: OK

! building  S_matrix  and  H'= S_inv * H ...
CALL Build_Hprime( system , basis , H_prime , S_matrix )

N = size(basis)

allocate( C_Psi    (N , order ) , source=C_zero )
allocate( C_k      (order     ) , source=C_zero )
allocate( Psi_temp (N         ) , source=C_zero )
allocate( Dual_bra (N         ) , source=C_zero )
allocate( Dual_ket (N         ) , source=C_zero )

norm_ref = real( dot_product(Psi_t,matmul(S_matrix,Psi_t)) )

! constants of evolution ...
inv = ( two * h_bar ) / E_range
tau = merge( E_range * delta_t / (two*h_bar) , save_tau * 1.15d0 , it == 2 )

! first convergence: best tau-parameter for k_ref ...
do
    CALL Convergence( Psi_t , k_ref , tau , H_prime , S_matrix , norm_ref , OK )

    if( OK ) then
        save_tau = tau
        exit
    else            
        tau = tau * 0.8d0
    end if
end do

! proceed evolution with best tau ...
C_k = coefficient(tau,order)
do
    C_Psi(:,1) = Psi_t(:)

    C_Psi(:,2) = matmul(H_prime,C_Psi(:,1))

    do j = 3 , k_ref
        C_Psi(:,j) = two * matmul(H_prime,C_Psi(:,j-1)) - C_Psi(:,j-2)
    end do

    forall(j=1:N) Psi_temp(j) = sum( C_k(:) * C_Psi(j,:) )

!   convergence criteria ...
    norm_test = dot_product( Psi_temp(:) , matmul(S_matrix , Psi_temp(:)) )
    if( abs( norm_test - norm_ref ) < norm_error ) then
        Psi_t(:) = Psi_temp(:)
    else
        do
            tau = tau * 0.975d0
            print*, "rescaling tau" , tau 
            CALL Convergence( Psi_t , k_ref , tau , H_prime , S_matrix , norm_ref , OK )
            if( OK ) exit
        end do
    end if

    t = t + (tau * inv)

    if( t >= MD_dt*frame_step*(it-1)  ) exit

end do

! prepare DUAL basis for local properties ...
DUAL_bra(:) = conjg(Psi_t)
DUAL_ket(:) = matmul(S_matrix,Psi_t)

! save populations(time) ...
QDyn%dyn(it,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

! clean and exit ...
it = it  + 1

deallocate( C_k , C_Psi , H_prime , S_matrix , DUAL_bra , DUAL_ket )

Print 186, t

include 'formats.h'

end subroutine Chebyshev
!
!
!
!===============================================================================
subroutine Convergence( Psi , k_ref , tau , H_prime , S_matrix , norm_ref , OK )
!===============================================================================
implicit none
complex*16  , intent(inout) :: Psi(:)
integer     , intent(inout) :: k_ref
real*8      , intent(in)    :: tau
real*8      , intent(in)    :: H_prime(:,:)
real*8      , intent(in)    :: S_matrix(:,:)
real*8      , intent(in)    :: norm_ref
logical     , intent(inout) :: OK

! local variables...
integer                     :: j , l , k , N  
real*8                      :: norm_temp
complex*16  , allocatable   :: C_Psi(:,:) , Psi_temp(:,:) , C_k(:)

N = size(Psi)

allocate( C_Psi    ( N , order ) , source=C_zero )
allocate( C_k      ( order     ) , source=C_zero )
allocate( Psi_temp ( N , 2     ) , source=C_zero )

OK = .false.

! get C_k coefficients ...
C_k = coefficient(tau,order)

C_Psi(:,1) = Psi(:)
C_Psi(:,2) = matmul( H_prime , C_Psi(:,1) )
forall( j=1:N ) Psi_temp(j,1) = sum( C_k(1:2) * C_Psi(j,1:2) )

k = 3
do
    C_Psi(:,k) = two * matmul( H_prime , C_Psi(:,k-1) ) - C_Psi(:,k-2)

    forall( j=1:N ) Psi_temp(j,2) = Psi_temp(j,1) + C_k(k)*C_Psi(j,k)

!   convergence criteria...
    l = count( abs(Psi_temp(:,2)-Psi_temp(:,1)) <= error )
    if( l == N ) then

        norm_temp = dot_product( Psi_temp(:,2) , matmul(S_matrix , Psi_temp(:,2)) )
        if( abs( norm_temp - norm_ref ) < norm_error ) then
            Psi(:) = Psi_temp(:,2)
            OK = .true.
            exit
        end if

    end if

    Psi_temp(:,1) =  Psi_temp(:,2)
    k             =  k + 1

    if( k == order+1 ) exit
end do

k_ref = k
 
deallocate( C_Psi , Psi_temp , C_k )

end subroutine Convergence
!
!
!
!======================================================
subroutine Build_Hprime( system , basis , H_prime , S )
!======================================================
implicit none
type(structure)                 , intent(in)  :: system
type(STO_basis)                 , intent(in)  :: basis(:)
real*8          , allocatable   , intent(out) :: H_prime(:,:)
real*8          , allocatable   , intent(out) :: S(:,:)

! local variables...
real*8  , allocatable   :: Hamiltonian(:,:)
real*8  , allocatable   :: S_inv(:,:)
integer                 :: i , j

CALL Overlap_Matrix( system , basis , S )

allocate( Hamiltonian ( size(basis) , size(basis) ) )

do j = 1 , size(basis)
    do i = 1 , j

        Hamiltonian(i,j) = huckel(i,j,S(i,j),basis)
        Hamiltonian(j,i) = Hamiltonian(i,j)

    end do
end do

! compute S_inverse...
CALL Invertion_Matrix( S , S_inv , size(basis) )

! allocate and compute H' = S_inv * H ...
allocate( H_prime ( size(basis) , size(basis) ) )

H_prime = matmul(S_inv,Hamiltonian)

deallocate( S_inv , Hamiltonian )

end subroutine Build_Hprime
!
!
!
!=====================================================
subroutine Invertion_Matrix( matrix , matrix_inv , N )
!=====================================================
implicit none
real*8                  , intent(in)  :: matrix(:,:)
real*8  , allocatable   , intent(out) :: matrix_inv(:,:)
integer                 , intent(in)  :: N

! local variables...
real*8  , allocatable   :: work(:)
integer , allocatable   :: ipiv(:)
real*8                  :: n_null
integer                 :: i , j , info , sparse

! compute inverse of S_matrix...
allocate( ipiv       ( N     ) )
allocate( work       ( N     ) )
allocate( matrix_inv ( N , N ) )

matrix_inv = matrix

CALL dsytrf( 'u' , N , matrix_inv , N , ipiv , work , N , info )
if ( info /= 0 ) then
    write(*,*) 'info = ',info,' in DSYTRF '
    stop
end if

CALL dsytri( 'u' , N , matrix_inv , N , ipiv , work , info )
if ( info /= 0 ) then
    write(*,*) 'info = ',info,' in DSYTRI '
    stop
end if

deallocate( ipiv , work )

do i = 2 , N
    do j = 1 , i - 1
        matrix_inv(i,j) = matrix_inv(j,i)
    end do
end do

end subroutine Invertion_Matrix
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

coefficient(1) = dbesjn(0,tau)
do k = 2 , k_max
   coefficient(k) = two * zi_k(mod(k,4)) * dbesjn(k-1,tau)
end do

end function coefficient
!
!
!
end module Chebyshev_m

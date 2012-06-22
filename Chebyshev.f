module Chebyshev_m

    use type_m              , g_time => f_time  
    use constants_m
    use parameters_m        , only : t_i , t_f , n_t , frame_step , DP_Field_                     
    use mkl95_blas
    use mkl95_lapack
    use ifport
    use Babel_m             , only : MD_dt
    use Overlap_Builder     , only : Overlap_Matrix
    use FMO_m               , only : FMO_analysis                   
    use QCmodel_Huckel      , only : Huckel ,                       &
                                     Huckel_with_FIELDS
    use Data_Output         , only : Populations

    public  :: Chebyshev , preprocess_Chebyshev

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-12
    real*8      , parameter :: norm_error   = 1.0d-12

! module variables ...
    real*8  , save :: save_tau 

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
Psi(li:li+N-1) = cmplx( wv_FMO(:) )
deallocate( wv_FMO )

! prepare DUAL basis for local properties ...
allocate( Dual_bra (size(basis)), source=C_zero )
allocate( Dual_ket (size(basis)), source=C_zero )
CALL Overlap_Matrix( system , basis , S_matrix )
DUAL_bra(:) = conjg( Psi )
DUAL_ket(:) = matmul( S_matrix , Psi )

allocate( Psi_bra (size(basis)), source=C_zero )
allocate( Psi_ket (size(basis)), source=C_zero )
Psi_bra = matmul( Psi , S_matrix )
Psi_ket = Psi

! save populations(time=t_i) ...
QDyn%dyn(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

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
real*8      , allocatable   :: H_prime(:,:)
real*8                      :: delta_t , tau , tau_max , norm_ref , norm_test
integer                     :: j , k_ref , N
logical                     :: OK

! building  S_matrix  and  H'= S_inv * H ...
CALL Build_Hprime( system , basis , H_prime )

N = size(basis)

allocate( C_Psi_bra   (N , order ) , source=C_zero )
allocate( C_Psi_ket   (N , order ) , source=C_zero )
allocate( C_k         (order     ) , source=C_zero )
allocate( Psi_tmp_bra (N         ) , source=C_zero )
allocate( Psi_tmp_ket (N         ) , source=C_zero )
allocate( Dual_bra    (N         ) , source=C_zero )
allocate( Dual_ket    (N         ) , source=C_zero )

norm_ref = dot_product( Psi_t_bra , Psi_t_ket )
print*, norm_ref

delta_t = merge( (t_f) / float(n_t) , MD_dt * frame_step , MD_dt == epsilon(1.0) )
tau_max = delta_t / h_bar

! constants of evolution ...

! trying to adapt time step for efficient propagation ...
tau = merge( delta_t / h_bar , save_tau * 1.15d0 , it == 2 )
! but tau should be never bigger than tau_max ...
tau = merge( tau_max , tau , tau > tau_max )

! first convergence: best tau-parameter for k_ref ...
do
    CALL Convergence( Psi_t_bra , Psi_t_ket , C_k , k_ref , tau , H_prime , norm_ref , OK )

    if( OK ) then
        t = t + tau * h_bar 
        save_tau = tau
        exit
    else            
        tau = tau * 0.9d0
    end if
end do

if( MD_dt*frame_step*(it-1)-t < tau*h_bar ) then
    tau = ( MD_dt*frame_step*(it-1)-t ) / h_bar
    C_k = coefficient(tau,order)
end if

! proceed evolution with best tau ...
do

    if( t >= MD_dt*frame_step*(it-1) ) exit

    C_Psi_bra(:,1) = Psi_t_bra(:)                                               ;   C_Psi_ket(:,1) = Psi_t_ket(:)

    C_Psi_bra(:,2) = matmul(C_Psi_bra(:,1),H_prime)                             ;   C_Psi_ket(:,2) = matmul(H_prime,C_Psi_ket(:,1))

    do j = 3 , k_ref
        C_Psi_bra(:,j) = two*matmul(C_Psi_bra(:,j-1),H_prime)-C_Psi_bra(:,j-2)  ;   C_Psi_ket(:,j) = two*matmul(H_prime,C_Psi_ket(:,j-1))-C_Psi_ket(:,j-2)
    end do

    forall(j=1:N)
        Psi_tmp_bra(j) = sum( C_k(:) * C_Psi_bra(j,:) )                         ;   Psi_tmp_ket(j) = sum( C_k(:) * C_Psi_ket(j,:) )
    end forall

!   convergence criteria ...
    norm_test = dot_product( Psi_tmp_bra , Psi_tmp_ket )
    if( abs( norm_test - norm_ref ) < norm_error ) then
        Psi_t_bra(:) = Psi_tmp_bra(:)                                           ;   Psi_t_ket(:) = Psi_tmp_ket(:)
    else
        do
            tau = tau * 0.975d0
            print*, "rescaling tau" , tau 
            CALL Convergence( Psi_t_bra , Psi_t_ket , C_k , k_ref , tau , H_prime , norm_ref , OK )
            if( OK ) exit
        end do
    end if

    t = t + (tau * h_bar)

    if( MD_dt*frame_step*(it-1)-t < tau*h_bar ) then
        tau = ( MD_dt*frame_step*(it-1)-t ) / h_bar
        C_k = coefficient(tau,order)
    end if

end do

! prepare DUAL basis for local properties ...
DUAL_bra = conjg(Psi_t_ket)
DUAL_ket = Psi_t_bra

! save populations(time) ...
QDyn%dyn(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

! clean and exit ...
deallocate( C_k , C_Psi_bra , C_Psi_ket , Psi_tmp_bra , Psi_tmp_ket , H_prime , DUAL_bra , DUAL_ket )

Print 186, t

include 'formats.h'

end subroutine Chebyshev
!
!
!
!===================================================================================================
subroutine Convergence( Psi_bra , Psi_ket , C_k , k_ref , tau , H_prime , norm_ref , OK )
!===================================================================================================
implicit none
complex*16  , intent(inout) :: Psi_bra(:)
complex*16  , intent(inout) :: Psi_ket(:)
complex*16  , intent(out)   :: C_k(:)
integer     , intent(inout) :: k_ref
real*8      , intent(in)    :: tau
real*8      , intent(in)    :: H_prime(:,:)
real*8      , intent(in)    :: norm_ref
logical     , intent(inout) :: OK

! local variables...
integer                     :: j , l , k , N  
real*8                      :: norm_tmp
complex*16  , allocatable   :: C_Psi_bra(:,:) , C_Psi_ket(:,:) , Psi_tmp_bra(:,:) , Psi_tmp_ket(:,:)

N = size(Psi_bra)

allocate( C_Psi_bra    ( N , order ) , source=C_zero )
allocate( C_Psi_ket    ( N , order ) , source=C_zero )
allocate( Psi_tmp_bra  ( N , 2     ) , source=C_zero )
allocate( Psi_tmp_ket  ( N , 2     ) , source=C_zero )

OK = .false.

! get C_k coefficients ...
C_k = coefficient(tau,order)

C_Psi_bra(:,1) = Psi_bra(:)                                                 ;   C_Psi_ket(:,1) = Psi_ket(:)
C_Psi_bra(:,2) = matmul(C_Psi_bra(:,1),H_prime)                             ;   C_Psi_ket(:,2) = matmul(H_prime,C_Psi_ket(:,1))

forall( j=1:N )
    Psi_tmp_bra(j,1) = sum( C_k(1:2) * C_Psi_bra(j,1:2) )                   ;   Psi_tmp_ket(j,1) = sum( C_k(1:2) * C_Psi_ket(j,1:2) )
end forall

k = 3
do

    C_Psi_bra(:,k) = two*matmul(C_Psi_bra(:,k-1),H_prime)-C_Psi_bra(:,k-2)  ;   C_Psi_ket(:,k) = two*matmul(H_prime,C_Psi_ket(:,k-1))-C_Psi_ket(:,k-2)

    forall( j=1:N )
        Psi_tmp_bra(j,2) = Psi_tmp_bra(j,1) + C_k(k)*C_Psi_bra(j,k)         ;   Psi_tmp_ket(j,2) = Psi_tmp_ket(j,1) + C_k(k)*C_Psi_ket(j,k)
    end forall

!   convergence criteria...
    l = count( abs( Psi_tmp_bra(:,2) - Psi_tmp_bra(:,1) ) <= error ) + &
        count( abs( Psi_tmp_ket(:,2) - Psi_tmp_ket(:,1) ) <= error )

    if( l == 2*N ) then
        norm_tmp = dot_product( Psi_tmp_bra(:,2) , Psi_tmp_ket(:,2) )
        if( abs( norm_tmp - norm_ref ) < norm_error ) then
            Psi_bra(:) = Psi_tmp_bra(:,2)                                   ;   Psi_ket(:) = Psi_tmp_ket(:,2)
            OK = .true.
            exit
        end if
    end if

    Psi_tmp_bra(:,1) =  Psi_tmp_bra(:,2)                                    ;   Psi_tmp_ket(:,1) =  Psi_tmp_ket(:,2)

    k = k + 1

    if( k == order+1 ) exit

end do

k_ref = k

deallocate( C_Psi_bra , C_Psi_ket , Psi_tmp_bra , Psi_tmp_ket )

end subroutine Convergence
!
!
!
!==================================================
subroutine Build_Hprime( system , basis , H_prime )
!==================================================
implicit none
type(structure)                 , intent(in)  :: system
type(STO_basis)                 , intent(in)  :: basis(:)
real*8          , allocatable   , intent(out) :: H_prime(:,:)

! local variables...
real*8  , allocatable   :: Hamiltonian(:,:)
real*8  , allocatable   :: S(:,:) , S_inv(:,:)
integer                 :: i , j

CALL Overlap_Matrix( system , basis , S )

allocate( Hamiltonian ( size(basis) , size(basis) ) )

If( DP_field_ ) then
 
    do j = 1 , size(basis)
        do i = 1 , j
     
            Hamiltonian(i,j) = huckel_with_FIELDS(i,j,S(i,j),basis)
            Hamiltonian(j,i) = Hamiltonian(i,j)

        end do
    end do  

else

    do j = 1 , size(basis)
        do i = 1 , j
     
            Hamiltonian(i,j) = huckel(i,j,S(i,j),basis)
            Hamiltonian(j,i) = Hamiltonian(i,j)

        end do
    end do  

end If

! compute S_inverse...
CALL Invertion_Matrix( S , S_inv , size(basis) )

! allocate and compute H' = S_inv * H ...
allocate( H_prime ( size(basis) , size(basis) ) )

H_prime = matmul(S_inv,Hamiltonian)

deallocate( S , S_inv , Hamiltonian )

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

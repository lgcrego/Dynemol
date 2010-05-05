module Chebyshev_m

    use type_m
    use constants_m
    use Overlap_Builder
    use QCmodel_Huckel
    use mkl95_blas
    use mkl95_lapack
    use ifport

    public  :: Chebyshev

    private

    integer     , parameter :: order     = 15
    real*8      , parameter :: delta_t   = 2.5d0
    real*8      , parameter :: E_range   = 2.0d0
    real*8      , parameter :: error     = 5.0d-6
    real*8      , parameter :: norm      = 1.0d-4
    integer     , parameter :: step_max  = 2300
    complex*16  , parameter :: z0        = ( 0.0d0 , 0.0d0 )

contains
!
!=================================================================================
subroutine Chebyshev( system , basis , MO , ind , time_init , tau_init , P1 , P2 )
!=================================================================================
implicit none
type(structure)                 , intent(in)    :: system
type(STO_basis)                 , intent(in)    :: basis(:)
complex*16                      , intent(inout) :: MO(:)
integer                         , intent(inout) :: ind
real*8                          , intent(inout) :: time_init
real*8                          , intent(inout) :: tau_init
integer                         , intent(inout) :: P1
integer                         , intent(inout) :: P2

! local variables...
complex*16  , allocatable   :: C_Psi(:,:)
complex*16  , allocatable   :: Psi_temp(:)
complex*16  , allocatable   :: Psi_T(:)
complex*16  , allocatable   :: C_k(:)
complex*16  , allocatable   :: A_r(:)
real*8      , allocatable   :: Hamiltonian(:,:)
real*8      , allocatable   :: S_matrix(:,:)
real*8                      :: tau , inv , ab , time_t , norm_parameter
integer                     :: i   , j   , K  , info   , counter        , N

! calcula of S_matrix, H, H' and the state to be propagate...
CALL INIT( system , basis , Hamiltonian , S_matrix )

N = size(ExCell_basis)

allocate( C_Psi    ( N , order ) )
allocate( Psi_temp ( N         ) )
allocate( Psi_T    ( N         ) )
allocate( C_k      ( order     ) )
allocate( A_r      ( N         ) )

if( ind == 1 ) time_init = 0.0d0

! defined some constants of evolution...
inv    = ( 2.0d0 * h_bar ) / E_range
Psi_T  = MO
time_t = time_init

norm_parameter = real( dot_product(Psi_T,matmul(S_matrix,Psi_T)) )

! file of data of the propagation...
if( ind == 1 ) then
    open(unit=11, file='population.dat', action='write', status='unknown')
    write(11,*), 0.0d0 , 1.0d0
else
    open(unit=11, file='population.dat', action='write', status='old', position='append')
end if

if( ind == 1 ) then
    tau = E_range * delta_t / (2.0d0*h_bar)
    do

        C_k(1) = dbesjn(0,tau)
        do j = 2 , order
            C_k(j) = 2.0d0 * (-zi)**(j-1) * dbesjn(j-1,tau)
        end do

        CALL Convergence( Psi_T , K , C_k , Hamiltonian , S_matrix , N , norm_parameter , info )

        if( info == 0 ) then
            A_r = matmul(S_matrix,Psi_T)
            write(11,*), time_t + tau*inv , real( dot_product( Psi_T(P1:P2) , A_r(P1:P2) ) )
            tau_init = tau
            exit
        end if

        tau = tau * 0.9825d0

        if( tau == 0.0d0 ) then
            print'("***** tau = 0.0d0 *****")'
            stop
        end if

    end do
else
    tau = tau_init * 1.15d0
    do

        C_k(1) = dbesjn(0,tau)
        do j = 2 , order
            C_k(j) = 2.0d0 * (-zi)**(j-1) * dbesjn(j-1,tau)
        end do

        CALL Convergence( Psi_T , K , C_k , Hamiltonian , S_matrix , N , norm_parameter , info )

        if( info == 0 ) then
            A_r = matmul(S_matrix,Psi_T)
            write(11,*), time_t + tau*inv , real( dot_product( Psi_T(P1:P2) , A_r(P1:P2) ) )
            tau_init = tau
            exit
        end if

        tau = tau * 0.9825d0

        if( tau == 0.0d0 ) then
            print'("***** tau = 0.0d0 *****")'
            stop
        end if

    end do
end if

counter  = 1

do

    if( counter >= step_max .AND. tau /= tau_init ) then
        counter = 0
        tau = tau_init
        C_k(1) = dbesjn(0,tau)
        do j = 2 , order
            C_k(j) = 2.0d0 * (-zi)**(j-1) * dbesjn(j-1,tau)
        end do
    end if
    
    C_Psi = z0

    C_Psi(:,1) = Psi_T(:)

    C_Psi(:,2) = matmul(Hamiltonian,C_Psi(:,1))

    do j = 3 , K
        C_Psi(:,j) = 2.0d0 * matmul(Hamiltonian,C_Psi(:,j-1)) - C_Psi(:,j-2)
    end do

    forall(j=1:N) Psi_temp(j) = sum( C_k(:) * C_Psi(j,:) )

!   convergence criteria...
    ab = sqrt( (sqrt( abs(dot_product(Psi_temp(:),matmul(S_matrix,Psi_temp(:)))**2) ) - norm_parameter)**2 )

    if( ab < norm ) then

        Psi_T(:) = Psi_temp(:)
        A_r = matmul(S_matrix,Psi_T)
        write(11,*), time_t + tau*inv , real( dot_product( Psi_T(P1:P2) , A_r(P1:P2) ) )
        Psi_temp = 0.0d0
        goto 30

    else

        do
            tau = tau * 0.975d0
            print*, "reallocate tau"
            if( tau == 0.0d0 ) then
                print*, "***** tau = 0.0d0 ****"
                stop
            end if
            C_k(1) = dbesjn(0,tau)
            do j = 2 , order
                C_k(j) = 2.0d0 * (-zi)**(j-1) * dbesjn(j-1,tau)
            end do
            CALL Convergence( Psi_T , K , C_k , Hamiltonian , S_matrix , N , norm_parameter , info )
            if( info == 0 ) then
                A_r = matmul(S_matrix,Psi_T)
                write(11,*), time_t + tau*inv , real( dot_product( Psi_T(P1:P2) , A_r(P1:P2) ) )
                goto 30
            end if
        end do

    end if

30  time_t = time_t + tau * inv

    if( time_t >= MD_dt*frame_step*ind ) exit
    
    counter = counter + 1

end do

close(11)

print*, "Evolution time = ", time_t*1.0d3 , "fs"

if( time_t >= t_f ) stop

MO        = Psi_T
ind       = ind + 1
time_init = time_t

deallocate( Psi_T , C_k , C_Psi , A_r , Hamiltonian , S_matrix )

end subroutine Chebyshev
!
!
!
!======================================================================================
subroutine Convergence( Psi , K , C_k , Hamiltonian , S_m , N , norm_parameter , info )
!======================================================================================
implicit none
complex*16  , intent(inout) :: Psi(:)
integer     , intent(inout) :: K
complex*16  , intent(in)    :: C_k(:)
real*8      , intent(in)    :: Hamiltonian(:,:)
real*8      , intent(in)    :: S_m(:,:)
integer     , intent(in)    :: N
real*8      , intent(in)    :: norm_parameter
integer     , intent(inout) :: info

! local variables...
complex*16  , allocatable   :: C_Psi(:,:) , Psi_temp(:,:)
complex*16  , allocatable   :: Temp(:)
real*8                      :: aa , bc
integer                     :: j , l

allocate( C_Psi    ( N , order ) )
allocate( Psi_temp ( N , 2     ) )
allocate( Temp     ( N         ) )

info = 1

C_Psi(:,1) = Psi(:)
C_Psi(:,2) = matmul(Hamiltonian,C_Psi(:,1))

forall(j=1:N) Psi_temp(j,1) = sum( C_k(:) * C_Psi(j,:) )

K = 2
K = K + 1

do

    C_Psi(:,K) = 2.0d0*matmul(Hamiltonian,C_Psi(:,j-1)) - C_Psi(:,K-2)

    forall(j=1:N) Psi_temp(j,2) = Psi_temp(j,1) + C_k(K) * C_Psi(j,K)

    l = 0
    do j = 1 , N
        aa = abs(sqrt( (Psi_temp(j,2) - Psi_temp(j,1))**2 )**2)
        if(  aa <= error ) l = l + 1
    end do

!   convergence criteria...
    if( l == N ) then
        bc = sqrt( (sqrt( abs(dot_product(Psi_temp(:,2),matmul(S_m,Psi_temp(:,2)))**2) ) - norm_parameter)**2 )
        if( bc < norm ) then
            Psi(:) = Psi_temp(:,2)
            info = 0
            exit
        end if

    end if

    K = K + 1

    if( K == order+1 ) then
        exit
    end if

    Psi_temp(:,1) = Psi_temp(:,2)

end do

deallocate( C_Psi , Psi_temp , Temp )

end subroutine Convergence
!
!
!
!========================================
subroutine INIT( system , basis , H , S )
!========================================
implicit none
type(structure)                 , intent(in)    :: system
type(STO_basis)                 , intent(in)    :: basis(:)
real*8          , allocatable   , intent(inout) :: H(:,:)
real*8          , allocatable   , intent(inout) :: S(:,:)

! local variables...
real*8  , allocatable   :: Hamiltonian(:,:)
real*8  , allocatable   :: S_inv(:,:)
integer                 :: i , j

CALL Overlap_Matrix( system , basis , S )

allocate( Hamiltonian ( size(basis) , size(basis) ) )

do j = 1 , size(basis)
    do i = 1 , j
        Hamiltonian(i,j) = huckel(i,j,S(i,j),basis)
    end do
end do

do j = 1 , size(basis) - 1
    do i = j + 1 , size(basis)
        Hamiltonian(i,j) = Hamiltonian(j,i)
    end do
end do

! compute S_inverse...
CALL Invertion_Matrix( S , S_inv , size(basis) )

! compute H'...
allocate( H ( size(basis) , size(basis) ) )

H = matmul(S_inv,Hamiltonian)

deallocate( S_inv )

end subroutine INIT
!
!
!
!=====================================================
subroutine Invertion_Matrix( matrix , matrix_inv , N )
!=====================================================
implicit none
real*8                  , intent(in)    :: matrix(:,:)
real*8  , allocatable   , intent(inout) :: matrix_inv(:,:)
integer                 , intent(in)    :: N

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

end module Chebyshev_m

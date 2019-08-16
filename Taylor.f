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
    use Data_Output         , only : Populations
    use Matrix_Math

    public  :: Propagation, dump_Qdyn

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

    If( eh_tag(n) == "XX" ) cycle

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


end module Taylor_m




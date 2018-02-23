#ifdef USE_GPU
module gpu_ElHl_Chebyshev_m

    use MPI
    use type_m              , g_time => f_time
    use blas95
    use lapack95
    use constants_m
    use ifport
    use MPI_definitions_m         , only : myCheby, ChebyCrew, ChebyComm, KernelComm, ChebyKernelComm, master
    use parameters_m              , only : t_i, frame_step, n_part, restart, QMMM
    use Structure_Builder         , only : Unit_Cell
    use Overlap_Builder           , only : Overlap_Matrix
    use FMO_m                     , only : FMO_analysis, eh_tag
    use Data_Output               , only : Populations
    use Chebyshev_m               , only : dump_Qdyn
    use Matrix_Math
    use execution_time_m

    public  :: gpu_ElHl_Chebyshev, gpu_preprocess_ElHl_Chebyshev

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-12
    real*8      , parameter :: norm_error   = 1.0d-12

! module variables ...
    real*8      ,   save          :: save_tau
    logical     ,   save          :: ready = .false.
    real*8      ,   allocatable   :: h0(:,:), H_prime(:,:)
    real*8      ,   allocatable   :: S_matrix(:,:)
    complex*16  ,   allocatable   :: Psi_t_bra(:,:) , Psi_t_ket(:,:)

contains
!
!
!
!==============================================================================================================
 subroutine gpu_preprocess_ElHl_Chebyshev( system , basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , it )
!==============================================================================================================
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
integer                         :: li , n , n_FMO
real*8          , allocatable   :: wv_FMO(:)
complex*16      , allocatable   :: ElHl_Psi(:,:)
type(R_eigen)                   :: FMO


n = size(basis)

! compute S before FMO_analysis (important) ...
call Overlap_Matrix( system , basis , S_matrix )

! after instantiating Overlap_matrix, processes wait outside ...
if (.not. ChebyCrew) return

! calculate Hamiltonian
allocate( h0(n,n), H_prime(n,n) )
call Huckel( basis , S_matrix , h0 )

! pin (for faster CPU <-> GPU transfers)
call GPU_Pin( S_matrix, n*n*8 )
call GPU_Pin( h0,       n*n*8 )
call GPU_Pin( H_prime,  n*n*8 )

!========================================================================
! prepare electron state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="E" )

! place the electron state in Structure's hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%El )
n_FMO = size(wv_FMO)
allocate( ElHl_Psi( n, n_part ) , source=C_zero )
ElHl_Psi(li:li+n_FMO-1,1) = dcmplx( wv_FMO(:) )
deallocate( wv_FMO )

!========================================================================
! prepare hole state ...
CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="H" )

! place the hole state in Structure's hilbert space ...
li = minloc( basis%indx , DIM = 1 , MASK = basis%Hl )
n_FMO = size(wv_FMO)
ElHl_Psi(li:li+n_FMO-1,2) = dcmplx( wv_FMO(:) )
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
allocate( Psi_t_bra( n, n_part ) )
allocate( Psi_t_ket( n, n_part ) )
call bra_x_op( Psi_t_bra, ElHl_Psi , S_matrix )
Psi_t_ket = ElHl_Psi
!==============================================

!==============================================
AO_bra = ElHl_Psi
AO_ket = ElHl_Psi
CALL QuasiParticleEnergies(AO_bra, AO_ket, H0)
!==============================================

If( .not. restart ) then
    ! save populations(time=t_i) ...
    QDyn%dyn(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )
    if (master) call dump_Qdyn( Qdyn , it )
end If

end subroutine gpu_preprocess_ElHl_Chebyshev
!
!
!
!=================================================================================================================
 subroutine gpu_ElHl_Chebyshev( system , basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , t , delta_t , it )
!=================================================================================================================
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
logical, save :: first_call = .true.
logical :: do_electron, do_hole
integer :: j, N, my_part, err, mpi_status(mpi_status_size), request_H_prime, t_stuff_sent, coords_sent
integer :: mpi_D_R = mpi_double_precision
integer :: mpi_D_C = mpi_double_complex
real*8  :: t_init , t_max , tau_max , tau , t_stuff(2)


N = size(basis)

t_init = t
! max time inside slice ...
t_max = delta_t*frame_step*(it-1)
! constants of evolution ...
tau_max = delta_t / h_bar


!----------------------------------
! Organize which process does what:
! (master) <-> (myCheby==0) process: electron
 do_electron = master
! (myCheby==1) process: hole
 do_hole = merge( .true., .false., (myCheby==1) )
! check if process doesn't belong here
 if (.not. (do_electron .or. do_hole)) return
! define index to address electron (1) or hole (2) in Psi_t_bra/ket
 if (do_electron) my_part = 1
 if (do_hole)     my_part = 2
!----------------------------------


! begin
100 continue

! compute overlap and hamiltonian
if(.not. first_call) then    ! already calculated in preprocess

    if (do_electron) then
        call MPI_Isend( system%coord, system%atoms*3, mpi_D_R, 1, 0, ChebyComm, coords_sent, err )
        call MPI_Request_Free( coords_sent, err )
    else
        call MPI_Recv( system%coord, system%atoms*3, mpi_D_R, 0, 0, ChebyComm, mpi_status, err )
    end if

    call Overlap_Matrix( system, basis, S_matrix )
    call Huckel( basis, S_matrix, h0 )

end if


! master bcasts QMMM
call MPI_Bcast( QMMM, 1, mpi_logical, 0, ChebyComm, err )


! master/slave sends/receives time info
if (do_electron) then
    t_stuff = [t_init, t_max]
    call MPI_Isend( t_stuff, 2, mpi_D_R, 1, 0, ChebyComm, t_stuff_sent, err )          !  ──┐
    call MPI_Request_Free( t_stuff_sent, err )                                         !    │
else if (do_hole) then                                                                 !    │
    call MPI_Recv( t_stuff, 2, mpi_D_R, 0, 0, ChebyComm, mpi_status, err )             ! <──┘
    t_init = t_stuff(1)
    t_max  = t_stuff(2)
end if


! trying to adapt time step for efficient propagation ...
tau = merge( tau_max, save_tau * 1.15d0, first_call )
! but tau should be never bigger than tau_max ...
tau = merge( tau_max, tau, tau > tau_max )

first_call = .false.


! send AO_bra/ket to Diabatic_Ehrenfest
if (do_electron .and. QMMM) then
    call MPI_Bcast( AO_bra, 2*N, mpi_D_C, 0, KernelComm, err )
    call MPI_Bcast( AO_ket, 2*N, mpi_D_C, 0, KernelComm, err )
end if


! propagate electron or hole
call PropagationElHl_gpucaller( N, S_matrix, h0, H_prime,                   &
                                AO_bra   (:,my_part), AO_ket   (:,my_part), &
                                Psi_t_bra(:,my_part), Psi_t_ket(:,my_part), &
                                t_init, t_max, tau, save_tau )

! remark: on return, AO_bra = S_inv * Psi_t_bra. AO_bra still lacks conjugation


! send/receive hole
if (do_hole) then
    call MPI_Send( Psi_t_bra(:,2), N, mpi_D_C, 0, 0, ChebyComm, err )
    call MPI_Send( Psi_t_ket(:,2), N, mpi_D_C, 0, 0, ChebyComm, err )
    call MPI_Send( AO_bra(:,2),    N, mpi_D_C, 0, 0, ChebyComm, err )
else if (do_electron) then
    call MPI_Recv( Psi_t_bra(:,2), N, mpi_D_C, 1, 0, ChebyComm, mpi_status, err )
    call MPI_Recv( Psi_t_ket(:,2), N, mpi_D_C, 1, 0, ChebyComm, mpi_status, err )
    call MPI_Recv( AO_bra(:,2),    N, mpi_D_C, 1, 0, ChebyComm, mpi_status, err )
end if


! (myCheb==1) process keeps going back forever
if (do_hole) goto 100 !begin


! only master (do_electron) reaches this part of the code
! send H_prime to process (myChebyKernel==2)
if (QMMM) call MPI_Isend( H_prime, N*N, mpi_D_R, 2, 0, ChebyKernelComm, request_H_prime, err )

t = t_init + (delta_t*frame_step)

! prepare DUAL basis for local properties ...
DUAL_bra = dconjg(Psi_t_ket)
DUAL_ket = Psi_t_bra

! prepare Slater basis for FORCE properties ...
AO_bra = dconjg(AO_bra)   ! conjugation of AO_bra is not done on the GPU
AO_ket = Psi_t_ket

call QuasiParticleEnergies( AO_bra, AO_ket, h0 )

! save populations(time) ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments, basis, DUAL_bra, DUAL_ket, t )
call dump_Qdyn( Qdyn , it )

call MPI_Wait( request_H_prime, mpi_status, err )
call stop_clock("cheby_gpu")

end subroutine gpu_ElHl_Chebyshev
!
!
!
!=========================================
subroutine Huckel( basis , S_matrix , h )
!=========================================
implicit none
type(STO_basis), intent(in)  :: basis(:)
real*8,          intent(in)  :: S_matrix(:,:)
real*8,          intent(out) :: h(:,:)

! local variables ...
real*8  :: k_eff, k_WH, c1, c2, c3, basis_j_IP, basis_j_k_WH
integer :: i, j, n

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

n = size(basis)

!$omp parallel do private(i,j,basis_j_IP,basis_j_k_WH,c1,c2,c3,k_WH,k_eff) default(shared) schedule(dynamic,1)
do j = 1, n

    basis_j_IP   = basis(j)%IP
    basis_j_k_WH = basis(j)%k_WH

    do i = 1, j-1

        c1 = basis(i)%IP - basis_j_IP
        c2 = basis(i)%IP + basis_j_IP
        c3 = (c1/c2)**2

        k_WH = (basis(i)%k_WH + basis_j_k_WH) / two

        k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

        h(i,j) = k_eff * S_matrix(i,j) * c2 / two

    end do

    h(j,j) = basis_j_IP

end do
!$omp end parallel do

call Matrix_Symmetrize( h, 'U' )     

end subroutine Huckel
!
!
!
!=======================================================
 subroutine QuasiParticleEnergies( AO_bra , AO_ket , H )
!=======================================================
implicit none
complex*16, intent(in) :: AO_bra(:,:)
complex*16, intent(in) :: AO_ket(:,:)
real*8,     intent(in) :: H(:,:)

!local variables ...
integer    :: i, j, n
complex*16 :: ket(2), erg(2)


n = size(AO_bra, 1)

erg = (0.d0, 0.d0)

!$omp parallel do private(i,j,ket) default(shared) reduction(+: erg)
do j = 1, n
    ket = AO_ket(j,:)
    do i = 1, n
        erg = erg + AO_bra(i,:)*H(i,j)*ket
    end do
end do
!$omp end parallel do

Unit_Cell% QM_wp_erg = erg

! QM_erg = E_occ - E_empty ; to be used in MM_dynamics energy balance ...
Unit_Cell% QM_erg = erg(1) - erg(2)

end subroutine QuasiParticleEnergies
!
!
!
!
end module gpu_ElHl_Chebyshev_m
#endif

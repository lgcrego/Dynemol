#ifdef USE_GPU
module gpu_ElHl_Chebyshev_m

    use MPI
    use type_m             , g_time => f_time  
    use blas95
    use lapack95
    use constants_m
    use ifport
    use Matrix_Math
    use MPI_definitions_m  , only : master , EnvCrew ,  myCheby ,   &       
                                    ChebyComm , ChebyCrew ,         &  
                                    ForceComm , ForceCrew ,         &
                                    KernelComm , ChebyKernelComm   
    use parameters_m       , only : t_i , frame_step , Coulomb_ ,   &
                                    EnvField_ , n_part, driver ,    &
                                    QMMM, CT_dump_step , HFP_Forces
    use Structure_Builder  , only : Unit_Cell 
    use Overlap_Builder    , only : Overlap_Matrix
    use FMO_m              , only : FMO_analysis , eh_tag    
    use Data_Output        , only : Populations 
    use Hamiltonians       , only : X_ij , even_more_extended_Huckel
    use Taylor_m           , only : Propagation, dump_Qdyn
    use ElHl_Chebyshev_m   , only : Build_Huckel ,                  &
                                    QuasiParticleEnergies ,         &
                                    preprocess_from_restart

    public  :: gpu_ElHl_Chebyshev, gpu_preprocess_ElHl_Chebyshev

    private

! module parameters ...
    integer     , parameter :: order        = 25
    real*8      , parameter :: error        = 1.0d-12
    real*8      , parameter :: norm_error   = 1.0d-12

! module variables ...
    logical     ,   save          :: first_call_ = .true.
    real*8      ,   save          :: save_tau
    real*8      ,   allocatable   :: h0(:,:)
    real*8      ,   allocatable   :: H_prime(:,:)
    real*8      ,   allocatable   :: S_matrix(:,:)
    complex*16  ,   allocatable   :: Psi_t_bra(:,:) , Psi_t_ket(:,:)

    interface gpu_preprocess_ElHl_Chebyshev
        module procedure gpu_preprocess_ElHl_Chebyshev
        module procedure preprocess_from_restart
    end interface

contains
!
!
!
!==============================================================================================================
 subroutine gpu_preprocess_ElHl_Chebyshev( system , basis , AO_bra , AO_ket , Dual_bra , Dual_ket , QDyn , it )
! used for normal start (see interface) ...
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
integer                         :: li , M , N , err 
integer                         :: mpi_D_R = mpi_double_precision
real*8          , allocatable   :: wv_FMO(:)
complex*16      , allocatable   :: ElHl_Psi(:,:)
type(R_eigen)                   :: FMO


! MUST compute S_matrix before FMO analysis ...
CALL Overlap_Matrix( system , basis , S_matrix )

N = size(basis)

allocate( h0(N,N), H_prime(N,N) )

!------------------------------------------------------------------------
If( master .OR. EnvCrew ) then 

   If( EnvField_ ) then
      ! EnCrew stay in even_more_extended_Huckel ...
      h0(:,:) = even_more_extended_Huckel( system , basis , S_matrix )
   Else
      h0(:,:) = Build_Huckel( basis , S_matrix )
   end If

   CALL MPI_BCAST( h0 , N*N , mpi_D_R , 0 , ForceComm , err )
   CALL MPI_BCAST( h0 , N*N , mpi_D_R , 0 , ChebyComm , err )

End If

If( ForceCrew ) then
   ! After instantiating S_matrix AND h0, ForceCrew leave to EhrenfestForce ...
   CALL MPI_BCAST( h0 , N*N , mpi_D_R , 0 , ForceComm , err )
   return
End If

! Another ChebyCrew  mate ...
If( myCheby == 1 ) CALL MPI_BCAST( h0 , N*N , mpi_D_R , 0 , ChebyComm , err )
!------------------------------------------------------------------------

!========================================================================

! pin (for faster CPU <-> GPU transfers)
call GPU_Pin( S_matrix, n*n*8 )
call GPU_Pin( h0,       n*n*8 )
call GPU_Pin( H_prime,  n*n*8 )

allocate( ElHl_Psi( N , n_part ) , source=C_zero )
!========================================================================
! prepare electron state ...
  CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="E" )

  ! place the electron state in Structure's Hilbert space ...
  li = minloc( basis%indx , DIM = 1 , MASK = basis%El )
  M  = size(wv_FMO)
  ElHl_Psi(li:li+M-1,1) = merge( dcmplx(wv_FMO(:)) , C_zero , eh_tag(1) == "el")
  deallocate( wv_FMO )
!========================================================================
! prepare hole state ...
  CALL FMO_analysis( system , basis, FMO=FMO , MO=wv_FMO , instance="H" )

  ! place the hole state in Structure's Hilbert space ...
  li = minloc( basis%indx , DIM = 1 , MASK = basis%Hl )
  M  = size(wv_FMO)
  ElHl_Psi(li:li+M-1,2) = merge( dcmplx(wv_FMO(:)) , C_zero , eh_tag(2) == "hl")
  deallocate( wv_FMO )
!========================================================================

!==============================================
! prepare DUAL basis for local properties ...
! DUAL_bra = (C*)^T    ;    DUAL_ket = S*C ...
  DUAL_bra = dconjg( ElHl_Psi )
  call op_x_ket( DUAL_ket, S_matrix , ElHl_Psi )
!==============================================

!==============================================
! vector states to be propagated ...
! Psi_bra = C^T*S       ;      Psi_ket = C ...
  allocate( Psi_t_bra(N,n_part) )
  allocate( Psi_t_ket(N,n_part) )
  call bra_x_op( Psi_t_bra, ElHl_Psi , S_matrix ) 
  Psi_t_ket = ElHl_Psi
!==============================================

!==============================================
! preprocess stuff for EhrenfestForce ...
  AO_bra = ElHl_Psi 
  AO_ket = ElHl_Psi 
  CALL QuasiParticleEnergies(AO_bra, AO_ket, H0)
!==============================================

! save populations(time=t_i) ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )
if (master) CALL dump_Qdyn( Qdyn , it )

! leaving S_matrix allocated

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
integer :: N, my_part, err, mpi_status(mpi_status_size), request_H_prime, t_stuff_sent, coords_sent, it_sync
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
10 continue

! compute overlap and hamiltonian
if(.not. first_call) then  

    if (do_electron) then
        call MPI_Isend( system%coord, system%atoms*3, mpi_D_R, 1, 0, ChebyComm, coords_sent, err )
        call MPI_Request_Free( coords_sent, err )
    else
        call MPI_Recv( system%coord, system%atoms*3, mpi_D_R, 0, 0, ChebyComm, mpi_status, err )
    end if

    call Overlap_Matrix( system, basis, S_matrix )

    If( Envfield_ ) then
        it_sync = it-1  ! <== for synchronizing EnvSetUp call therein ...
        h0(:,:) = even_more_extended_Huckel( system , basis , S_matrix , it_sync)
    else
        h0(:,:) = Build_Huckel( basis , S_matrix )
    end If

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
if (do_hole) goto 10 !begin


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

if( mod(it,CT_dump_step) == 0 ) CALL dump_Qdyn( Qdyn , it )

call MPI_Wait( request_H_prime, mpi_status, err )

include 'formats.h'

end subroutine gpu_ElHl_Chebyshev
!
!
!
end module gpu_ElHl_Chebyshev_m
#endif

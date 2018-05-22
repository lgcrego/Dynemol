! Program for computing Ehrenfest forces from Huckel Hamiltonian
module DiabaticEhrenfest_Builder

    use MPI
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only  : verbose 
    use MPI_definitions_m       , only  : myForce , master , npForce ,           &
                                          KernelComm , ForceComm , ForceCrew ,   &
                                          myKernel , ChebyKernelComm 
    use Overlap_Builder         , only  : Overlap_Matrix

    public :: EhrenfestForce 

    private

    interface EhrenfestForce 
        module procedure Diabatic_Ehrenfest
    end interface EhrenfestForce 

    !module variables ...
    integer     , allocatable   :: BasisPointer(:) , DOS(:)
    real*8      , allocatable   :: rho_eh(:,:) , A_ad_nd(:,:) , B_ad_nd(:,:) , Kernel(:,:) 
    real*8      , allocatable   :: grad_S(:,:) , F_vec(:) , F_snd(:,:,:) , F_rcv(:,:,:) , H_prime(:,:)
    logical     , allocatable   :: mask(:,:)

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]
    real*8  , parameter :: delta      = 1.d-8
    
contains
!
!
!=================================================================
 subroutine Diabatic_Ehrenfest( system , basis , AO_bra , AO_ket )
!=================================================================
 use MD_read_m   , only  : atom
 implicit none
 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 complex*16      , optional , intent(in)    :: AO_bra(:,:)
 complex*16      , optional , intent(in)    :: AO_ket(:,:)

! local variables ... 
 integer :: i , j , N , xyz , err , request_H_prime, Kernel_bcast
 integer :: mpi_status(mpi_status_size)
 integer :: mpi_D_R = mpi_double_precision
 integer :: mpi_D_C = mpi_double_complex

 real*8  , allocatable :: X_ij(:,:) 

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 logical , parameter :: T_ = .true. , F_ = .false.
 

!================================================================================================
! some preprocessing ...
!================================================================================================

N = size(basis)

If( .NOT. allocated(Kernel) ) then
        allocate( Kernel (N,N)                         )
        allocate( grad_S (N,10)                        )
        allocate( F_snd  (system%atoms,system%atoms,3) )
        allocate( F_rcv  (system%atoms,system%atoms,3) )
        allocate( F_vec  (system%atoms)                ) 
end If

! ForceCrew in stand-by ...
99 If( .not. master ) CALL MPI_BCAST( system%coord , system%atoms*3 , mpi_D_R , 0 , ForceComm, err )

! preprocess overlap matrix for Pulay calculations, all ForceCrew must do it ...
 CALL Overlap_Matrix( system , basis )
 CALL preprocess    ( system )

If( myKernel == 1 ) then

    CALL MPI_BCAST( AO_bra, 2*N , mpi_D_C , 0 , KernelComm , err ) 
    CALL MPI_BCAST( AO_ket, 2*N , mpi_D_C , 0 , KernelComm , err ) 

    if( .not. allocated(rho_eh) ) then
        allocate( rho_eh  (N,N) )
        allocate( A_ad_nd (N,N) )
        allocate( H_prime (N,N) )
#ifdef USE_GPU
        allocate( X_ij    (N,N) )
        call GPU_Pin( X_ij,    8*N**2 )
        call GPU_Pin( Kernel,  8*N**2 )
        call GPU_Pin( H_prime, 8*N**2 )
        call GPU_Pin( A_ad_nd, 8*N**2 )
#else
        allocate( B_ad_nd (N,N) )
#endif
    end If


#ifdef USE_GPU
    call MPI_Irecv( H_prime, N*N, mpi_D_R, 0, 0, ChebyKernelComm, request_H_prime, err)
#else
    CALL MPI_IBCAST( H_prime, N*N, mpi_D_R, 0, ChebyKernelComm, request_H_prime, err )
#endif

    ! build up electron-hole density matrix ...
    forall( i=1:N , j=1:N ) rho_eh(i,j) = real( AO_ket(j,1)*AO_bra(i,1) ) - real( AO_ket(j,2)*AO_bra(i,2) )
    
    A_ad_nd = ( rho_eh + transpose(rho_eh) ) / two

    CALL Huckel_stuff( basis , X_ij ) 

    call MPI_Wait( request_H_prime, mpi_status, err )

#ifdef USE_GPU
    call EhrenfestKernel_gpu( N, H_prime, A_ad_nd, X_ij, Kernel )
#else
    CALL gemm( H_prime , A_ad_nd , B_ad_nd )

    Kernel = X_ij * A_ad_nd - B_ad_nd
#endif

end If

call MPI_BCAST( Kernel , N*N , mpi_D_R , 1 , ForceComm , err )

! Run, Forrest, Run ...
!================================================================================================
! set all forces to zero beforehand ...

forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero
F_rcv = D_zero ; F_snd = D_zero ; F_vec = D_zero 

do i = myForce+1 , system% atoms , npForce

    If( system%QMMM(i) == "MM" .OR. system%flex(i) == F_ ) cycle
    CALL Ehrenfest( system, basis, i ) 

end do


call MPI_reduce( F_snd , F_rcv , 3*system%atoms**2 , MPI_double_precision , mpi_SUM , 0 , ForceComm , err )

if( master ) then
    !$omp parallel do private(i,xyz) default(shared)
    do i = 1 , system% atoms
    do xyz = 1 , 3
       atom(i)% Ehrenfest(xyz) = two * sum( F_rcv(:,i,xyz) ) * eVAngs_2_Newton 
    end do
    end do
    !$omp end parallel do
end if

deallocate( mask )

! return Extended ForceCrew to stand-by ...
If( ForceCrew ) goto 99

include 'formats.h'

end subroutine Diabatic_Ehrenfest
!
!
!
!===========================================
 subroutine Ehrenfest( system, basis, site ) 
!===========================================
use Semi_empirical_parms , only: atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
integer          , intent(in)    :: site 

! local variables ...
integer :: i , j , xyz , size_basis , jL , L , indx
integer :: k , ik , DOS_atom_k , BasisPointer_k 

! local arrays ...
integer , allocatable :: pairs(:)
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) 
real*8                :: tmp_coord(3) , delta_b(3)


verbose    = .false.
size_basis = size(basis)
grad_S     = D_zero

!force on atom site ...
k = site 
DOS_atom_k     =  atom( system% AtNo(k) )% DOS
BasisPointer_k =  system% BasisPointer(k) 

allocate( pairs , source = pack([( L , L=1,system% atoms )] , mask(:,K)) )

! save coordinate ...
tmp_coord = system% coord(k,:)


do xyz = 1 , 3

       delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )
   
       system% coord (k,:) = tmp_coord + delta_b
       CALL Overlap_Matrix( system , basis , S_fwd , purpose = "Pulay" , site = K )
       
       system% coord (k,:) = tmp_coord - delta_b
       CALL Overlap_Matrix( system , basis , S_bck , purpose = "Pulay" , site = K )

       forall( j=1:DOS_Atom_K ) grad_S(:,j) = ( S_fwd( : , BasisPointer_K+j ) - S_bck( : , BasisPointer_K+j ) ) / (TWO*delta)

       !==============================================================================================
       F_vec = D_zero

       !$OMP parallel do schedule(dynamic) private(indx,iK,jL,i,j,L) default(shared) reduction(+:F_vec)
       do indx = 1, size(pairs)
         
         L = pairs(indx)
         do jL = 1 , DOS(L)
            j = BasisPointer(L) + jL

            do iK = 1 , DOS_atom_K
              i = BasisPointer_K + iK

              ! adiabatic and non-adiabatic components of the Force ...
              F_vec(L) = F_vec(L) - grad_S(j,iK) * Kernel(i,j)

            end do

         end do

       end do
       !$OMP end parallel do
       !==============================================================================================
       
       ! anti-symmetric F_snd (action-reaction) ...
       F_snd(K,K,xyz) = D_zero
       do L = K+1, system% atoms
          F_snd(L,K,xyz) =   F_vec(L)
          F_snd(K,L,xyz) = - F_vec(L) ! F_snd(L,K,xyz) 
       end do
 
end do 

! recover original system ...
system% coord (K,:) = tmp_coord

deallocate(pairs)

end subroutine Ehrenfest
!
!
!
!=====================================
 subroutine Huckel_stuff( basis , Xi ) 
!=====================================
implicit none
type(STO_basis) , intent(in) :: basis(:)
real*8          , allocatable , intent(out) :: Xi(:,:)

!local variables ...
integer :: i, j, n
real*8  :: k_eff, k_WH, c1, c2, c3
real*8  :: basis_j_IP, basis_j_k_WH

n = size(basis)

if (.not. allocated(Xi)) allocate( Xi(n,n) )

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

!$omp parallel private(i,j,basis_j_IP,basis_j_k_WH,c1,c2,c3,k_WH,k_eff) default(shared)
!$omp do schedule(dynamic,1)
do j = 1, n

    basis_j_IP   = basis(j)%IP
    basis_j_k_WH = basis(j)%k_WH

    do i = 1, j-1

        c1 = basis(i)%IP - basis_j_IP
        c2 = basis(i)%IP + basis_j_IP
        c3 = (c1/c2)**2

        k_WH = (basis(i)%k_WH + basis_j_k_WH) * half

        k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

        Xi(i,j) = k_eff * c2 * half

    end do

    Xi(j,j) = basis_j_IP

end do
!$omp end do

!$omp do
do i = 1, n
    Xi( i+1:n, i ) = Xi( i, i+1:n )
end do
!$omp end do nowait
!$omp end parallel

end subroutine Huckel_stuff
!
!
!
!
!============================
 subroutine Preprocess( sys ) 
!============================
use Semi_empirical_parms , only: atom
implicit none
type(structure) , intent(in) :: sys

!local variables ...
real*8                :: R_LK
integer               :: K , L
logical               :: flag1 , flag2 , flag3

If( .NOT. allocated(BasisPointer) ) allocate( BasisPointer(sys%atoms) , DOS(sys%atoms) )

Allocate( mask(sys%atoms,sys%atoms) , source = .false. )

do K = 1   , sys% atoms
   do L = K+1 , sys% atoms
   
       R_LK = sqrt(sum( (sys%coord(K,:)-sys%coord(L,:))**2 ) )
   
       flag1 = R_LK < cutoff_Angs  
        
       flag2 = sys% flex(K) .AND. sys% flex(L)
   
       flag3 = (sys% QMMM(L) == "QM") .AND. (sys% QMMM(K) == "QM")
   
       mask(L,K) = flag1 .AND. flag2 .AND. flag3

   end do
   BasisPointer(K) = sys% BasisPointer(K) 
   DOS(K)          = atom( sys% AtNo(K) )% DOS
end do    

end subroutine Preprocess
!
!
end module DiabaticEhrenfest_Builder

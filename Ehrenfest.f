! Program for computing Ehrenfest forces from Huckel Hamiltonian
module Ehrenfest_Builder

    use MPI
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use MPI_definitions_m       , only  : myForce , master , npForce , myKernel , KernelComm , ForceComm , KernelCrew , ForceCrew 
    use parameters_m            , only  : driver , verbose , n_part , QMMM
    use Structure_Builder       , only  : Unit_Cell 
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: EhrenfestForce 

    private

    !module variables ...
    integer     , allocatable   :: BasisPointer(:) , DOS(:)
    real*8      , allocatable   :: A_ad_nd(:,:) , B_ad_nd(:,:) , Kernel(:,:) , rho_eh(:,:) , tool(:,:) , X_ij(:,:)
    real*8      , allocatable   :: grad_S(:,:) , F_vec(:) , F_snd(:,:,:) , F_rcv(:,:,:) 
    logical     , allocatable   :: mask(:,:)
    logical :: first_time = .true.

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]
    real*8  , parameter :: delta      = 1.d-8

contains
!
!
!
!==================================================================
 subroutine EhrenfestForce( system , basis , QM , MO_bra , MO_ket )
!==================================================================
 use MD_read_m              , only  : atom
 implicit none
 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 type(R_eigen)   , optional , intent(in)    :: QM
 complex*16      , optional , intent(inout) :: MO_bra(:,:)
 complex*16      , optional , intent(inout) :: MO_ket(:,:)

! local variables ... 
 integer :: i , j , N , xyz , err  
 integer :: mpi_status(mpi_status_size) , request
 integer :: mpi_D_R = mpi_double_precision
 integer :: mpi_D_C = mpi_double_complex

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 logical , parameter :: T_ = .true. , F_ = .false.

!================================================================================================
! some preprocessing ...
!================================================================================================

N = size(basis)

If( .NOT. allocated(Kernel) ) then
        allocate( Kernel (N,N)                         )
        allocate( F_snd  (system%atoms,system%atoms,3) )
        allocate( F_rcv  (system%atoms,system%atoms,3) )
        allocate( F_vec  (system%atoms)                ) 
end If

! ForceCrew in stand-by ...
99 If( .not. master ) CALL MPI_BCAST( system%coord , system%atoms*3 , mpi_D_R , 0 , ForceComm, err )

! preprocess overlap matrix for Pulay calculations, all ForceCrew must do it ...
 CALL Overlap_Matrix( system , basis )
 CALL preprocess    ( system )

if( KernelCrew ) then

    ! KernelCrew in stand-by to receive data from master ...
    CALL MPI_BCAST( QM%erg , N   , mpi_D_R , 0 , KernelComm , err )
    CALL MPI_BCAST( QM%L   , N*N , mpi_D_R , 0 , KernelComm , err )
    CALL MPI_BCAST( MO_ket , N*2 , mpi_D_C , 0 , KernelComm , err )

    If( .NOT. allocated(rho_eh) ) then
        allocate( rho_eh  (N,N) )
        allocate( B_ad_nd (N,N) )
        allocate( tool    (N,N) )
        If( myKernel == 1 ) then
             allocate( A_ad_nd (N,N) )
             CALL Huckel_stuff( basis )
        end If
    end If

    select case (myKernel)

        case (1)

           MO_bra = conjg( MO_ket )

           ! build up electron-hole density matrix ...
           forall( i=1:N , j=1:N ) rho_eh(i,j) = real( MO_ket(j,1)*MO_bra(i,1) - MO_ket(j,2)*MO_bra(i,2) )

           CALL MPI_ISend( rho_eh , N*N , mpi_D_R , 2 , 0 , KernelComm , request , err )
           CALL MPI_Request_Free( request , err )

           CALL symm( rho_eh , QM%L , tool )
           CALL gemm( QM%L , tool , A_ad_nd , 'T' , 'N' )

           CALL MPI_Recv( B_ad_nd , N*N , mpi_D_R , 2 , mpi_any_tag , KernelComm , mpi_status , err )

           Kernel = X_ij * A_ad_nd - B_ad_nd

        case (2)  ! <== firstmate of KernelCrew ...

           CALL MPI_Recv( rho_eh , N*N , mpi_D_R , 1 , mpi_any_tag , KernelComm , mpi_status , err )

           forall( j=1:N ) rho_eh(:,j) = QM%erg(j) * rho_eh(:,j) 

           CALL gemm( rho_eh , QM%L , tool )
           CALL gemm( tool , QM%L  , B_ad_nd , 'T' , 'N' )

           CALL MPI_ISend( B_ad_nd , N*N , mpi_D_R , 1 , 0 , KernelComm , request , err )
           CALL MPI_Request_Free( request , err )

    end select

end if

call MPI_BCAST( Kernel , N*N , mpi_D_R , 1 , ForceComm , err )

! Run, Forrest, Run ...
!================================================================================================
! set all forces to zero beforehand ...

forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero
F_rcv = D_zero ; F_snd = D_zero ; F_vec = D_zero 

select case( driver )

    case( "slice_AO" )

        do i = myForce+1 , system% atoms , npForce

            If( system%QMMM(i) == "MM" .OR. system%flex(i) == F_ ) cycle
            CALL Ehrenfest_AO( system, basis, i ) 

        end do

        call MPI_reduce( F_snd , F_rcv , 3*system%atoms**2 , MPI_double_precision , mpi_SUM , 0 , ForceComm , err )

        if( master ) then
            do i = 1 , system% atoms
            do xyz = 1 , 3
               atom(i)% Ehrenfest(xyz) = two * sum( F_rcv(:,i,xyz) ) * eVAngs_2_Newton 
            end do
            end do
        end if

    case( "slice_ElHl" )

        
end select

deallocate( mask )

! return Extended ForceCrew to stand-by ...
If( ForceCrew .OR. KernelCrew ) goto 99

include 'formats.h'

end subroutine EhrenfestForce
!
!
!
!==============================================
 subroutine Ehrenfest_AO( system, basis, site ) 
!==============================================
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

verbose = .false.
size_basis = size(basis)
If( .NOT. allocated(grad_S) ) allocate( grad_S( size_basis , 10 ) )
grad_S = D_zero

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

       !$OMP parallel do schedule(dynamic,3) private(iK,jL,i,j,L) default(shared) reduction(+:F_vec)
       do indx = 1 , size(pairs)
         
         L = pairs(indx)
         do jL = 1 , DOS(L)
            j = BasisPointer(L) + jL

           do iK = 1 , DOS_atom_K
              i = BasisPointer_K + iK

              ! adiabatic and non-adiabatic components of the Force ...
              F_vec(L) = F_vec(L) -  grad_S(j,iK) * Kernel(i,j)

           end do   
         end do

       end do
       !$OMP end parallel do
       !==============================================================================================
 
       ! anti-symmetric F_snd (action-reaction) ...
       do L = K+1, system% atoms
          F_snd(L,K,xyz) =   F_vec(L)
          F_snd(K,L,xyz) = - F_snd(L,K,xyz) 
       end do
       F_snd(K,K,xyz) = D_zero
 
end do 

! recover original system ...
system% coord (K,:) = tmp_coord

deallocate(pairs)

end subroutine Ehrenfest_AO
!
!
!
!
!================================
 subroutine Huckel_stuff( basis ) 
!================================
use Hamiltonians , only : X => X_ij
implicit none
type(STO_basis) , intent(in) :: basis(:)

!local variables ...
integer :: i , j , N

N = size(basis)

allocate ( X_ij(N,N) )

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

do j = 1 , N
do i = j , N

         X_ij(i,j) = X( i , j , basis )

         X_ij(j,i) = X_ij(i,j)

end do
end do

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
!
!
end module Ehrenfest_Builder

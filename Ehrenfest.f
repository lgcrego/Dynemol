! Program for computing Ehrenfest forces from Huckel Hamiltonian
module Ehrenfest_Builder

    use MPI
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use MPI_definitions_m       , only  : myForce , master , npForce , myKernel , KernelComm , ForceComm , KernelCrew , ForceCrew , world
    use parameters_m            , only  : driver , verbose , n_part , QMMM
    use Structure_Builder       , only  : Unit_Cell 
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: EhrenfestForce 

    private

    !module variables ...
    integer , allocatable :: BasisPointer(:) , DOS(:) , scheduler(:,:)
    real*8  , allocatable :: A_ad_nd(:,:) , B_ad_nd(:,:) , Kernel(:,:) , rho_eh(:,:) , tool(:,:) , X_ij(:,:)
    real*8  , allocatable :: grad_S(:,:) , F_vec(:) , F_snd(:,:,:) , F_rcv(:,:,:) 
    logical , allocatable :: mask(:,:)

    logical :: done_npForceSchedule = .false.

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
 integer :: i , N , xyz , err , mytasks , myatom
 integer :: mpi_D_R = mpi_double_precision
 logical :: job_done , QMMM_done , job_status(2)

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 logical , parameter :: T_ = .true. , F_ = .false.

!======================================================================================
! some preprocessing ...
!======================================================================================

N = size(basis)

If( .NOT. allocated(Kernel) ) then
        allocate( Kernel (N,N)                         )
        allocate( F_snd  (system%atoms,system%atoms,3) )
        allocate( F_rcv  (system%atoms,system%atoms,3) )
        allocate( F_vec  (system%atoms)                ) 
end If

! Force+Kernel Crew in stand-by ...
99 If( .not. master ) then

       CALL MPI_BCAST( system%coord , system%atoms*3 , mpi_D_R , 0 , ForceComm, err )

       CALL MPI_BCAST( job_status , 2 , mpi_logical , 0 , world , err ) 
       job_done  = job_status(1)
       QMMM_done = job_status(2)
       If( QMMM_done ) then ! <== Force+Kernel Crew pack and wait ...
           call FINALIZE
       elseif( job_done ) then
           call packing
           call MPI_FINALIZE(err)
           STOP
       end if

end if

! preprocess overlap matrix for Pulay calculations, all Force+Kernel Crew must do it ...
 CALL Overlap_Matrix( system , basis )
 CALL preprocess    ( system , mytasks)

if( KernelCrew ) then
    call get_Kernel( basis , QM , MO_bra , MO_ket )
    end if
call MPI_BCAST( Kernel , N*N , mpi_D_R , 1 , ForceComm , err )

! set all forces to zero beforehand ...
forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero
F_rcv = D_zero ; F_snd = D_zero ; F_vec = D_zero 

select case( driver )

    case( "slice_AO" )

        ! Run, Forrest, Run ...
        do i = 1 , mytasks

            myatom = scheduler(i,myForce+1)

            If( system%QMMM(myatom) == "MM" .OR. system%flex(myatom) == F_ ) cycle
            CALL Ehrenfest_AO( system, basis, myatom ) 

        end do

        call MPI_reduce( F_snd , F_rcv , 3*system%atoms**2 , MPI_double_precision , mpi_SUM , 0 , ForceComm , err )

        if( master ) then
            do i = 1 , system% atoms
            do xyz = 1 , 3
               atom(i)% Ehrenfest(xyz) = two * sum( F_rcv(:,i,xyz) ) * eVAngs_2_Newton 
            end do
            end do
        end if
        
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
!=====================================================
 subroutine get_Kernel( basis , QM , MO_bra , MO_ket ) 
!=====================================================
implicit none
type(STO_basis)          , intent(in)    :: basis(:)
type(R_eigen) , optional , intent(in)    :: QM
complex*16    , optional , intent(inout) :: MO_bra(:,:)
complex*16    , optional , intent(inout) :: MO_ket(:,:)

! local variables ... 
integer :: i , j , N , err 
integer :: mpi_status(mpi_status_size) , request
integer :: mpi_D_R = mpi_double_precision
integer :: mpi_D_C = mpi_double_complex

N = size(basis)

! KernelCrew in stand-by to receive data from master ...
CALL MPI_BCAST( QM%erg , N   , mpi_D_R , 0 , KernelComm , err )
CALL MPI_BCAST( QM%L   , N*N , mpi_D_R , 0 , KernelComm , err )
CALL MPI_BCAST( MO_bra , N*2 , mpi_D_C , 0 , KernelComm , err )
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

       ! build up electron-hole density matrix ...
       forall( i=1:N , j=1:N ) rho_eh(i,j) = real( MO_ket(j,1)*MO_bra(i,1) - MO_ket(j,2)*MO_bra(i,2) )
       tool   = transpose(rho_eh)
       rho_eh = ( rho_eh + tool ) / two

       CALL MPI_ISend( rho_eh , N*N , mpi_D_R , 2 , 0 , KernelComm , request , err )
       CALL MPI_Request_Free( request , err )

       CALL symm( rho_eh , QM%L , tool )
       CALL gemm( QM%L , tool , A_ad_nd , 'T' , 'N' )

       CALL MPI_Recv( B_ad_nd , N*N , mpi_D_R , 2 , mpi_any_tag , KernelComm , mpi_status , err )

       Kernel = X_ij * A_ad_nd - B_ad_nd  ! <== all this to calculate Kernel ...

    case (2) 

       CALL MPI_Recv( rho_eh , N*N , mpi_D_R , 1 , mpi_any_tag , KernelComm , mpi_status , err )

       forall( j=1:N ) rho_eh(:,j) = QM%erg(j) * rho_eh(:,j) 

       CALL gemm( rho_eh , QM%L , tool )
       CALL gemm( tool , QM%L  , B_ad_nd , 'T' , 'N' )

       CALL MPI_ISend( B_ad_nd , N*N , mpi_D_R , 1 , 0 , KernelComm , request , err )
       CALL MPI_Request_Free( request , err )

end select

end subroutine get_Kernel
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
!======================================
 subroutine Preprocess( sys , mytasks ) 
!======================================
use Semi_empirical_parms , only: atom
implicit none
type(structure) , intent(in)  :: sys
integer         , intent(out) :: mytasks

!local variables ...
real*8  :: R_LK
integer :: K , L
logical :: flag1 , flag2 , flag3

integer :: NFold , remainder 
integer :: i , j , j1 , j2 , step_j , aux

If( .NOT. allocated(BasisPointer) ) allocate( BasisPointer(sys%atoms) , DOS(sys%atoms) )

! define (update) matrix of atom pairs for HFP calculations ...
!-------------------------------------------------------------------
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
!-------------------------------------------------------------------

! define Load Balance for MPI HFP calculations in EhrenfestForce ...
!-------------------------------------------------------------------
If( .not. done_npForceSchedule ) then

   NFold = int(sys% atoms/npForce) + merge( 1 , 0 , mod(sys% atoms,npForce)/=0 )
   remainder = mod(sys% atoms,npForce)
   
   allocate( scheduler(nFold,npForce) , source = -1 )
   
   k = 0; j1 = 1 ; j2 = npForce ; step_j = 1
   do i = 1 , NFold
      do j = j1 , j2 , step_j
   
         k = k + 1
         scheduler(i,j) = merge( k , 0 , k <= sys% atoms)
   
      end do
   
      aux = j1
      j1  = j - step_j
      j2  = aux
      step_j = step_j * (-1)
   
   end do

   done_npForceSchedule = .true.

End If

mytasks = count( scheduler(:,myForce+1) /= 0 )
!-------------------------------------------------------------------

end subroutine Preprocess
!
!
!
!
!==================
 subroutine packing
!==================
implicit none

deallocate( Kernel , F_snd , F_rcv , F_vec )

if( KernelCrew ) then
     deallocate( rho_eh , B_ad_nd , tool )
     If( myKernel == 1 ) deallocate( A_ad_nd )
end if

end subroutine packing
!
!
!
!
!===================
 subroutine FINALIZE
!===================
implicit none

!local variables ...
integer :: err
logical :: job_status(2) , job_done

98 CALL MPI_BCAST( job_status , 2 , mpi_logical , 0 , world , err )
job_done = job_status(1)

If( job_done ) then
    call packing
    call MPI_FINALIZE(err)
    STOP
else
    goto 98
end if

end subroutine FINALIZE
!
!
!
!
end module Ehrenfest_Builder

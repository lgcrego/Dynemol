! Program for computing Ehrenfest forces from Huckel Hamiltonian with Coherent-Switch-Decay-of-Mixing
module CSDM_workers

    use MPI
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use MD_read_m         , only: atom
    use MPI_definitions_m , only: myKernel, KernelComm, ForceComm, KernelCrew, ForceCrew, ForceCrewComm , world, myforce
    use MPI_definitions_m , only: myAxisComm , AxisCrew , myaxis_rank , np_per_axis
    use Overlap_Builder   , only: Overlap_Matrix

    public :: Ehrenfest_workers 

    private 

    !module variables ...
    integer                                :: dim_E , dim_N , PST(2) , mytasks
    integer , allocatable , dimension(:)   :: BasisPointer, DOS , list
    integer , allocatable , dimension(:,:) :: my_axis_todo_list
    real*8  , allocatable , dimension(:)   :: erg, F_vec
    real*8  , allocatable , dimension(:,:) :: Kernel, grad_S , QL, Phi, d_NA, Xij, tmp_El, tmp_Hl
    logical , allocatable , dimension(:,:) :: mask

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

contains
!
!
!
!==============================================
 subroutine Ehrenfest_workers( system , basis )
!==============================================
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(in)    :: basis(:)

! local variables ... 
real*8  , allocatable :: Force(:)
real*8  :: real_dumb
integer :: i , j , xyz , N , err , int_dumb
integer :: mpi_D_R = mpi_double_precision
logical :: job_done

dim_E = size(basis)
N = dim_E 

! Force+Kernel Crew in stand-by ...
99 CALL MPI_BCAST( system%coord , system%atoms*3 , mpi_D_R , 0 , ForceComm , err )

! preprocess overlap matrix for Pulay calculations, all Force+Kernel Crew must do it ...                                                                                
CALL Overlap_Matrix( system , basis )

CALL MPI_BCAST( job_done , 1 , mpi_logical , 0 , world , err ) 
If( job_done ) then ! <== Force+Kernel Crew pack and stop here ...
    call MPI_FINALIZE(err)
    STOP
    end if

CALL preprocess( system , basis )

if( KernelCrew ) then
   CALL get_Kernel( basis )
   end if
CALL MPI_BCAST( Kernel , N*N , mpi_D_R , 0 , ForceCrewComm , err )

xyz = mod(myforce,3)+1
CALL EhrenfestForce( system , basis , Force , xyz )

if(KernelCrew) then

    CALL MPI_GatherV( Force , system%atoms , mpi_D_R , real_dumb , int_dumb , int_dumb , mpi_D_R , 0 , KernelComm , err ) 

    CALL MPI_GatherV( tmp_El , dim_N*dim_E , mpi_D_R , real_dumb , int_dumb , int_dumb , mpi_D_R , 0 , KernelComm , err ) 
    CALL MPI_GatherV( tmp_Hl , dim_N*dim_E , mpi_D_R , real_dumb , int_dumb , int_dumb , mpi_D_R , 0 , KernelComm , err ) 

end if

deallocate( mask , Force , tmp_El , tmp_Hl )

! return Extended ForceCrew to stand-by ...
If( ForceCrew ) goto 99

include 'formats.h'

end subroutine Ehrenfest_workers
!
!
!
!
!=========================================================
 subroutine EhrenfestForce( system , basis , Force , xyz )
!=========================================================
use Semi_empirical_parms , only: ChemAtom => atom
implicit none
type(structure)       , intent(inout) :: system
type(STO_basis)       , intent(in)    :: basis(:)
real*8  , allocatable , intent(out)   :: Force(:)
integer               , intent(in)    :: xyz

! local parameters ...
integer , parameter :: xyz_key(3) = [1,2,3]
real*8  , parameter :: delta = 1.d-8
real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 

! local variables ...
integer :: i , j , jL , L , indx , err
integer :: k , ik , DOS_k , BP_k , task , AllAtoms , address
integer :: mpi_D_R = mpi_double_precision

! local arrays ...
integer , allocatable :: pairs(:)
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:), F_mtx(:,:) 
real*8                :: tmp_coord(3) , delta_b(3) 

AllAtoms = system%atoms

allocate( F_mtx(AllAtoms,AllAtoms) )
allocate( Force(AllAtoms) )
F_mtx = d_zero ; Force  = d_zero

grad_S = d_zero

do task = 1 , mytasks

    k = my_axis_todo_list( task , myAxis_rank+1 )

    !force on atom site ...
    DOS_k = ChemAtom( system% AtNo(k) )% DOS
    BP_k  = system% BasisPointer(k) 

    allocate( pairs , source = pack([( L , L=1,AllAtoms )] , mask(:,K)) )
    
    ! save coordinate ...
    tmp_coord = system% coord(k,:)

    delta_b = delta * merge(D_one , d_zero , xyz_key == xyz )
 
    system% coord (k,:) = tmp_coord + delta_b
    CALL Overlap_Matrix( system , basis , S_fwd , purpose = "Pulay" , site = K )

    system% coord (k,:) = tmp_coord - delta_b
    CALL Overlap_Matrix( system , basis , S_bck , purpose = "Pulay" , site = K )

    ! grad_S is an anti-symmetric matrix 
    do j = 1 , DOS_k
       grad_S( BP_k+1: , BP_k+j  ) = ( S_fwd( BP_k+1: , BP_k+j ) - S_bck( BP_k+1: , BP_k+j ) ) / (TWO*delta) 
       grad_S( BP_k+j  , BP_k+1: ) = -grad_S( BP_k+1:,BP_k+j )
    end do

    !==============================================================================================
    F_vec = d_zero

    !$OMP parallel do schedule(dynamic,3) private(iK,jL,i,j,L) default(shared) reduction(+:F_vec)
    do indx = 1 , size(pairs)
      
       L = pairs(indx)
       do jL = 1 , DOS(L)
          j = BasisPointer(L) + jL

          do iK = 1 , DOS_K
             i = BP_K + iK

             ! adiabatic and non-adiabatic components of the Force ...
             F_vec(L) = F_vec(L) -  grad_S(j,i) * Kernel(i,j)

          end do   
       end do

    end do
    !$OMP end parallel do
    !==============================================================================================
     
    ! anti-symmetric F_mtx (action-reaction) ...
    do L = K+1, AllAtoms
       F_mtx(K,L) =   F_vec(L)
       F_mtx(L,K) = - F_mtx(K,L) 
    end do
    F_mtx(K,K) = d_zero

    Force(K) = two * sum( F_mtx(K,:) ) * eVAngs_2_Newton

    ! calculation of d_NA ...
    d_NA = NAcoupling( grad_S( : , BP_K+1 : BP_K+DOS_k) , DOS_k , BP_K )  ! <== units = 1/Angs

    address = findloc( list , value=k , dim=1 )
    do concurrent (j=1:dim_E)
       tmp_El(address,j) = d_NA(j,1)
       tmp_Hl(address,j) = d_NA(j,2)
       enddo

    ! recover original system ...
    system% coord (K,:) = tmp_coord
            
    deallocate(pairs)
    ! ready for next atom in system

end do 
deallocate( F_mtx )

if( myAxis_rank == 0 ) then
    call MPI_reduce( MPI_in_Place , Force , AllAtoms , mpi_D_R , mpi_SUM , 0 , myAxisComm , err )
    call MPI_reduce( MPI_in_Place , tmp_El , dim_N*dim_E , mpi_D_R , mpi_SUM , 0 , myAxisComm , err ) 
    call MPI_reduce( MPI_in_Place , tmp_Hl , dim_N*dim_E , mpi_D_R , mpi_SUM , 0 , myAxisComm , err ) 
else
    call MPI_reduce( Force , Force , AllAtoms , mpi_D_R , mpi_SUM , 0 , myAxisComm , err )
    call MPI_reduce( tmp_El , tmp_El , dim_N*dim_E , mpi_D_R , mpi_SUM , 0 , myAxisComm , err ) 
    call MPI_reduce( tmp_Hl , tmp_Hl , dim_N*dim_E , mpi_D_R , mpi_SUM , 0 , myAxisComm , err ) 
endif

end  subroutine EhrenfestForce
!
!
!
!
!==============================
 subroutine get_Kernel( basis ) 
!==============================
implicit none
type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
integer :: i , j , N , err 
integer :: mpi_status(mpi_status_size) , request
integer :: mpi_D_R = mpi_double_precision
integer :: mpi_D_C = mpi_double_complex
real*8     , allocatable :: A_ad_nd(:,:) , B_ad_nd(:,:) , rho_eh(:,:) , tool(:,:)
complex*16 , allocatable :: MObra(:,:) , MOket(:,:)

N = dim_E

allocate( rho_eh  (N,N) )
allocate( B_ad_nd (N,N) )
allocate( tool    (N,N) )
allocate( MObra   (N,2) )
allocate( MOket   (N,2) )

! KernelCrew in stand-by to receive data from master ...
CALL MPI_BCAST( MObra  , N*2 , mpi_D_C , 0 , KernelComm , err )
CALL MPI_BCAST( MOket  , N*2 , mpi_D_C , 0 , KernelComm , err )

select case (myKernel)

    case (1)

          ! build up electron-hole density matrix ...
          forall( i=1:N , j=1:N ) rho_eh(i,j) = real( MOket(j,1)*MObra(i,1) - MOket(j,2)*MObra(i,2) )
          tool   = transpose(rho_eh)
          rho_eh = ( rho_eh + tool ) / two

          CALL MPI_ISend( rho_eh , N*N , mpi_D_R , 2 , 0 , KernelComm , request , err )
          CALL MPI_Request_Free( request , err )

          allocate( A_ad_nd (N,N) )
          CALL symm( rho_eh , QL , tool )
          CALL gemm( QL , tool , A_ad_nd , 'T' , 'N' )

          CALL MPI_Recv( B_ad_nd , N*N , mpi_D_R , 2 , mpi_any_tag , KernelComm , mpi_status , err )

          Kernel = Xij * A_ad_nd - B_ad_nd  ! <== all this to calculate Kernel ...

          deallocate( A_ad_nd )

    case (2) 

          CALL MPI_Recv( rho_eh , N*N , mpi_D_R , 1 , mpi_any_tag , KernelComm , mpi_status , err )

          forall( j=1:N ) rho_eh(:,j) = erg(j) * rho_eh(:,j) 

          CALL gemm( rho_eh , QL , tool )
          CALL gemm( tool , QL  , B_ad_nd , 'T' , 'N' )

          CALL MPI_Send( B_ad_nd , N*N , mpi_D_R , 1 , 0 , KernelComm , err )

end select

deallocate( rho_eh , B_ad_nd , tool , MObra , MOket )

end subroutine get_Kernel
!
!
!
!
!================================================
 function NAcoupling( grad_Slice , DOSk , BPk ) &
 result(d_NA)
!================================================
implicit none
real*8  , intent(in)  :: grad_Slice(:,:)
integer , intent(in)  :: DOSk
integer , intent(in)  :: BPk
! result ...
real*8  , allocatable :: d_NA(:,:)

! local variables ... 
integer               :: j , j1 , j2 , dima , dimb
real*8  , allocatable :: Mat1(:,:) , A(:,:) , R1(:,:) , R2(:,:)
real*8  , allocatable :: Mat2(:,:) , B(:,:) 

j1 = BPk + 1
j2 = BPk + DOSk

dima = size(grad_Slice(:,1))
dimb = size(grad_Slice(1,:))

! temporary arrays ...
allocate( A(dima,2) , R1(dima,2) , R2(dima,2) , d_NA(dima,2) , Mat2(dima,dima) )

Phi(:,1) = QL(PST(1),:)
Phi(:,2) = QL(PST(2),:)

do concurrent (j=1:dima) shared(QL,erg,Mat2)
   Mat2(:,j) = QL(:,j)*erg(:)
   end do

allocate( B(dimb,2) , Mat1(dima,dimb) )
Mat1 = grad_Slice * Xij(:,j1:j2)

!===============================================
CALL gemm( Mat1 , Phi(j1:j2,:) , A )
CALL gemm( QL , A , R1 )

CALL gemm( grad_Slice , Phi(j1:j2,:) , A ) 
CALL gemm( Mat2 , A , R2 )

d_NA = R1 - R2
!===============================================

!===============================================
CALL gemm( Mat1 , Phi , B , transa = 'T' ) 
CALL gemm( QL(:,j1:j2) , B , R1 )

do concurrent ( j=1:2 ) shared(erg,PST,Phi,A)
   A(:,j) = Phi(:,j) * erg(PST(j))
   end do

CALL gemm( grad_Slice , A , B , transa = 'T' ) 
CALL gemm( QL(:,j1:j2) , B , R2 )

d_NA = d_NA + (R1-R2)
!===============================================

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! checklist
if( abs( d_NA(PST(2),1)-d_NA(PST(1),2) > high_prec ) ) then
    Print*, "WARNING: failed high precision test in NAcoupling"
    end if
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

d_NA( PST(1) , 1 ) = d_zero
d_NA( PST(2) , 2 ) = d_zero

deallocate( Mat1 , Mat2 , A , B , R1 , R2 )

end function NAcoupling
!
!
!
!====================================
 subroutine preprocess( sys , basis )
!====================================
use Semi_empirical_parms , only: ChemAtom => atom
implicit none
type(structure) , intent(in) :: sys
type(STO_basis) , intent(in) :: basis(:)

! local variables ...
real*8  :: R_LK
integer :: i , j , K , L , err
integer :: mpi_D_R = mpi_double_precision 
logical :: flag1 , flag2 , flag3
logical , save :: first_time = .true.                                                                                                                           

if( first_time ) then
    call setup_Module( sys , basis )
    first_time = F_
    end if

allocate( tmp_El (dim_N,dim_E) , source = d_zero )
allocate( tmp_Hl (dim_N,dim_E) , source = d_zero )

CALL MPI_BCAST( erg , dim_E       , mpi_D_R     , 0 , ForceComm ,  err )
CALL MPI_BCAST( QL  , dim_E*dim_E , mpi_D_R     , 0 , ForceComm ,  err )
CALL MPI_BCAST( PST , 2           , mpi_Integer , 0 , ForceComm ,  err )

!-------------------------------------------------------------------
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
   DOS(K)          = ChemAtom( sys% AtNo(K) )% DOS
end do    
!-------------------------------------------------------------------

end subroutine Preprocess
!
!
!
!
!======================================
 subroutine setup_Module( sys , basis )
!======================================
implicit none
type(structure) , intent(in) :: sys
type(STO_basis) , intent(in) :: basis(:)

! local variables ...
integer :: NFold , remainder 
integer :: i , j , k , j1 , j2 , step_j , aux , L

allocate( grad_S  (dim_E,dim_E) )
allocate( Kernel  (dim_E,dim_E) )
allocate( QL      (dim_E,dim_E) )
allocate( erg     (dim_E)       )
allocate( Phi     (dim_E, 2)    )
allocate( d_NA    (dim_E, 2)    )

allocate( F_vec(sys%atoms) )

CALL Huckel_stuff( basis , Xij )

! list of atoms subject to Ehrenfest force ...
allocate( list , &
source = pack( [( L , L=1,sys%atoms )] , sys%QMMM(:) == "QM" .AND. sys%flex(:) == T_ ) &
)

dim_N = count( sys%QMMM == "QM" .AND. sys%flex == T_ )

!-------------------------------------------------------------------
! define Load Balance for MPI HFP calculations in EhrenfestForce ...
!-------------------------------------------------------------------
NFold = int(dim_N/np_per_axis) + merge( 1 , 0 , mod(dim_N,np_per_axis)/=0 )
remainder = mod(dim_N,np_per_axis)

allocate( my_axis_todo_list(nFold,np_per_axis) , source = 0 )

k = 0; j1 = 1 ; j2 = np_per_axis ; step_j = 1
do i = 1 , NFold
   do j = j1 , j2 , step_j

      k = k + 1
      if( k > size(list) ) exit

      my_axis_todo_list(i,j) = list(k)

   end do

   aux = j1
   j1  = j - step_j
   j2  = aux
   step_j = step_j * (-1)

end do
mytasks = count( my_axis_todo_list(:,myAxis_rank+1) /= 0 )
!-------------------------------------------------------------------

end subroutine setup_Module
!
!
!
!
!======================================
 subroutine Huckel_stuff( basis , Xij ) 
!======================================
use Hamiltonians , only : X_ij
implicit none
type(STO_basis)  , intent(in) :: basis(:)
real*8           , allocatable , intent(out) :: Xij(:,:)

!local variables ...
integer :: i , j

allocate ( Xij(dim_E,dim_E) )

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

do j = 1 , dim_E 
do i = j , dim_E 

         Xij(i,j) = X_ij( i , j , basis )
         Xij(j,i) = Xij(i,j) 

         end do
         end do

end subroutine Huckel_stuff
!
!
!
!
end module CSDM_workers

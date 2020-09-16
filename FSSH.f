! Program for computing Surface Hopping forces and transition probabilities from eHuckel Hamiltonian

module Surface_Hopping

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only  : driver , verbose , n_part , QMMM
    use Structure_Builder       , only  : Unit_Cell 
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: SH_Force 

    private

    !module variables ...
    integer                     :: mm , PES(2)
    integer     , allocatable   :: BasisPointer(:) , DOS(:) 
    real*8      , allocatable   :: Kernel(:,:) , X_(:,:) , QL(:,:) , Phi(:,:) , erg(:)
    real*8      , allocatable   :: grad_S(:,:) , F_vec(:) , F_mtx(:,:,:) , Rxd_NA(:,:)
    logical     , allocatable   :: mask(:,:)

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]
    real*8  , parameter :: delta      = 1.d-8

contains
!
!
!
!==========================================
 subroutine SH_Force( system , basis , QM )
!==========================================
 use MD_read_m     , only  : atom
 use parameters_m  , only  : electron_state , hole_state
 implicit none
 type(structure)   , intent(inout) :: system
 type(STO_basis)   , intent(in)    :: basis(:)
 type(R_eigen)     , intent(in)    :: QM

! local variables ... 
 integer :: i , j , nn 

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 logical , parameter :: T_ = .true. , F_ = .false.

!================================================================================================
! some preprocessing ...
!================================================================================================
PES(1) = electron_state
PES(2) = hole_state

mm = size(basis)
nn = n_part

If( .NOT. allocated(F_mtx)  ) allocate( F_mtx   (system%atoms,system%atoms,3) , source=D_zero )
If( .NOT. allocated(F_vec)  ) allocate( F_vec   (system%atoms)                , source=D_zero )
If( .NOT. allocated(Kernel) ) then
    allocate( Kernel (mm,mm) , source = D_zero )
    allocate( QL     (mm,mm) , source = QM%L   )
    allocate( erg    (mm)    , source = QM%erg )
    allocate( Phi    (mm, 2) )
    Phi(:,1) = QM%L(PES(1),:)
    Phi(:,2) = QM%L(PES(2),:)
    allocate( grad_S (mm,9 ) )
end if

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis )
CALL preprocess    ( system )

CALL Huckel_stuff( basis , X_ )

do j = 1 , mm
   kernel(:,j) = ( X_(:,j) - QM%erg(PES(1)) ) * Phi(:,1) * Phi(j,1)   &      ! <== electron part 
               - ( X_(:,j) - QM%erg(PES(2)) ) * Phi(:,2) * Phi(j,2)          ! <== hole part 
end do

!
!================================================================================================
! set all forces to zero beforehand ...
! using %Ehrenfest(:) to store the SH force ...
forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero

! setup new nonadiabtic coupling vector <Psi/dPhi/dt> ...
allocate( Rxd_NA(mm,2) , source = D_zero )

! Run, Forrest, Run ...

do i = 1 , system% atoms
    If( system%QMMM(i) == "MM" .OR. system%flex(i) == F_ ) cycle
    atom(i)% Ehrenfest = SHForce( system, basis, i ) * eVAngs_2_Newton 
end do

deallocate( mask , X_ , F_vec , F_mtx , QL , Phi , erg , Rxd_NA , Kernel , grad_S )

include 'formats.h'

end subroutine SH_Force
!
!
!
!=======================================================
 function SHForce( system, basis, site ) result(Force)
!=======================================================
use Semi_empirical_parms , only: atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
integer          , intent(in)    :: site 

! local variables ...
integer :: i , j , xyz , jL , L , indx
integer :: k , ik , DOS_atom_k , BasisPointer_k 

! local arrays ...
integer , allocatable :: pairs(:)
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) 
real*8                :: Force(3) , tmp_coord(3) , delta_b(3) 

verbose    = .false.
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
     
       ! anti-symmetric F_mtx (action-reaction) ...
       do L = K+1, system% atoms
          F_mtx(K,L,xyz) =   F_vec(L)
          F_mtx(L,K,xyz) = - F_mtx(K,L,xyz) 
       end do
       F_mtx(K,K,xyz) = D_zero

       Force(xyz) = two * sum( F_mtx(K,:,xyz) )

       Rxd_NA = Rxd_NA &
              + TransitionMatrix( grad_S(:,:DOS_atom_k) , k , DOS_atom_k , BasisPointer_k , xyz )

end do 

! recover original system ...
system% coord (K,:) = tmp_coord
        
deallocate(pairs)

end function SHForce
!
!
!
!=============================================================================
 function TransitionMatrix( grad_S , k , DOS_atom_k , BasisPointer_k , xyz ) &
 result(R)
!=============================================================================
use MD_read_m , only: atom  
implicit none
real*8  , intent(in) :: grad_S(:,:)
integer , intent(in) :: k
integer , intent(in) :: DOS_atom_k
integer , intent(in) :: BasisPointer_k
integer , intent(in) :: xyz

! local variables ... 
integer               :: i , j , j1 , j2 , dima , dimb
real*8  , allocatable :: Mat1(:,:) , A(:,:) , R1(:,:) , R2(:,:)
real*8  , allocatable :: Mat2(:,:) , B(:,:) , R(:,:)

j1 = BasisPointer_K + 1
j2 = BasisPointer_K + DOS_atom_K

dima = size(grad_S(:,1))
dimb = size(grad_S(1,:))

! temporary arrays ...
allocate( A(dima,2) , R1(dima,2) , R2(dima,2) , R(dima,2) , Mat2(dima,dima) )
forall(j=1:dima) Mat2(:,j) = QL(:,j)*erg(:)

allocate( B(dimb,2) , Mat1(dima,dimb) )
Mat1 = grad_S * X_(:,j1:j2)

!===============================================
CALL gemm( Mat1 , Phi(j1:j2,:) , A )
CALL gemm( QL , A , R1 )

CALL gemm( grad_S , Phi(j1:j2,:) , A ) 
CALL gemm( Mat2 , A , R2 )

R = R1 - R2
!===============================================

!===============================================
CALL gemm( Mat1 , Phi , B , transa = 'T' ) 
CALL gemm( QL(:,j1:j2) , B , R1 )

forall( j=1:2 ) A(:,j) = Phi(:,j) * erg(PES(j))
CALL gemm( grad_S , A , B , transa = 'T' ) 
CALL gemm( QL(:,j1:j2) , B , R2 )

R = R + (R1-R2)
!===============================================

forall( i=1:dima , j=1:2 , i/=PES(j) ) R(i,j) = atom(k)%vel(xyz) * R(i,j) / ( erg(i) - erg(PES(j)) )

deallocate( Mat1 , Mat2 , A , B , R1 , R2 )

end function TransitionMatrix
!
!
!
!
!=====================================
 subroutine Huckel_stuff( basis , Xi ) 
!=====================================
use Hamiltonians , only : X_ij
implicit none
type(STO_basis)  , intent(in) :: basis(:)
real*8           , allocatable , intent(out) :: Xi(:,:)

!local variables ...
integer :: i , j

allocate ( Xi(mm,mm) )

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

do j = 1 , mm
do i = j , mm

         Xi(i,j) = X_ij( i , j , basis )

         Xi(j,i) = Xi(i,j) 

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
end module Surface_Hopping

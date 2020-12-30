! Program for computing Ehrenfest forces from Huckel Hamiltonian
module Ehrenfest_Builder

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only  : driver , verbose , n_part , QMMM
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: EhrenfestForce , store_Hprime

    private

    interface EhrenfestForce 
        module procedure Ehrenfest_Adiabatic
        module procedure Ehrenfest_Diabatic
    end interface EhrenfestForce 

    !module variables ...
    integer                     :: mm 
    integer     , allocatable   :: BasisPointer(:) , DOS(:)
    real*8      , allocatable   :: A_ad_nd(:,:) , B_ad_nd(:,:) , Kernel(:,:) , rho_eh(:,:) , aux(:,:) 
    real*8      , allocatable   :: grad_S(:,:) , F_vec(:) , F_mtx(:,:,:) , H_prime(:,:)
    logical     , allocatable   :: mask(:,:)

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]
    real*8  , parameter :: delta      = 1.d-8

contains
!
!
!
!========================================================================================
 subroutine Ehrenfest_Adiabatic( system , basis , MO_bra , MO_ket , QM , representation )
!========================================================================================
 use MD_read_m              , only  : atom
 implicit none
 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 complex*16                 , intent(in)    :: MO_bra(:,:)
 complex*16                 , intent(in)    :: MO_ket(:,:)
 type(R_eigen)              , intent(in)    :: QM
 character(len=*)           , intent(in)    :: representation

! local variables ... 
 integer                  :: i , j , nn
 real*8     , allocatable :: X_ij(:,:) 
 complex*16 , allocatable :: AO_bra(:,:) , AO_ket(:,:)

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 logical , parameter :: T_ = .true. , F_ = .false.

!================================================================================================
! some preprocessing ...
!================================================================================================
mm = size(basis)
nn = n_part

If( .NOT. allocated(F_mtx ) ) allocate( F_mtx   (system%atoms,system%atoms,3) , source=D_zero )
If( .NOT. allocated(F_vec ) ) allocate( F_vec   (system%atoms)                , source=D_zero )
If( .NOT. allocated(rho_eh) ) then
    allocate( grad_S  (mm,10) )
    allocate( rho_eh  (mm,mm) )
    allocate( A_ad_nd (mm,mm) )
    allocate( B_ad_nd (mm,mm) )
    allocate( aux     (mm,mm) )
    allocate( Kernel  (mm,mm) )
end if

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis )
CALL preprocess    ( system )

select case( representation )

    case( "MO" )
         ! build up electron-hole density matrix ...
         forall( i=1:mm , j=1:mm ) rho_eh(i,j) = real( MO_ket(j,1)*MO_bra(i,1) - MO_ket(j,2)*MO_bra(i,2) )
         aux   = transpose(rho_eh)
         rho_eh = ( rho_eh + aux ) / two 
        
         CALL symm( rho_eh , QM%L , aux )
         CALL gemm( QM%L , aux , A_ad_nd , 'T' , 'N' )
        
         forall( j=1:mm ) rho_eh(:,j) = QM%erg(j) * rho_eh(:,j) 
         CALL gemm( rho_eh , QM%L , aux )
         CALL gemm( aux , QM%L  , B_ad_nd , 'T' , 'N' )

         CALL Huckel_stuff( basis , X_ij )

    case( "AO" )
         If( .NOT. allocated(AO_bra) ) then
              allocate( AO_bra  (mm,n_part) )
              allocate( AO_ket  (mm,n_part) )
         end If
         ! coefs of <k(t)| in AO basis 
         CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , QM%L , mm , MO_bra , mm , C_zero , AO_bra , mm )
         ! coefs of |k(t)> in AO basis 
         CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , QM%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

         ! build up electron-hole density matrix ...
         forall( i=1:mm , j=1:mm ) rho_eh(i,j) = real( AO_ket(j,1)*AO_bra(i,1) - AO_ket(j,2)*AO_bra(i,2) )
         aux     = transpose(rho_eh)
         A_ad_nd = ( rho_eh + aux ) / two 

         CALL Huckel_stuff( basis , X_ij )
         CALL S_invH( system , basis , X_ij )

         CALL gemm( H_prime , A_ad_nd , B_ad_nd )

end select

Kernel = X_ij * A_ad_nd - B_ad_nd
!
!================================================================================================

! set all forces to zero beforehand ...
forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero

! Run, Forrest, Run ...
select case( driver )
    case( "slice_AO" , "slice_FSSH")

        do i = 1 , system% atoms
            If( system%QMMM(i) == "MM" .OR. system%flex(i) == F_ ) cycle
            atom(i)% Ehrenfest = Ehrenfest( system, basis, i ) * eVAngs_2_Newton 
        end do

end select

deallocate( mask , X_ij , F_vec , F_mtx )

include 'formats.h'

end subroutine Ehrenfest_Adiabatic
!
!
!
!
!=================================================================
 subroutine Ehrenfest_Diabatic( system , basis , AO_bra , AO_ket )
!=================================================================
 use MD_read_m   , only  : atom
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(in)    :: basis(:)
 complex*16      , intent(in)    :: AO_bra(:,:)
 complex*16      , intent(in)    :: AO_ket(:,:)

! local variables ... 
 integer               :: i , j , nn
 real*8  , allocatable :: X_ij(:,:) 

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
 logical , parameter :: T_ = .true. , F_ = .false.

!================================================================================================
! some preprocessing ...
!================================================================================================

mm = size(basis)
nn = n_part

If( .NOT. allocated(F_mtx ) ) allocate( F_mtx   (system%atoms,system%atoms,3) , source=D_zero )
If( .NOT. allocated(F_vec ) ) allocate( F_vec   (system%atoms)                , source=D_zero )
If( .NOT. allocated(rho_eh) ) then
    allocate( grad_S  (mm,10) )
    allocate( rho_eh  (mm,mm) )
    allocate( A_ad_nd (mm,mm) )
    allocate( B_ad_nd (mm,mm) )
    allocate( Kernel  (mm,mm) )
end if

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis )
CALL preprocess    ( system )

! build up electron-hole density matrix ...
forall( i=1:mm , j=1:mm ) rho_eh(i,j) = real( AO_ket(j,1)*AO_bra(i,1) - AO_ket(j,2)*AO_bra(i,2) )

CALL Huckel_stuff( basis , X_ij )

A_ad_nd = ( rho_eh + transpose(rho_eh) ) / two 

CALL gemm( H_prime , A_ad_nd , B_ad_nd )

Kernel = X_ij * A_ad_nd - B_ad_nd
!
!================================================================================================

! set all forces to zero beforehand ...
forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero

! Run, Forrest, Run ...
do i = 1 , system% atoms
    If( system%QMMM(i) == "MM" .OR. system%flex(i) == F_ ) cycle
    atom(i)% Ehrenfest = Ehrenfest( system, basis, i ) * eVAngs_2_Newton 
end do

deallocate( mask , X_ij , F_vec , F_mtx )

include 'formats.h'

end subroutine Ehrenfest_Diabatic
!
!
!
!=======================================================
 function Ehrenfest( system, basis, site ) result(Force)
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

end do 

! recover original system ...
system% coord (K,:) = tmp_coord
        
deallocate(pairs)

end function Ehrenfest
!
!
!
!
!==========================================
 subroutine S_invH( system , basis , X_ij )
!==========================================
implicit none
type(structure)     , intent(in)  :: system
type(STO_basis)     , intent(in)  :: basis(:)
real*8              , intent(in)  :: X_ij(:,:)

! local variables...
real*8 , allocatable :: H(:,:) , S(:,:) , S_inv(:,:) 

CALL Overlap_Matrix( system , basis , S )

allocate( H(mm,mm) , source = X_ij*S )

! compute S_inverse...
call Inversion_Matrix( S , S_inv )
deallocate( S )

! allocate and compute H' = S_inv * H ...
If( .not. allocated(H_prime) ) allocate( H_prime (mm,mm) )

CALL symm( S_inv , H , H_prime )

deallocate( S_inv , H )

end subroutine S_invH
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
!
!=================================================
subroutine Inversion_Matrix( matrix , matrix_inv )
!=================================================
implicit none
real*8                  , intent(in)  :: matrix(:,:)
real*8  , allocatable   , intent(out) :: matrix_inv(:,:)

! local variables...
real*8  , allocatable   :: work(:)
integer , allocatable   :: ipiv(:)
integer                 :: i , j , info , N

N = size( matrix(:,1) )

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

end subroutine Inversion_Matrix
!
!
!
!===========================================
 subroutine store_Hprime( N , fetch_Hprime )
!===========================================
implicit none
integer , intent(in) :: N
real*8  , intent(in) :: fetch_Hprime(:,:)

If( .not. allocated(H_prime) ) allocate( H_prime (N,N) )
H_prime = fetch_Hprime

end subroutine store_Hprime
!
!
end module Ehrenfest_Builder

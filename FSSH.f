! Program for computing Surface Hopping forces and transition probabilities from eHuckel Hamiltonian

module Surface_Hopping

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only  : driver, verbose, n_part, QMMM, electron_state, hole_state
    use Structure_Builder       , only  : Unit_Cell 
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: SH_Force 

    private

    !module parameters ...
    integer , parameter   :: xyz_key(3) = [1,2,3]
    real*8  , parameter   :: delta      = 1.d-8
    logical , parameter   :: T_ = .true. , F_ = .false.

    character(len=7), parameter :: method = "Tully"

    !module variables ...
    integer               :: mm , PES(2)
    integer , allocatable :: PB(:) , DOS(:) 
    real*8  , allocatable :: Kernel(:,:) , X_(:,:) , QL(:,:) , Phi(:,:) , erg(:) 
    real*8  , allocatable :: grad_S(:,:) , F_vec(:) , F_mtx(:,:,:) , Rxd_NA(:,:) , pastQR(:,:) 
    logical , allocatable :: mask(:,:)

contains
!
!
!
!=====================================================================
 subroutine SH_Force( system , basis , MO_bra , MO_ket , QM , t_rate )
!=====================================================================
 use MD_read_m     , only  : atom
! args
 implicit none
 type(structure)   , intent(inout) :: system
 type(STO_basis)   , intent(in)    :: basis(:)
 type(R_eigen)     , intent(in)    :: QM
 complex*16        , intent(in)    :: MO_bra(:,:)
 complex*16        , intent(in)    :: MO_ket(:,:)
 real*8            , intent(in)    :: t_rate

! local variables ... 
 integer                   :: i , j , nn , xyz
 real*8, allocatable, save :: rho_eh(:,:) , g_switch(:,:)

! local parameters ...
 real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 

!==================================================================
! some preprocessing ...
!==================================================================
mm = size(basis)
nn = n_part

call setup_Module( system , basis , QM , g_switch , rho_eh )

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis )
CALL preprocess    ( system )

Phi(:,1) = QM%L(PES(1),:)
Phi(:,2) = QM%L(PES(2),:)

do concurrent (j = 1:mm) shared(kernel,X_,Phi,QM)
   kernel(:,j) = ( X_(:,j) - QM%erg(PES(1)) ) * Phi(:,1) * Phi(j,1)   &      ! <== electron part 
               - ( X_(:,j) - QM%erg(PES(2)) ) * Phi(:,2) * Phi(j,2)          ! <== hole part 
   end do

!================================================================================================
! set all forces to zero beforehand ...
do concurrent( i=1:system% atoms ) 
   atom(i)% Ehrenfest(:) = D_zero
   end do

! setup nonadiabtic coupling vector <Psi/dPhi/dt>, used in Tully's method ...
allocate( Rxd_NA(mm,2) , source = D_zero )

! using %Ehrenfest(:) to store the SH force ...
do xyz = 1 , 3
    atom(:)% Ehrenfest(xyz) = SHForce( system , basis , xyz ) * eVAngs_2_Newton 
    end do

! perform FSSH step ...
call FSSH( QM%R , MO_bra , MO_ket , t_rate , rho_eh , g_switch ) 


deallocate( mask , F_vec , F_mtx , QL , erg , Rxd_NA )

include 'formats.h'

end subroutine SH_Force
!
!
!
!====================================================
 function SHForce( system, basis, xyz ) result(Force)
!====================================================
use Semi_empirical_parms , only: atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
integer          , intent(in)    :: xyz 
real*8, dimension(system%atoms)  :: Force

! local variables ...
integer :: i , j , jL , L , indx
integer :: k , ik , DOSk , BPk 

! local arrays ...
integer , allocatable :: pairs(:)
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) 
real*8                :: tmp_coord(3) , delta_b(3) 

 verbose    = .false.

 grad_S = D_zero

 do K = 1 , system% atoms
 
     If( system%QMMM(k) == "MM" .OR. system%flex(k) == F_ ) then
        cycle
     endif
 
     !force on atom site ...
     DOSk = atom( system% AtNo(k) )% DOS
     BPk  = system% BasisPointer(k) 
     
     allocate( pairs , source = pack([( L , L=1,system% atoms )] , mask(:,K)) )
     
     ! save coordinate ...
     tmp_coord = system% coord(k,:)
 
     delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )
  
     system% coord (k,:) = tmp_coord + delta_b
     CALL Overlap_Matrix( system , basis , S_fwd , purpose = "Pulay" , site = K )
 
     system% coord (k,:) = tmp_coord - delta_b
     CALL Overlap_Matrix( system , basis , S_bck , purpose = "Pulay" , site = K )

     ! grad_S is an anti-symmetric matrix 
     do j = 1 , DOSk
        grad_S( BPk+1: , BPk+j  ) = ( S_fwd( BPk+1: , BPk+j ) - S_bck( BPk+1: , BPk+j ) ) / (TWO*delta) 
        grad_S( BPk+j  , BPk+1: ) = -grad_S( BPk+1:,BPk+j )
     end do
 
     !==============================================================================================
     F_vec = D_zero
 
     !$OMP parallel do schedule(dynamic,3) private(iK,jL,i,j,L) default(shared) reduction(+:F_vec)
     do indx = 1 , size(pairs)
       
       L = pairs(indx)
       do jL = 1 , DOS(L)
          j = PB(L) + jL
 
         do iK = 1 , DOSk
            i = BPk + iK
 
            ! adiabatic and non-adiabatic components of the Force ...
            F_vec(L) = F_vec(L) -  grad_S(j,i) * Kernel(i,j)
 
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

     Force(K) = two * sum( F_mtx(K,:,xyz) )

     ! Rxd_NA = dot_product(velocity,force_NA) for El abd Hl 
     ! summing over system%atoms (internal loop) and xyz (external loop)
     If( method == "Tully" ) then
         Rxd_NA = Rxd_NA + TransitionMatrix( grad_S( : , BPk+1:BPk+DOSk ) , k , DOSk , BPk , xyz )
     end If
 
     ! recover original system ...
     system% coord (K,:) = tmp_coord
             
     deallocate(pairs)
     ! ready for next atom in system
end do 
   
end function SHForce
!
!
!
!================================================================
 function TransitionMatrix( grad_Slice , k , DOSk , BPk , xyz ) &
 result(R)
!================================================================
use MD_read_m , only: atom  
implicit none
real*8  , intent(in) :: grad_Slice(:,:)
integer , intent(in) :: k
integer , intent(in) :: DOSk
integer , intent(in) :: BPk
integer , intent(in) :: xyz

! local paranters ...
! convertion factor for nuclear velocity: m/s (MM) to Ang/ps (QM)
real*8  :: V_factor = 1.d-2  

! local variables ... 
integer               :: i , j , j1 , j2 , dima , dimb
real*8  , allocatable :: Mat1(:,:) , A(:,:) , R1(:,:) , R2(:,:)
real*8  , allocatable :: Mat2(:,:) , B(:,:) , R(:,:)

j1 = BPk + 1
j2 = BPk + DOSk

dima = size(grad_Slice(:,1))
dimb = size(grad_Slice(1,:))

! temporary arrays ...
allocate( A(dima,2) , R1(dima,2) , R2(dima,2) , R(dima,2) , Mat2(dima,dima) )

do concurrent (j=1:dima) shared(QL,erg,Mat2)
   Mat2(:,j) = QL(:,j)*erg(:)
   end do

allocate( B(dimb,2) , Mat1(dima,dimb) )
Mat1 = grad_Slice * X_(:,j1:j2)

!===============================================
CALL gemm( Mat1 , Phi(j1:j2,:) , A )
CALL gemm( QL , A , R1 )

CALL gemm( grad_Slice , Phi(j1:j2,:) , A ) 
CALL gemm( Mat2 , A , R2 )

R = R1 - R2
!===============================================

!===============================================
CALL gemm( Mat1 , Phi , B , transa = 'T' ) 
CALL gemm( QL(:,j1:j2) , B , R1 )

do concurrent ( j=1:2 ) shared(erg,PES,Phi,A)
   A(:,j) = Phi(:,j) * erg(PES(j))
   end do

CALL gemm( grad_Slice , A , B , transa = 'T' ) 
CALL gemm( QL(:,j1:j2) , B , R2 )

R = R + (R1-R2)
!===============================================

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! the minus sign guarantees consistency with the Force
R = -R

! checklist
if( abs( R(PES(2),1)-R(PES(1),2) > high_prec ) ) then
    Print*, "WARNING: failed high precision test in TransitionMatrix"
    end if
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

do concurrent ( i=1:dima , j=1:2 , i/=PES(j) ) shared(R,erg,atom)
   R(i,j) = (atom(k)%vel(xyz)*V_factor) * R(i,j) / ( erg(i) - erg(PES(j)) )
   end do

R( PES(1) , 1 ) = 0
R( PES(2) , 2 ) = 0

deallocate( Mat1 , Mat2 , A , B , R1 , R2 )

end function TransitionMatrix
!
!
!
!====================
 function Omega( QR ) 
!====================
implicit none
real*8  , intent(in) :: QR(:,:)

! local variables ...
integer :: i
real*8  , allocatable :: newQR(:,:) , Omega(:,:)
logical               :: flip 
logical , save        :: done = F_

allocate( newQR (mm,2 ) )
allocate( Omega (mm,2 ) )

if( .not. done ) then
    ! setup environment ...
    allocate( pastQR (mm,mm) , source=QR )
    done = T_
else
    ! calculate g_switch via Scattering Matrix (Omega): DynEMol method ...
    do concurrent (i=1:mm) shared(QR,pastQR) local(flip)
       flip = dot_product( QR(:,i) , pastQR(:,i) ) < 0
       if(flip) pastQR(:,i) = -pastQR(:,i)
       end do

    newQR = QR(:,PES(:))

    call gemm( pastQR , newQR , Omega , 'T' )    

    do i=1,2
       Omega(PES(i),i) = d_zero
       end do

    pastQR = QR

end if

deallocate( newQR )

end function Omega
!
!
!
!===========================================================
 subroutine FSSH( QR , MO_bra , MO_ket , t_rate , rho_eh , g_switch )
!===========================================================
implicit none
! args
real*8    , intent(in)   :: QR       (:,:)
complex*16, intent(in)   :: MO_bra   (:,:)
complex*16, intent(in)   :: MO_ket   (:,:)
real*8    , intent(in)   :: t_rate
real*8    , intent(inout):: rho_eh   (:,:)
real*8    , intent(inout):: g_switch (:,:)

! local variables
integer              :: i , j
real*8               :: rn
real*8, allocatable  :: base(:,:)

! this loop: Re(rho_ij)/rho_ii, j=1(el), 2(hl)
do j = 1 , 2 
   rho_eh(:,j) = real( MO_ket(:,j) * MO_bra(PES(j),j) ) 
   rho_eh(:,j) = rho_eh(:,j) / rho_eh( PES(j) , j )
   end do

select case ( method )
    
       case( "Tully" ) 
       g_switch = two * t_rate * rho_eh * Rxd_NA

       case( "Dynemol" )
       g_switch = two * rho_eh * Omega(QR)

       case default
       stop "wrong FSSH method"

       end select


allocate( base(0:mm,2) , source=D_zero )

call random_number(rn)

base(0,:) = D_zero
do j = 1 , 2
do i = 1 , mm
   base(i,j) = base(i-1,j) + max(d_Zero,g_switch(i,j)) 
   if( rn > base(i-1,j) .AND. rn <= base(i,j) ) then
       PES(j) = i     
       cycle
       end if
   end do
   end do









deallocate( base ) 

end subroutine FSSH
!
!
!
!==================================================================
 subroutine setup_Module( system , basis , QM , g_switch , rho_eh )
!==================================================================
implicit none
! args
type(structure)     , intent(in)    :: system
type(R_eigen)       , intent(in)    :: QM
type(STO_basis)     , intent(in)    :: basis(:)
real*8, allocatable , intent(inout) :: g_switch(:,:)
real*8, allocatable , intent(inout) :: rho_eh(:,:)

allocate( F_mtx  (system%atoms,system%atoms,3) , source=D_zero   )
allocate( F_vec  (system%atoms)                , source=D_zero   )
allocate( QL     (mm,mm)                       , source = QM%L   )
allocate( erg    (mm)                          , source = QM%erg )

If( .NOT. allocated(grad_S) ) then

    PES(1) = electron_state
    PES(2) = hole_state

    call init_random_seed()

    allocate( Kernel   (mm,mm) )
    allocate( grad_S   (mm,mm) )
    allocate( Phi      (mm, 2) )
    allocate( rho_eh   (mm, 2) )
    allocate( g_switch (mm, 2) )

    CALL Huckel_stuff( basis , X_ )

    end if

end subroutine setup_Module
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

If( .NOT. allocated(PB) ) allocate( PB(sys%atoms) , DOS(sys%atoms) )

Allocate( mask(sys%atoms,sys%atoms) , source = .false. )

do K = 1   , sys% atoms
   do L = K+1 , sys% atoms
   
       R_LK = sqrt(sum( (sys%coord(K,:)-sys%coord(L,:))**2 ) )
   
       flag1 = R_LK < cutoff_Angs  
        
       flag2 = sys% flex(K) .AND. sys% flex(L)
   
       flag3 = (sys% QMMM(L) == "QM") .AND. (sys% QMMM(K) == "QM")
   
       mask(L,K) = flag1 .AND. flag2 .AND. flag3

   end do
   PB(K) = sys% BasisPointer(K) 
   DOS(K)          = atom( sys% AtNo(K) )% DOS
end do    

end subroutine Preprocess
!
!
!
!==============================
 subroutine init_random_seed ()
!==============================
implicit none

!local variables ...
integer :: seed(5)

seed = [10051965,27092004,2092002,22021967,-76571]

call random_seed(put=seed(1:5))
    
end subroutine init_random_seed
!
!
end module Surface_Hopping

! Program for computing Ehrenfest forces from Huckel Hamiltonian with Coherent-Switch-Decay-of-Mixing
module Ehrenfest_CSDM

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m    , only: verbose, n_part,electron_state, hole_state
    use Overlap_Builder , only: Overlap_Matrix

    public :: Ehrenfest , PST , d_NA_El , d_NA_Hl

    private

    !module variables ...
    integer                                  :: space , PST(2) , Fermi , dim_3N
    integer , allocatable , dimension(:)     :: BasisPointer, DOS
    real*8  , allocatable , dimension(:)     :: erg, F_vec
    real*8  , allocatable , dimension(:,:)   :: Xij, Kernel, grad_S , QL, Phi, d_NA, d_NA_El, d_NA_Hl              
    real*8  , allocatable , dimension(:,:,:) :: F_mtx
    logical , allocatable , dimension(:,:)   :: mask

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

contains
!
!
!
!=========================================================================================
 subroutine Ehrenfest( system , basis , MO_bra , MO_ket , MO_TDSE_bra , MO_TDSE_ket , QM )
!=========================================================================================
 implicit none
 type(structure)  , intent(inout) :: system
 type(STO_basis)  , intent(in)    :: basis(:)
 complex*16       , intent(in)    :: MO_bra(:,:)
 complex*16       , intent(in)    :: MO_ket(:,:)
 complex*16       , intent(in)    :: MO_TDSE_bra(:,:)
 complex*16       , intent(in)    :: MO_TDSE_ket(:,:)
 type(R_eigen)    , intent(in)    :: QM

! local variables ... 
 integer              :: i , j , nn
 real*8 , allocatable :: A_ad_nd(:,:) , B_ad_nd(:,:) , aux(:,:) , rho_EH(:,:) 
 integer              :: GetNewPST(2)
    integer :: oldPST(2)
    logical :: jump
!================================================================================================
! some preprocessing ...
!================================================================================================
space = size(basis)
nn = n_part

CALL setup_Module( system , basis , QM , A_ad_nd , B_ad_nd , rho_EH , aux )

! build up electron-hole density matrix ...
forall( i=1:space , j=1:space ) rho_EH(i,j) = real( MO_ket(j,1)*MO_bra(i,1) - MO_ket(j,2)*MO_bra(i,2) )
aux = transpose(rho_EH)
rho_EH = ( rho_EH + aux ) / two 

CALL symm( rho_EH , QM%L , aux )
CALL gemm( QM%L , aux , A_ad_nd , 'T' , 'N' )

forall( j=1:space ) rho_EH(:,j) = erg(j) * rho_EH(:,j) 
CALL gemm( rho_EH , QM%L , aux )
CALL gemm( aux , QM%L  , B_ad_nd , 'T' , 'N' )

Kernel = Xij * A_ad_nd - B_ad_nd

deallocate( A_ad_nd , B_ad_nd , rho_EH , aux )
!================================================================================================

oldPST=PST
print*, oldPST

GetNewPST = getPointerState( QM%R , MO_TDSE_bra , MO_TDSE_ket )   
PST = GetNewPST

jump = merge( T_ , F_ , any(oldPST /= PST) )
print*, "jump = ", jump
print*, PST

if(jump) pause

CALL EhrenfestForce( system , basis )

deallocate( mask , F_vec , F_mtx , QL , erg , d_NA , Phi )

include 'formats.h'

end subroutine Ehrenfest
!
!
!
!
!===========================================
 subroutine EhrenfestForce( system , basis )
!===========================================
use Semi_empirical_parms , only: ChemAtom => atom
use MD_read_m            , only: atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)

! local parameters ...
integer , parameter :: xyz_key(3) = [1,2,3]
real*8  , parameter :: delta = 1.d-8
real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 

! local variables ...
integer :: i , j , xyz , jL , L , indx , atm_counter
integer :: k , ik , DOS_k , BP_k 

! local arrays ...
integer , allocatable :: pairs(:)
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) , Force(:)
real*8                :: tmp_coord(3) , delta_b(3) 

verbose = .false.

allocate( Force(system%atoms) )

! set atomic forces to zero beforehand ...
forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = d_zero

! Run, Forrest, Run ...
do xyz = 1 , 3

        grad_S = d_zero
        Force  = d_zero 
       
        atm_counter = 0
  
        do k = 1 , system% atoms
        
            If( system%QMMM(k) == "MM" .OR. system%flex(k) == F_ ) cycle

            atm_counter = atm_counter + 1
        
            !force on atom site ...
            DOS_k = ChemAtom( system% AtNo(k) )% DOS
            BP_k  = system% BasisPointer(k) 
        
            allocate( pairs , source = pack([( L , L=1,system% atoms )] , mask(:,K)) )
            
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
            do L = K+1, system% atoms
               F_mtx(K,L,xyz) =   F_vec(L)
               F_mtx(L,K,xyz) = - F_mtx(K,L,xyz) 
            end do
            F_mtx(K,K,xyz) = d_zero
        
            Force(K) = two * sum( F_mtx(K,:,xyz) )
        
            ! calculation of d_NA ...
            d_NA = NAcoupling( grad_S( : , BP_K+1 : BP_K+DOS_k) , DOS_k , BP_K )  ! <== units = 1/Angs
       
            indx = (atm_counter-1)*3 + xyz 
            do concurrent (j=1:space)
               d_NA_El(indx,j) = d_NA(j,1)
               d_NA_Hl(indx,j) = d_NA(j,2)
               enddo
        
            ! recover original system ...
            system% coord (K,:) = tmp_coord
                    
            deallocate(pairs)
            ! ready for next atom in system
        
        end do 

        atom(:)% Ehrenfest(xyz) = Force(:) * eVAngs_2_Newton

end do

end  subroutine EhrenfestForce
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
!
!===============================================================
 function getPointerState( QR , MO_TDSE_bra , MO_TDSE_ket ) &
 result(newPST)
!===============================================================
implicit none
! args
real*8                  , intent(in)  :: QR(:,:)
complex*16              , intent(in)  :: MO_TDSE_bra(:,:)
complex*16              , intent(in)  :: MO_TDSE_ket(:,:)

! local variables ...
integer              :: i , j , newPST(2)
real*8               :: rn
real*8, allocatable  :: rho(:,:) , base(:,:) , g_switch(:,:)

allocate( rho(space, 2) )
! this loop: Symm. Re(rho_ij)/rho_ii, j=1(el), 2(hl)
do j = 1 , 2
   rho(:,j) = real( MO_TDSE_ket(:,j)*MO_TDSE_bra(PST(j),j) + MO_TDSE_ket(PST(j),j)*MO_TDSE_bra(:,j) ) / TWO
   rho(:,j) = rho(:,j) / rho( PST(j) , j )
   end do

allocate(g_switch(space,2))
g_switch(:,:) = two * rho * Omega(QR)

deallocate( rho )

call random_number(rn)

newPST = PST
allocate( base(0:space,2) , source=d_zero )
base(0,:) = d_zero
do j = 1 , 2 
   do i = 1 , space
      base(i,j) = base(i-1,j) + max(d_Zero,g_switch(i,j)) 
      if( rn > base(i-1,j) .AND. rn <= base(i,j) ) then
          newPST(j) = i     
          cycle
          end if
      end do
      end do

if( newPST(1) > Fermi .AND. newPST(2) <= Fermi ) then
         ! do nothing, transitions are allowed
   elseif( newPST(1) == newPST(2) ) then
         ! electron/hole annihilation
         ! system returns to GS
         newPST(1:2) = Fermi 
   elseif( (newPST(1) == PST(2)) .AND. (newPST(2) == PST(1)) ) then
         ! electron/hole exchange transition
         ! system returns to GS
         newPST(1:2) = Fermi 
   else
         ! transitions not allowed
         newPST = PST  
   end if

If( newPST(1) < newPST(2) ) then
    CALL system("sed '11i >>> ATTENTION: electron below hole state <<<' warning.signal |cat")
    stop 
    end If

! to be used in NACoupling ...
Phi(:,1) = QL(newPST(1),:)
Phi(:,2) = QL(newPST(2),:)

deallocate( base , g_switch ) 

end function getPointerState
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

!local saved variables ...
logical               , save :: done = F_
real*8  , allocatable , save :: pastQR(:,:)

allocate(newQR(space,2))
allocate(Omega(space,2))

if( .not. done ) then
    ! setup environment ...
    allocate( pastQR (space,space) , source=QR )
    done = T_
else
    ! used to calculate g_switch via Scattering Matrix (Omega): DynEMol method ...
    do concurrent (i=1:space) shared(QR,pastQR) local(flip)
       flip = dot_product( QR(:,i) , pastQR(:,i) ) < 0
       if(flip) pastQR(:,i) = -pastQR(:,i)
       end do

    newQR = QR(:,PST(:))

    call gemm( pastQR , newQR , Omega , 'T' )    

    !change sign for hole wvpckt ...
    Omega(:,2) = -Omega(:,2)

    do i=1,2
       Omega(PST(i),i) = d_zero
       end do

    pastQR = QR

end if

deallocate( newQR )

end function Omega
!
!
!
!=================================================================================
 subroutine setup_Module( system , basis , QM , A_ad_nd , B_ad_nd , rho_EH , aux )
!=================================================================================
implicit none
type(structure)               , intent(in)    :: system
type(R_eigen)                 , intent(in)    :: QM
type(STO_basis)               , intent(in)    :: basis(:)
real*8          , allocatable , intent(inout) :: A_ad_nd(:,:)
real*8          , allocatable , intent(inout) :: B_ad_nd(:,:) 
real*8          , allocatable , intent(inout) :: rho_EH(:,:) 
real*8          , allocatable , intent(inout) :: aux(:,:) 

! local variables ...

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis ) 
CALL preprocess( system )

allocate( F_mtx(system%atoms,system%atoms,3) , source=d_zero )
allocate( F_vec(system%atoms)                , source=d_zero )

allocate( A_ad_nd (space,space) )
allocate( B_ad_nd (space,space) )
allocate( rho_EH  (space,space) )
allocate( aux     (space,space) )
allocate( Phi     (space, 2) )
allocate( erg     (space)       , source = QM%erg )
allocate( QL      (space,space) , source = QM%L   )
allocate( d_NA    (space, 2)    , source = d_zero )

If( .NOT. allocated(grad_S) ) then
    allocate( grad_S  (space,space) )
    allocate( Kernel  (space,space) )

    dim_3N = 3*count( system%QMMM == "QM" .AND. system%flex == T_ )
    allocate( d_NA_El (dim_3N,space) , source = d_zero )
    allocate( d_NA_Hl (dim_3N,space) , source = d_zero )

    PST(1) = electron_state
    PST(2) = hole_state

    CALL init_random_seed()

    CALL Huckel_stuff( basis , Xij )

    Fermi = QM%Fermi_state
end if

end subroutine setup_Module
!
!
!
!
!=====================================
 subroutine Huckel_stuff( basis , Xij ) 
!=====================================
use Hamiltonians , only : X_ij
implicit none
type(STO_basis)  , intent(in) :: basis(:)
real*8           , allocatable , intent(out) :: Xij(:,:)

!local variables ...
integer :: i , j

allocate ( Xij(space,space) )

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

do j = 1 , space 
do i = j , space 

         Xij(i,j) = X_ij( i , j , basis )
         Xij(j,i) = Xij(i,j) 

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
!
!
end module Ehrenfest_CSDM

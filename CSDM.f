! Program for computing Ehrenfest forces from Huckel Hamiltonian with Coherent-Switch-Decay-of-Mixing
module Ehrenfest_CSDM

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use FMO_m            , only: PointerState
    use Structure_Builder, only: Unit_Cell
    use parameters_m     , only: verbose, n_part,electron_state, hole_state
    use Overlap_Builder  , only: Overlap_Matrix

    public :: Ehrenfest , PST , dNA_El , dNA_Hl , NewPointerState

    private

    !module variables ...
    integer                                            :: dim_E , PST(2) , Fermi , dim_N
    integer           , allocatable , dimension(:)     :: BasisPointer, DOS
    real*8            , allocatable , dimension(:)     :: erg, F_vec
    real*8            , allocatable , dimension(:,:)   :: Xij, Kernel, grad_S , QL, Phi, d_NA
    type(d_NA_vector) , allocatable , dimension(:,:)   :: dNA_El, dNA_Hl              
    real*8            , allocatable , dimension(:,:,:) :: F_mtx
    logical           , allocatable , dimension(:,:)   :: mask

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

contains
!
!
!
!=============================================================
 subroutine Ehrenfest( system , basis , MO_bra , MO_ket , QM )
!=============================================================
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(in)    :: basis(:)
 complex*16      , intent(in)    :: MO_bra(:,:)
 complex*16      , intent(in)    :: MO_ket(:,:)
 type(R_eigen)   , intent(in)    :: QM

! local variables ... 
 integer              :: i , j , nn
 real*8 , allocatable :: A_ad_nd(:,:) , B_ad_nd(:,:) , aux(:,:) , rho_EH(:,:) 

!============================================================
! some preprocessing ...
!============================================================
dim_E = size(basis)
nn = n_part

CALL setup_Module( system , basis , QM , A_ad_nd , B_ad_nd , rho_EH , aux )

! build up electron-hole density matrix ...
forall( i=1:dim_E , j=1:dim_E ) rho_EH(i,j) = real( MO_ket(j,1)*MO_bra(i,1) - MO_ket(j,2)*MO_bra(i,2) )
aux = transpose(rho_EH)
rho_EH = ( rho_EH + aux ) / two 

CALL symm( rho_EH , QM%L , aux )
CALL gemm( QM%L , aux , A_ad_nd , 'T' , 'N' )

forall( j=1:dim_E ) rho_EH(:,j) = erg(j) * rho_EH(:,j) 
CALL gemm( rho_EH , QM%L , aux )
CALL gemm( aux , QM%L  , B_ad_nd , 'T' , 'N' )

Kernel = Xij * A_ad_nd - B_ad_nd

deallocate( A_ad_nd , B_ad_nd , rho_EH , aux )
!============================================================

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
       
            do concurrent (j=1:dim_E)
               dNA_El(atm_counter,j)% vec(xyz) = d_NA(j,1)
               dNA_Hl(atm_counter,j)% vec(xyz) = d_NA(j,2)
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
!
!==============================================================
 subroutine NewPointerState( system , bra , ket , QM , t_rate )
!==============================================================
implicit none
! args
type(structure), intent(in):: system
complex*16     , intent(in):: bra(:,:)
complex*16     , intent(in):: ket(:,:)
type(R_eigen)  , intent(in):: QM
real*8         , intent(in):: t_rate

! local variables ...
integer             :: i , j , newPST(2)
real*8              :: rn , EH_jump
real*8, allocatable :: rho(:,:) , base(:,:) , P_switch(:,:) , B_kl(:,:) , Omega(:,:)
character(len=7)    :: method

allocate( rho(dim_E, 2) )

! this loop: Symm. Re(rho_ij)/rho_ii, j=1(el), 2(hl) ...
do j = 1 , 2
   rho(:,j) = real( ket(:,j)*bra(PST(j),j) + ket(PST(j),j)*bra(:,j) ) / TWO
   rho(:,j) = rho(:,j) / rho( PST(j) , j )
   enddo

!============================================
! choose method = "Dynemol" or "Tully" ...
  method = "Dynemol"
!============================================
! both methods are equivalent ...
allocate(P_switch(dim_E,2))
if ( method == "Dynemol" ) then

   call Dynemol_way(QM,Omega)
   P_switch(:,:) = two * rho * Omega
   deallocate(Omega)

else

   call Tully_way( system , rho , B_kl )
   forall( j=1:2 ) P_switch(:,j) = t_rate * B_kl(:,j) 
   deallocate(B_kl)

end if
!============================================

call random_number(rn)

newPST = PST
allocate( base(0:dim_E,2) , source=d_zero )
base(0,:) = d_zero
do j = 1 , 2 
   do i = 1 , dim_E
      base(i,j) = base(i-1,j) + max(d_Zero,P_switch(i,j)) 
      if( rn > base(i-1,j) .AND. rn <= base(i,j) ) then
          newPST(j) = i     
          exit
          end if
      end do
      end do

EH_jump = (QM%erg(newPST(1)) - QM%erg(PST(1))) - (QM%erg(newPST(2)) - QM%erg(PST(2)))

if( EH_jump > Unit_Cell% MD_Kin ) then
   !transitions are not allowed ; energy forbidden 
   return
   endif

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
   endif

If( newPST(1) < newPST(2) ) then
    CALL warning("ATTENTION: electron below hole state")
    stop 
    end If

deallocate( rho , base , P_switch ) 

PST = newPST

end subroutine NewPointerState
!
!
!
!====================================
 subroutine Dynemol_way( QM , Omega ) 
!====================================
implicit none
type(R_eigen)        , intent(in)  :: QM
real*8 , allocatable , intent(out) :: Omega(:,:)

! local variables ...
integer               :: i , j
real*8  , allocatable :: newQ(:,:) 
logical               :: flip 

!local saved variables ...
logical               , save :: done = F_
real*8  , allocatable , save :: pastQ(:,:)

allocate( newQ  (dim_E,2) )
allocate( Omega (dim_E,2) )

if( .not. done ) then
    ! setup environment ...
    allocate( pastQ (dim_E,dim_E) )
    pastQ = QM%R
    done = T_
else
    ! used to calculate P_switch via Scattering Matrix (Omega): DynEMol method ...
    do i = 1 , dim_E 
       flip = dot_product( QM%L(i,:) , pastQ(:,i) ) < 0
       if(flip) pastQ(:,i) = -pastQ(:,i)
       end do

    forall(j=1:2) newQ(:,j) = QM%L(PST(j),:)

    call gemm( pastQ , newQ , Omega , 'T' )    

    !change sign for hole wvpckt ...
    Omega(:,2) = -Omega(:,2)

    do j=1,2
       Omega(PST(j),j) = d_zero
       end do

    pastQ = QM%R
end if

deallocate( newQ ) 

end subroutine Dynemol_way
!
!
!
!===========================================
 subroutine Tully_way( system , rho , B_kl ) 
!===========================================
use MD_read_m , only: atom
implicit none
type(structure)               , intent(in)  :: system
real*8          , allocatable , intent(in)  :: rho(:,:)
real*8          , allocatable , intent(out) :: B_kl(:,:)

! local parameters ...
real*8, parameter :: V_factor  = 1.d-2   ! <== converts nuclear velocity: m/s (MM) to Ang/ps (QM)

!local variables ...
integer             :: i , j , n
real*8, allocatable :: v_x_dNA(:,:)

!local saved variables ...
logical               , save :: done = F_
real*8  , allocatable , save :: past_rho(:,:) 

allocate( v_x_dNA (dim_E,2) , source = d_zero )
allocate( B_kl    (dim_E,2) )

if( .not. done ) then
    ! setup environment ...
    allocate( past_rho (dim_E,2) )
    past_rho = rho
    done = T_
else
    do i = 1 , dim_E
         do n = 1 , system%atoms
            If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
            v_x_dNA(i,1) = v_x_dNA(i,1) + dot_product( atom(n)% vel(:) , dNA_El(n,i)% vec(:) )
            v_x_dNA(i,2) = v_x_dNA(i,2) - dot_product( atom(n)% vel(:) , dNA_Hl(n,i)% vec(:) )
            end do
            enddo
    v_x_dNA = v_x_dNA * V_factor

    forall( j=1:2 ) B_kl(:,j) = - two * past_rho(:,j) * v_x_dNA(:,j)

    do j=1,2
       B_kl(PST(j),j) = d_zero
       end do

    past_rho = rho
end if

deallocate( v_x_dNA ) 

end subroutine Tully_way
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
integer :: i , j

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis ) 
CALL preprocess( system )

allocate( F_mtx(system%atoms,system%atoms,3) , source=d_zero )
allocate( F_vec(system%atoms)                , source=d_zero )

allocate( A_ad_nd (dim_E,dim_E) )
allocate( B_ad_nd (dim_E,dim_E) )
allocate( rho_EH  (dim_E,dim_E) )
allocate( aux     (dim_E,dim_E) )
allocate( Phi     (dim_E, 2) )
allocate( erg     (dim_E)       , source = QM%erg )
allocate( QL      (dim_E,dim_E) , source = QM%L   )
allocate( d_NA    (dim_E, 2)    , source = d_zero )

If( .NOT. allocated(grad_S) ) then
    allocate( grad_S  (dim_E,dim_E) )
    allocate( Kernel  (dim_E,dim_E) )

    dim_N = count( system%QMMM == "QM" .AND. system%flex == T_ )
    allocate( dNA_El (dim_N,dim_E) )
    allocate( dNA_Hl (dim_N,dim_E) )

    do concurrent( i=1:dim_N , j=1:dim_E )
       dNA_El (i,j)% vec(:) = d_zero
       dNA_Hl (i,j)% vec(:) = d_zero
       enddo

    PST(1) = PointerState(1)
    PST(2) = PointerState(2)

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
!call random_seed()
    
end subroutine init_random_seed
!
!
!
!
end module Ehrenfest_CSDM

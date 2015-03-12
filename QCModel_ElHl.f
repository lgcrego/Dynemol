#include "GPU.h"
 
 module QCModel_Huckel_ElHl

    use type_m
    use omp_lib
    use constants_m
    use parameters_m                , only : DP_Field_  ,       &
                                             Induced_ ,         &
                                             Coulomb_ ,         &
                                             driver ,           &
                                             verbose
    use f95_precision
    use blas95
    use lapack95
    use Overlap_Builder             , only : Overlap_Matrix
    use DP_potential_m              , only : DP_phi
    use Coulomb_SMILES_m            , only : Build_Coulomb_potential
    use DP_main_m                   , only : DP_matrix_AO
    use Polarizability_m            , only : Induced_DP_phi

    public :: EigenSystem_ElHl , Huckel 

    private

    ! module variables ...
    complex*16  , allocatable   :: V_coul(:,:)
    real*8      , allocatable   :: H_DP(:,:)
    logical                     :: PTheory

contains
!
!
!
!==============================================================================================
subroutine EigenSystem_ElHl( system , basis , AO_bra , AO_ket , QM_el , QM_hl , flag1 , flag2 )
!==============================================================================================
implicit none
type(structure)                            , intent(in)    :: system
type(STO_basis)                            , intent(in)    :: basis(:)
complex*16      , optional                 , intent(in)    :: AO_bra(:,:)
complex*16      , optional                 , intent(in)    :: AO_ket(:,:)
type(R_eigen)                              , intent(inout) :: QM_el
type(R_eigen)                              , intent(inout) :: QM_hl
integer         , optional                 , intent(inout) :: flag1
integer         , optional                 , intent(in)    :: flag2

! local variables ...
real*8  , ALLOCATABLE   :: V_coul_El(:) , V_coul_Hl(:) 
real*8  , ALLOCATABLE   :: h(:,:) , S_matrix(:,:) 
integer                 :: i , j 
 
PTheory = present(AO_bra) 

CALL Overlap_Matrix( system , basis , S_matrix )

If ( Coulomb_ ) then

    CALL Build_Coulomb_Potential( system , basis , AO_bra , AO_ket , V_coul , V_coul_El , V_coul_Hl )

else
    ! does not calculate V_coul further ...
    allocate( V_coul_El (size(basis)) , source = D_zero )
    allocate( V_coul_Hl (size(basis)) , source = D_zero )

end If

!-----------------------------------------------------------------------
!           Electron Hamiltonian : upper triangle of V_coul ...

ALLOCATE( h (size(basis),size(basis)) , source = Huckel( basis , S_matrix ) )

if( DP_field_ .OR. Induced_ ) then

    CALL H_DP_Builder( basis , S_matrix )

    h = h + H_DP
    if( PTheory ) then
        forall( j=1:size(basis) ) h(j,j) = h(j,j) + ( AO_bra(j,1)*AO_ket(j,1) * H_DP(j,j) ) + V_coul_El(j)
    else
        forall( j=1:size(basis) ) h(j,j) = h(j,j) + H_DP(j,j) + V_coul_El(j)
    end if

else

    forall( j=1:size(basis) ) h(j,j) = h(j,j) + V_coul_El(j)

end if

! eigensystem for ELECTRON wavepacket ...
CALL Build_MO_basis( h , S_matrix , QM_el , AO_bra , AO_ket , flag1 , flag2 , instance="el" )

!-----------------------------------------------------------------------
!            Hole Hamiltonian : lower triangle of V_coul ...

! re-initialize the hamiltonian ...
h = Huckel( basis , S_matrix )

if( DP_field_ .OR. Induced_ ) then

    CALL H_DP_Builder( basis , S_matrix )

    h = h + transpose(H_DP)
    if( PTheory ) then
        forall( j=1:size(basis) ) h(j,j) = h(j,j) - ( AO_bra(j,2)*AO_ket(j,2) * H_DP(j,j) )  + V_coul_Hl(j)
    else
        forall( j=1:size(basis) ) h(j,j) = h(j,j) + H_DP(j,j) + V_coul_Hl(j)
    end if

else

    forall( j=1:size(basis) ) h(j,j) = h(j,j) + V_coul_Hl(j)

end if

deallocate( V_coul_El , V_coul_Hl )

! eigensystem for HOLE wavepacket ...
CALL Build_MO_basis( h , S_matrix , QM_hl , AO_bra , AO_ket , flag1 , flag2 , instance="hl" )

If( allocated(V_coul) ) deallocate( V_coul )

end subroutine EigenSystem_ElHl
!
!
!
!=================================================================================================
subroutine Build_MO_basis( H_matrix , S_matrix , QM , AO_bra , AO_ket , flag1 , flag2 , instance )
!=================================================================================================
implicit none
real*8                      ,  allocatable  , intent(inout) :: H_matrix(:,:)
real*8                      ,  allocatable  , intent(inout) :: S_matrix(:,:)
type(R_eigen)                               , intent(inout) :: QM
complex*16      , optional                  , intent(in)    :: AO_bra(:,:)
complex*16      , optional                  , intent(in)    :: AO_ket(:,:)
integer         , optional                  , intent(inout) :: flag1
integer         , optional                  , intent(in)    :: flag2
character(*)                                , intent(in)    :: instance

! local variables ...
real*8         , ALLOCATABLE   :: Lv(:,:) , Rv(:,:) 
real*8         , ALLOCATABLE   :: dumb_s(:,:) 
integer                        :: i , info , basis_size
character(1)                   :: uplo

uplo = merge( 'U' , 'L' , instance == "el" )

basis_size = size( H_matrix(:,1) )

ALLOCATE( dumb_s(basis_size,basis_size) )

! clone S_matrix because SYGVD will destroy it ... 
dumb_s = S_matrix

If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(basis_size)) 

CALL SYGVD( H_matrix , dumb_s , QM%erg , 1 , 'V' , uplo , info )

If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '
If ( present(flag1) ) flag1 = info

DEALLOCATE(dumb_s)

!     ---------------------------------------------------
!   ROTATES THE HAMILTONIAN:  H --> H*S_inv 
!
!   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S.|C> 
!
!   Rv = <AO|MO> coefficients
!     ---------------------------------------------------

ALLOCATE( Lv(basis_size,basis_size) )

Lv = H_matrix

If( instance == "hl" ) DEALLOCATE(H_matrix)

! garantees continuity between basis:  Lv(old)  and  Lv(new) ...
If( (driver == "slice_MOt") .AND. (flag2 > 1) ) CALL phase_locking( Lv , QM%R , QM%erg )

ALLOCATE( Rv(basis_size,basis_size) )

CALL gemm(S_matrix,Lv,Rv,'N','N',D_one,D_zero)

If( instance == "hl" ) DEALLOCATE( S_matrix )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

If( .NOT. allocated(QM%L) ) ALLOCATE( QM%L(basis_size,basis_size) ) 
! eigenvectors in the rows of QM%L
QM%L = transpose(Lv)
DEALLOCATE( Lv )

If( .NOT. ALLOCATED(QM%R) ) ALLOCATE( QM%R(basis_size,basis_size) )
! eigenvectors in the columns of QM%R
QM%R = Rv
DEALLOCATE( Rv )

!  the order of storage is the ascending order of eigenvalues
!----------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! Perturbation terms ...
! dipole potential interaction ...
If( PTheory .and. DP_Field_ ) CALL V_DP_off_Diagonal( QM , AO_bra , AO_ket , instance )
If( allocated(H_DP) ) deallocate( H_DP )

! electron-hole interaction ...
If( PTheory .and. Coulomb_  ) CALL V_Coul_off_Diagonal( QM , instance )


! save energies of the TOTAL system 
If( instance == "hl") then
    OPEN(unit=9,file='hl_UNI-ergs.dat',status='unknown')
else
    OPEN(unit=9,file='el_UNI-ergs.dat',status='unknown')
end IF
do i = 1 , basis_size
    write(9,*) i , QM%erg(i)
end do
CLOSE(9)  

If( verbose ) Print*, '>> EigenSystem done <<'

end subroutine Build_MO_basis
!
!
!
!========================================
 pure function Huckel( basis , S_matrix )
!========================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

real*8  , allocatable   :: Huckel(:,:)

! local variables ... 
real*8  :: k_eff , k_WH , c1 , c2 , c3
integer :: i , j

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

ALLOCATE( Huckel(size(basis),size(basis)) , source = D_zero )

do j = 1 , size(basis)

    do i = 1 , j - 1

        c1 = basis(i)%IP - basis(j)%IP
        c2 = basis(i)%IP + basis(j)%IP

        c3 = (c1/c2)*(c1/c2)

        k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

        k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

        Huckel(i,j) = k_eff * S_matrix(i,j) * (basis(i)%IP + basis(j)%IP) / two

        Huckel(j,i) = Huckel(i,j)

    end do

    Huckel(j,j) = basis(j)%IP

end do

end function Huckel
!
!
!
!==========================================
subroutine H_DP_Builder( basis , S_matrix )
!==========================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

! local variables ... 
real*8  :: DP(4) = D_zero
real*8  :: vector(3)
integer :: i , j
logical :: flag

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

ALLOCATE( H_DP(size(basis),size(basis)) , source = D_zero )

!$OMP PARALLEL DO private ( i , j , flag , DP , vector ) schedule( GUIDED , 10 )
do j = 1 , size(basis)
    do i = 1 , j

        flag = ( abs(S_matrix(i,j)) > mid_prec )

        if( flag ) then

            If( Induced_  ) DP = Induced_DP_phi( i , j , basis )

            If( DP_field_ ) DP = DP + DP_phi( i , j , basis )

            vector = DEBYE_inv * DP_matrix_AO(i,j,:) + S_matrix(i,j) * ( Ri(basis(i)) - r0(basis(i),basis(j)) )

            H_DP(i,j) = S_matrix(i,j) * DP(1) + dot_product( vector(1:3) , DP(2:4) )

        end if

    end do
end do  
!$OMP END PARALLEL DO

end subroutine H_DP_Builder
!
!
!
!==================================
pure function r0( basisi , basisj )
!==================================
implicit none
type(STO_basis) , intent(in)    :: basisi
type(STO_basis) , intent(in)    :: basisj

! local variables ...
real*8  :: r0(3)

r0(1) = ( basisi%x + basisj%x ) / two
r0(2) = ( basisi%y + basisj%y ) / two
r0(3) = ( basisi%z + basisj%z ) / two

end function r0
!
!
!
!========================
pure function Ri( basis )
!========================
implicit none
type(STO_basis) , intent(in)    :: basis

! local variables ...
real*8  :: Ri(3)

Ri(1) = basis%x
Ri(2) = basis%y
Ri(3) = basis%z

end function Ri
!
!
!
!=========================================
 subroutine phase_locking( Lv , CR , Erg )
!=========================================
implicit none
real*8      , intent(inout) :: Lv(:,:)
real*8      , intent(in)    :: CR(:,:)
real*8      , intent(inout) :: Erg(:)

! local variables ...
real*8      , allocatable  :: temp_Lv(:,:) , Energies(:) , MO_ovlp(:,:) , old_Rv(:,:)
integer     , allocatable  :: ind(:)
real*8                     :: val
integer                    :: N , i , j , pos

N = size( CR(:,1) )

allocate( old_Rv   ( N , N ) )
allocate( temp_Lv  ( N , N ) )
allocate( MO_ovlp  ( N , N ) )
allocate( Energies ( N     ) )
allocate( ind      ( N     ) )

old_Rv = CR

! MO overlap
CALL gemm(Lv,old_Rv,MO_ovlp,'T','N',D_one,D_zero)

! correction of crossing states ...
temp_Lv  = Lv
Energies = Erg

ind = maxloc( abs(transpose(MO_ovlp)) , dim=1 )

forall(i=1:N)
    Lv(:,i) = temp_Lv(:,ind(i))
    Erg(i)  = Energies(ind(i))
end forall

deallocate( temp_Lv , Energies , MO_ovlp , ind )

! correction of the phases ...
do i = 1 , N
    if( dot_product(Lv(:,i),old_Rv(:,i)) < D_zero ) then
        Lv(:,i) = - Lv(:,i)
    end if
end do

deallocate( old_Rv )

end subroutine phase_locking
!
!
!
!==============================================================
subroutine V_DP_off_diagonal( QM , AO_bra , AO_ket , instance )
!==============================================================
implicit none
type(R_eigen)   , intent(inout) :: QM
complex*16      , intent(in)    :: AO_bra(:,:)
complex*16      , intent(in)    :: AO_ket(:,:)
character*2     , intent(in)    :: instance

! local variables ...
integer                   :: i , j , n_basis
complex*16  , allocatable :: V(:,:) , A(:,:)

n_basis = size(QM%erg)

allocate( V(n_basis,n_basis) , source=C_zero )

select case( instance )

    case( "el" )
        do j = 1 , n_basis
            do i = 1 , j - 1
                V(i,j) = AO_bra(i,1) * AO_ket(j,1) * H_DP(i,j)
                V(j,i) = conjg( V(i,j) )
            end do
            V(j,j) = C_zero
        end do

    case( "hl" )
        do j = 1 , n_basis
            do i = 1 , j - 1
                V(i,j) = - AO_bra(i,2) * AO_ket(j,2) * H_DP(i,j)
                V(j,i) =   conjg( V(i,j) )
            end do
            V(j,j) = C_zero
        end do

end select

deallocate( H_DP )

!===============================================
!         FIRST ORDER Perturbation
!===============================================
! energy correction ...

allocate( A(n_basis,n_basis) , source=C_zero )

CALL DZgemm( 'N' , 'N' , n_basis , n_basis , n_basis , C_one , QM%L , n_basis , V , n_basis , C_zero , A , n_basis )

do i = 1 , n_basis
    QM%erg(i) = QM%erg(i) + real( sum( A(i,:) * QM%L(i,:) ) )
end do

deallocate( A , V )

end subroutine V_DP_off_diagonal
!
!
!
!==============================================
subroutine V_coul_off_diagonal( QM , instance )
!==============================================
implicit none
type(R_eigen)   , intent(inout) :: QM
character*2     , intent(in)    :: instance

! local variables ...
integer                   :: i , j , n_basis
complex*16  , allocatable :: V(:,:) , A(:,:)

n_basis = size(QM%erg)

allocate( V(n_basis,n_basis) , source=V_coul )

select case( instance )

    case( "el" )
        do j = 1 , n_basis
            do i = j+1 , n_basis 
                V(i,j) = conjg( V_coul(j,i) )
            end do
            V(j,j) = C_zero
        end do

    case( "hl" )
        do j = 1 , n_basis
            do i = 1 , j-1
                V(i,j) = conjg( V_coul(j,i) )
            end do
            V(j,j) = C_zero
        end do

end select

!===============================================
!         FIRST ORDER Perturbation
!===============================================
! energy correction ...

allocate( A(n_basis,n_basis) , source=C_zero )

CALL DZgemm( 'N' , 'N' , n_basis , n_basis , n_basis , C_one , QM%L , n_basis , V , n_basis , C_zero , A , n_basis )

do i = 1 , n_basis
    QM%erg(i) = QM%erg(i) + real( sum( A(i,:) * QM%L(i,:) ) )
end do

deallocate( A , V )

end subroutine V_coul_off_diagonal
!
!
!
end module QCModel_Huckel_ElHl

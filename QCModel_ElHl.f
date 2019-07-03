#include "GPU.h"
 
 module QCModel_Huckel_ElHl

    use type_m
    use omp_lib
    use f95_precision
    use blas95
    use lapack95
    use constants_m
    use parameters_m                , only : DP_Field_  ,       &
                                             Induced_ ,         &
                                             Coulomb_ ,         &
                                             driver ,           &
                                             verbose
    use Overlap_Builder             , only : Overlap_Matrix
    use DP_potential_m              , only : DP_phi
    use Coulomb_SMILES_m            , only : Build_Coulomb_potential
    use DP_main_m                   , only : DP_matrix_AO
    use Polarizability_m            , only : Induced_DP_phi
    use QCModel_Huckel              , only : h0 => Huckel,      &
                                             even_more_extended_Huckel

    public :: EigenSystem_ElHl 

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
integer                 :: j , N

N = size(basis)
 
PTheory = present(AO_bra) 

CALL Overlap_Matrix( system , basis , S_matrix )

If ( Coulomb_ ) then

    CALL Build_Coulomb_Potential( system , basis , AO_bra , AO_ket , V_coul , V_coul_El , V_coul_Hl )

else
    ! does not calculate V_coul further ...
    allocate( V_coul_El (N) , source = D_zero )
    allocate( V_coul_Hl (N) , source = D_zero )

end If

!-----------------------------------------------------------------------
!           Electron Hamiltonian : upper triangle of V_coul ...

ALLOCATE( h (N,N) , source = Huckel( basis , S_matrix ) )

if( DP_field_ .OR. Induced_ ) then

    allocate( H_DP(N,N) , source = D_zero )
    H_DP = even_more_extended_Huckel( system , basis , S_matrix )

    h = h + transpose(H_DP)
    if( PTheory ) then
        forall( j=1:N ) h(j,j) = h(j,j) + ( AO_bra(j,1)*AO_ket(j,1) * H_DP(j,j) ) + V_coul_El(j)
    else
        forall( j=1:N ) h(j,j) = h(j,j) + H_DP(j,j) + V_coul_El(j)
    end if

else

    forall( j=1:N ) h(j,j) = h(j,j) + V_coul_El(j)

end if

! eigensystem for ELECTRON wavepacket ...
CALL Build_MO_basis( h , S_matrix , QM_el , AO_bra , AO_ket , flag1 , flag2 , instance="el" )

!-----------------------------------------------------------------------
!            Hole Hamiltonian : lower triangle of V_coul ...

! re-initialize the hamiltonian ...
h = Huckel( basis , S_matrix )

if( DP_field_ .OR. Induced_ ) then

    allocate( H_DP(N,N) , source = D_zero )
    H_DP = even_more_extended_Huckel( system , basis , S_matrix )

    h = h + H_DP
    if( PTheory ) then
        forall( j=1:N ) h(j,j) = h(j,j) - ( AO_bra(j,2)*AO_ket(j,2) * H_DP(j,j) )  + V_coul_Hl(j)
    else
        forall( j=1:N ) h(j,j) = h(j,j) + H_DP(j,j) + V_coul_Hl(j)
    end if

else

    forall( j=1:N ) h(j,j) = h(j,j) + V_coul_Hl(j)

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
integer                        :: i , info , N
character(1)                   :: uplo

uplo = merge( 'U' , 'L' , instance == "el" )

N = size( H_matrix(:,1) )

ALLOCATE( dumb_s(N,N) )

! clone S_matrix because SYGVD will destroy it ... 
dumb_s = S_matrix

If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(N)) 

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

ALLOCATE( Lv(N,N) )

Lv = H_matrix

If( instance == "hl" ) DEALLOCATE(H_matrix)

ALLOCATE( Rv(N,N) )

CALL gemm(S_matrix,Lv,Rv,'N','N',D_one,D_zero)

If( instance == "hl" ) DEALLOCATE( S_matrix )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

If( .NOT. allocated(QM%L) ) ALLOCATE( QM%L(N,N) ) 
! eigenvectors in the rows of QM%L
QM%L = transpose(Lv)
DEALLOCATE( Lv )

If( .NOT. ALLOCATED(QM%R) ) ALLOCATE( QM%R(N,N) )
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
do i = 1 , N
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

! local variables ... 
integer :: i , j , N
real*8  , allocatable   :: Huckel(:,:)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

N = size(basis)
ALLOCATE( Huckel(N,N) , source = D_zero )

do j = 1 , N
    do i = 1 , j

        Huckel(i,j) = h0( i , j , S_matrix(i,j) , basis )

        Huckel(j,i) = Huckel(i,j)

    end do
end do

end function Huckel
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

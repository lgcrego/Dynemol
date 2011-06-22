 module QCModel_Huckel_ElHl

    use type_m
    use omp_lib
    use constants_m
    use parameters_m                , only : DP_Field_  ,       &
                                             driver ,           &
                                             verbose
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Overlap_Builder             , only : Overlap_Matrix
    use DP_potential_m              , only : DP_phi
    use Coulomb_m                   , only : Build_Coulomb_potential

    public :: EigenSystem_ElHl

    private

 contains
!
!
!
!=============================================================================
 subroutine EigenSystem_ElHl( system , basis , QM_el , QM_hl , flag1 , flag2 )
!=============================================================================
 implicit none
 type(structure)                             , intent(in)    :: system
 type(STO_basis)                             , intent(in)    :: basis(:)
 type(R_eigen)                               , intent(inout) :: QM_el
 type(R_eigen)                               , intent(inout) :: QM_hl
 integer          , optional                 , intent(inout) :: flag1
 integer          , optional                 , intent(in)    :: flag2

! local variables ...
 complex*16 , ALLOCATABLE   :: V_coul(:,:) , V_coul_El(:) , V_coul_Hl(:) 
 real*8     , ALLOCATABLE   :: h(:,:) , S_matrix(:,:)
 integer                    :: i , j 


CALL Overlap_Matrix( system , basis , S_matrix )

CALL Build_Coulomb_Potential( S_matrix , basis , V_coul , V_coul_El , V_coul_Hl )

ALLOCATE( h(size(basis),size(basis)) )

!-----------------------------------------------------------------------

If( DP_field_ ) then
    !$OMP PARALLEL DO schedule( GUIDED , 10 )
    do j = 1 , size(basis)
        do i = 1 , j
     
            h(i,j) = huckel_with_FIELDS(i,j,S_matrix(i,j),basis)

        end do
    end do  
    !$OMP END PARALLEL DO
end If

!-----------------------------------------------------------------------
!           Electron Hamiltonian : upper triangle of V_coul ...

do j = 1 , size(basis)
    do i = 1 , j-1
     
        h(i,j) = huckel(i,j,S_matrix(i,j),basis) + real(V_coul(i,j))

    end do

    h(j,j) = huckel(j,j,S_matrix(j,j),basis) + real(V_coul_El(j))

end do  

! eigensystem for electron-wavepacket ...
CALL Build_MO_basis( h , S_matrix , QM_el , flag1 , flag2 , instance="el" )

print*, maxval(QM_el%erg) , minval(real(V_coul)) , minval(real(V_coul_El))

!-----------------------------------------------------------------------
!            Hole Hamiltonian : lower triangle of V_coul ...

h(:,:) = D_zero

do j = 1 , size(basis)
    do i = j+1 , size(basis)
     
        h(i,j) = huckel(i,j,S_matrix(i,j),basis) + real(V_coul(i,j))

    end do

    h(j,j) = huckel(j,j,S_matrix(j,j),basis) + real(V_coul_Hl(j))

end do  

CALL Build_MO_basis( h , S_matrix , QM_hl , flag1 , flag2 , instance="hl" )

print*, maxval(QM_hl%erg) , minval(real(V_coul)) , minval(real(V_coul_Hl))

!-----------------------------------------------------------------------

deallocate( V_coul , V_coul_El , V_coul_Hl )

end subroutine EigenSystem_ElHl
!
!
!
!================================================================================
 subroutine Build_MO_basis( H_matrix , S_matrix , QM , flag1 , flag2 , instance )
!================================================================================
 implicit none
 real*8                      ,  allocatable  , intent(inout) :: H_matrix(:,:)
 real*8                      ,  allocatable  , intent(inout) :: S_matrix(:,:)
 type(R_eigen)                               , intent(inout) :: QM
 integer          , optional                 , intent(inout) :: flag1
 integer          , optional                 , intent(in)    :: flag2
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

! save energies of the TOTAL system 
 If( instance == "hl") then
    OPEN(unit=9,file='hl-UNI-ergs.dat',status='unknown')
 else
    OPEN(unit=9,file='el-UNI-ergs.dat',status='unknown')
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
!====================================
 pure function Huckel(i,j,S_ij,basis)
!====================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel
 real*8  :: k_eff , k_WH , c1 , c2 , c3

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

 k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

 huckel = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / two

 IF(i == j) huckel = basis(i)%IP

 end function Huckel
!
!
!================================================
 pure function Huckel_with_FIELDS(i,j,S_ij,basis)
!================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel_with_FIELDS
 real*8  :: k_eff , k_WH , c1 , c2 , c3
 logical :: flag

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
    
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

 k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

 huckel_with_FIELDS = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / two

 IF( i == j ) huckel_with_FIELDS = basis(i)%IP 

 flag = ( abs(S_ij) > mid_prec ) 

 IF( flag )  huckel_with_FIELDS = huckel_with_FIELDS + S_ij*DP_phi(i,j,basis)

end function Huckel_with_FIELDS
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
end module QCModel_Huckel_ElHl

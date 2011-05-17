 module QCModel_Huckel

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

    public :: EigenSystem , Huckel , Huckel_with_FIELDS

    private

    interface EigenSystem
        module procedure EigenSystem
        module procedure EigenSystem_just_erg
    end interface

 contains
!
!
!
!=============================================================
 subroutine EigenSystem( system , basis , QM , flag1 , flag2 )
!=============================================================
 implicit none
 type(structure)                             , intent(in)    :: system
 type(STO_basis)                             , intent(in)    :: basis(:)
 type(C_eigen)                               , intent(inout) :: QM
 integer          , optional                 , intent(inout) :: flag1
 integer          , optional                 , intent(in)    :: flag2

! local variables ...
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) 
 real*8  , ALLOCATABLE :: h(:,:) , dumb_s(:,:) , S_matrix(:,:)
 integer               :: i , j , info

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 CALL Overlap_Matrix(system,basis,S_matrix)

 If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(size(basis))) 

 ALLOCATE(h(size(basis),size(basis)),dumb_s(size(basis),size(basis)))

 CALL start_clock
! clone S_matrix because SYGVD will destroy it ... 
 dumb_s = S_matrix

 If( DP_field_ ) then

!$OMP PARALLEL DO schedule( GUIDED , 10 )
    do j = 1 , size(basis)
        do i = 1 , j
     
            h(i,j) = huckel_with_FIELDS(i,j,S_matrix(i,j),basis)

        end do
    end do  
!$OMP END PARALLEL DO

 else
 
    do j = 1 , size(basis)
        do i = 1 , j
     
            h(i,j) = huckel(i,j,S_matrix(i,j),basis)

        end do
    end do  

 end If
 CALL stop_clock("Hamiltonian")

 CALL SYGVD(h,dumb_s,QM%erg,1,'V','U',info)

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

 ALLOCATE(Lv(size(basis),size(basis)))

 Lv = h

 DEALLOCATE(h)

 ! garantees continuity between basis:  Lv(old)  and  Lv(new) ...
 If( (driver == "slice_MOt") .AND. (flag2 > 1) ) CALL phase_locking( Lv , QM%R , QM%erg )

 ALLOCATE(Rv(size(basis),size(basis)))

 CALL gemm(S_matrix,Lv,Rv,'N','N',D_one,D_zero)

 DEALLOCATE( S_matrix )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

 If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(size(basis),size(basis))) 
! eigenvectors in the rows of QM%L
 QM%L = cmplx( transpose(Lv) )
 DEALLOCATE( Lv )

 If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(size(basis),size(basis)))
! eigenvectors in the columns of QM%R
 QM%R = cmplx( Rv )
 DEALLOCATE( Rv )

!  the order of storage is the ascending order of eigenvalues
!----------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! save energies of the TOTAL system 
 OPEN(unit=9,file='system-ergs.dat',status='unknown')
    do i = 1 , size(basis)
        write(9,*) i , QM%erg(i)
    end do
 CLOSE(9)  

 If( verbose ) Print*, '>> EigenSystem done <<'

 end subroutine EigenSystem
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
complex*16  , intent(in)    :: CR(:,:)
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

old_Rv = real(CR)

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
!=======================================================
 subroutine EigenSystem_just_erg( system , basis , erg )
!=======================================================
 implicit none
 type(structure)  , intent(in)    :: system
 type(STO_basis)  , intent(in)    :: basis(:)
 real*8           , intent(out)   :: erg( size(basis) )

! local variables ...
 real*8  , ALLOCATABLE :: h(:,:) , S_matrix(:,:)
 integer               :: i , j , info , N

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 N = size(basis)

 CALL Overlap_Matrix(system,basis,S_matrix)

 ALLOCATE( h(N,N) )

 If( DP_field_ ) then

    do j = 1 , N
        do i = 1 , j
     
            h(i,j) = huckel_with_FIELDS(i,j,S_matrix(i,j),basis)

        end do
    end do  

 else

    do j = 1 , N
        do i = 1 , j
     
            h(i,j) = huckel(i,j,S_matrix(i,j),basis)

        end do
    end do  

 end If

 CALL SYGVD(h,S_matrix,erg,1,'N','U',info)

 If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '

 DEALLOCATE( h , S_matrix )

 end subroutine EigenSystem_just_erg
!
!
!
end module QCModel_Huckel

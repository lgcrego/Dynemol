#include "GPU.h"

 module QCModel_Huckel

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use omp_lib
    use constants_m
    use parameters_m                , only : DP_Field_  ,       &
                                             Induced_ ,         &
                                             driver ,           &
                                             verbose
    use Overlap_Builder             , only : Overlap_Matrix
    use DP_potential_m              , only : DP_phi
    use DP_main_m                   , only : DP_matrix_AO
    use polarizability_m            , only : Induced_DP_phi

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
!=============================================
 subroutine EigenSystem( system , basis , QM )
!=============================================
use Matrix_math
implicit none
type(structure)                             , intent(in)    :: system
type(STO_basis)                             , intent(in)    :: basis(:)
type(R_eigen)                               , intent(inout) :: QM

! local variables ...
real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:)  
real*8  , ALLOCATABLE :: h(:,:) , S_matrix(:,:)
real*8  , ALLOCATABLE :: dumb_S(:,:) , tool(:,:) , S_eigen(:) 
integer               :: i , j , N , info 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

N = size(basis)

CALL Overlap_Matrix( system , basis , S_matrix )

If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(N))

Allocate(      h(N,N) )
Allocate( dumb_S(N,N) )

! clone S_matrix because SYGVD will destroy it ...
dumb_s = S_matrix

If( DP_field_ .OR. Induced_ ) then

!$OMP PARALLEL DO schedule( GUIDED , 10 )
    do j = 1 , N
        do i = j , N
    
            h(i,j) = huckel_with_FIELDS(i,j,S_matrix(i,j),basis)

        end do
    end do  
!$OMP END PARALLEL DO

else
    do j = 1 , N
        do i = j , N

            h(i,j) = huckel(i,j,S_matrix(i,j),basis) 

        end do
    end do
end If

CALL SYGVD( h , dumb_S , QM%erg , 1 , 'V' , 'L' , info )

If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '

select case ( driver ) 

    case ( "slice_AO" )

          !---------------------------------------------------
          !   ROTATES THE HAMILTONIAN:  H --> H*S_inv 
          !
          !   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S.|C> 
          !
          !   Rv = <AO|MO> coefficients
          !---------------------------------------------------

          Allocate( Lv(size(basis),size(basis)) )

          Lv = h

          Deallocate(h)

          Allocate( Rv(size(basis), size(basis)) )

          !CALL gemm(S_matrix,Lv,Rv,'N','N',D_one,D_zero)
          call Multiply( S_matrix, Lv, Rv )

          DEALLOCATE( S_matrix )

          !----------------------------------------------------------
          !  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

          If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(size(basis),size(basis))) 
          ! eigenvectors in the rows of QM%L
          QM%L = transpose(Lv) 
          Deallocate( Lv )

          If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(size(basis),size(basis)))
          ! eigenvectors in the columns of QM%R
          QM%R = Rv
          Deallocate( Rv )

    case ("slice_FSSH" )    

          !--------------------------------------------------------
          ! Overlap Matrix Factorization: S^(1/2) ...
          Allocate( S_eigen(N) )

          dumb_s = S_matrix

          CALL SYEVD(dumb_S , S_eigen , 'V' , 'L' , info)

          Allocate( tool(N,N) , source = transpose(dumb_S) )

          forall( i=1:N ) tool(:,i) = sqrt(S_eigen) * tool(:,i)

          !now S_matrix = S^(1/2) Lowdin Orthogonalization matrix ...
          CALL gemm(dumb_S , tool , S_matrix , 'N' , 'N')

          DEALLOCATE( S_eigen , dumb_S , tool )

          !---------------------------------------------------
          !RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S^(1/2).|C> 
          !
          !normalizes the L&R eigenvectors as < L(i) | R(i) > = 1
          !---------------------------------------------------

          Allocate( Lv(N,N) )
          Allocate( Rv(N,N) )

          Lv = h
          Deallocate( h )

          If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(N,N)) 
          ! eigenvectors in the rows of QM%L
          QM%L = transpose(Lv) 

          ! Rv = S^(1/2) * Lv ...
          CALL symm( S_matrix , Lv , Rv )

          If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(N,N))
          ! eigenvectors in the columns of QM%R
          QM%R = Rv

          Deallocate( Lv , Rv , S_matrix )

end select

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! save energies of the TOTAL system ...
OPEN(unit=9,file='system-ergs.dat',status='unknown')
    do i = 1 , N
        write(9,*) i , QM%erg(i)
    end do
CLOSE(9)  

If( verbose ) Print*, '>> EigenSystem done <<'

end subroutine EigenSystem
!
!
!
!============================================
 pure function Huckel( i , j , S_ij , basis )
!============================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel
 real*8  :: k_eff , k_WH , c1 , c2 , c3 , c4

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 if (i == j) then
    huckel = basis(i)%IP + basis(i)%V_shift
 else
    c1 = basis(i)%IP - basis(j)%IP
    c2 = basis(i)%IP + basis(j)%IP

    c3 = (c1/c2)*(c1/c2)

    c4 = (basis(i)%V_shift + basis(j)%V_shift)*HALF

    k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

    k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

    huckel = k_eff*S_ij*c2/two + c4*S_ij
 endif

 end function Huckel
!
!
!
!========================================================
 pure function Huckel_with_FIELDS( i , j , S_ij , basis )
!========================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: DP(4)
 real*8  :: r0(3) , Ri(3) , vector(3)
 real*8  :: Huckel_with_FIELDS
 real*8  :: k_eff , k_WH , c1 , c2 , c3
 logical :: flag

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
   
 DP = D_zero

 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

 k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

 huckel_with_FIELDS = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / two

 IF( i == j ) huckel_with_FIELDS = basis(i)%IP 

 flag = ( abs(S_ij) > mid_prec ) 

 IF( flag ) then
    
     r0(1) = ( basis(i)%x + basis(j)%x ) / two
     r0(2) = ( basis(i)%y + basis(j)%y ) / two
     r0(3) = ( basis(i)%z + basis(j)%z ) / two

     Ri(1) = basis(i)%x
     Ri(2) = basis(i)%y
     Ri(3) = basis(i)%z

     vector = DEBYE_inv * DP_matrix_AO(i,j,:) + S_ij * ( Ri - r0 )    ! <== in Angs

     If( Induced_  ) DP = Induced_DP_phi( i , j , basis )

     If( DP_field_ ) DP = DP + DP_phi( i , j , basis )

     huckel_with_FIELDS = huckel_with_FIELDS + S_ij*DP(1) + dot_product( vector(1:3) , DP(2:4) )

 end if

end function Huckel_with_FIELDS
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

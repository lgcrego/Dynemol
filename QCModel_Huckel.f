#include "GPU.h"

 module QCModel_Huckel

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use omp_lib
    use constants_m
    use parameters_m     , only : EnvField_ , Induced_ , driver , verbose , restart , SO_coupl , extmagfield
    use Overlap_Builder  , only : Overlap_Matrix
    use Hamiltonians     , only : X_ij , even_more_extended_Huckel , spin_orbit_h , MF_interaction
    use Matrix_Math

    public :: EigenSystem , S_root_inv 

    private

    interface EigenSystem
        module procedure EigenSystem
        module procedure EigenSystem_just_erg
    end interface

    ! module variables ...
    real*8 , allocatable :: S_root_inv(:,:) 

 contains
!
!
!
!==================================================
 subroutine EigenSystem( system , basis , QM , it )
!==================================================
use Matrix_math
implicit none
type(structure)             , intent(in)    :: system
type(STO_basis)             , intent(in)    :: basis(:)
type(R_eigen)               , intent(inout) :: QM
integer          , optional , intent(in)    :: it


! local variables ...
real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) 
real*8  , ALLOCATABLE :: h(:,:) , h_SO(:,:) , h_MF(:,:) , S_matrix(:,:) , S_root(:,:)
real*8  , ALLOCATABLE :: dumb_S(:,:) , tool(:,:) , S_eigen(:)
real*8                :: suml , sums , sumj
integer               :: i , j , k , l1 , l2 , N , info 
logical , save        :: first_call_ = .true.

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

N = size(basis)

CALL Overlap_Matrix( system , basis , S_matrix )

If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(N))

Allocate(      h(N,N) )
Allocate( dumb_S(N,N) )

! clone S_matrix because SYGVD will destroy it ...
dumb_s = S_matrix

If( EnvField_ .OR. Induced_ ) then
    h(:,:) = even_more_extended_Huckel( system , basis , S_matrix , it ) 
else
    h(:,:) = Build_Huckel( basis , S_matrix ) 
end If

if( SO_coupl ) then
    CALL spin_orbit_h( basis , h_SO , S_matrix )
    h = h + h_SO
end if

if( extmagfield ) then
    CALL MF_interaction( basis , h_MF , S_matrix )
    h = h + h_MF
end if

CALL SYGVD( h , dumb_S , QM%erg , 1 , 'V' , 'L' , info )
if ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '

select case ( driver ) 

    case default

          !---------------------------------------------------
          !   ROTATES THE HAMILTONIAN:  H --> H*S_inv 
          !
          !   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S.|C> 
          !
          !   normalizes the L&R eigenvectors as < L(i) | R(i) > = 1 
          !---------------------------------------------------

          Allocate( Lv(N,N) )
          Allocate( Rv(N,N) )

          Lv = h
          Deallocate(h)

          If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(N,N)) 
          ! eigenvectors in the rows of QM%L
          QM%L = transpose(Lv) 

          ! Rv = S * Lv ...
          call Multiply( S_matrix, Lv, Rv )

          DEALLOCATE( S_matrix )

          If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(N,N))
          ! eigenvectors in the columns of QM%R
          QM%R = Rv

          Deallocate( Lv , Rv )

    case ("slice_FSSH" )    

          !--------------------------------------------------------
          ! Overlap Matrix Factorization: S^(1/2) ...

          dumb_s = S_matrix

          Allocate( S_eigen(N)  )

          CALL SYEVD(dumb_S , S_eigen , 'V' , 'L' , info)

          Allocate( tool(N,N) , source = transpose(dumb_S) )

          forall( i=1:N ) tool(:,i) = sqrt(S_eigen) * tool(:,i)

          allocate( S_root(N,N) )
          CALL gemm(dumb_S , tool , S_root , 'N' , 'N')

          !now S_root   = S^(1/2) Lowdin Orthogonalization matrix ...
          !now S_matrix = S ...

          DEALLOCATE( S_eigen , tool )
!          DEALLOCATE( S_eigen , dumb_S , tool )

          dumb_S = S_matrix

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
          ! keeping the nonorthogonal representation of %L for future use ...
          QM%L = transpose(Lv) 

          If( first_call_ .AND. (.NOT. restart) ) then

              ! Rv = S * Lv ...
              call symm( S_matrix, Lv, Rv )
              call invert( S_root )
              first_call_ = .false.
          else

              ! Rv = S^(1/2) * Lv ...
              ! Lowding representation ...
              CALL symm( S_root , Lv , Rv )
          end If

          If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(N,N))
          ! eigenvectors in the columns of QM%R
          QM%R = Rv
          
          Deallocate( Lv , Rv , S_matrix , S_root )

end select

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! save energies of the TOTAL system ...
OPEN(unit=9,file='system-ergs.dat',status='unknown')
    do i = 1 , N
        suml = 0.0d0
        sums = 0.0d0
        sumj = 0.0d0
        do j = 1 , N
            suml = suml + QM % L( i , j ) * QM % R( j , i ) * dfloat( basis( j ) % l )
            sums = sums + QM % L( i , j ) * QM % R( j , i ) * HALF * dfloat( basis( j ) % s )
            sumj = sumj + QM % L( i , j ) * QM % R( j , i ) * basis( j ) % j
        end do
!        if( mod(i,2) == 0 ) then
!            write(9,fmt='(i5,4f14.8,i5,a5,i5,f14.8)') i, QM%erg(i), suml, sums, sumj, i, "-->", i - 1, QM % erg( i ) - QM % erg( i - 1 )
!        else
            write(9,fmt='(i5,4f14.8,i5,a5,i5,f14.8)') i, QM%erg(i), suml, sums, sumj
!        end if
    end do
CLOSE(9)  

If( verbose ) Print*, '>> EigenSystem done <<'

end subroutine EigenSystem
!
!
!
!===================================================
 function Build_Huckel( basis , S_matrix ) result(h)
!===================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

! local variables ... 
integer :: i , j , N , N2
real*8  , allocatable   :: h(:,:)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
N = size(basis)

if( SO_coupl ) then
    N2 = 2 * size(basis)
    ALLOCATE( h(N2,N2) , source = D_zero )
else
    ALLOCATE( h(N,N) , source = D_zero )
end if

do j = 1 , N
    do i = j , N

        h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j) 

    end do
end do

if( SO_coupl ) then
 
    do j = N + 1 , N2
        do i = j , N2

            h( i , j ) = h( i - N , j - N )

        end do
    end do

end if

end function Build_Huckel
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
 integer               :: i , j , N , info 

 N = size(basis)

 CALL Overlap_Matrix(system,basis,S_matrix)

 ALLOCATE( h(N,N) )

 If( EnvField_ ) then

    h(:,:) = even_more_extended_Huckel( system , basis , S_matrix ) 

 else

    do j = 1 , N
        do i = 1 , j
     
            h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j)

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
!===========================
 subroutine invert( matrix )
!===========================
implicit none
real*8  , intent(inout) :: matrix(:,:) 

!local variables ...
integer :: N

N = size(matrix(:,1))

CALL syInvert( matrix, return_full ) ! <== matrix content is destroyed and matrix_inv is returned
#define matrix_inv matrix
allocate( S_root_inv(N,N) , source = matrix_inv )
#undef matrix_inv

end subroutine invert
!
!
!
end module QCModel_Huckel

#include "GPU.h"

 module QCModel_Huckel

    use MPI
    use omp_lib
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use MPI_definitions_m           , only : ForceCrew , master 
    use parameters_m                , only : EnvField_ , Induced_ , verbose 
    use Overlap_Builder             , only : Overlap_Matrix
    use Hamiltonians                , only : X_ij , even_more_extended_Huckel

    public :: EigenSystem 

    private

    interface EigenSystem
        module procedure EigenSystem
        module procedure EigenSystem_just_erg
    end interface

 contains
!
!
!
!==================================================
 subroutine EigenSystem( system , basis , QM , it )
!==================================================
use Matrix_math
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)
type(R_eigen)    , intent(inout) :: QM
integer          , optional , intent(in) :: it

! local variables ...
integer               :: i , N , info 
real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) , dumb_S(:,:) , h(:,:) , S_matrix(:,:)

CALL Overlap_Matrix( system , basis , S_matrix )

! After instantiating Overlap_matrix, processes wait outside ...
If( .not. master ) return

N = size(basis)

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

CALL SYGVD( h , dumb_S , QM%erg , 1 , 'V' , 'L' , info )
If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '

!     ---------------------------------------------------
!   RIGHT EIGENVECTOR ALSO CHANGE: C^R = SC 
!
!   normalizes the L&R eigenvectors as < L(i) | R(i) > = 1
!     ---------------------------------------------------

Allocate( Lv(N,N) )
Allocate( Rv(N,N) )

Lv = h
Deallocate( h , dumb_S )

If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(N,N)) 
! eigenvectors in the rows of QM%L
QM%L = transpose(Lv) 

! Rv = S * Lv ...
CALL symm( S_matrix , Lv , Rv )

If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(N,N))
! eigenvectors in the columns of QM%R
QM%R = Rv

Deallocate( Lv , Rv , S_matrix )

! save energies of the TOTAL system 
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
!===================================================
 function Build_Huckel( basis , S_matrix ) result(h)
!===================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

! local variables ... 
integer :: i , j , N
real*8  , allocatable   :: h(:,:)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

N = size(basis)
ALLOCATE( h(N,N) , source = D_zero )

do j = 1 , N
  do i = j , N

        h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j) 

    end do
end do

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
end module QCModel_Huckel

#include "GPU.h"

 module QCModel_Huckel

    use MPI
    use omp_lib
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use MPI_definitions_m           , only : ForceCrew , myEigen , EigenComm , EigenCrew , master 
    use parameters_m                , only : DP_Field_ , Induced_ , verbose 
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
integer               :: mpi_D_R = mpi_double_precision
real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:)  
real*8  , ALLOCATABLE :: h(:,:) , S_matrix(:,:)
real*8  , ALLOCATABLE :: dumb_S(:,:) , tool(:,:) , S_eigen(:) 
integer               :: i , N , info , err , mpi_status(mpi_status_size)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

CALL Overlap_Matrix( system , basis , S_matrix )

! After instantiating Overlap_matrix, processes wait outside ...
If( .not. EigenCrew ) return

N = size(basis)

If( master ) then

     ! Send S_matrix to EigenComm mate ...     
     CALL MPI_Send( S_matrix , N*N , mpi_D_R , 1 , 0 , EigenComm , err )

     If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(N))

     Allocate( h(N,N) )

     If( DP_field_ .OR. Induced_ ) then
         h(:,:) = even_more_extended_Huckel( system , basis , S_matrix , it )
     else
         h(:,:) = Build_Huckel( basis , S_matrix )
     end If

     CALL SYGVD( h , S_matrix , QM%erg , 1 , 'V' , 'L' , info )
     If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '

     ! save energies of the TOTAL system 
     OPEN(unit=9,file='system-ergs.dat',status='unknown')
         do i = 1 , N
             write(9,*) i , QM%erg(i)
         end do
     CLOSE(9)  

else If( myEigen == 1 ) then 

     do  ! <== myEigen = 1 dwells in here Forever ...
         !--------------------------------------------------------
         ! Overlap Matrix Factorization: S^(1/2) ...

         CALL MPI_Recv( S_matrix, N*N, mpi_D_R, 0, mpi_any_tag, EigenComm, mpi_status, err )

         ! clone S_matrix because SYGVD and SYEV will destroy it ... 
         Allocate( dumb_S(N,N) , source = S_matrix )

         Allocate( S_eigen(N) )

         CALL SYEVD(dumb_S , S_eigen , 'V' , 'L' , info)

         Allocate( tool(N,N) , source = transpose(dumb_S) )

         forall( i=1:N ) tool(:,i) = sqrt(S_eigen) * tool(:,i)

         ! now S_matrix = S^(1/2) Lowdin Orthogonalization matrix ...
         CALL gemm(dumb_S , tool , S_matrix , 'N' , 'N')

         DEALLOCATE( S_eigen  )
         DEALLOCATE( dumb_S   )
         DEALLOCATE( tool     )

         CALL MPI_Send( S_matrix , N*N , mpi_D_R , 0 , 0 , EigenComm , err )

     end do  !<== Return to Forever ...

end If

!     ---------------------------------------------------
!   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S^(1/2).|C> 
!
!   normalizes the L&R eigenvectors as < L(i) | R(i) > = 1
!     ---------------------------------------------------

Allocate( Lv(N,N) )
Allocate( Rv(N,N) )

Lv = h
Deallocate( h )

If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(N,N)) 
! eigenvectors in the rows of QM%L
QM%L = transpose(Lv) 

CALL MPI_Recv( S_matrix , N*N , mpi_D_R , 1 , mpi_any_tag , EigenComm , mpi_status , err )

! Rv = S^(1/2) * Lv ...
CALL symm( S_matrix , Lv , Rv )

If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(N,N))
! eigenvectors in the columns of QM%R
QM%R = Rv

Deallocate( Lv , Rv , S_matrix )

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

 If( DP_field_ ) then

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

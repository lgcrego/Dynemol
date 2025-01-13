module Alignment_routines

use f95_precision
use blas95
use lapack95
use types_m
use Constants_m
use util_m         , only : det , Frobenius_norm
use RW_routines    , only : Read_from_XYZ
use GMX_routines   , only : Read_GROMACS , Dump_pdb
use EDT_util_m     , only : parse_this

public :: alignment 

private

! module variables ...
real*8  , allocatable :: Q_all(:,:), P_all(:,:)
integer , allocatable :: indx_list_Q(:) ,  indx_list_P(:)
character(len=30) :: file_1 , file_2

contains

!========================
subroutine Alignment(  )
!========================
implicit none 

! local variables ...
type(universe)  :: sys_Q
type(universe)  :: sys_P
character(len=1)  :: y_or_n

CALL system("clear")

CALL Read_Structures( sys_Q , sys_P )

! assign landmarks ...
write(*,'(/a)',advance='no') 'need to assign landmarks? (y,n=default):  '
read (*,*) y_or_n
Print *, ""

if( y_or_n == "y" ) &
then
    CALL Iteractive_Closest_Point( sys_Q , sys_P )
else
    CALL Procrustes_Analysis( sys_Q , sys_P )
end if

end subroutine Alignment
!
!
!
!
!===================================================
subroutine Iteractive_Closest_Point( sys_Q , sys_P )
!===================================================
implicit none 
integer :: i
type(universe) , intent(inout) :: sys_Q
type(universe) , intent(inout) :: sys_P

!local variables ...
real*8 , allocatable :: Q(:,:) , P(:,:) , U(:,:) , aux(:,:) , diff(:,:)

CALL assign_landmarks( sys_Q , sys_P , Q , P )

allocate( diff , source=Q-P)
write(*,10) " Frobenius distance between ", size(indx_list_Q) , "landmark points = ", Frobenius_norm(diff)

allocate( aux , mold=Q_all )

do i = 1 , 10
       !---------------------------------------------------------------------------
       CALL Procrustes_Rotation_Matrix_U( Q , P , U )
       
       aux = Q_all 
       CALL gemm( U , aux , Q_all , 'N' , 'N' )
       
       Q = Q_all(:,indx_list_Q(:))
       
       !translate to centroid ...
       P_all = translate_to_centroid( P_all , P )
       P     = translate_to_centroid(  P )
       Q_all = translate_to_centroid( Q_all , Q )
       Q     = translate_to_centroid(  Q )
       
       diff = Q - P
       write(*,10) " Frobenius distance between ", size(indx_list_Q) , " landmark points = ", Frobenius_norm(diff)
       !---------------------------------------------------------------------------
       
       CALL get_distance_mtx( sys_Q , sys_P , Q_all , P_all )
       
       CALL Associate_Cloud_Points( Q_all , P_all , Q , P )
end do

CALL Save_Structures( sys_Q , sys_P , Q_all , P_all )

! mission accomplished ...
stop

10 format(a28,I4,a19,F10.6)

end subroutine Iteractive_Closest_Point
!
!
!
!
!==============================================
subroutine Procrustes_Analysis( sys_Q , sys_P )
!==============================================
implicit none 
type(universe) , intent(inout) :: sys_Q
type(universe) , intent(inout) :: sys_P

! local variables ...
integer :: j , N_of_atoms
real*8  , allocatable :: U(:,:) , P(:,:) , Q(:,:) , diff(:,:) , aux(:,:)

! build Q and P matrices ...
N_of_atoms = sys_Q%N_of_atoms
    
allocate( Q(3,N_of_atoms) )
allocate( P(3,N_of_atoms) )
do j=1,N_of_atoms
   Q(:,j) = sys_Q%atom(j)%xyz(:)
   P(:,j) = sys_P%atom(j)%xyz(:)
   end do

allocate( diff , source=Q-P)
Print *, "Frobenius distance of input files = ", Frobenius_norm(diff)

!translate to centroid ...
Q = translate_to_centroid( Q )
P = translate_to_centroid( P )

diff = Q - P
print*, "Frobenius distance of input files at centroid positions = ",  Frobenius_norm(diff)

CALL Procrustes_Rotation_Matrix_U(Q,P,U)

! apply rotation ...
allocate( aux , source=Q )
CALL gemm( U , aux , Q , 'N' , 'N' )

diff = Q - P
print*, "Frobenius distance after Precustes Rotation = ",  Frobenius_norm(diff)

CALL Save_Structures( sys_Q , sys_P , Q , P )

! mission accomplished ...
stop

end subroutine Procrustes_Analysis
!
!
!
!==========================================
subroutine Procrustes_Rotation_Matrix_U( Q , P , U )
!==========================================
implicit none 
real*8                , intent(inout) :: Q(:,:)
real*8                , intent(inout) :: P(:,:)
real*8  , allocatable , intent(out)   :: U(:,:)

!========================================================
!
!     EXPLANATION
!
!========================================================

! local variables ...
integer :: m , n , info
real*8 :: det_V , det_W , det_VW
real*8 , allocatable :: A(:,:) , V(:,:) , W_T(:,:) , S(:) , aux(:,:)

m = size(Q(:,1))
n = size(Q(1,:))

! define matrix A = Q*P^T ...
allocate( A(m,m) )
CALL gemm( Q , P  , A , 'N' , 'T' )

allocate( V  (m,m) )
allocate( W_T(m,m) )
allocate( S  (m)   )
allocate( aux(m,m) )

aux = A
CALL gesvd( aux , S , V , W_T , job='N' , info=info )

det_V = det(V)
det_W = det(transpose(W_T))

det_VW = det_V*det_W

allocate( U (m,m) )

if( det_VW > 0.0 ) then
    aux = transpose(W_T)
else
    aux = transpose(W_T)
    aux(:,m) = -aux(:,m)
end if

CALL gemm( aux , V  , U , 'N' , 'T' )

deallocate( A , V , W_T , S , aux )

end subroutine Procrustes_Rotation_Matrix_U
!
!
!
!===================================================
subroutine get_distance_mtx( sys_Q , sys_P , Q , P )
!===================================================
implicit none 
type(universe) , intent(inout) :: sys_Q
type(universe) , intent(inout) :: sys_P
real*8         , intent(in)    :: Q(:,:)
real*8         , intent(in)    :: P(:,:)

! local variables ...
integer :: i , j , k , n , m
integer , allocatable :: tmp_Q(:) , tmp_P(:)
real*8  , allocatable :: dist_mtx(:,:)

n = size(Q(1,:))
m = size(P(1,:))

allocate( tmp_Q(n) , source=0 )
allocate( tmp_P(m) , source=0 )

allocate(dist_mtx(m,n)) 

k = 0

do j = 1 , n
do i = 1 , m

   if( sys_Q% atom(j)% MMSymbol .NE. sys_P% atom(i)% MMSymbol ) &
   then
       dist_mtx(i,j) = infty
   else
      dist_mtx(i,j) = sum( (Q(:,i) - P(:,j))**2 )
   end if
   dist_mtx(j,i) = dist_mtx(i,j)

   if( dist_mtx(i,j) <= 1.0 ) &
   then
       k = k + 1
       tmp_Q(k) = j
       tmp_P(k) = i 
   end if

end do
end do

deallocate( indx_list_Q , indx_list_P )

allocate( indx_list_Q , source=Pack(tmp_Q,tmp_Q/=0) )
allocate( indx_list_P , source=Pack(tmp_P,tmp_P/=0) )

deallocate( tmp_Q , tmp_P )

end subroutine get_distance_mtx
!
!
!
!=========================================================
subroutine Associate_Cloud_Points( Q_all , P_all , Q , P )
!=========================================================
implicit none 
real*8, intent(in)                  :: Q_all(:,:)
real*8, intent(in)                  :: P_all(:,:)
real*8, allocatable , intent(inout) :: Q(:,:)
real*8, allocatable , intent(inout) :: P(:,:)

! local variables ...
integer :: j , n

n = size(indx_list_Q)

if( allocated(Q) ) deallocate( Q , P )
allocate  ( Q(3,n) , P(3,n) )

! define matrix of landmarks Q and P ...
do j = 1 , n
   Q(:,j) = Q_all(:,indx_list_Q(j))
   P(:,j) = P_all(:,indx_list_P(j))
   end do

end subroutine Associate_Cloud_Points
!
!
!
!================================================
function translate_to_centroid( A , B ) result(C)
!================================================
implicit none 
real*8           , intent(in) :: A(:,:)
real*8, optional , intent(in) :: B(:,:)

!local variables ...
integer :: i , n
real*8  :: centroid(3)
real*8  , allocatable :: C(:,:)

n = size(A(1,:))

! calculate centroids  
if( present(B) ) &
then
    centroid = get_centroid( B )
else
    centroid = get_centroid( A )
end if

allocate( C , mold=A )
!translate to centroid ...
do i=1,n
   C(:,i) = A(:,i) - centroid(:)
   end do

end function translate_to_centroid
!
!
!
!
!===================================================
subroutine Assign_Landmarks( sys_Q , sys_P , Q , P )
!===================================================
implicit none 
type(universe)        , intent(in) :: sys_Q
type(universe)        , intent(in) :: sys_P
real*8  , allocatable , intent(out):: Q(:,:)
real*8  , allocatable , intent(out):: P(:,:)

! local variables ...
integer               :: j , N_of_atoms
logical               :: done
character(len=80)     :: line 

write(*,'(/a)') '>> list atom indices separated by space, press Enter at the end'

done = .false.
do while( .NOT. done ) 
     ! landmarks for file_1 ...
     write(*,'(/3a)') 'landmark atoms for  ', trim(file_1),' :  '
     read (*,'(a)') line
     indx_list_Q = parse_this(line)
     
     ! landmarks for file_2 ...
     write(*,'(/3a)') 'landmark atoms for  ', trim(file_2),' :  '
     read (*,'(a)') line
     indx_list_P = parse_this(line)

     if( size(indx_list_Q) /= size(indx_list_P) ) &
     then
		 write(*,'(/a)') '>> quantity of landmarks does not agree; repeat <<<'
     else 
         done = .true.
     end if
end do

! define matrix of atomic positions for the entire structure ...
N_of_atoms = sys_Q%N_of_atoms
allocate( Q_all(3,N_of_atoms) )
allocate( P_all(3,N_of_atoms) )
do j=1,N_of_atoms
   Q_all(:,j) = sys_Q%atom(j)%xyz(:)
   P_all(:,j) = sys_P%atom(j)%xyz(:)
   end do

! define matrix of landmarks ...
CALL Associate_Cloud_Points( Q_all , P_all , Q , P )

end subroutine Assign_Landmarks
!
!
!
!==========================================
subroutine Read_Structures( sys_Q , sys_P )
!==========================================
implicit none 
type(universe) , intent(out) :: sys_Q
type(universe) , intent(out) :: sys_P

! local variables ...
integer           :: option
character(len=30) :: f_name 

! read structures to align ...
write(*,'(a)',advance='no') '> format of the input files to be aligned: (1) pdb or (2) xyz : '
read(*,*) option

! read file_1 ...
write(*,'(a)',advance='no') 'name of first file (without extension):  '
read (*,*) f_name
file_1 = merge( adjustr(trim(f_name))//".pdb" , adjustr(trim(f_name))//".xyz" , option==1 )
select case (option)
   case(1) ! <== pdb format
            CALL Read_GROMACS( sys_Q , file_name=file_1)
   case(2) ! <== xyz format
            CALL Read_from_XYZ( sys_Q , file_name=file_1)
end select

! read file_2 ...
write(*,'(a)',advance='no') 'name of second file (without extension):  '
read (*,*) f_name
file_2 = merge( adjustr(trim(f_name))//".pdb" , adjustr(trim(f_name))//".xyz" , option==1 )
select case (option)
   case(1) ! <== pdb format
            CALL Read_GROMACS( sys_P , file_name=file_2)
   case(2) ! <== xyz format
            CALL Read_from_XYZ( sys_P , file_name=file_2)
end select

end subroutine Read_Structures
!
!
!
!===================================================
subroutine Save_Structures( sys_Q , sys_P , Q , P )
!===================================================
implicit none 
type(universe) , intent(inout) :: sys_Q
type(universe) , intent(inout) :: sys_P
real*8         , intent(in)    :: Q(:,:)
real*8         , intent(in)    :: P(:,:)

! local variables ...
integer :: j

do j = 1 , sys_Q % N_of_atoms
   sys_Q%atom(j)%xyz(:) = Q(:,j) 
   sys_P%atom(j)%xyz(:) = P(:,j) 
   end do

write(*,'(/a,a,a)') '> ',trim(file_1),' saved to seed-1.pdb'
call Dump_pdb( sys_Q , "seed-1.pdb" )

write(*,'(/a,a,a)') '> ',trim(file_2),' saved to seed-2.pdb'
call Dump_pdb( sys_P , "seed-2.pdb" )

end subroutine Save_Structures
!
!
!
!
!===============================================
pure function get_centroid( A ) result(centroid)
!===============================================
implicit none 
real*8 , intent(in) :: A(:,:)

!local variables ...
integer :: i , n
real*8  :: centroid(3)

n = size(A(1,:))

! calculate centroids  
do i = 1 , 3
   centroid(i) = sum( A(i,:) ) / n
   end do

end function get_centroid
!
!
!
!
end module Alignment_routines

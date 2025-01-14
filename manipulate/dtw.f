module DTW_routines

use types_m
use Constants_m
use util_m      , only : read_general_file , count_lines , read_CSV_file

public :: dtw_driver , dtw_medoid , dba

private

       ! module types ...
       type dtw
            real*8  , allocatable :: seq1(:)
            real*8  , allocatable :: seq2(:)
            real*8  , allocatable :: cost_mtx(:,:) 
            integer , allocatable :: path(:,:)  !<== lowest cost path
            real*8  :: alignment_cost
            integer :: path_size
       endtype dtw
       
       ! module variables ...
       type(dtw) , allocatable :: dtw_pair(:,:)
       
contains
!
!
!
!
!========================
subroutine dba(  )
!========================
implicit none 


! local variables ...
integer :: i , j , k , i1 , i2 , nl , nseqs , iostat , medoid , anomalous , subsequence_size , f_unit, nrows
real*8  :: dumb , aux
real*8    , allocatable :: data_mtx(:,:) , cost_vec(:) , buffer(:)
type(dtw) , allocatable :: dtw_vec(:)
character(len=3)  :: string
character(len=16) :: f_name
character(len=9) , allocatable :: dados(:,:)

!##########################################################3
!nseqs = 100 ! <== number of sequences to be analyzed ...
!nl = count_lines("quantum_E.dat001")
!
!allocate( data_mtx(nl,0:nseqs) )
!
!do j = 1 , nseqs
!
!     write(string,'(i3.3)') j
!     f_name = "quantum_E.dat"//string
!     print*, f_name
!
!     OPEN( unit=3 , file=f_name , status='old' , action="read" )
!
!     do i = 1 , nl
!          read(3,*,IOSTAT=iostat) dumb , data_mtx(i,j) , dumb
!          if (iostat < 0) stop "Unexpected EoF"
!     end do
!
!     close(3)
!
!end do
!##########################################################3
call read_CSV_file( "Trace.dat", dados )

nrows = size( dados(:,1) )
nl    = size( dados(1,:) ) - 1

nseqs = count(dados(:,1)=="1")

allocate( data_mtx(nl,0:nseqs) )

k=0
do j = 1 , nrows
     if( dados(j,1) == "1" ) then
         k = k + 1
         do i = 1 , nl
              read(dados(j,i+1),*) data_mtx(i,k)
              end do
     end if
end do
!###################################################################

call dtw_medoid( data_mtx(:,1:nseqs), medoid, anomalous, cost_vec, instance="light" )
print*, medoid
!data_mtx(:,0) = data_mtx(:,medoid)
data_mtx(:,0) = data_mtx(:,4)
!do i = 1 , nl
!     data_mtx(i,0) = sum(data_mtx(i,1:nseqs))/nseqs
!end do

OPEN(file="traces.dat",status='unknown',newunit=f_unit)
      do i = 1 , nl
         write(f_unit,*) i , data_mtx(i,0)  , data_mtx(i,1) , data_mtx(i,2) , data_mtx(i,3) , data_mtx(i,4)
      end do
close(f_unit)


allocate( dtw_vec(nseqs) )
do j = 1 , nseqs

     allocate(dtw_vec(j)%seq1(nl))
     allocate(dtw_vec(j)%seq2(nl))

     dtw_vec(j)%seq1 = data_mtx(:,0)
     dtw_vec(j)%seq2 = data_mtx(:,j)

     call dynamic_time_warping( dtw_vec(j) , instance="light" )

end do


!OPEN(unit=14,file="dtw1.dat",status='unknown')
!OPEN(unit=24,file="dtw2.dat",status='unknown')
OPEN(unit=34,file="path.dat",status='unknown')
      do i = 1 , dtw_vec(1)%path_size
         write(34,*) dtw_vec(1)%path(i,1) , dtw_vec(1)%path(i,2) , data_mtx(dtw_vec(1)%path(i,2),1)
!         write(14,*) i , dtw_vec(2)%seq1(dtw_vec(2)%path(i,1)) 
!         write(24,*) i , dtw_vec(2)%seq2(dtw_vec(2)%path(i,2))
      end do
!close(14)
!close(24)
close(34)



!================================
! DBA self-consistency loop ...
!================================
OPEN(file="cost.dat",status='unknown',newunit=f_unit)
do k = 1 , 100

       do i = 1 , nl
              aux = 0.0
              subsequence_size = 0
              do j = 1 , nseqs

                     i1 = findloc( dtw_vec(j)%path(:,1) , value=i , dim=1 )
                     i2 = findloc( dtw_vec(j)%path(:,1) , value=i , dim=1 , back=.true. )

                     subsequence_size = subsequence_size + i2-i1+1

                     i1 = dtw_vec(j)%path(i1,2)
                     i2 = dtw_vec(j)%path(i2,2)

!if(i==118)&
!               write(*,'(4I5)') i , i1, i2 , subsequence_size

                     aux = aux + sum( data_mtx(i1:i2,j) )

!if(i==118)&
!              print*, i, data_mtx(i1:i2,j)

              end do
              data_mtx(i,0) = aux / subsequence_size

!if(i==118)&
!              print*, i, data_mtx(i,0), data_mtx(i,0)



       end do

       do j = 1 , nseqs
              dtw_vec(j)%seq1 = data_mtx(:,0)
              dtw_vec(j)%path = 0
              call dynamic_time_warping( dtw_vec(j) , instance="light" )
       end do

       write(f_unit,*) k , sum( dtw_vec(:)%alignment_cost )

end do
close(f_unit)




OPEN(file="avrg.dat",status='unknown',newunit=f_unit)
      do i = 1 , nl
         write(f_unit,*) i , data_mtx(i,0) 
      end do
close(f_unit)



!================================

stop
end subroutine dba
!
!
!
!========================
subroutine dtw_driver(  )
!========================
implicit none 

! local variables ...
integer :: i 
real*8 , allocatable :: sequence1(:,:) , sequence2(:,:)
type(dtw) :: test

call read_general_file( sequence1 )
allocate( test%seq1 , source = sequence1(:,2) ) 

call read_general_file( sequence2 )
allocate( test%seq2 , source = sequence2(:,2) ) 

call dynamic_time_warping( test )

print*, test%alignment_cost

OPEN(unit=14,file="dtw1.dat",status='unknown')
OPEN(unit=24,file="dtw2.dat",status='unknown')
OPEN(unit=34,file="path.dat",status='unknown')
OPEN(unit=44,file="path_cost.dat",status='unknown')
      do i = 1 , test%path_size
         write(34,*) test%path(i,1) , test%path(i,2) 
         write(14,*) i , test%seq1(test%path(i,1)) 
         write(24,*) i , test%seq2(test%path(i,2))
         write(44,*) i , test%cost_mtx( test%path(i,1),test%path(i,2) )
      end do
close(14)
close(24)
close(34)

stop

end subroutine dtw_driver
!
!
!
!
!===========================================================================
subroutine dtw_medoid( data_mtx , medoid , anomalous , cost_vec , instance )
!===========================================================================
implicit none 
real*8                           , intent(in)  :: data_mtx(:,:) 
integer                          , intent(out) :: medoid
integer                          , intent(out) :: anomalous
real*8  , allocatable            , intent(out) :: cost_vec(:)
character(*)          , optional , intent(in)  :: instance

! local variables ...
integer :: i , j , k , nl , nseqs , npairs
integer , allocatable :: indx_table(:,:)

nl     = size(data_mtx(:,1))
nseqs  = size(data_mtx(1,:))
npairs = nseqs*(nseqs-1)/2

! create vector of sequence pairs ...
allocate( dtw_pair(nseqs,nseqs) )
allocate( indx_table(nseqs,nseqs) )

k=0
do i = 1 , nseqs
     do j = i+1 , nseqs

          k = k + 1
          allocate( dtw_pair(i,j)%seq1(nl) , source=data_mtx(:,i))
          allocate( dtw_pair(i,j)%seq2(nl) , source=data_mtx(:,j))

          call dynamic_time_warping( dtw_pair(i,j) , instance="light" )

          dtw_pair(j,i)%alignment_cost = dtw_pair(i,j)%alignment_cost

          indx_table(i,j) = k

     end do
     dtw_pair(i,i)%alignment_cost = 0.0
end do

allocate( cost_vec(nseqs) )
do j = 1 , nseqs
     cost_vec(j) = sum( dtw_pair(:,j)%alignment_cost )
end do

medoid    = minloc(cost_vec,dim=1)
anomalous = maxloc(cost_vec,dim=1)

if( present(instance) ) then
    if( instance == "light" ) deallocate( dtw_pair )
    end if

end subroutine dtw_medoid
!
!
!
!
!==============================================
subroutine dynamic_time_warping( a , instance )
!==============================================
implicit none 
type(dtw)               , intent(inout) :: a
character(*) , optional , intent(in)    :: instance


! local variables ...
logical :: flag
integer :: i , j , k , m , n , maxsize
integer , allocatable :: tmp(:,:)

m = size(a%seq1)
n = size(a%seq2)

allocate( a%cost_mtx( 0:m+1 , 0:n+1 ) , source=infty)
a%cost_mtx(0,0) = 0.d0

do j = 1 , n
   do i = 1 , m 
      a%cost_mtx(i,j) = dist(a%seq1(i),a%seq2(j)) + minval( [ a%cost_mtx(i-1,j-1) , a%cost_mtx(i-1,j) , a%cost_mtx(i,j-1) ] ) 
end do
end do

a%alignment_cost = a%cost_mtx(m,n)/(m+n)

! find lowest-cost path by backtracking the path ...
maxsize = m + n - 1
allocate( tmp(maxsize,2) )

tmp(1,1) = m
tmp(1,2) = n
do k = 2 , maxsize
     i = tmp(k-1,1) - 1
     j = tmp(k-1,2) - 1

     if( a%cost_mtx(i,j) < a%cost_mtx(i+1,j) ) then
         tmp(k,1) = i
     else
         tmp(k,1) = i+1
     endif
     tmp(k,2) = j

     if( a%cost_mtx(tmp(k,1),tmp(k,2)) > a%cost_mtx(i,j+1) ) then
         tmp(k,1) = i
         tmp(k,2) = j+1
     endif

     flag = (tmp(k,1) == 1) .and.  (tmp(k,2) == 1)
     if( flag ) exit
end do

a%path_size = k
a%path = reverse( tmp(1:a%path_size,:) )

deallocate( tmp )

if( present(instance) ) then
    if( instance == "light" ) deallocate(a%cost_mtx)    
    end if



!tmp(1,1) = 1
!tmp(1,2) = 1
!do k = 2 , maxsize
!
!     i = tmp(k-1,1) + 1
!     j = tmp(k-1,2) + 1
!
!     if( a%cost_mtx(i,j) < a%cost_mtx(i-1,j) ) then
!         tmp(k,1) = i
!     else
!         tmp(k,1) = i-1
!     endif
!     tmp(k,2) = j
!
!     if( a%cost_mtx(tmp(k,1),tmp(k,2)) > a%cost_mtx(i,j-1) ) then
!         tmp(k,1) = i
!         tmp(k,2) = j-1
!     endif
!
!     flag = (tmp(k,1) == size_seq1) .and.  (tmp(k,2) == size_seq2)
!     if( flag ) exit
!
!end do
!a%path_size = k
!allocate( a%path , source=tmp(1:a%path_size,1:2) )
!deallocate(tmp)

end subroutine dynamic_time_warping
!
!
!
!
!======================================
function dist( a , b ) result(distance)
!======================================
implicit none 
real*8 , intent(in) :: a , b

!local variables ...
real*8 :: distance

distance = dabs(a-b)

end function dist
!
!
!
!
!==============================
function reverse( a ) result(b)
!==============================
implicit none 
integer, intent(in) :: a(:,:)

!local variables ...
integer :: i , indx , size_a
integer , allocatable :: b(:,:)

allocate( b, mold=a)

size_a = size(a(:,1))

do i = 1 , size_a
     indx = size_a - i + 1
     b(i,:) = a(indx,:)
end do 

end function reverse
!
!
!
!
end module DTW_routines

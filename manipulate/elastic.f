module Elastic_Band

use types_m
use Read_Parms
use RW_routines         , only : Sort_Chemical_Elements
use FUNCTION_routines
use GMX_routines        , only : Dump_PDB

public :: Build_Configurations , Pick_Configurations

private

contains
!
!
!
!===============================================
subroutine Pick_Configurations( trj , file_type)
!===============================================
implicit none
type(universe)  , allocatable   , intent(inout) :: trj(:)
integer                         , intent(in)    :: file_type

! local varibles ...
integer                         :: n , N_of_frames
integer         , allocatable   :: frames(:) , N_of_elements(:) 
character(21)                   :: string
character(4)    , parameter     :: file_format(2)=["vasp","pdb"] 

CALL system( "clear" )

! read frame indices ...
write(*,'(/a)',advance='no') "Number of frames (if > size(trj) select all): "
read (*,'(i3.3)'           ) N_of_frames

! if N_of_frames > size(trj) select 100 snapshots equally spaced ...
N_of_frames = min( N_of_frames,size(trj) )
allocate( frames(N_of_frames) )

If( N_of_frames < size(trj) ) then

    write(*,'(/a)') "Frame indices : "
    do n = 1 , N_of_frames
        read(*,'(i2.2)') frames(n)
    end do

else

    frames = [ ( n , n=1,size(trj) ) ]

end if

! Pre-process system & environment ...
CALL Pre_process( trj(1) , N_of_elements )

! save frames in elastic band format ...
do n = 1 , N_of_frames

    write(string,'(A14,i7.7)')  'frame index = ',frames(n)

    CALL save_Confgs( trj(frames(n)) , N_of_elements , string , n-1 , file_format(file_type) )   

end do
 
stop

end subroutine Pick_Configurations
!
!
!
!=================================================
subroutine Build_Configurations( sys , file_type )
!=================================================
implicit none 
type(universe)      , target , intent(inout) :: sys
integer                      , intent(in)    :: file_type

! local variables ...
integer                             :: i , i1 , i2 , n , N_Confgs 
integer             , allocatable   :: N_of_elements(:) 
real*8                              :: versor(3) , Trans(3)
real*8                              :: bond_stretch , bond_length , norm , distance , step
type(real_interval)                 :: bond
character(len=96)                   :: string
character(4)        , parameter     :: file_format(2)=["vasp","pdb"] 

! input data ...
write(*,'(2/a)'  ) '>>> REMEMBER TO DEFINE FRAGMENT <<< '
write(*,'(3/a)') '> Number of Configurations ? '
write(*,'(a)',advance='no') '# Confgs = '
read (*,'(I3)' ) N_Confgs
!
write(*,'(a)') '> index of the first atom' 
read (*,'(I3)') i1
write(*,'(a)') '> index of the second atom (>> IN THE FRAGMENT <<)'
read (*,'(I3)') i2
!
write(*,'(a)') '> enter minimum bond length (as a Real number) '
read (*,'(F8.4)') bond%inicio
!
write(*,'(a)') '> enter maximum bond length (as a Real number) '
read (*,'(F8.4)') bond%fim

! Pre-process system & environment ...
CALL Pre_process( sys , N_of_elements )

! define translation vector ...
norm = dsqrt( sum( (sys%atom(i2)%xyz(:)-sys%atom(i1)%xyz(:))**2 ) )
versor(:) = ( sys % atom(i2) % xyz(:) - sys % atom(i1) % xyz(:) ) / norm 

! move the fragment to the initial position ...
Trans(:) = versor(:) * (norm-bond%inicio)
forall( i=1:sys%N_of_atoms , sys%atom(i)%fragment=='F' ) sys%atom(i)%xyz(:) = sys%atom(i)%xyz(:) - Trans(:)

! generate configurations ...
bond_length = bond%inicio
write(string,'(A14,F7.4,A2,A72)') 'bond length = ',bond_length,'   ',sys%Surface_Characteristics
CALL save_Confgs( sys , N_of_elements , string , 0 , file_format(file_type) )

bond_stretch = bond%fim - bond%inicio
step = bond_stretch / N_Confgs
do n = 1 , N_Confgs

    forall( i=1:sys%N_of_atoms , sys%atom(i)%fragment=='F' ) sys%atom(i)%xyz(:) = sys%atom(i)%xyz(:) + versor(:)*step

    bond_length = dsqrt( sum( (sys%atom(i2)%xyz(:)-sys%atom(i1)%xyz(:))**2 ) )
   
    write(string,'(A14,F7.4,A2,A72)') 'bond length = ',bond_length,'   ',sys%Surface_Characteristics

    CALL save_Confgs( sys , N_of_elements , string , n , file_format(file_type) )   

end do

stop

end subroutine Build_Configurations
!
!
!
!===============================================================================
subroutine save_Confgs( structure , N_of_elements , title , frame , file_format)
!===============================================================================
implicit none 
type(universe) , intent(inout) :: structure
integer        , intent(in) :: N_of_elements(:)
character(*)   , intent(in) :: title
integer        , intent(in) :: frame
character(*)   , intent(in) :: file_format

!	local variables ...
integer             :: i , j
character(len=3)    :: string
character(len=11)   :: dir
real*8  , parameter :: bottom_spacer = 3.0 !(Angs)
type(universe)      :: sys

write(string,'(i3.3)') frame
dir = "images/"//string//"/"
CALL system("mkdir "//dir)


select case( file_format )

    case('vasp') !-----------------------------------------------------------

        allocate( sys%atom(structure%N_of_atoms) )
        sys = structure

!       sort by chemical element ...
        CALL Sort_Chemical_Elements( sys )

!       rescale coordinates by the unit cell vectors.
        sys%atom%xyz(1) = ( sys%atom%xyz(1) - minval(sys%atom%xyz(1))                 ) / sys%box(1)
        sys%atom%xyz(2) = ( sys%atom%xyz(2) - minval(sys%atom%xyz(2))                 ) / sys%box(2)
        sys%atom%xyz(3) = ( sys%atom%xyz(3) - minval(sys%atom%xyz(3)) + bottom_spacer ) / sys%box(3)

!       ----------------------------------------------
!             generate input file for VASP
!       ----------------------------------------------

        OPEN(unit=4,file=dir//'POSCAR',status='unknown')

        write(4,'(A96)'   ) title
        write(4,'(A2)'    ) '1.' 
        write(4,'(3F12.5)') sys%box(1)     , 0.d0        , 0.d0
        write(4,'(3F12.5)') 0.d0           , sys%box(2)  , 0.d0 
        write(4,'(3F12.5)') 0.d0           , 0.d0        , sys%box(3)  
        write(4,*         ) ( N_of_elements(i),i=1,size(N_of_elements) )
        write(4,'(A18)'   ) 'Selective Dynamics'
        write(4,'(A6)'    ) 'Direct'

        do i = 1 , size(sys%atom)
            write(4,100) (sys%atom(i)%xyz(j),j=1,3) , 'F' , 'F' , 'F' 
        end do

        deallocate( sys%atom )

    case('pdb') !-----------------------------------------------------------

        OPEN(unit=4,file=dir//'seed.pdb',status='unknown')

        CALL Dump_pdb( structure , title )

end select

CLOSE(4)

100   FORMAT(F10.5,F10.5,F10.5,3A3)

end subroutine save_Confgs
!
!
!
!============================================
subroutine Pre_process( sys , N_of_elements )
!============================================
implicit none 
type(universe)                  , intent(inout) :: sys
integer         , allocatable   , intent(out)   :: N_of_elements(:)

! local variables ...
integer                         :: indx
integer     , allocatable       :: temp(:)
type(universe)                  :: dumb
character(len=3)                :: chemical

CALL system( "clear" )

! system information ...
print*, 'Ti = ', count(sys%atom%symbol=='Ti')
print*, 'O  = ', count(sys%atom%symbol=='O')
print*, 'C  = ', count(sys%atom%symbol=='C')
print*, 'H  = ', count(sys%atom%symbol=='H')
print*, 'N  = ', count(sys%atom%symbol=='N')
print*, 'S  = ', count(sys%atom%symbol=='S')
print*, 'Al = ', count(sys%atom%symbol=='Al')

! need to sort by chemical element ...
allocate( dumb%atom(sys%N_of_atoms) )
dumb = sys
CALL Sort_Chemical_Elements( dumb )

print*, dumb%atom%symbol

! find abundance of chemical elements
allocate( temp(20) )
chemical = dumb%atom(1)%symbol
temp(1) = count(dumb%atom%symbol==chemical)
indx = 1
do 
    if( sum(temp(1:indx)) == dumb%N_of_atoms ) EXIT
    chemical = dumb % atom(sum(temp(1:indx))+1) % symbol
    indx = indx + 1
    temp(indx) = count(dumb%atom%symbol==chemical)
end do
allocate( N_of_elements(indx) , source=temp )
deallocate( dumb%atom , temp )

! preparing the environment ...
CALL system("rm -r -f images")
CALL system("mkdir images")

end subroutine pre_process
!
!
!
!--------------------------
subroutine  sort2(ra,ira)
!--------------------------
real*8  , intent(inout) :: ra(:)
integer , intent(inout) :: ira(:)

! local variables
real    :: rra
integer :: irra, l, n, ir, i, j

!----------------------------------------------------------
!  SORT A(I) , SO THAT THE ELEMENTS IA(I) FOLLOW TOGETHER
!----------------------------------------------------------

      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         rra  = ra(l)
         irra = ira(l)
      else
         rra  = ra(ir)
         irra = ira(ir)
         ra(ir)  = ra(1)
         ira(ir) = ira(1)
         ir = ir - 1
         if(ir .eq. 1) then
             ra(1)  = rra
             ira(1) = irra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(ra(j) .lt. ra(j+1)) j = j + 1
        endif
      if(rra .lt. ra(j)) then
        ra(i)  = ra(j)
        ira(i) = ira(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      ra(i)  = rra
      ira(i) = irra
      goto 10

end subroutine sort2
!
!
!
end module Elastic_Band

module crystal_routines

use types_m
use Read_Parms
use Constants_m
use GMX_routines
use RW_routines         , only : Initialize_System

contains
!
!
!
!====================================
 subroutine buildup_crystal( crystal )
!====================================
implicit none 
type(universe)              , intent(inout) :: crystal

! local variables ...
real*8           , allocatable :: basis_coords(:,:)
real*8                         :: Bravais(3,3) , ref(3)
integer                        :: i , j , k , L , m , n(3)
integer                        :: cell_size = 0
character(len=1)               :: string , again
character(len=3) , allocatable :: basis_element(:)
character(len=3)               :: element

CALL system( "clear" )

! system characteristics ...
write(*,'(/a)',advance='no') '> enter System Characteristics: '
read (*,*) crystal% System_Characteristics

do
! define Periodic vectors and store in Bravais matrix ...
    do i = 1 , 3

       write(string,'(i1)') i
       write(*,'(/a)') '> enter A'//string//' crystal axis in cartesian system:'
       write(*,'(a)',advance='no') 'A'//string//'_x = '
       read (*,*) Bravais(1,i)
       write(*,'(a)',advance='no') 'A'//string//'_y = '
       read (*,*) Bravais(2,i)
       write(*,'(a)',advance='no') 'A'//string//'_z = '
       read (*,*) Bravais(3,i)

    end do

    write(*,'(/a)') '> OK ? (y/n)'
    read (*,'(a)') again
    if (again /= 'n') EXIT 
end do

do
! define periodic coordinates ...
    write(*,'(/a)') '> enter maximum coordinates (n1,n2,n3) of Bravais lattice:'
    write(*,'(a)',advance='no') 'n1 = '
    read (*,*) n(1)
    write(*,'(a)',advance='no') 'n2 = '
    read (*,*) n(2)
    write(*,'(a)',advance='no') 'n3 = '
    read (*,*) n(3)

    write(*,'(/a)') '> OK ? (y/n)'
    read (*,'(a)') again
    if (again /= 'n') EXIT 
end do

! choose element of periodic lattice ...
write(*,'(/a)',advance='no') '> Symbol of element at the periodic lattice: '
read (*,*) element

! define crystal basis ...
write(*,'(/a)',advance='no') '> enter size of crystal basis (0 or unit cell size): '
read (*,*) cell_size
allocate( basis_coords (3,cell_size) )
allocate( basis_element  (cell_size) )

do
    do m = 1 , cell_size

       write(*,'(/a)',advance='no') '>  Symbol of element at the basis : '
       read (*,*) basis_element(m)

       write(string,'(i1)') m
       write(*,'(a)') 'atom '//string//' (coords in cartesian system) :'
       write(*,'(a)',advance='no') 'x = '
       read (*,*) basis_coords(1,m)
       write(*,'(a)',advance='no') 'y = '
       read (*,*) basis_coords(2,m)
       write(*,'(a)',advance='no') 'z = '
       read (*,*) basis_coords(3,m)

    end do

    write(*,'(/a)') '> OK ? (y/n)'
    read (*,'(a)') again
    if (again /= 'n') EXIT 
end do

! Building crystal ...
write(*,'(/a)') 'Building the Crystal ...'

crystal% N_of_atoms = (n(1)+1) * (n(2)+1) * (n(3)+1) * m
allocate( crystal% atom( crystal% N_of_atoms) )

CALL Initialize_System( crystal )

crystal% atom% Symbol   = element
crystal% atom% MMSymbol = crystal% atom% Symbol 

L = 0
do i = 0 , n(1)
   do j = 0 , n(2)
      do k = 0 , n(3)

         L = L + 1
         crystal% atom(L)% xyz(1) = Bravais(1,1)*i + Bravais(1,2)*j + Bravais(1,3)*k
         crystal% atom(L)% xyz(2) = Bravais(2,1)*i + Bravais(2,2)*j + Bravais(2,3)*k
         crystal% atom(L)% xyz(3) = Bravais(3,1)*i + Bravais(3,2)*j + Bravais(3,3)*k

         ref = crystal% atom(L)% xyz

         do m = 1 , cell_size
            L = L + 1
            crystal% atom(L)% xyz(:) = ref(:) + basis_coords(:,m)
         end do

      end do
   end do
end do

! define crystal box ...
crystal% box(:) = max( Bravais(:,1)*(n(1)+1) , Bravais(:,2)*(n(2)+1) , Bravais(:,3)*(n(3)+1) )

end subroutine buildup_crystal
!
!
!
end module crystal_routines

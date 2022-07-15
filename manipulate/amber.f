module amber_routines

use types_m
use constants_m

public :: amber_stuff

private

contains
!
!
!========================================
 subroutine amber_stuff
!========================================
implicit none

! local varibles ...
character(len=1)                :: wait
integer                         :: choice

CALL system( "clear" )

! select from options menu ...
do
    write(*,'(/a)') ' (1)  = velcity from input.ncrst file (ASCII format)'
!   write(*,'(/a)') ' (2)  = Translation Operation'
!   write(*,'(/a)') ' (3)  = select frame '      
!   write(*,'(/a)') ' (4)  = save PDB trajectory '
!   write(*,'(/a)') ' (5)  = re-GROUP molecules ( may need to use AD-HOC & Translation before; first frame must be united )'
!   write(*,'(/a)') ' (6)  = DELETE fragment ( uses AD-HOC )'
!   write(*,'(/a)') ' (7)  = RMSD of frames'
!   write(*,'(/a)') ' (8)  = produce trajectory from single PDB frame'
!   write(*,'(/a)') ' (9)  = Interpolate between frames'
!   write(*,'(/a)') ' (10) = Reverse time direction in trajectory'
!   write(*,'(/a)') ' (11) = Replicate structure'
    write(*,'(/a)') ' (0)  = DONE                '
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(I)') choice 

    select case( choice )

        case( 1 ) 

            CALL velocity_from_ncrst

       case default
            exit

    end select

end do    

write(*,'(/a)',advance='no') 'press ENTER '
read (*,'(a)') wait

end subroutine amber_stuff
!
!
!==============================
 subroutine velocity_from_ncrst
!==============================
implicit none

! local variables ...
integer              :: i , j , n , ioerr , N_of_atoms 
real*8               :: dummy_real
character(1)         :: dummy_remark
real*8 , allocatable :: v_xyz(:,:)

OPEN(unit=3,file='input.ncrst',status='old',iostat=ioerr,err=11)

read(3,'(a)') dummy_remark
read(3,*)     N_of_atoms

! going through coordinates ...
do i = 1 , ceiling(N_of_atoms/TWO) 
    read(3,*) dummy_real
end do

! now reading the velocities ...
allocate( v_xyz(N_of_atoms,3) )

do i = 1 , ceiling(N_of_atoms/TWO)-1  ! <== always stop before the last line ...
    read(3,*,iostat=ioerr) ( (v_xyz((i-1)*2+n ,j) , j=1,3) , n=1,2 )
end do
read(3,*,iostat=ioerr) ( (v_xyz((i-1)*2+n ,j) , j=1,3) , n=1,merge(2,mod(N_of_atoms,2),mod(N_of_atoms,2)==0) )

close(3)

! converting from amber units to dynemol units (Angstrom/ps) ...
v_xyz = v_xyz / 20.455d0

!saving velocities ...
OPEN(unit=4,file='velocity_MM.inpt',status='new')
do i = 1 , N_of_atoms
    write(4,*) (v_xyz(i,j) , j=1,3)
end do
write(*,'(/a)') '>>>  Saving velocity_MM.inpt  <<<'

11 if( ioerr > 0 ) stop "input.ncrst file not found; terminating execution"

end subroutine velocity_from_ncrst
!
!
!
end module amber_routines

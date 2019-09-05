module RW_routines

use Read_Parms          , only : atom , atomic_mass , Symbol_2_AtNo , AtNo_2_Symbol
use types_m             , only : universe
use FUNCTION_routines
use diagnosis_m    

public :: Read_from_Poscar ,        &
          Read_from_XYZ ,           &
          Read_from_MD_Urht ,       &
          view_XYZ ,                &
          view_YAEHMOP ,            &
          save_POSCAR ,             &
          save_MD_Urht ,            &
          Sort_Chemical_Elements ,  &
          Initialize_System

private

! module types ...

type date_time
    character(len=8)  :: date
    character(len=10) :: time
end type date_time

type chemicals
    integer           :: N_of_elements
    character(len=3)  :: type_of_elements
end type chemicals

contains
!
!
!
!===================================
subroutine Read_from_poscar(system) 
!===================================
implicit none
type(universe) , intent(out) :: system

! local variables
integer          , allocatable  :: atom_No(:)
integer                         :: i , j , ioerr , InputStatus , indx , N_of_elements
real*8                          :: x0 , y0 , z0 , factor
character(len=3)                :: temp(20) , dumb
character(len=3), allocatable   :: element(:) 

! read System Characteristics and list of chemical elements ...
OPEN(unit=3,file='poscar.dat',form="formatted",access="sequential",status='old',iostat=ioerr,err=10)
read(3 , '(A72)') system%Surface_Characteristics

! read multiplication factor for coordinates ...
read(3 , *) factor

! reads the unit cell vectors ...
read(3,*) x0 , y0 , z0
system%box(1) = dsqrt(x0*x0 + y0*y0 + z0*z0)
read(3,*) x0 , y0 , z0
system%box(2) = dsqrt(x0*x0 + y0*y0 + z0*z0)
read(3,*) x0 , y0 , z0
system%box(3) = dsqrt(x0*x0 + y0*y0 + z0*z0)

system%box = system%box * factor

! reads the number of atoms of each species ...
write(*,'(/a)',advance='no') 'Number of distinct chemical elements (in the order that appears in poscar.dat): '
read(*,*) N_of_elements
allocate( atom_No(N_of_elements) )
allocate( element(N_of_elements) )

read ( 3 , * ) ( element(i), i=1,N_of_elements)
read ( 3 , * ) ( atom_No(i), i=1,N_of_elements)

system%N_of_atoms = sum(atom_No)
allocate( system%atom(system%N_of_atoms) )

CALL Initialize_System( system )

! reads the coordinates from vasp-CONTCAR
do i = 1 , system%N_of_atoms

    read(3,*,iostat=ioerr) (system%atom(i)%xyz(j),j=1,3) , (system%atom(i)%TorF(j),j=1,3)
    print * , i

    do j = 1 , size(element)
        if( (sum(atom_No(:j-1))+1 <= i) .AND. (i <= sum(atom_No(:j))) ) system%atom(i)%symbol = element(j)
    end do

end do
close(3)

! get the Atomic_Number ...
CALL Symbol_2_AtNo(system%atom)

! get the Atomic_Masses ...
system%atom%mass = Atomic_Mass(system%atom%AtNo)

! rescale coordinates by the unit cell vectors'
system%atom%xyz(1) = system%atom%xyz(1)*system%box(1)
system%atom%xyz(2) = system%atom%xyz(2)*system%box(2)
system%atom%xyz(3) = system%atom%xyz(3)*system%box(3)

! fixing atoms to the original unit-cell
!where( system%atom%xyz(3) > 0.8*system%box(3) ) system%atom%xyz(3) = system%atom%xyz(3) - system%box(3)
!system%atom%xyz(3) = system%atom%xyz(3) - minval( system%atom%xyz(3) )

! find the surface_position + add a spacer of 1.0 Angs
!system%Surface_Boundary = maxval(system%atom%xyz(3) , system%atom%symbol=='Ti') + 1.0

10 if( ioerr > 0 )then
    stop 'poscar.dat file not found; terminating execution'
end if

end subroutine Read_from_poscar
!
!
!
!====================================
 subroutine Read_from_MD_Urht(system)
!====================================
implicit none
type(universe) , intent(out) :: system

!	local variables
integer :: i , j , ioerr , N_of_atoms , dummy_int
real*8  :: dummy_real

OPEN(unit=3,file='config.inpt',status='old',iostat=ioerr,err=10)
10 if( ioerr > 0 )then
    stop '"config.inpt" file not found; terminating execution'
end if

! read the unit cell vectors ...
read(3,*) system%box(1) , system%box(2) , system%box(3)

! find the number of atoms N_of_atoms ...
i = 1
do
    read(3,fmt=100,iostat = ioerr) dummy_int, dummy_real, dummy_real, dummy_real, dummy_real
    if(ioerr /= 0) EXIT
    
    i = i + 1
end do

system%N_of_atoms = i-1

! rewind and read configuration ...
rewind 3
read(3,*) dummy_real , dummy_real , dummy_real

allocate( system%atom(system%N_of_atoms) )
CALL Initialize_System( system )

! read the coordinates ...
do i = 1 , system%N_of_atoms 

    read(3,fmt=100,iostat = ioerr) system%atom(i)%AtNo , (system%atom(i)%xyz(j),j=1,3) , system%atom(i)%charge

end do
close(3)

! get chemical symbol ...
CALL AtNo_2_Symbol(system%atom)

system%atom%MMSymbol = system%atom%Symbol

100 format(t4,I6,t13,F11.5,t24,F10.5,t35,F10.5,t46,F10.6)

end subroutine Read_from_MD_Urht
!
!
!
!================================
 subroutine Read_from_XYZ(system)
!================================
implicit none
type(universe) , intent(out) :: system

!	local variables
integer :: i , j , ioerr , N_of_atoms

OPEN(unit=3,file='input.xyz',status='old',iostat=ioerr,err=10)
read(3,*) system%N_of_atoms
read(3,*) system%Surface_Characteristics

! reads the unit cell vectors for Direct coordinate mode
read(3,*) system%box(1)
read(3,*) system%box(2)
read(3,*) system%box(3)

allocate( system%atom(system%N_of_atoms) )
CALL Initialize_System( system )

! reads the coordinates 
do i = 1 , system%N_of_atoms 

    read(3,*,iostat=ioerr) system%atom(i)%symbol , (system%atom(i)%xyz(j),j=1,3) 

    if(ioerr < 0) EXIT
    print*,i

end do
close(3)

! get the Atomic_Number ...
CALL Symbol_2_AtNo(system%atom)

! get the Atomic_Masses ...
system%atom%mass = Atomic_Mass(system%atom%AtNo)

10 if( ioerr > 0 )then
    stop 'input.xyz file not found; terminating execution'
end if

end subroutine Read_from_XYZ
!
!
!
!===========================
 subroutine view_XYZ(system)
!===========================
implicit none 
type(universe)  , intent(inout) :: system

!	local variables
integer             :: i , j
real*8  , parameter :: bottom_spacer = 0.0 !(Angs)

!----------------------------------------------
!     generate input file for XYZ
!----------------------------------------------

!system%atom%xyz(3) = system%atom%xyz(3)-minval(system%atom%xyz(3)) + bottom_spacer

OPEN(unit=9,file='seed.xyz',status='unknown')

write(9,*) size(system%atom)
write(9,10) system%Surface_Characteristics , (system%box(j),j=1,3)
DO i = 1 , size(system%atom)
    write(9,100) system%atom(i)%symbol , (system%atom(i)%xyz(j),j=1,3)
END DO

print*, system%atom%fragment

10    FORMAT(a72,3f10.5)
100   FORMAT(a4,f10.5,f10.5,f10.5)

end subroutine view_XYZ
!
!
!
!===============================
 subroutine view_Yaehmop(system)
!===============================
implicit none 
type(universe) , intent(inout) :: system

!	local variables
integer             :: i , j
real*8  , parameter :: bottom_spacer = 0.0 !(Angs)

!----------------------------------------------
!     generate input file for YAEHMOP
!----------------------------------------------

system%atom%xyz(3) = system%atom%xyz(3)-minval(system%atom%xyz(3)) + bottom_spacer

OPEN(unit=9,file='seed',status='unknown')

write(9,*) system%Surface_Characteristics
write(9,*) " "
write(9,*) "Molecular"      
write(9,*) " "
write(9,*) "Geometry"
write(9,*) " "
write(9,*) size(system%atom)
write(9,*) " " 
DO i = 1 , size(system%atom)
    write(9,101) i, system%atom(i)%symbol , (system%atom(i)%xyz(j),j=1,3)
END DO
write(9,*) " "
write(9,*) "just geom"
CLOSE(9)

print*, system%atom%fragment

101   FORMAT(I5,A4,F10.5,F10.5,F10.5)

end subroutine view_Yaehmop
!
!
!
!=============================
subroutine save_POSCAR(system)
!=============================
implicit none 
type(universe) , intent(inout) :: system

! local variables ...
integer                                     :: i , j , indx 
integer           , allocatable             :: N_of_elements(:) 
character(len=1)                            :: choice
character(len=3)                            :: chemical
character(len=3)  , allocatable             :: type_of_elements(:) 
type(chemicals)   , allocatable             :: structure(:)
real*8                          , parameter :: bottom_spacer = 0.d0 !(Angs)

write(*,'(/a)',advance='no') 'Import T or F flags for selective dynamics from file TorF.dat (n) ? (y/n)  '
read (*,'(a)') choice

if( choice == 'y' ) CALL Import_TAGs(system)

! sort by chemical elemen ...
CALL Sort_Chemical_Elements(system)

! find abundance of chemical elements ...
allocate( structure(20) )
chemical  = system%atom(1)%symbol
structure(1) % N_of_elements    = count(system%atom%symbol==chemical)
structure(1) % type_of_elements = chemical
indx = 1
do 
    if( sum(structure(1:indx)%N_of_elements) == system%N_of_atoms ) EXIT
    chemical = system % atom(sum(structure(1:indx)%N_of_elements)+1) % symbol
    indx = indx + 1
    structure(indx) % N_of_elements    = count(system%atom%symbol==chemical)
    structure(indx) % type_of_elements = chemical
end do

allocate   ( N_of_elements(indx)    , source=structure(:)%N_of_elements    )
allocate   ( type_of_elements(indx) , source=structure(:)%type_of_elements )
deallocate ( structure )

! rescale coordinates by the unit cell vectors ...
system%atom%xyz(1) = ( system%atom%xyz(1) - minval(system%atom%xyz(1))                 ) / system%box(1)
system%atom%xyz(2) = ( system%atom%xyz(2) - minval(system%atom%xyz(2))                 ) / system%box(2)
system%atom%xyz(3) = ( system%atom%xyz(3) - minval(system%atom%xyz(3)) + bottom_spacer ) / system%box(3)

CALL diagnosis( system )

!----------------------------------------------
!     generate input file for VASP
!----------------------------------------------
OPEN(unit=4,file='poscar.xyz',status='unknown')

write(4,10        ) system%Surface_Characteristics , " # " , type_of_elements , "#  "
write(4,'(A2)'    ) '1.' 
write(4,'(3F12.5)') system%box(1)     , 0.d0           , 0.d0
write(4,'(3F12.5)') 0.d0              , system%box(2)  , 0.d0 
write(4,'(3F12.5)') 0.d0              , 0.d0           , system%box(3)  
write(4,*         ) ( N_of_elements(i),i=1,size(N_of_elements) )
write(4,'(A18)'   ) 'Selective Dynamics'
write(4,'(A6)'    ) 'Direct'

do i = 1 , size(system%atom)
    write(4,100) (system%atom(i)%xyz(j),j=1,3) , (system%atom(i)%TorF(j),j=1,3) 
end do

CLOSE(4)

write(*,119) system%atom%symbol

! total number of atoms of given type ...
Print 120 , type_of_elements(1) , N_of_elements(1)
do i = 2 , size(N_of_elements)
    Print 121 ,  type_of_elements(i) , N_of_elements(i)
end do

Print*, ""
Print*, "saved in poscar.xyz"

Print 122 , 'number of atoms allowed to move = ', count(system%atom%TorF(1)=='T')


10    FORMAT(t1,A72,t73,A3,20A3)
100   FORMAT(F10.5,F10.5,F10.5,3A3)
119   FORMAT(/,30a3)
120   FORMAT(/,1x,A2,' atoms  = ',I5)
121   FORMAT(1x,A2,' atoms  = ',I5)
122   FORMAT(/,t2,A34,I5)

end subroutine save_POSCAR
!
!
!
!============================
 subroutine save_MD_Urht(sys)
!============================
implicit none
type(universe) , intent(in) :: sys

! local variables ...
integer :: i , j 

OPEN(unit=3,file='seed.inpt',status='unknown')

! write unit cell vectors ...
write(3,*) (sys%box(i) , i=1,3)

! write the coordinates ...
do i = 1 , sys%N_of_atoms 

    write(3,fmt=100) sys%atom(i)%AtNo , (sys%atom(i)%xyz(j),j=1,3) , sys%atom(i)%charge

end do
close(3)

100 format(t4,I6,t13,F11.5,t24,F10.5,t35,F10.5,t46,F10.6)

end subroutine save_MD_Urht
!
!
!
!==============================
 subroutine Import_TAGs(system)
!==============================
implicit none 
type(universe) , intent(inout) :: system

! local variables ...
integer      :: i , j , N_of_tags, io_err, ioerr
character(20) :: dumb
logical      :: exist

OPEN(unit=33,file='TorF.dat',status='old',iostat=ioerr,err=10)

N_of_tags = 0
do
    if ( io_err < 0 ) exit
    read( unit=33 , fmt=11 , iostat=io_err ) dumb 
    N_of_tags = N_of_tags + 1
end do
N_of_tags = N_of_tags - 1

If( N_of_tags /= system%N_of_atoms ) stop ' >> number of TorF tags differs from N_of_atoms << ' 

rewind(33)
do i = 1 , system % N_of_atoms
    read(33,*) ( system % atom(i) % TorF(j) , j=1,3 ) 
end do

close(33)

10 if( ioerr > 0 ) stop ' >> file TorF.dat not found ! << '
11 format(A)

end subroutine Import_TAGs
!
!
!
!========================================
subroutine Sort_Chemical_Elements(system)
!========================================
implicit none
type(universe) , intent(inout) :: system

!	local variables
integer        :: i , j , N_of_elements , N_of_atoms ,  iptr
type(universe) :: temp

allocate( temp%atom(1) )

N_of_atoms = size(system%atom)

do i = 1 , N_of_atoms-1

    iptr = i

    do j = i+1 , N_of_atoms
        if( LGT( system%atom(j)%symbol , system%atom(iptr)%symbol ) ) then
            iptr = j
        end if
    end do

    if( i /= iptr ) then
        temp%atom(1)      = system%atom(i)
        system%atom(i)    = system%atom(iptr)
        system%atom(iptr) = temp%atom(1)
    end if

end do

end subroutine Sort_Chemical_Elements
!
!
!
!================================
subroutine Initialize_System( a )
!================================
implicit none
type(universe)  :: a

! local variables ...
integer :: i

a % time = 0.d0

forall( i=1:3 ) 
    a % atom % xyz(i)  = 0.d0
    a % atom % TorF(i) = "F"
end forall

a % atom % mass      = 0.d0
a % atom % charge    = 0.d0
a % atom % AtNo      = 0
a % atom % nrcg      = 0
a % atom % nresid    = 0
a % atom % resid     = "XXX"
a % atom % Symbol    = "XXX"
a % atom % MMsymbol  = "XXX"
a % atom % fragment  = "X"
a % atom % copy      = .false.
a % atom % delete    = .false.
a % atom % translate = .false.
a % atom % rotate    = .false.
a % atom % group     = .false.

a % N_of_Surface_Atoms      = 0
a % N_of_Solvent_Atoms      = 0
a % N_of_Solvent_Molecules  = 0

end subroutine Initialize_System
!
!
!
end module RW_routines

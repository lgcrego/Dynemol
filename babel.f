 module Babel_m

    use type_m
    use Allocation_m

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

    character(len=72) , public :: System_Characteristics

    PUBLIC :: Read_from_XYZ , Read_from_Poscar , seed_Yaehmop , seed_VASP 

    PRIVATE

 contains
!
!
!
 subroutine Read_from_XYZ(unit_cell)

 type(structure) , intent(out) :: unit_cell

 character(len=2) :: dumb_char 
 character(len=1) :: fragment
 real*8           :: k_WH
 integer          :: i , j , j1 , j2 , n_kWH , n_fragments , N_of_Configurations

 OPEN(unit=3,file='xyz.dat',status='old')   

! start reading the structure characteristics
 read(3,*) System_Characteristics
 read(3,*) N_of_Configurations      ! <== No of configurations 
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell

! allocating arrays
 CALL Allocate_UnitCell(unit_cell)

! Acceptor (A) , Donor (D) or Molecule (M)
! fragment atoms must be packed together
 read(3,*) n_fragments
 do i = 1 , n_fragments
    read(3,*) j1 , j2 , fragment
    forall(j=j1:j2) unit_cell%fragment(j) = fragment
 end do

! defining the k_WH parameter for the atom 
! No of k_WH to be used 
 read(3,*) n_kWH                    
 do i = 1 , n_kWH
    read(3,*) j1 , j2 , k_WH
    forall(j=j1:j2) unit_cell%k_WH(j) = k_WH
 end do

! reading coordinates
 read(3,*) dumb_char
 DO j = 1 , unit_cell%atoms
     read(3,*) unit_cell%symbol(j), (unit_cell%coord(j,i), i=1,3)
 END DO   

 CLOSE(3)

! unit_cell dimensions
 unit_cell%T_xyz(1) = T_x
 unit_cell%T_xyz(2) = T_y
 unit_cell%T_xyz(3) = T_z

 CALL Symbol_2_AtNo(unit_cell)
 
 print 70, System_Characteristics

 include 'formats.h'

 end subroutine Read_from_XYZ
!
!
!
 subroutine Read_from_Poscar(unit_cell)

 type(structure) , intent(out) :: unit_cell

 real*8            :: x0 , y0 , z0
 real*8            :: a , b , c
 real*8            :: kWH
 integer           :: n_kWH , n_fragments , i , j , j1 , j2 , n , boundary_atom
 integer           :: N_of_atoms , N_of_elements , N_of_cluster_atoms , N_of_Configurations
 character(len=1)  :: TorF , fragment
 logical           :: flag

 real*8           , allocatable :: xyz(:,:)
 integer          , allocatable :: atom_No(:)
 character(len=2) , allocatable :: element(:) , symbol(:)

 OPEN(unit=3,file='poscar.in',status='old')

! start reading the structure characteristics
 read(3,*) System_Characteristics
 read(3,*) N_of_Configurations      ! <== No of configurations 
 read(3,*) unit_cell%atoms          ! <== No of atoms in the ORIGINAL supercell
 
! allocating arrays
 CALL Allocate_UnitCell(unit_cell)

! Acceptor (A) , Donor (D) or Molecule (M)
! fragment atoms must be packed together
 read(3,*) n_fragments
 do i = 1 , n_fragments
    read(3,*) j1 , j2 , fragment
    forall(j=j1:j2) unit_cell%fragment(j) = fragment
 end do

! defining the k_WH parameter for the atom 
! No of k_WH to be used 
 read(3,*) n_kWH                    
 do i = 1 , n_kWH
    read(3,*) j1 , j2 , k_WH
    forall(j=j1:j2) unit_cell%k_WH(j) = k_WH
 end do

! reads the unit cell vectors
 read(3,*) x0 , y0 , z0
 a = dsqrt(x0*x0 + y0*y0 + z0*z0)
 read(3,*) x0 , y0 , z0
 b = dsqrt(x0*x0 + y0*y0 + z0*z0)
 read(3,*) x0 , y0 , z0
 c = dsqrt(x0*x0 + y0*y0 + z0*z0)

! reads the chemical elements
 read(3,*) N_of_elements 
 allocate(element(N_of_elements),atom_No(N_of_elements))
 read(3,*) (element(i),i=1,N_of_elements)

! reads the number of atoms for each species
 read(3,*) (atom_No(i),i=1,N_of_elements)

 N_of_atoms = sum(atom_No(1:N_of_elements))
 if( N_of_atoms /= unit_cell%atoms) pause 'Number of atoms do not match !!!'

 allocate(xyz(N_of_atoms,3),symbol(N_of_atoms))

! reads the atom that sets the interface 
 read(3,*) boundary_atom

! reads the coordinates from vasp-POSCAR
 do i = 1 , N_of_atoms

    read(3,*) xyz(i,1) , xyz(i,2) , xyz(i,3) , TorF , TorF , TorF
 
    do j = 1 , N_of_elements
        if( (sum(atom_No(:j-1))+1 <= i) .AND. (i <= sum(atom_No(:j))) ) symbol(i) = element(j)
    end do

 end do
 close(3)

! rescale coordinates by the unit cell vectors
 xyz(:,1) = xyz(:,1)*a
 xyz(:,2) = xyz(:,2)*b
 xyz(:,3) = xyz(:,3)*c

! separate cluster from molecule atoms
 n = 1
 do i = 1 , N_of_atoms
    flag = xyz(i,3) > xyz(boundary_atom,3) ! <== cluster atoms first
    if (flag) then
        unit_cell%symbol(n)  = symbol(i)
        unit_cell%coord(n,:) = xyz(i,:)
        n = n + 1
    end if
 end do
 N_of_cluster_atoms = n - 1

 do i = 1 , N_of_atoms
    flag = xyz(i,3) <= xyz(boundary_atom,3) ! <== molecule atoms at the end
    if (flag) then
        unit_cell%symbol(n)  = symbol(i)
        unit_cell%coord(n,:) = xyz(i,:)
        n = n + 1
    end if
 end do

 CALL Symbol_2_AtNo(unit_cell)

 deallocate(xyz,symbol,element,atom_No)

 include 'formats.h'

 end subroutine Read_from_Poscar
!
!
!
!
 subroutine Symbol_2_AtNo(a)

 type(structure) :: a

 DO i = 1 , a%atoms

    select case(a%symbol(i))
        case( 'H') 
            a%AtNo(i) = 1 
        case( 'C') 
            a%AtNo(i) = 6 
        case( 'N') 
            a%AtNo(i) = 7 
        case( 'O') 
            a%AtNo(i) = 8 
        case( 'F') 
            a%AtNo(i) = 9 
        case( 'AL','Al') 
            a%AtNo(i) = 13 
        case( 'S','s') 
            a%AtNo(i) = 16 
        case( 'CL','Cl') 
            a%AtNo(i) = 17 
        case( 'TI','Ti') 
            a%AtNo(i) = 22 
        case( 'MN','Mn') 
            a%AtNo(i) = 25 
        case default
            print*, ' >> unkown atom found ; execution terminated <<' 
            stop
    end select

 END DO

 end subroutine Symbol_2_AtNo
!
!
!
!
 subroutine seed_Yaehmop(system)

 type(structure)   , intent(in) :: system

!----------------------------------------------
!     generate input file for   YAEHMOP
!----------------------------------------------

      OPEN(unit=9,file='yaehmop',status='unknown')

      write(9,*) System_Characteristics
      write(9,*) " "
      write(9,*) "Molecular"      
      write(9,*) " "
      write(9,*) "Geometry"
      write(9,*) " "
      write(9,*) system%atoms
      write(9,*) " "

      DO j = 1 , system%atoms

        write(9,100)j, system%symbol(j) , (system%coord(j,i),i=1,3)

      END DO

      write(9,*) " "
      write(9,*) "just geom"

      CLOSE(9)

!------------------------------------------------------------

 include 'formats.h'

 end subroutine seed_Yaehmop
!
!
!
 subroutine seed_VASP(system)

 use EHT_parameters

 type(structure)   , intent(in) :: system

 real*8 :: x_min , y_min , z_min
 real*8 , dimension(system%atoms) :: x , y , z
 integer                          :: i , j
 integer , dimension(size(atom))  :: list
 integer , dimension(30)          :: short_list
 character(len=1) , parameter     :: TorF = 'T' 

!----------------------------------------------
!     generate input file for  VASP 
!----------------------------------------------

 x_min = minval(system%coord(:,1))
 y_min = minval(system%coord(:,2))
 z_min = minval(system%coord(:,3))

! rescale coordinates by the unit cell vectors'
 x = (system%coord(:,1)-x_min)/T_x
 y = (system%coord(:,2)-y_min)/T_y
 z = (system%coord(:,3)-z_min)/T_z

 do i = 1 , system%atoms
    list(system%AtNo(i)) = list(system%AtNo(i)) + 1
 end do

 j = 0
 do i = size(list) , 1 , -1
    if( list(i) /= 0) then
        j = j + 1
        short_list(j) = list(i)
    end if
 end do

 OPEN(unit=9,file='poscar.out',status='unknown')

 write(9,*) System_Characteristics
 write(9,*) '1.0'
 write(9,'(F12.5,F12.5,F12.5)') T_x , 0.d0 , 0.d0
 write(9,'(F12.5,F12.5,F12.5)') 0.d0 , T_y , 0.d0 
 write(9,'(F12.5,F12.5,F12.5)') 0.d0 , 0.d0, T_z 
 write(9,*) short_list(1:j)
 write(9,*) 'Selective dynamics'
 write(9,*) 'Direct'

 do i = size(list) , 1 , -1
    if( list(i) /= 0 ) then
        do j = 1 , system%atoms

            if( system%AtNo(j) == i ) write(9,101) x(j), y(j), z(j) , TorF , TorF , TorF

        end do
    end if
 end do

 CLOSE(9)
     
!------------------------------------------------------------

 include 'formats.h'

 end subroutine seed_VASP
!
!
!
end module Babel_m



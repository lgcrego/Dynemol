module solvent_routines

use EDIT_routines , only : atomic , molecular , universe

contains 
!
!
!
!===============================
subroutine Read_Molecule(dye_mol)
!===============================
implicit none
type(molecular) , intent(inout) :: dye_mol

!   local variables
integer :: i , j , ioerr 

! read solvent molecule
OPEN(unit=3,file='molecule.dat',status='old',iostat=ioerr,err=10)

read(3,*) dye_mol%N_of_atoms    
allocate( dye_mol%atom(dye_mol%N_of_atoms) )

read(3,*) dye_mol%Dye_Characteristics  
Print*  , dye_mol%Dye_Characteristics

do i = 1 , dye_mol%N_of_atoms
    read(3,*,iostat=ioerr) dye_mol%atom(i)%symbol , (dye_mol%atom(i)%xyz(j),j=1,3) 
end do

close(3)
Print* , dye_mol%atom%symbol

! find the geo_center of dye_molecule
forall(i=1:3) 
    dye_mol%atom%xyz(i) = dye_mol%atom%xyz(i) - sum( dye_mol%atom%xyz(i) ) / dye_mol%N_of_atoms
end forall

! finding the radius size of the dye-molecule
dye_mol%radius = maxval( dsqrt(dye_mol%atom%xyz(1)**2 + dye_mol%atom%xyz(2)**2 + dye_mol%atom%xyz(3)**2) )


10 if( ioerr > 0 )then
    stop 'molecule.dat file not found; terminating execution'
end if

end subroutine Read_Solvent
!
!
!
!=================================================
subroutine Place_Solvent_Molecules(system,dye_mol)
!=================================================
implicit none
type(universe)  , intent(inout) :: system
type(molecular) , intent(inout) :: dye_mol

! local variables
integer        :: i , j , k , counter , indx , N_x , N_y , Old_No_of_atoms , New_No_of_atoms
real*8         :: x0 , y0 , delta_x , delta_y , geo_center_x , geo_center_y
type(universe) :: temp
!
!
! input data
write(*,'(/a,f8.5,a,f8.5)') '> area of the surface : ',system%box(1),' x ',system%box(2)
write(*,'(/a,f8.6)') '> radius of solvent molecule : ',dye_mol%radius
write(*,'(/a)') '> distribution of solvent Molecules in the solvent layer (N_x,N_y) ?'
write(*,'(/a)', advance='no') 'N_x = '
read (*,'(i2)') N_x 
write(*,'(/a)', advance='no') 'N_y = '
read (*,'(i2)') N_y 

! prepare temporary system
allocate( temp%atom( system%N_of_atoms ) )
temp = system
Old_No_of_atoms = system%N_of_atoms 
New_No_of_atoms = Old_No_of_atoms + (N_x*N_y)*dye_mol%N_of_atoms

! redefine system 
system%N_of_atoms             = New_No_of_atoms
system%N_of_Solvent_atoms     = N_x*N_y * dye_mol%N_of_atoms
system%N_of_Solvent_molecules = N_x*N_y

deallocate( system%atom )
allocate( system%atom(New_No_of_atoms) )
system%atom( 1:Old_No_of_atoms ) = temp%atom( 1:Old_No_of_atoms )
deallocate( temp%atom )
allocate( system%solvent(system%N_of_Solvent_molecules) )

! introduce solvent molecules on the first solvent layer
x0 = minval(system%atom%xyz(1))
y0 = minval(system%atom%xyz(2))
delta_x = system%box(1)/N_x
delta_y = system%box(2)/N_y

counter = 0
do i = 0 , N_x-1
    do j = 0 , N_y-1

        geo_center_x = x0 + (i+0.5)*delta_x
        geo_center_y = y0 + (j+0.5)*delta_y

        do k = 1 , dye_mol%N_of_atoms

            indx = k + (counter * dye_mol%N_of_atoms)

            system%atom(Old_No_of_atoms+indx)%xyz(1) = geo_center_x + dye_mol%atom(k)%xyz(1) 
            system%atom(Old_No_of_atoms+indx)%xyz(2) = geo_center_y + dye_mol%atom(k)%xyz(2)  
            system%atom(Old_No_of_atoms+indx)%xyz(3) = (system%Surface_Boundary+3) + dye_mol%atom(k)%xyz(3) 

            system%atom(Old_No_of_atoms+indx)%symbol = dye_mol%atom(k)%symbol  

        end do

        counter = counter + 1

        system%solvent(counter)%xyz(1) = geo_center_x
        system%solvent(counter)%xyz(2) = geo_center_y
        system%solvent(counter)%xyz(3) = system%Surface_Boundary + 3

    end do
end do

forall(i = Old_No_of_atoms:New_No_of_atoms)
    system%atom(i)%fragment = 'S'
    system%atom(i)%TorF     = 'T'
end forall

end subroutine Place_Solvent_Molecules
!
end module solvent_routines

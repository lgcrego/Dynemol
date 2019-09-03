module solvent_routines

    use constants_m
    use Read_Parms      , only : atomic , atomic_mass , Symbol_2_AtNo
    use types_m         , only : molecular , universe
    use RW_routines     , only : Read_from_XYZ
    use GMX_routines    , only : Dump_pdb
    use diagnosis_m

    public  ::  Include_Solvent

private

    type mol_PTR
        integer :: Nx
        integer :: Ny
        integer :: layer
    end type mol_PTR

contains 
!
!
!
!==========================================
subroutine Include_Solvent(structure)
!==========================================
implicit none
type(universe) , intent(inout) :: structure

! local variables ...
character(len=1) :: choice
type(molecular)  :: sol_mol


CALL system( "clear" )

write(*,'(/a)') ' (w) = wet the surface '
write(*,'(/a)') ' (c) = read solvent box '
write(*,'(/a)') ' (b) = build solvent box - I '
write(*,'(/a)') ' (x) = build solvent box - Standard '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice

select case( choice )

    case( 'w' ) 
        CALL Read_Solvent_Molecule   ( sol_mol )
        CALL Place_solvent_Molecules ( structure,sol_mol )
        
    case( 'x' ) 
        CALL Solid_Liquid_Interface ( structure )
        
    case( 'c' )
        CALL Read_Solvent_Cube( structure )
        CALL dump_pdb( structure )
        write(*,'(/a)') '>>>  Saving seed.pdb  <<<'

    case( 'b' )
        CALL Build_Solvent_Cube( structure )
        CALL dump_pdb( structure )
        write(*,'(/a)') '>>>  Saving seed.pdb  <<<'

end select

end subroutine Include_Solvent
!
!
!
!====================================
subroutine Build_Solvent_Cube(system)
!====================================
implicit none
type(universe)  , intent(out)   :: system

! local variables
integer                         :: i , j , k , n , n_grid , indx , nresid , mol_indx
real*8                          :: v_BOX , v_cell , a_cell , density , solvent_mass 
real*8                          :: R_vector(3)
real*8          , parameter     :: cm_2_Angs_factor = 1.d24
character(len=5)                :: string
type(molecular)                 :: sol_mol

!-----------------------------------------------------------------
! provide solvent molecule file ...
print*, "be prepared to provide file (solvent.dat) with solvent molecule info: "
print*, " "
print*, "N_of_atoms"
print*, "residue  &  Solvent Characteristics"
print*, "symbol xyz MMSymbol charge (free format)"
print*, " .      .     .        . "
print*, " .      .     .        . "
print*, " .      .     .        . "
!-----------------------------------------------------------------
 
! information input ...
write(*,'(/a)',advance='no' ) '> Number of solvent molecules (N^3 = 8 , 27 , 64 , 125 , 216 , 343 , 512 , 729 , 1000) =  ?  '
read (*,'(i4)'              ) system%N_of_Solvent_Molecules 
write(*,'(/a)', advance='no') 'Density of the solvent (g/cm^3) = '
read (*,'(f10.5)'           ) density 


CALL Read_Solvent_Molecule( sol_mol )

! define cube parameters ...
solvent_mass = sum( sol_mol % atom % mass ) * u_mass * system % N_of_Solvent_Molecules

V_BOX = (solvent_mass / density) * cm_2_Angs_factor
V_cell = V_BOX /  system%N_of_Solvent_Molecules
a_cell = V_cell**third

system%N_of_atoms = system%N_of_Solvent_Molecules * sol_mol%N_of_atoms
system%box = V_BOX**third

! build solvent cube ...
allocate( system%atom(system%N_of_atoms) )

n_grid = system%N_of_Solvent_Molecules**third
 
mol_indx = 0
do k = 0 , n_grid-1     ;   R_vector(3) = k * a_cell
do j = 0 , n_grid-1     ;   R_vector(2) = j * a_cell
do i = 0 , n_grid-1     ;   R_vector(1) = i * a_cell

            do n = 1 , sol_mol%N_of_atoms

                system % atom(mol_indx+n) % xyz      =  (R_vector + a_cell / 2.d0) + sol_mol % atom(n) % xyz

                system % atom(mol_indx+n) % symbol   =  sol_mol % atom(n) % symbol  
                system % atom(mol_indx+n) % MMSymbol =  sol_mol % atom(n) % MMSymbol  
                system % atom(mol_indx+n) % charge   =  sol_mol % atom(n) % charge  

            end do

            mol_indx = mol_indx + sol_mol%N_of_atoms
end do
end do
end do

system%atom%resid    = sol_mol%resid
system%atom%fragment = 'S'

! get the Atomic_Number ...
CALL Symbol_2_AtNo(system%atom)

! get the Atomic_Masses ...
system%atom%mass = Atomic_Mass(system%atom%AtNo)

! include solvent molecule information in universe ...
allocate( system%solvent%atom(sol_mol%N_of_atoms) )
system%solvent = sol_mol
system%solvent%N_of_atoms = sol_mol%N_of_atoms

! MMSymbols for the box ...
indx    = 0
nresid  = 0
do  
    nresid = nresid + 1
    do j = 1 , sol_mol%N_of_atoms
        system%atom(indx+j)%MMSymbol = sol_mol%atom(j)%MMsymbol
        system%atom(indx+j)%nresid   = nresid 
    end do
    indx = indx + sol_mol%N_of_atoms

    if( indx >= system%N_of_atoms ) EXIT
end do

! update system characteristics ...
write(string,'(i5.5)') system % N_of_Solvent_molecules
system%Surface_Characteristics = trim(sol_mol%Solvent_Characteristics)//'-'//string//'solvent molecules'

end subroutine Build_Solvent_Cube
!
!
!
!
!=====================================
subroutine Read_Solvent_Cube( system )
!=====================================
implicit none
type(universe)  , intent(inout) :: system

! local variables ...
integer                         :: i , j , indx , nresid
real*8                          :: distance
real*8                          :: xyz_min(3) , GC(3) , box_center(3)
character(len=1)                :: shrink
logical         , allocatable   :: mask(:)
type(atomic)    , allocatable   :: temp(:)
type(molecular)                 :: sol_mol


CALL Read_from_XYZ( system )

CALL Read_Solvent_Molecule( sol_mol )

system%atom%resid    = sol_mol%resid
system%atom%fragment = 'S'
system%solvent%N_of_atoms = sol_mol%N_of_atoms

! determine box side ...
forall( i=1:3 ) system%box(i) = maxval(system%atom%xyz(i)) - minval(system%atom%xyz(i))

box_center = system%box / two

Print 100, system%box(1) , system%box(2) , system%box(3) 

! place the solvent in a first quadrant
forall( i=1:3 ) xyz_min(i) = minval( system%atom%xyz(i) )
forall( i=1:3 ) system%atom%xyz(i) = system%atom%xyz(i) - xyz_min(i)

! change the side of the box ...
write(*,'(/a)',advance='no') '> Change the side of the box (y/n) ? '
read (*,'(a)') shrink
if( shrink == 'y' ) then

    allocate( mask(system%N_of_atoms) )
    mask = .TRUE.

    write(*,'(/a)'             ) ' new size of the box : '
    write(*,'(/a)',advance='no') ' L_x = '
    read (*,'(f8.3)'           ) system%box(1)
    write(*,'(/a)',advance='no') ' L_y = '
    read (*,'(f8.3)'           ) system%box(2)
    write(*,'(/a)',advance='no') ' L_z = '
    read (*,'(f8.3)'           ) system%box(3)

    indx    = 0
    do  
        forall( i=1:3 ) GC(i) = sum( system % atom (indx+1 : indx+sol_mol%N_of_atoms) % xyz(i) ) / sol_mol%N_of_atoms 
        
        do i = 1 , 3
            distance = GC(i) - box_center(i)
            if( dabs(distance) > (system%box(i)/two - 2.d0*sol_mol%radius) ) then
                mask( indx+1 : indx+sol_mol%N_of_atoms ) = .FALSE. 
!                system % atom( indx+1 : indx+sol_mol%N_of_atoms ) % symbol = 'H'
                EXIT
            end if
        end do

        indx = indx + sol_mol%N_of_atoms 

        if( indx >= system%N_of_atoms ) EXIT
    end do

    ! redine system ...
    system%N_of_atoms = count( mask )
    allocate( temp(system%N_of_atoms) )
    temp = pack( system%atom , mask , temp )
    CALL move_alloc( from=temp , to=system%atom )

end if

! fix the MMSymbols for the box ...
indx    = 0
nresid  = 0
do  
    nresid = nresid + 1
    do j = 1 , sol_mol%N_of_atoms
        system%atom(indx+j)%MMSymbol = sol_mol%atom(j)%MMsymbol
        system%atom(indx+j)%nresid   = nresid 
    end do
    indx = indx + sol_mol%N_of_atoms

    if( indx >= system%N_of_atoms ) EXIT
end do

system%N_of_Solvent_Molecules = nresid

write(*,'(/a36,i5)') '>>>  Number of solvent molecules = ' , nresid

100 FORMAT(1x,'box = (',f7.3,','f7.3,','f7.3,')')

end subroutine Read_Solvent_Cube
!
!
!
!========================================
subroutine Read_Solvent_Molecule(sol_mol)
!========================================
implicit none
type(molecular) , intent(inout) :: sol_mol

! local variables ...
integer :: i , j , ioerr 

! read solvent molecule and MM paramenters ...
!-----------------------------------------------------------------
OPEN(unit=3,file='solvent.dat',status='old',iostat=ioerr,err=10)

read(3,*) sol_mol%N_of_atoms    
allocate( sol_mol%atom(sol_mol%N_of_atoms) )

read(3,'(a3,1x,a60)') sol_mol%resid , sol_mol%Solvent_Characteristics  
write(*,'(/,1x,a/)'  ) sol_mol%Solvent_Characteristics

do i = 1 , sol_mol%N_of_atoms
    read(3,*,iostat=ioerr) sol_mol%atom(i)%symbol , (sol_mol%atom(i)%xyz(j),j=1,3) , sol_mol%atom(i)%MMSymbol , sol_mol%atom(i)%charge
end do

close(3)
!-----------------------------------------------------------------
Print*, sol_mol%atom%MMSymbol
write(*, '(/)') 

! translate solvent coordinates to the geometric center ...
forall(i=1:3) 
    sol_mol%atom%xyz(i) = sol_mol%atom%xyz(i) - sum( sol_mol%atom%xyz(i) ) / sol_mol%N_of_atoms
end forall

! finding the radius size of the solvent ...
sol_mol%radius = maxval( dsqrt(sol_mol%atom%xyz(1)**2 + sol_mol%atom%xyz(2)**2 + sol_mol%atom%xyz(3)**2) )

! get the Atomic_Number ...
CALL Symbol_2_AtNo(sol_mol%atom)

! get the Atomic_Masses ...
sol_mol%atom%mass = Atomic_Mass(sol_mol%atom%AtNo)

10 if( ioerr > 0 )then
    stop 'solvent.dat file not found; terminating execution'
end if

end subroutine Read_Solvent_Molecule
!
!
!
!=================================================
subroutine Place_Solvent_Molecules(system,sol_mol)
!=================================================
implicit none
type(universe)  , intent(inout) :: system
type(molecular) , intent(inout) :: sol_mol

! local variables
integer                         :: i , j , k , n
integer                         :: molecule_counter , indx , Nx , Ny , layer , N_layers , Old_No_of_atoms , New_No_of_atoms , N_erase , n_void
real*8                          :: x0 , y0 , delta_x , delta_y , geo_center_x , geo_center_y
type(mol_PTR)   , allocatable   :: void(:)
type(universe)                  :: temp
character(len=1)                :: Y_or_N
character(len=3)                :: string
 
! information input ...
write(*,'(/a,f8.5,a,f8.5)') '> area of the surface : ',system%box(1),' x ',system%box(2)
write(*,'(/a,f8.4)') '> radius of solvent molecule : ',sol_mol%radius
write(*,'(/a)', advance='no') '> number of solvent layers : '
read (*,'(i2)') N_layers 
write(*,'(/a)') '> distribution of solvent Molecules in the solvent layer (Nx,Ny) ?'
write(*,'(/a)', advance='no') 'Nx = '
read (*,'(i2)') Nx 
write(*,'(/a)', advance='no') 'Ny = '
read (*,'(i2)') Ny 
! eliminate solvent molecules ...
write(*,'(/a)') '> Will eliminate any solvent molecules (Y/N) ?'
read (*,'(a)' ) Y_or_N
if( (Y_or_N=='Y') .OR. (y_or_N=='y') ) then
    write(*,'(/a)', advance='no') '> How many ? '
    read (*,'(i2)') N_erase
    allocate( void(n_erase) )
    do n = 1 , N_erase
        write(*,'(/a)') '> solvent molecule (Ix,Iy) ?'
        write(*,'(/a)', advance='no') 'Ix = ' 
        read (*,'(i2)') void(n)%Nx                  ;   void(n)%Nx = void(n)%Nx - 1
        write(*,'(/a)', advance='no') 'Iy = '
        read (*,'(i2)') void(n)%Ny                  ;   void(n)%Ny = void(n)%Ny - 1
        write(*,'(/a)', advance='no') 'Layer = '
        read (*,'(i2)') void(n)%layer               ;   void(n)%layer = void(n)%layer - 1
    end do
else
    N_erase = 0
    allocate( void(1) )
    void(1)%Nx    = - 1
    void(1)%Ny    = - 1
    void(1)%layer = - 1
end if

! prepare temporary system ...
New_No_of_atoms = system%N_of_atoms + ( Nx*Ny*N_layers - N_erase ) * sol_mol%N_of_atoms 
allocate( temp%atom( New_No_of_atoms ) , source=system%atom )

Old_No_of_atoms = system%N_of_atoms 

! redefine system ... 
system % N_of_atoms             =  New_No_of_atoms
system % N_of_Solvent_molecules =  Nx*Ny*N_layers - N_erase 
system % N_of_Solvent_atoms     =  system % N_of_Solvent_molecules * sol_mol%N_of_atoms

! find the interface position + add a spacer of 2 Angs
system%Surface_Boundary = maxval(system%atom%xyz(3) , system%atom%symbol=='Ti') + 2.0

! introduce solvent molecules ...
x0 = minval(system%atom%xyz(1))
y0 = minval(system%atom%xyz(2))
delta_x = system%box(1)/Nx
delta_y = system%box(2)/Ny

n_void           = 1
molecule_counter = 0

do layer = 0 , N_layers-1
    do i = 0 , Nx-1
        do j = 0 , Ny-1

            geo_center_x = x0 + (i+0.5)*delta_x
            geo_center_y = y0 + (j+0.5)*delta_y

            if( (i==void(n_void)%Nx) .AND. (j==void(n_void)%Ny) .AND. (layer==void(n_void)%layer) ) then

                n_void = n_void + 1

            else

                do k = 1 , sol_mol%N_of_atoms

                    indx = k + (molecule_counter * sol_mol%N_of_atoms)

                    temp % atom(Old_No_of_atoms+indx) % xyz(1)   =  geo_center_x + sol_mol % atom(k) % xyz(1) 
                    temp % atom(Old_No_of_atoms+indx) % xyz(2)   =  geo_center_y + sol_mol % atom(k) % xyz(2)  
                    temp % atom(Old_No_of_atoms+indx) % xyz(3)   =  (system%Surface_Boundary+1.0) + layer*(2.5*sol_mol%radius) + sol_mol%atom(k)%xyz(3) 

                    temp % atom(Old_No_of_atoms+indx) % symbol   =  sol_mol % atom(k) % symbol  
                    temp % atom(Old_No_of_atoms+indx) % MMSymbol =  sol_mol % atom(k) % MMSymbol  
                    temp % atom(Old_No_of_atoms+indx) % charge   =  sol_mol % atom(k) % charge  

                    temp % atom(Old_No_of_atoms+indx) % nrcg     =  Old_No_of_atoms + indx - k + 1

                end do

                molecule_counter = molecule_counter + 1

            end if

        end do
    end do
end do

deallocate(void)

forall( i = Old_No_of_atoms+1:New_No_of_atoms )
    temp % atom(i) % fragment = 'S'
    temp % atom(i) % TorF     = 'T'
end forall

CALL move_alloc(from=temp%atom,to=system%atom)

! get the Atomic_Number ...
CALL Symbol_2_AtNo(system%atom)

! get the Atomic_Masses ...
system%atom%mass = Atomic_Mass(system%atom%AtNo)

! include solvent molecule information in universe ...
allocate( system%solvent%atom(sol_mol%N_of_atoms) )
system%solvent = sol_mol

! update system characteristics ...
write(string,'(i3.3)') system % N_of_Solvent_molecules
system%Surface_Characteristics = trim(system%Surface_Characteristics)//'-'//string//'solvent molecules'

end subroutine Place_Solvent_Molecules
!
!
!
!==========================================
subroutine Solid_Liquid_interface( system )
!==========================================
implicit none
type(universe)  , intent(inout)   :: system

! local variables
integer                         :: Nx , Ny , N_layers , New_No_of_atoms , Old_No_of_atoms
integer                         :: i , j , k , n , counter , layer , nresid
real                            :: random
real*8                          :: x0 , y0 , delta_x , delta_y , delta_z , geo_center_x , geo_center_y , geo_center_z
real*8                          :: cell_height , cell_area , density , solvent_mass 
real*8                          :: sol_mol_mass , R_wigner_mol
real*8          , parameter     :: cm_2_Angs_factor = 1.d8
character(1)                    :: YorN
character(3)                    :: solv_resid , surf_resid , string
type(molecular)                 :: sol_mol
type(atomic)    , allocatable   :: temp(:)

! input parameters ...
write(*,'(/a)',advance='no') 'Place Solvent on top of residue : '
read (*,'(a)'              ) surf_resid
write(*,'(/a)',advance='no') 'Type of Solvent residue : '
read (*,'(a)'              ) solv_resid
write(*,'(/a)',advance='no') 'Density of the solvent (g/cm^3) = '
read (*,'(f10.5)'          ) density 

CALL Read_Solvent_Molecule( sol_mol )

! Wigner radius of the solvent molecule ...
sol_mol_mass = sum( sol_mol % atom % mass ) * u_mass 
R_wigner_mol = ( sol_mol_mass/(four*PI*third*density) )**third * cm_2_Angs_factor 

write(*,fmt=100) 'Wigner radius of the molecule = ', R_wigner_mol ,'(',sol_mol%radius,')'
write(*,fmt=101) 'Area of the interface = ', system % box(1) ,'x', system % box(2)

! information input ...
write(*,'(/a)') '> distribution of solvent Molecules in the solvent layer (Nx,Ny) ?'
write(*,'(/a)', advance='no') 'Nx = '
read (*,'(i2)') Nx 
write(*,'(/a)', advance='no') 'Ny = '
read (*,'(i2)') Ny 
write(*,'(/a)',advance='no' ) 'Number of layers = '
read (*,'(i4)') N_layers

! dimensions of the solvent box ...
solvent_mass = sol_mol_mass * Nx * Ny * N_layers
cell_area    = system % box(1) * system % box(2)
cell_height  = solvent_mass / (cell_area * density) * cm_2_Angs_factor**3

! find the interface position ...
system % atom % xyz(3) = system % atom % xyz(3) - minval( system % atom % xyz(3) )
system % Surface_Boundary = maxval( system % atom % xyz(3) , system % atom % resid == surf_resid ) 
system % box(3) = system % Surface_Boundary + cell_height

! prepare temporary system ...
New_No_of_atoms = system%N_of_atoms + ( Nx*Ny*N_layers ) * sol_mol%N_of_atoms 
Old_No_of_atoms = system%N_of_atoms 
allocate( temp( New_No_of_atoms ) )
temp(1:Old_No_of_atoms) = system%atom

! introduce solvent molecules ...
x0 = minval(system % atom % xyz(1))
y0 = minval(system % atom % xyz(2))
delta_x = system%box(1) / Nx
delta_y = system%box(2) / Ny
delta_z = cell_height   / N_layers

counter = Old_No_of_atoms
!nresid  = system % atom(Old_No_of_atoms) % nresid + 2
nresid  = maxval( system % atom % nresid ) + 1

CALL random_seed

do layer = 0 , N_layers-1
    do i = 0 , Nx-1
        do j = 0 , Ny-1

            CALL random_number( random )

            geo_center_x = x0 + (i+0.5)*delta_x
            geo_center_y = y0 + (j+0.5)*delta_y
            geo_center_z = system % Surface_Boundary + (layer+0.5)*delta_z + (two*random-1.0)*delta_z/four + (-1)**(i+j)*delta_z/five 

            do k = 1 , sol_mol%N_of_atoms

                temp (counter+k) % xyz(1)   =  geo_center_x + sol_mol % atom(k) % xyz(1) 
                temp (counter+k) % xyz(2)   =  geo_center_y + sol_mol % atom(k) % xyz(2)  
                temp (counter+k) % xyz(3)   =  geo_center_z + sol_mol % atom(k) % xyz(3) 

                temp (counter+k) % symbol   =  sol_mol % atom(k) % symbol  
                temp (counter+k) % MMSymbol =  sol_mol % atom(k) % MMSymbol  
                temp (counter+k) % charge   =  sol_mol % atom(k) % charge  

                temp (counter+k) % nresid   =  nresid

            end do

            counter = counter + sol_mol%N_of_atoms
            nresid  = nresid   + 1

        end do
    end do
end do

forall( i = Old_No_of_atoms+1:New_No_of_atoms )
    temp (i) % fragment = 'S'
    temp (i) % TorF     = 'T'
    temp (i) % resid    = solv_resid
end forall

! redefine system ... 
CALL move_alloc(from=temp,to=system%atom)

system % N_of_atoms             =  New_No_of_atoms
system % N_of_Solvent_molecules =  Nx*Ny*N_layers 
system % N_of_Solvent_atoms     =  system % N_of_Solvent_molecules * sol_mol%N_of_atoms

! get the Atomic_Number ...
CALL Symbol_2_AtNo(system%atom)

! get the Atomic_Masses ...
system%atom%mass = Atomic_Mass(system%atom%AtNo)

! include solvent molecule information in universe ...
allocate( system%solvent%atom(sol_mol%N_of_atoms) )
system%solvent = sol_mol

! mark S molecules to be deleted ? ...
write(*,'(/a)',advance='no') '>>> Solvent molecules to be deleted ? (y/n) '
read (*,'(a)'              ) YorN
if( YorN == "y" ) CALL Mark_S_molecules( system )

! update system characteristics ...
write(string,'(i3.3)') system % N_of_Solvent_molecules
system%Surface_Characteristics = trim(system%Surface_Characteristics)//'-'//string//'solvent molecules'

! Formats ...
100 format(a32,f7.4,a3,f7.4,a3,/)
101 format(a24,f10.4,a2,f8.4,/) 

end subroutine Solid_Liquid_interface
!
!
!
!====================================
subroutine Mark_S_molecules( system )
!====================================
implicit none
type(universe)  , intent(inout)   :: system

! local variables
integer                         :: i , nresid , begin_nresid , end_nresid 
real*8                          :: CG(3) , distance_sqr , distance

begin_nresid =  minval( system%atom%nresid , system%atom%fragment == "S" )
end_nresid   =  maxval( system%atom%nresid , system%atom%fragment == "S" )

do nresid = begin_nresid , end_nresid

    forall( i=1:3 ) CG(i) = sum( system%atom%xyz(i) , system%atom%nresid == nresid ) / system%solvent%N_of_atoms

    distance_sqr = sum( (CG - system%atom(289)%xyz)**2 )

    distance = sqrt( distance_sqr )

    If( distance <= 6.7 ) where( system%atom%nresid == nresid ) system%atom%fragment = "D" 

end do

Write(*,'(/a,i3)') ">>> Molecules marked for deletion = ", count(system%atom%fragment == "D")/system%solvent%N_of_atoms

end subroutine Mark_S_molecules
!
!
!
end module solvent_routines

module EDIT_routines

use types_m
use Read_Parms
use Constants_m
use GMX_routines

interface ReGroup
    module procedure ReGroup_Molecule
    module procedure ReGroup_Surface
end interface

contains
!
!
!
!========================
subroutine Copy( system )
!========================
implicit none 
type(universe)  , intent(inout) :: system

! local variables ...
integer                         :: i , indx , nresid , N_copies , New_N_of_atoms , Old_N_of_atoms
integer         , allocatable   :: copies(:)
type(universe)                  :: temp

N_copies = count( system%atom%copy )    
allocate( copies(N_copies) )

indx = 1
do i = 1 , system%N_of_Atoms
    If( system%atom(i)%copy ) then
        copies(indx) = i
        indx = indx + 1
    end if
end do
        
! prepare temp system ...
Old_N_of_atoms  = system%N_of_atoms 
New_N_of_atoms  = system%N_of_atoms + N_copies
allocate( temp%atom(New_N_of_atoms) )
temp%atom(1:Old_N_of_atoms) = system%atom

forall( indx=1:N_copies ) temp%atom( indx + Old_N_of_atoms ) = system%atom( copies(indx) )

! finish copy ...
CALL move_alloc( from=temp%atom , to=system%atom )

print*, '>>> copy done <<<'

! translate copy ...
copies = [( i , i = Old_N_of_atoms+1 , New_N_of_atoms )]
CALL Translation( system , copies )

! input information ...
write(*,'(2/a)',advance='no') "residue # of the reference (use zero (0) to keep nresid invariant) :  "
read(*,*) nresid

! finalize copy ...
system % atom(copies) % fragment = "F"
system % atom(copies) % nresid   = system % atom(copies) % nresid + nresid
system % N_of_atoms              = New_N_of_atoms

deallocate( copies )

end subroutine Copy
!
!
!
!========================================
subroutine Translation( system , copies )
!========================================
implicit none 
type(universe)              , intent(inout) :: system
integer        , optional   , intent(in)    :: copies(:)

! local variables ...
real*8       :: T_vector(3) ,T_versor(3) , distance
integer      :: i , option , at2 , at1

! define translation vector
write(*,'(a)',advance='no') '> Use cartesian axis(1) or ad-hoc vector(2)? : '
read(*,*) option

select case (option)

   case(1)
        ! use cartesian vectors 
        write(*,'(/a)') '> enter translation vector (T_x,T_y,T_z) as a Real number:'
        write(*,'(a)',advance='no') 'T_x = '
        read (*,'(f8.4)') T_vector(1)
        write(*,'(a)',advance='no') 'T_y = '
        read (*,'(f8.4)') T_vector(2)
        write(*,'(a)',advance='no') 'T_z = '
        read (*,'(f8.4)') T_vector(3)

   case(2) 
        ! define translation vector 
        write(*,'(/a)') '> define vector: at1 ======> at2'
        write(*,'(a)',advance='no') 'index of atom 1 = '
        read(*,*) at1
        write(*,'(a)',advance='no') 'index of atom 2 = '
        read(*,*) at2
        write(*,'(/a)') '> translation distance (Real number):'
        read (*,'(f8.4)') distance

        T_versor = (system% atom(at2)% xyz - system% atom(at1)% xyz) / sqrt(sum( (system% atom(at2)% xyz-system% atom(at1)% xyz)**2) )
        T_vector = distance * T_versor

end select

if( present(copies) ) then
    forall( i=1:size(copies) ) system%atom(copies(i))%xyz = system%atom(copies(i))%xyz + T_vector
else
    forall( i=1:system%N_of_atoms , system%atom(i)%translate ) system%atom(i)%xyz = system%atom(i)%xyz + T_vector
end if

print*, '>>> translation done <<<'

end subroutine translation
!
!
!
!============================
subroutine Rotation( system )
!============================
implicit none 
type(universe) , intent(inout) :: system

!	local variables
type(universe)  :: temp
real*8          :: R_x(3,3) , R_y(3,3) , R_z(3,3) , R_v(3,3) , pivot(3) , v(3) , angle
integer         :: i , j , N_of_atoms , pivot_atom, option , at1, at2, at3
character(1)    :: axis 

N_of_atoms = size(system%atom)

! define rotation axis
write(*,'(a)',advance='no') '> Rotation around: (1)-cartesian axis , (2)-ad-hoc vector ,  (3)-normal vector?   '
read(*,*) option

select case (option)

   case(1)
        ! use cartesian vectors 
        ! define rotation ...
        write(*,'(a)',advance='no') '> enter angle for clockwise rotation (degrees)  = '
        read (*,'(f8.3)') angle
        write(*,'(\a)',advance='no') '> enter axis  = '
        read (*,*) axis
        write(*,'(\a)',advance='no') '> enter pivot atom  = '
        read (*,*) pivot_atom

   case(2) 
        ! define translation vector 
        write(*,'(/a)') '> define vector: at1 ======> at2'
        write(*,'(a)',advance='no') 'index of atom 1 = '
        read(*,*) at1
        write(*,'(a)',advance='no') 'index of atom 2 = '
        read(*,*) at2

        ! the versor ...
        v = (system% atom(at2)% xyz - system% atom(at1)% xyz) / sqrt(sum( (system% atom(at2)% xyz-system% atom(at1)% xyz)**2) )

        axis = 'v'
        write(*,'(a)',advance='no') '> enter angle for clockwise rotation (degrees)  = '
        read (*,'(f8.3)') angle
        write(*,'(\a)',advance='no') '> enter pivot atom  = '
        read (*,*) pivot_atom

   case(3) 
        ! define rotation vector 
        write(*,'(/a)') '> define vector perpendicular to the plane: '
        write(*,'(/a)') '            at1    at3                      '
        write(*,'(a)')  '              \    /                        '
        write(*,'(a)')  '               \  /                         '
        write(*,'(a)')  '                \/                          '
        write(*,'(a)')  '                at2                         '
        write(*,'(/a)',advance='no') 'index of atom 1 = '
        read(*,*) at1
        write(*,'(a)',advance='no') 'index of atom 2 = '
        read(*,*) at2
        write(*,'(a)',advance='no') 'index of atom 3 = '
        read(*,*) at3

        ! the versor ...
        v = vector_product(system,at1,at2,at3)

        axis = 'v'
        write(*,'(a)',advance='no') '> enter angle for clockwise rotation (degrees)  = '
        read (*,'(f8.3)') angle
        write(*,'(\a)',advance='no') '> enter pivot atom  = '
        read (*,*) pivot_atom

end select
   

allocate( temp%atom(N_of_atoms) )

temp = system

angle = angle * (PI/180.d0)
pivot = system%atom(pivot_atom)%xyz 

select case (axis)

    case( "x" )
        !------------------------
        R_x      =  0.d0
        R_x(1,1) =  1.d0
        R_x(2,2) =  dcos(angle)
        R_x(2,3) = -dsin(angle)
        R_x(3,2) =  dsin(angle)
        R_x(3,3) =  dcos(angle)
        !------------------------
        forall( i=1:N_of_atoms , j=1:3 , system%atom(i)%rotate ) 
            system%atom(i)%xyz(j) = sum( R_x(j,:) * (temp%atom(i)%xyz(:)-pivot(:)) )  + pivot(j)
        end forall

    case( "y" )
        !------------------------
        R_y      =  0.d0
        R_y(2,2) =  1.d0
        R_y(1,1) =  dcos(angle)
        R_y(1,3) =  dsin(angle)
        R_y(3,1) = -dsin(angle)
        R_y(3,3) =  dcos(angle)
        !------------------------
        forall( i=1:N_of_atoms , j=1:3 , system%atom(i)%rotate ) 
            system%atom(i)%xyz(j) = sum( R_y(j,:) * (temp%atom(i)%xyz(:)-pivot(:)) )  + pivot(j)
        end forall

    case( "z" )
        !------------------------
        R_z      =  0.d0
        R_z(3,3) =  1.d0
        R_z(1,1) =  dcos(angle)
        R_z(1,2) = -dsin(angle)
        R_z(2,1) =  dsin(angle)
        R_z(2,2) =  dcos(angle)
        !------------------------
        forall( i=1:N_of_atoms , j=1:3 , system%atom(i)%rotate ) 
            system%atom(i)%xyz(j) = sum( R_z(j,:) * (temp%atom(i)%xyz(:)-pivot(:)) )  + pivot(j)
        end forall

    case( "v" )
        !------------------------
        R_v(1,1) = v(1)**2 + (v(2)**2+v(3)**2)*cos(angle)
        R_v(1,2) = v(1)*v(2)*(1-cos(angle)) - v(3)*sin(angle)
        R_v(1,3) = v(1)*v(3)*(1-cos(angle)) + v(2)*sin(angle)
        R_v(2,1) = v(1)*v(2)*(1-cos(angle)) + v(3)*sin(angle)
        R_v(2,2) = v(2)**2 + (v(1)**2+v(3)**2)*cos(angle)
        R_v(2,3) = v(2)*v(3)*(1-cos(angle)) - v(1)*sin(angle)
        R_v(3,1) = v(1)*v(3)*(1-cos(angle)) - v(2)*sin(angle)
        R_v(3,2) = v(2)*v(3)*(1-cos(angle)) + v(1)*sin(angle)
        R_v(3,3) = v(3)**2 + (v(1)**2+v(2)**2)*cos(angle)
        !------------------------
        forall( i=1:N_of_atoms , j=1:3 , system%atom(i)%rotate ) 
            system%atom(i)%xyz(j) = sum( R_v(j,:) * (temp%atom(i)%xyz(:)-pivot(:)) ) + pivot(j)
        end forall

end select 

print*, '>>> rotation done <<<'

end subroutine rotation
!
!
!
!============================
subroutine Reflection(system)
!============================
implicit none
type(universe) , intent(inout) :: system

! local variables ...
integer          :: xyz
character(len=1) :: axis

write(*,'(/a)') ' Choose the axis (x,y,z) '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') axis

select case( axis )
    case( 'x' ) 
        xyz = 1
    case( 'y' )
        xyz = 2
    case( 'z' )
        xyz = 3 
end select

system%atom%xyz(xyz) = -system%atom%xyz(xyz)
system%atom%xyz(xyz) =  system%atom%xyz(xyz) - minval(system%atom%xyz(xyz))

end subroutine Reflection
!
!
!
!====================================
subroutine Eliminate_Fragment(system)
!====================================
implicit none
type(universe) , intent(inout) :: system

!local variables
integer          :: New_No_of_atoms
character(len=1) :: choice
type(universe)   :: temp

New_No_of_atoms = count( .NOT. system%atom%delete )
allocate( temp%atom( New_No_of_atoms ) , source=system%atom )

temp%atom = pack( system%atom, .NOT. system%atom%delete )

CALL move_alloc(from=temp%atom,to=system%atom)
system%N_of_atoms = New_No_of_atoms

end subroutine Eliminate_Fragment
!
!
!
!==================================
subroutine ReGroup_Molecule(system)
!==================================
implicit none
type(universe) , intent(inout) :: system

!local variables
integer                 :: nr , i , j , indx1 , indx2 
real*8                  :: delta(3)

do nr = minval(system%atom%nresid) , maxval(system%atom%nresid)

    ! atomic pointers of molecule with nresidue = nr
    indx1 = minval( [(i , i=1,size(system%atom))] , (system%atom%nresid == nr) )
    indx2 = maxval( [(i , i=1,size(system%atom))] , (system%atom%nresid == nr) )

    if( system%atom(indx1)%group == .false. ) cycle

    do i = indx1+1 , indx2

        delta = system%atom(i)%xyz - system%atom(indx1)%xyz

        do j = 1 , 3
             if( abs(delta(j)) > system%box(j)*HALF ) system%atom(i)%xyz(j) = system%atom(i)%xyz(j) - sign( system%box(j) , delta(j) )
        end do 

    end do

end do

end subroutine ReGroup_Molecule
!
!
!
!==============================
subroutine ReGroup_Surface(trj)
!==============================
type(universe)  , allocatable   , intent(inout)   :: trj(:)

! local variables ...
integer              :: i , j
real*8 , allocatable :: pm(:)

!local parameters ...
real*8 , parameter   :: ThreeFourth = three/four
real*8 , parameter   :: TwoThird = two/three

allocate( pm(trj(1)%N_of_atoms) )

! fixing atoms to the original unit-cell ; THIS IS ONLY FOR EXTENDED SURFACES 
do j = 1 , size(trj)
do i = 1 , 3
    where( dabs(trj(j)%atom%xyz(i)-trj(1)%atom%xyz(i)) > TwoThird*trj(1)%box(i) ) 

        pm = sign( 1.d0 , trj(1)%atom%xyz(i) - trj(j)%atom%xyz(i) )

        trj(j)%atom%xyz(i) = trj(j)%atom%xyz(i) + pm * trj(1)%box(i)

    end where
end do
end do

end subroutine ReGroup_Surface
!
!
!
!==================================
subroutine Include_Fragment(system)
!==================================
implicit none
type(universe)  , intent(inout) :: system

!local variables
integer                         :: i , j , indx , nresid
real*8          , parameter     :: border = 2.2d0
real*8                          :: GC(3) , new_GC(3) , translate(3) , radius
real*8          , allocatable   :: distance(:)
logical         , allocatable   :: mask(:)
character(1)                    :: choice , order , place , extension
character(12)                   :: file_name
type(universe)                  :: temp , fragment

! STDIN info ...
write(*,'(/a)') ' Type of residue (use capital letters) ?     '
write(*,'(/a)') ' C = cluster' 
write(*,'(/a)') ' F = fragment' 
write(*,'(/a)') ' S = solvent'
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice
write(*,'(/a)',advance='no') ' Include before (b) or after (a) ?    '
read (*,'(a)') order 
write(*,'(/a)',advance='no') ' Format of fragment input file :    pdb-1     /     gro-2'
read (*,'(i1)') extension  

! read fragment from input.gro ...
select case( extension )
    case("1")
        file_name = "fragment.pdb"
    case("2")
        file_name = "fragment.gro"
    end select
CALL Read_GROMACS( fragment, file_name )

! calculate the geometric center of the frgment ...
forall(i=1:3) GC(i) = sum( fragment%atom%xyz(i) ) / fragment%N_of_atoms

! finding the radius of the fragment ...
radius = maxval( dsqrt( (fragment%atom%xyz(1)-GC(1))**2 + (fragment%atom%xyz(2)-GC(2))**2 + (fragment%atom%xyz(3)-GC(3))**2) )

! re-place the GC ...
write(*,'(/a)',advance='no') ' Place the geometric center (y/n) ?    '
read (*,'(a)') place
if( place == "y" ) then
    write(*,'(/a)'             ) ' position to place the GC : '
    write(*,'(/a)',advance='no') ' GC_x = '
    read (*,'(f8.4)'           ) new_GC(1)
    write(*,'(/a)',advance='no') ' GC_y = '
    read (*,'(f8.4)'           ) new_GC(2)
    write(*,'(/a)',advance='no') ' GC_z = '
    read (*,'(f8.4)'           ) new_GC(3)
else
    new_GC = system%box / 2.d0    
end if

! translate fragment to the new geometric center ...
translate = new_GC - GC
forall(i=1:3) fragment%atom%xyz(i) = fragment%atom%xyz(i) + translate(i)

! clean the space around the fragment ...
allocate( distance  ( system%N_of_atoms ) )
allocate( mask      ( system%N_of_atoms ) )

forall( i=1:system%N_of_atoms ) distance(i) = dsqrt( (system%atom(i)%xyz(1) - new_GC(1))**2 + &
                                                     (system%atom(i)%xyz(2) - new_GC(2))**2 + &  
                                                     (system%atom(i)%xyz(3) - new_GC(3))**2 )  
mask = .TRUE. 
do i = 1 , system%N_of_atoms
    if( distance(i) <= radius + border ) then 
        where( system%atom%nresid == system%atom(i)%nresid) mask = .FALSE.
    end if
end do

allocate( temp%atom( count(mask) ) )
temp%atom = pack( system%atom , mask )
CALL move_alloc( from=temp%atom , to=system%atom )
system%N_of_atoms = count( mask )

! fix nresid for the solvent ...
indx    = 0
nresid  = 0
do  
    nresid = nresid + 1
    do j = 1 , system%solvent%N_of_atoms
        system%atom(indx+j)%nresid = nresid 
    end do
    indx = indx + system%solvent%N_of_atoms 
    if( indx >= system%N_of_atoms ) EXIT
end do

! merge the systems ...
temp%N_of_atoms = system%N_of_atoms + fragment%N_of_atoms 
allocate( temp%atom( temp%N_of_atoms ) )

select case( order )
    case( "b" )
        temp % atom( 1 : fragment%N_of_atoms )                  = fragment % atom(:)
        temp % atom( fragment%N_of_atoms+1 : temp%N_of_atoms )  = system % atom(:) 

        temp % atom( 1 : fragment%N_of_atoms ) % nresid = 1
        forall( i = fragment%N_of_atoms+1 : temp%N_of_atoms ) temp % atom(i) % nresid = temp % atom(i) % nresid + 1

    case( "a" )
        temp % atom( 1 : system%N_of_atoms )                    = system % atom(:)
        temp % atom( system%N_of_atoms+1 : temp%N_of_atoms )    = fragment % atom(:) 

        temp % atom( system%N_of_atoms+1 : temp%N_of_atoms ) % nresid = temp % atom(system%N_of_atoms) % nresid + 1
end select

CALL move_alloc( from=temp%atom , to=system%atom )
system%N_of_atoms = temp%N_of_atoms

write(*,'(/a)',advance='no') '>>> Dump seed.pdb (y/n) ?  '
read (*,'(a)') choice
if( choice == 'y' ) then
    CALL Dump_pdb(system)
    write(*,'(/a)') '>>>  Saving seed.pdb  <<<'
end if

deallocate(fragment%atom , distance , mask)

end subroutine Include_Fragment
!
!
!
!==============================
 subroutine Replicate( system )
!==============================
implicit none
type(universe) , intent(inout) :: system

! local variables ...
type(atomic)           , allocatable    :: temp(:)
type(integer_interval)                  :: n_x , n_y , n_z
integer                                 :: New_N_of_atoms , i , j , k , n , counter , S_counter , F_counter , Replication_Factor 
integer                                 :: max_nresid = 0 , P_counter 
character(len=3)                        :: string
logical                , save           :: done = .false.

! reading parameters ...
If( .NOT. done ) then
    write(*,'(a)') '> enter replication configuration parameters (integer)'
    write(*,'(a)',advance='no') 'n_x%inicio (<= 0) = '
    read (*,'(I200)') n_x%inicio 
    write(*,'(a)',advance='no') 'n_x%fim    (>= 0) = '
    read (*,'(I200)') n_x%fim
    write(*,'(a)',advance='no') 'n_y%inicio (<= 0) = '
    read (*,'(I200)') n_y%inicio
    write(*,'(a)',advance='no') 'n_y%fim    (>= 0) = '
    read (*,'(I200)') n_y%fim
    write(*,'(a)',advance='no') 'n_z%inicio (<= 0) = '
    read (*,'(I200)') n_z%inicio 
    write(*,'(a)',advance='no') 'n_z%fim =  (>= 0) = '
    read (*,'(I200)') n_z%fim

    done = .true.
end If

! creating the temp array ...
Replication_Factor  = ( n_x%fim - n_x%inicio + 1 ) * ( n_y%fim - n_y%inicio + 1 ) * ( n_z%fim - n_z%inicio + 1 ) 
New_N_of_atoms      = system%N_of_atoms * Replication_factor

allocate( temp( New_N_of_atoms ) )
temp(1:system%N_of_atoms) = system%atom

! replicating ...
forall( i = 1:Replication_Factor ) temp( system%N_of_atoms*(i-1)+1:system%N_of_atoms*i ) = system%atom

counter   = 0
S_counter = 0   ! <== solvent 
F_counter = 0   ! <== fragment
P_counter = 0   ! <== polymer

If( any(system % atom % fragment == "P") ) max_nresid = maxval( system % atom %  nresid )   ! <== polymer

do k = n_z%inicio , n_z%fim 
do j = n_y%inicio , n_y%fim 
do i = n_x%inicio , n_x%fim 

    do n = 1 , system%N_of_atoms

        temp (counter+n) % xyz(1) = system % atom(n) % xyz(1) + i*system%box(1)
        temp (counter+n) % xyz(2) = system % atom(n) % xyz(2) + j*system%box(2)
        temp (counter+n) % xyz(3) = system % atom(n) % xyz(3) + k*system%box(3)

        if( temp(counter+n) % fragment == 'S' ) temp(counter+n)%nresid = system%atom(n)%nresid + S_counter
        if( temp(counter+n) % fragment == 'F' ) temp(counter+n)%nresid = system%atom(n)%nresid + F_counter
        if( temp(counter+n) % fragment == 'P' ) temp(counter+n)%nresid = system%atom(n)%nresid + P_counter

        if( abs(i)+abs(j)+abs(k) == 0 ) temp(counter+n) % copy = .true.

    end do

    counter    =  counter   +  system % N_of_atoms
    S_counter  =  S_counter +  maxval( system%atom%nresid , system%atom%fragment == 'S' )       &
                            -  minval( system%atom%nresid , system%atom%fragment == 'S' ) + 1
    F_counter  =  F_counter +  maxval( system%atom%nresid , system%atom%fragment == 'F' )       &
                            -  minval( system%atom%nresid , system%atom%fragment == 'F' ) + 1
    P_counter  =  P_counter +  max_nresid

end do
end do
end do    

! final setting of the system array ...
CALL move_alloc(from=temp,to=system%atom)
system%N_of_atoms = New_N_of_atoms

system%box(1) = ( n_x%fim - n_x%inicio + 1 ) * system%box(1)
system%box(2) = ( n_y%fim - n_y%inicio + 1 ) * system%box(2)
system%box(3) = ( n_z%fim - n_z%inicio + 1 ) * system%box(3)

system%N_of_Surface_atoms      =  system%N_of_Surface_atoms      * Replication_factor
system%N_of_Solvent_atoms      =  system%N_of_Solvent_atoms      * Replication_factor
system%N_of_Solvent_molecules  =  system%N_of_Solvent_molecules  * Replication_factor

where(system%atom%copy == .false.) system%atom%delete = .true.
CALL Eliminate_fragment( system )    

write(string,'(i3.3)') Replication_Factor
system%Surface_Characteristics = trim(system%Surface_Characteristics)//'-Replicated:'//string

end subroutine Replicate
!
!
!
!=======================================
 subroutine Nonbonding_Topology(system)
!=======================================
implicit none
type(universe) , intent(in) :: system

! local variables ...
integer             :: ati , atj , N_of_pairs
real*8              :: cutoff , ScaleFactor , distance
logical             :: mask
character(len=3)    :: residue
character(len=10)   :: file_name

! reading parameters ...
write(*,'(a)',advance='no') 'cutoff radius (Real)  = '
read (*,'(F8.4)') cutoff
write(*,'(a)',advance='no') 'scale factor  (Real)  = '
read (*,'(F8.4)') ScaleFactor
write(*,'(a)',advance='no') 'residue name  (len=3 , use capital letters) = '
read (*,'(A3)') residue


! count number of nonbonding pairs ...
N_of_pairs = 0
do ati = 1 , system%N_of_atoms
    do atj = ati+1 , system%N_of_atoms

        distance = sqrt( sum( (system%atom(ati)%xyz(:) - system%atom(atj)%xyz(:))**2 ) )

        mask = distance >= cutoff

        If( mask ) N_of_pairs = N_of_pairs + 1

    end do
end do

! save nonbonding pairs ...
file_name = residue//".inpt14"
OPEN(unit=3,file=file_name,status='new')
write(3,'(I9)') N_of_pairs

do ati = 1 , system%N_of_atoms
    do atj = ati+1 , system%N_of_atoms

        distance = sqrt( sum( (system%atom(ati)%xyz(:) - system%atom(atj)%xyz(:))**2 ) )

        mask = distance >= cutoff

        If( mask ) write(unit=3,fmt=10) ati , atj , ScaleFactor

    end do
end do

close(3)

10 format(2I6,F9.4)

end subroutine Nonbonding_Topology
!
!
!
!
!======================================================
 function vector_product(sys,at1,at2,at3) result(w_vec)
!======================================================
implicit none
type(universe) , intent(in) :: sys
integer        , intent(in) :: at1
integer        , intent(in) :: at2
integer        , intent(in) :: at3

!local variables ...
real*8 :: u_vec(3) , v_vec(3) , w_vec(3) 

u_vec = (sys% atom(at1)% xyz - sys% atom(at2)% xyz) 
v_vec = (sys% atom(at3)% xyz - sys% atom(at2)% xyz) 

w_vec(1) = u_vec(3)*v_vec(2) - u_vec(2)*v_vec(3)
w_vec(2) = u_vec(1)*v_vec(3) - u_vec(3)*v_vec(1)
w_vec(3) = u_vec(2)*v_vec(1) - u_vec(1)*v_vec(2)

w_vec = w_vec / sqrt( dot_product(w_vec,w_vec) )

end function vector_product
!
!
!
end module EDIT_routines

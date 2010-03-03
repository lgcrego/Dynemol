! Subroutine from preparing a liquid environment

module Solvated_m

    use type_m
    use constants_m
    use Allocation_m            , only : Allocate_UnitCell
    use Babel_m                 , only : System_Characteristics , trj
    use Structure_Builder       , only : Unit_Cell
    use Semi_Empirical_Parms    , only : Define_EH_Parametrization

    public :: Prepare_Solvated_System , DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 

    private

contains
!
!
!=============================================================
 subroutine Prepare_Solvated_System( Solvated_System , frame )
!=============================================================
implicit none
type(universe)   , intent(inout)  :: Solvated_System
integer          , intent(in)     :: frame

! local variables ...
integer                                         :: i , system_PBC_size , solvent_PBC_size , system_size
real*8                                          :: solvation_radius , Molecule_CG(3)
real*8              , allocatable               :: distance(:)
logical             , allocatable               :: mask(:)
type(atomic)        , allocatable               :: system_PBC(:) , system(:)
type(molecular)     , allocatable   , target    :: solvent_PBC(:) 
type(int_pointer)   , allocatable               :: nres(:) 

! local parameters ; number of 3D PBC unit-cells ...
integer , parameter :: PBC_Factor = 27  

If( nnx+nny+mmx+mmy /= 0 ) Pause " >>> Using Replication in Solvated_M <<< "

! identify the CG of the fragment ...
forall( i=1:3 ) Molecule_CG(i) = sum( trj(frame)%atom%xyz(i) , trj(frame)%atom%fragment == "M" ) / count(trj(frame)%atom%fragment == "M")

! place origin at GC ...
forall( i=1:trj(frame)%N_of_Atoms             ) trj(frame) % atom(i) % xyz(:)   = trj(frame) % atom(i) % xyz(:)   - Molecule_CG(:)
forall( i=1:trj(frame)%N_of_Solvent_Molecules ) trj(frame) % solvent(i) % CG(:) = trj(frame) % solvent(i) % CG(:) - Molecule_CG(:) 

! define the PBC system ...
system_PBC_size  = trj(frame) % N_of_atoms             *  PBC_Factor
solvent_PBC_size = trj(frame) % N_of_Solvent_Molecules *  PBC_Factor

allocate( system_PBC  (system_PBC_size)  )
allocate( solvent_PBC (solvent_PBC_size) )

CALL Apply_PBC( trj(frame) , system_PBC , solvent_PBC , nres )

solvation_radius = minval( trj(frame)%box ) * two / three

allocate( distance(solvent_PBC_size) )

! distance of solvent CG from Molecule CG ...
forall( i=1:solvent_PBC_size ) distance(i) = sqrt( sum(solvent_PBC(i)%CG**2) )

where( distance > solvation_radius ) solvent_PBC%nresid = 0

forall( i=1:system_PBC_size ) system_PBC(i)%nresid = nres(i)%PTR

allocate( mask(system_PBC_size) )

mask = (system_PBC%nresid /= 0 .AND. system_PBC%fragment == "S") .OR. (system_PBC%copy_No == 0 .AND. system_PBC%fragment /= "S") 

! define the solvation cell ...
system_size = count(mask)
allocate( system(system_size) )
system = pack(system_PBC , mask , system)

! build Solvated_System
CALL move_alloc( from=system , to=Solvated_System%atom )

CALL Sort_Fragments( Solvated_System )

Solvated_System%N_of_Atoms              =  system_size
Solvated_System%N_of_Solvent_Molecules  =  count( solvent_PBC%nresid /= 0 )
Solvated_System%System_Characteristics  =  trj(frame)%System_Characteristics
Solvated_System%box                     =  trj(frame)%box

deallocate( system_PBC , solvent_PBC , distance , mask ) 

if( frame == 1 ) then
    CALL dump_pdb( Solvated_System )
    write(*,'(/a)') '>>>  Saving seed.pdb  <<<'
end if

end subroutine Prepare_Solvated_System
!
!
!
!=========================================================
subroutine Apply_PBC( trj , trj_PBC , solvent_PBC , nres )
!=========================================================
implicit none
type(universe)                      , intent(in)        :: trj
type(atomic)        , target        , intent(inout)     :: trj_PBC(:) 
type(molecular)     , target        , intent(inout)     :: solvent_PBC(:)
type(int_pointer)   , allocatable   , intent(out)       :: nres(:) 

! local variables ...
integer :: i , j , k , n , atom , molecule , nresid , trj_PBC_size , solvent_PBC_size

! local parameters ; number of 3D PBC unit-cells ...
integer ,   parameter   :: PBC_Factor = 27  
integer ,   parameter   :: tres = 3 , nove = 9

! replicating the central cell to the surrounding cells ...
forall( i = 1:PBC_Factor ) 

    trj_PBC    ( trj%N_of_Atoms*(i-1)+1             : trj%N_of_atoms*i             ) = trj%atom

    solvent_PBC( trj%N_of_Solvent_Molecules*(i-1)+1 : trj%N_of_Solvent_Molecules*i ) = trj%solvent

end forall

trj_PBC_size      =  size( trj_PBC ) 
solvent_PBC_size  =  size( solvent_PBC )

! defining the coordinates for the surrounding cells ...
atom     = 0
molecule = 0
nresid   = 0

do k = -1,+1 
do j = -1,+1
do i = -1,+1 

    forall( n = 1:trj%N_of_Solvent_Molecules )

        solvent_PBC(molecule+n) % CG(1)    =  trj % solvent(n) % CG(1) + i * trj % box(1)
        solvent_PBC(molecule+n) % CG(2)    =  trj % solvent(n) % CG(2) + j * trj % box(2)
        solvent_PBC(molecule+n) % CG(3)    =  trj % solvent(n) % CG(3) + k * trj % box(3)
 
        solvent_PBC(molecule+n) % nresid   =  trj % solvent(n) % nresid + nresid
        solvent_PBC(molecule+n) % copy_No  =  i + j*tres + k*nove  

    end forall

    forall( n = 1:trj%N_of_atoms )

        trj_PBC(atom+n) % xyz(1)   =  trj % atom(n) % xyz(1) + i * trj % box(1)
        trj_PBC(atom+n) % xyz(2)   =  trj % atom(n) % xyz(2) + j * trj % box(2)
        trj_PBC(atom+n) % xyz(3)   =  trj % atom(n) % xyz(3) + k * trj % box(3)

        trj_PBC(atom+n) % nresid   =  trj % atom(n) % nresid + nresid
        trj_PBC(atom+n) % copy_No  =  i + j*tres + k*nove  

    end forall

    atom     =  atom     + trj % N_of_atoms
    molecule =  molecule + trj % N_of_Solvent_Molecules
    nresid   =  nresid   + trj % atom(trj%N_of_Atoms)%nresid

end do
end do
end do

allocate( nres(trj_PBC_size) )

forall( i=1:trj_PBC_size )  nres(i)%PTR => trj_PBC(i)%nresid

do j = 1 , solvent_PBC_size

    forall( i=1:trj_PBC_size , trj_PBC(i)%nresid == solvent_PBC(j)%nresid ) nres(i)%PTR => solvent_PBC(j)%nresid

end do

end subroutine Apply_PBC
!
!
!
!=======================
subroutine Dump_pdb(sys)
!=======================
implicit none 
type(universe)      , intent(inout) ::  sys

! local variables ...
integer ::  i , j , k , nr 

!----------------------------------------------
!     generate pdb file for GROMACS
!----------------------------------------------

OPEN(unit=4,file='seed.pdb',status='unknown')

write(4,6) 'COMPND' , sys%System_Characteristics
write(4,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do i = 1 , sys%N_of_atoms

            write(4,2)  'HETATM'                        ,  &    ! <== non-standard atom
                        i                               ,  &    ! <== global number
                        sys%atom(i)%MMSymbol            ,  &    ! <== atom type
                        ' '                             ,  &    ! <== alternate location indicator
                        sys%atom(i)%residue             ,  &    ! <== residue name
                        ' '                             ,  &    ! <== chain identifier
                        sys%atom(i)%nresid              ,  &    ! <== residue sequence number
                        ' '                             ,  &    ! <== code for insertion of residues
                        ( sys%atom(i)%xyz(k) , k=1,3 )  ,  &    ! <== xyz coordinates 
                        1.00                            ,  &    ! <== occupancy
                        0.00                            ,  &    ! <== temperature factor
                        ' '                             ,  &    ! <== segment identifier
                        ' '                             ,  &    ! <== here only for tabulation purposes
                        sys%atom(i)%Symbol              ,  &    ! <== chemical element symbol
                        ' '                                     ! <== charge on the atom
end do

write(4,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , sys%N_of_atoms , 0 , sys%N_of_atoms , 0
write(4,*) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,a2)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a6,a72)

end subroutine Dump_pdb
!
!
!==================================
subroutine Sort_Fragments( system )
!==================================
implicit none
type(universe) , intent(inout) :: system

!	local variables
integer        :: i , j , N_of_elements , N_of_atoms ,  iptr
type(universe) :: temp

allocate( temp%atom(1) , temp%solvent(1) )

N_of_atoms = size(system%atom)

do i = 1 , N_of_atoms-1

    iptr = i

    do j = i+1 , N_of_atoms
        if( LGT( system%atom(j)%fragment , system%atom(iptr)%fragment ) ) then
            iptr = j
        end if
    end do

    if( i /= iptr ) then
        temp%atom(1)      = system%atom(i)
        system%atom(i)    = system%atom(iptr)
        system%atom(iptr) = temp%atom(1)
    end if

end do

deallocate( temp%atom , temp%solvent )

end subroutine Sort_Fragments
!
!
!
!=========================================
 subroutine DeAllocate_TDOS( TDOS , flag )
!=========================================
implicit none
type(f_grid)    , intent(inout) :: TDOS
character(*)    , intent(in)    :: flag

! local parameter ...
integer , parameter   :: npoints = 1500

select case( flag )

    case( "alloc" )
        allocate( TDOS%grid(npoints) )
        allocate( TDOS%func(npoints) )
        allocate( TDOS%average(npoints) , source = 0.d0)

    case( "dealloc" )
        deallocate( TDOS%grid , TDOS%func , TDOS%average )

end select

end subroutine DeAllocate_TDOS
!
!
!
!=========================================
 subroutine DeAllocate_PDOS( PDOS , flag )
!=========================================
implicit none
type(f_grid)    , allocatable   , intent(inout) :: PDOS(:)
character(*)                    , intent(in)    :: flag

! local parameter ...
integer , parameter   :: npoints = 1500

! local variable ...
integer :: i , N_of_residues

select case( flag )

    case( "alloc" )

        if( allocated(trj) ) then
            N_of_residues = size(trj(1)%list_of_residues) 
        else
            N_of_residues = size(unit_cell%list_of_residues)
        end if

        allocate( PDOS(N_of_residues) )

        do i = 1 , N_of_residues
            allocate( PDOS(i)%grid(npoints) )
            allocate( PDOS(i)%func(npoints) )
            allocate( PDOS(i)%average(npoints) , source = 0.d0)
        end do

    case( "dealloc" )
        deallocate(PDOS)

end select

end subroutine DeAllocate_PDOS
!
!
!
!=========================================
 subroutine DeAllocate_SPEC( SPEC , flag )
!=========================================
implicit none
type(f_grid)    , intent(inout) :: SPEC
character(*)    , intent(in)    :: flag

! local parameter ...
integer , parameter   :: npoints = 1500

select case ( flag )

    case( "alloc" )
        allocate( SPEC%grid   (npoints) )
        allocate( SPEC%func   (npoints) )
        allocate( SPEC%peaks  (npoints) )
        allocate( SPEC%average(npoints) , source = 0.d0)

    case( "dealloc")
        deallocate( SPEC%grid , SPEC%func , SPEC%peaks , SPEC%average )

end select

end subroutine DeAllocate_SPEC
!
!
end module Solvated_m

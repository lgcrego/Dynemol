! Subroutine from preparing a liquid environment

module Solvated_m

    use type_m
    use constants_m
    use parameters_m            , only : nnx , nny , PBC , T_ , F_
    use Allocation_m            , only : Allocate_UnitCell
    use Babel_m                 , only : System_Characteristics , trj
    use Structure_Builder       , only : Unit_Cell

    public :: Prepare_Solvated_System , DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 

    private

! module variables
integer                         , save  :: N_of_Solvent_Molecules 
integer         , allocatable   , save  :: data_B(:,:)
logical                         , save  :: first_time = .true.
logical                         , save  :: done_preprocess
type(universe)                  , save  :: Solvated_System_0

contains
!
!
!=============================================================
 subroutine Prepare_Solvated_System( Solvated_System , frame )
!=============================================================
implicit none
type(universe)  , intent(out)   :: Solvated_System
integer         , intent(in)    :: frame

! local variables ...
integer                                         :: system_PBC_size , solvent_PBC_size , system_size , N_of_insiders , i
real*8                                          :: solvation_radius , solute_CG(3)
real*8              , allocatable               :: distance(:)
logical             , allocatable               :: mask(:)
type(atomic)        , allocatable               :: system_PBC(:) , system(:)
type(molecular)     , allocatable   , target    :: solvent_PBC(:) 
type(int_pointer)   , allocatable               :: nres(:)

! local parameters ; number of 3D PBC unit-cells ...
integer , parameter :: PBC_Factor = 27

! check-list ...
if( .not. any(trj(frame)%atom%solute) ) Pause " >>> Solute is not tagged <<< " 

If( nnx+nny+sum(PBC) /= 0 ) Pause " >>> Using Replication in Solvated_M <<< "

! identify the CG of the solute ...
forall( i=1:3 ) solute_CG(i) = sum( trj(frame)%atom%xyz(i) , trj(frame)%atom%solute == .true. ) / count(trj(frame)%atom%solute)

! place origin at GC ...
forall( i=1:trj(frame)%N_of_Atoms             ) trj(frame) % atom(i) % xyz(:)   = trj(frame) % atom(i) % xyz(:)   - solute_CG(:)
forall( i=1:trj(frame)%N_of_Solvent_Molecules ) trj(frame) % solvent(i) % CC(:) = trj(frame) % solvent(i) % CC(:) - solute_CG(:) 

! define the PBC System ...
system_PBC_size  = trj(frame) % N_of_atoms             *  PBC_Factor
solvent_PBC_size = trj(frame) % N_of_Solvent_Molecules *  PBC_Factor

allocate( system_PBC  (system_PBC_size)  )
allocate( solvent_PBC (solvent_PBC_size) )

solvation_radius = minval( trj(frame)%box ) * two / four   !2.8d0  ! <== best value 

allocate( distance(solvent_PBC_size) )

CALL Apply_PBC( trj(frame) , system_PBC , solvent_PBC , nres , distance )

If( first_time ) then

    N_of_insiders = count( distance <= solvation_radius )

    N_of_Solvent_Molecules = N_of_insiders

    first_time = .false.

else

    N_of_insiders = N_of_Solvent_Molecules 

end If

solvent_PBC( N_of_insiders+1 : solvent_PBC_size ) % nr = 0

forall( i=1:system_PBC_size ) system_PBC(i)%nr = nres(i)%PTR

allocate( mask(system_PBC_size) )

mask = (system_PBC%nr /= 0 .AND. system_PBC%fragment == "S") .OR. (system_PBC%copy_No == 0 .AND. system_PBC%fragment /= "S") 

! define the solvation cell ...
system_size = count(mask)
allocate( system(system_size) )
System = pack(system_PBC , mask , system)

! build Solvated_System
CALL move_alloc( from=System , to=Solvated_System%atom )

CALL Sort_Fragments( Solvated_System )

Solvated_System%N_of_Atoms              =  system_size
Solvated_System%N_of_Solvent_Molecules  =  count( solvent_PBC%nr /= 0 )
Solvated_System%System_Characteristics  =  trj(frame)%System_Characteristics
Solvated_System%box                     =  trj(frame)%box

! Correction of the solvent positions ...
if( done_preprocess == T_ ) CALL Correction_of_Solvent_Positions( Solvated_System )

! Preprocess for routine Correction_of_Solvent_Positions ...
if( done_preprocess == F_ ) CALL preprocess_CSP( data_B , Solvated_System , done_preprocess )

deallocate( system_PBC , solvent_PBC , distance , mask ) 

if( frame == 1 ) then
    CALL dump_pdb( Solvated_System )
    write(*,'(/a)') '>>>  Saving seed.pdb  <<<'
end if

end subroutine Prepare_Solvated_System
!
!
!
!=============================================================
subroutine preprocess_CSP( data_B , System , done_preprocess )
!=============================================================
implicit none
integer         , allocatable   , intent(out)   :: data_B(:,:)
type(universe)                  , intent(in)    :: System
logical                         , intent(inout) :: done_preprocess

! local variable ...
integer :: i , size_sys

size_sys = size( System%atom )

allocate( Solvated_System_0%atom(size_sys) )

allocate( data_B ( count(System%atom%fragment == 'S') , 2 ) )

Solvated_System_0 % atom = System % atom

data_B(:,1) = pack( System%atom%nr                  , System%atom%fragment == 'S' )
data_B(:,2) = pack( [( i , i=1,System%N_of_atoms )] , System%atom%fragment == 'S' )

done_preprocess = T_

end subroutine preprocess_CSP
!
!
!
!===================================================
subroutine Correction_of_Solvent_Positions( System )
!===================================================
implicit none
type(universe)  , intent(inout) :: System

! local variables ...
integer                 :: i , j , l , var , counter , size_S_sys , N_of_atoms_in_S_mol , flux
integer , allocatable   :: data_A(:,:) , atoms_out(:) , atoms_in(:)
type(universe)          :: System_temp

size_S_sys          = count(System%atom%fragment == 'S')
N_of_atoms_in_S_mol = size_S_sys / system%N_of_Solvent_Molecules

! store indices of atoms and nr's in data_A ...
allocate( data_A ( size_S_sys , 2 ) )

data_A(:,1) = pack( System%atom%nr                  , System%atom%fragment == 'S' )
data_A(:,2) = pack( [( i , i=1,System%N_of_atoms )] , System%atom%fragment == 'S' )

! number of molecules that have entered/exited the sphere ...
counter = 0
do i = 1 , size(data_A(:,1)) , N_of_atoms_in_S_mol
    do j = 1 , size(data_B(:,1)) , N_of_atoms_in_S_mol

        if( data_A(i,1) == data_B(j,1) ) counter = counter + 1
       
    end do
end do

! fix the system if there is flux of solvent molecules in/out of droplet ...
flux = size_S_sys - (N_of_atoms_in_S_mol*counter)

if( flux /= I_zero ) then

    allocate( atoms_out( flux ) , atoms_in( flux ) )

    ! save atoms that exit the sphere ...
    var = 1
    do i = 1 , size_S_sys , N_of_atoms_in_S_mol

        If( count(data_A(:,1) == data_B(i,1)) == I_zero ) then

            do l = 0 , N_of_atoms_in_S_mol - 1 
                atoms_out(var) = data_B(i,2) + l
                var = var + 1
            end do

        end if

    end do

    ! save atoms that entered the sphere ...
    var = 1
    do i = 1 , size_S_sys , N_of_atoms_in_S_mol

        If( count(data_B(:,1) == data_A(i,1)) == I_zero) then

            do l = 0 , N_of_atoms_in_S_mol - 1
                atoms_in(var) = data_A(i,2) + l
                var = var + 1
            end do

        end if

    end do

    deallocate( data_A )

    ! exchange atoms that exited by atoms that entered ...
    allocate( System_temp%atom( size(atoms_in) ) )

    System_temp%atom(:) = System%atom(atoms_in(:))

    System%atom = Solvated_System_0%atom
        
    System%atom(atoms_out(:)) = System_temp%atom(:)

    deallocate( atoms_in , atoms_out ,  System_temp%atom )

end if

end subroutine Correction_of_Solvent_Positions
!
!
!
!====================================================================
subroutine Apply_PBC( trj , trj_PBC , solvent_PBC , nres , distance )
!====================================================================
implicit none
type(universe)                      , intent(in)        :: trj
type(atomic)        , target        , intent(inout)     :: trj_PBC(:) 
type(molecular)     , target        , intent(inout)     :: solvent_PBC(:)
type(int_pointer)   , allocatable   , intent(out)       :: nres(:) 
real*8                              , intent(inout)     :: distance(:)

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

        solvent_PBC(molecule+n) % CC(1)    =  trj % solvent(n) % CC(1) + i * trj % box(1)
        solvent_PBC(molecule+n) % CC(2)    =  trj % solvent(n) % CC(2) + j * trj % box(2)
        solvent_PBC(molecule+n) % CC(3)    =  trj % solvent(n) % CC(3) + k * trj % box(3)
 
        solvent_PBC(molecule+n) % nr       =  trj % solvent(n) % nr + nresid
        solvent_PBC(molecule+n) % copy_No  =  i + j*tres + k*nove  

    end forall

    forall( n = 1:trj%N_of_atoms )

        trj_PBC(atom+n) % xyz(1)   =  trj % atom(n) % xyz(1) + i * trj % box(1)
        trj_PBC(atom+n) % xyz(2)   =  trj % atom(n) % xyz(2) + j * trj % box(2)
        trj_PBC(atom+n) % xyz(3)   =  trj % atom(n) % xyz(3) + k * trj % box(3)

        trj_PBC(atom+n) % nr       =  trj % atom(n) % nr + nresid
        trj_PBC(atom+n) % copy_No  =  i + j*tres + k*nove  

    end forall

    atom     =  atom     + trj % N_of_atoms
    molecule =  molecule + trj % N_of_Solvent_Molecules
    nresid   =  nresid   + trj % atom(trj%N_of_Atoms)%nr

end do
end do
end do

! distance of [solvent molecule CC] with [solute CG] at the origin ...
forall( i=1:solvent_PBC_size ) distance(i) = sqrt( sum(solvent_PBC(i)%CC**2) )

! sort solvent_PBC molecules by increasing distance from the origin ...
CALL Sort2(solvent_PBC,distance)

! nresidue[atomic] => nresidue[molecule]
allocate( nres(trj_PBC_size) )

forall( i=1:trj_PBC_size )  nres(i)%PTR => trj_PBC(i)%nr

do j = 1 , solvent_PBC_size

    forall( i=1:trj_PBC_size , trj_PBC(i)%nr == solvent_PBC(j)%nr ) nres(i)%PTR => solvent_PBC(j)%nr

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
integer ::  i , k 

!----------------------------------------------
!     generate pdb file for GROMACS
!----------------------------------------------

OPEN(unit=4,file='seed.pdb',status='unknown')

write(4,6) sys%System_Characteristics
write(4,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do i = 1 , sys%N_of_atoms

            write(4,2)  'HETATM'                        ,  &    ! <== non-standard atom
                        i                               ,  &    ! <== global number
                        sys%atom(i)%MMSymbol            ,  &    ! <== atom type
                        ' '                             ,  &    ! <== alternate location indicator
                        sys%atom(i)%residue             ,  &    ! <== residue name
                        ' '                             ,  &    ! <== chain identifier
                        sys%atom(i)%nr                  ,  &    ! <== residue sequence number
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
subroutine Sort_Fragments( System )
!==================================
implicit none
type(universe) , intent(inout) :: system

!	local variables
integer        :: i , j , N_of_atoms ,  iptr
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
type(f_grid)    , optional  , intent(inout) :: TDOS
character(*)                , intent(in)    :: flag

! local parameter ...
integer , parameter   :: npoints = 1500

select case( flag )

    case( "alloc" )
        allocate( TDOS%grid       (npoints) , source = 0.d0 )
        allocate( TDOS%func2      (npoints,2) , source = 0.d0 )
        allocate( TDOS%peaks2     (npoints,2) , source = 0.d0 )
        allocate( TDOS%occupation (npoints) , source = 0.d0 )
        allocate( TDOS%average2    (npoints,2) , source = 0.d0 )

    case( "dealloc" )
        deallocate( TDOS%grid , TDOS%func2 , TDOS%peaks2 , TDOS%occupation , TDOS%average2 )

end select

end subroutine DeAllocate_TDOS
!
!
!
!=========================================
 subroutine DeAllocate_PDOS( PDOS , flag )
!=========================================
implicit none
type(f_grid)    , optional  , allocatable   , intent(inout) :: PDOS(:)
character(*)                                , intent(in)    :: flag

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
            allocate( PDOS(i)%grid       (npoints) , source = 0.d0 )
            allocate( PDOS(i)%func       (npoints) , source = 0.d0 )
            allocate( PDOS(i)%func2      (npoints,2) , source = 0.d0 )
            allocate( PDOS(i)%peaks      (npoints) , source = 0.d0 )
            allocate( PDOS(i)%peaks2     (npoints,2) , source = 0.d0 )
            allocate( PDOS(i)%occupation (npoints) , source = 0.d0 )
            allocate( PDOS(i)%average    (npoints) , source = 0.d0 )
            allocate( PDOS(i)%average2   (npoints,2) , source = 0.d0 )

            if( allocated(trj) ) then
                PDOS(i) % residue = trj(1) % list_of_residues(i) 
            else
                PDOS(i) % residue = unit_cell % list_of_residues(i)
            end if
            
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
type(f_grid)    , optional  ,  intent(inout) :: SPEC
character(*)                , intent(in)     :: flag

! local parameter ...
integer , parameter   :: npoints = 1500

select case ( flag )

    case( "alloc" )
        allocate( SPEC%grid   (npoints) , source = D_zero )
        allocate( SPEC%func   (npoints) , source = D_zero )
        allocate( SPEC%peaks  (npoints) , source = D_zero )
        allocate( SPEC%average(npoints) , source = D_zero )

    case( "dealloc")
        deallocate( SPEC%grid , SPEC%func , SPEC%peaks , SPEC%average )

end select

end subroutine DeAllocate_SPEC
!
!
!
!=========================
 subroutine  Sort2(ra,xra)
!=========================
implicit none
type(molecular) , intent(inout) :: ra(:)
real*8          , intent(inout) :: xra(:)

! local variables ...
integer         :: l, n, ir, i, j
real*8          :: xrra
type(molecular) :: rra

!-----------------------------------------------------------
!  SORT ra(I) , SO THAT THE ELEMENTS xra(I) FOLLOWS TOGETHER
!-----------------------------------------------------------
      n = size(ra)
      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         xrra = xra(l)
         rra  = ra(l)
      else
         xrra = xra(ir)
         rra  = ra(ir)
         xra(ir) = xra(1)
         ra(ir)  = ra(1)
         ir = ir - 1
         if(ir .eq. 1) then
             xra(1) = xrra
             ra(1)  = rra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(xra(j) .lt. xra(j+1)) j = j + 1
        endif
      if(xrra .lt. xra(j)) then
        xra(i) = xra(j)
        ra(i)  = ra(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      xra(i) = xrra
      ra(i)  = rra
      goto 10

end subroutine sort2
!
!
!
!
end module Solvated_m

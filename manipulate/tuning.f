module FUNCTION_routines

use types_m             , only : universe , atomic
use diagnosis_m
use Read_Parms          , only : atomic_mass , Symbol_2_AtNo

type residues
    integer                        :: N_of_frags
    integer          , allocatable :: N_of_frag_elements(:)
    character(len=1) , allocatable :: set_of_frags(:)
end type residues

type(residues) , public :: res

contains
!
!
!
!
!=========================================
subroutine ad_hoc_tuning( system , frame )
!=========================================
implicit none
type(universe)  , intent(inout) :: system
integer         , optional      :: frame

! local variables ...
real*8       :: delta_t = 0.d0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      system % atom(:) % etc
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!----------------------------------
!      define Delta_t, if not yet 
!----------------------------------

 delta_t = 2.5d-4

! If( system%time == 0.d0 .AND. present(frame) ) system % time = delta_t * (frame -1)

!----------------------------------
!      define SPECIAL atoms 
!----------------------------------

!----------------------------------
!      define MM atom types 
!----------------------------------

!----------------------------------
!      define fragment's
!----------------------------------

!----------------------------------
!      define operations: 
! copy, delete, translate, rotate, group
!----------------------------------
system % atom % group = .true.
!----------------------------------
!      define resid's
!----------------------------------

!----------------------------------
!      define nresid's
!----------------------------------

!where(system % atom % Symbol /= "Si") system % atom %  nresid = 1

!----------------------------------
!     Selective_Dynamics
!----------------------------------

! Move or not to Move ...

!----------------------------------
!       charge of the atoms 
!----------------------------------

!OPEN(unit=33,file='charge.dat',status='old',iostat=ioerr,err=11)
!close(33)

!11 if( ioerr > 0 )then
!    stop 'charge.dat file not found; terminating execution'
!end if

!----------------------------------
!CALL Information_from_ITP( system ) 
!----------------------------------


end subroutine ad_hoc_tuning
!
!
!
!==================================
subroutine Define_Fragments(system)
!==================================
implicit none
type(universe)  , intent(inout) :: system

! local variables ...
character(1)  , allocatable :: fragment(:)

! setting up residue structure ...
res%N_of_frags = 7
allocate( res%set_of_frags(res%N_of_frags) , res%N_of_frag_elements(res%N_of_frags) )

res%set_of_frags        =  ["C","M","F","S","B","I","@"] 
res%N_of_frag_elements  =  0

!====== fragment table =======
!   C = cluster
!   M = molecule
!   F = fragment (general)
!   S = solvent
!   P = polymer
!   B = bulk
!   I = interface
!   D = delete
!   @ = cloning
!=============================

! helper array ...
allocate( fragment(system%N_of_atoms) )

! define MM atom types ...
!include 'TiO2.f'

CALL Symbol_2_AtNo(system%atom)

deallocate( fragment )

end subroutine Define_Fragments
!
!
!
!========================================
subroutine Information_from_ITP( system )
!========================================
implicit none
type(universe)  , intent(inout) :: system

! local variables ...
integer                         :: i , ioerr , fragment_size , nr , counter
type(atomic)    , allocatable   :: fragment(:)
character(72)                   :: fragment_characteristics , dumb

! reads data from itp file ...
OPEN(unit=3,file='fragment.itp',status='old',iostat=ioerr,err=10)

    read(3,*) fragment_size
    read(3,*) fragment_characteristics
    read(3,*) dumb

    allocate( fragment(fragment_size) )

    do i = 1 , fragment_size

        read(3,20) nr                       , &
                 fragment(i) % MMSymbol     , &
                 fragment(i) % nresid       , &
                 fragment(i) % resid        , &
                 fragment(i) % Symbol       , &
                 fragment(i) % nrcg         , &
                 fragment(i) % charge       , &
                 fragment(i) % mass        

    end do

close(3)

counter = 0
do i = 1 , system%N_of_atoms

    if( system%atom(i)%fragment == "F" ) then
        counter = counter + 1
        system % atom(i) % MMSymbol = fragment(counter) % MMSymbol
        system % atom(i) % resid    = fragment(counter) % resid   
        system % atom(i) % Symbol   = fragment(counter) % Symbol
        system % atom(i) % nrcg     = fragment(counter) % nrcg    
        system % atom(i) % charge   = fragment(counter) % charge  
        system % atom(i) % mass     = fragment(counter) % mass    
    end if
end do

10 if( ioerr > 0 )then
    stop 'fragment.itp file not found; terminating execution'
end if

20 format(i6,a11,i7,a7,a6,i6,f12.3,f11.4)

end subroutine Information_from_ITP
!
!
!
!=====================================================
 subroutine Solvent_residue_groups( sys , N_of_atoms )
!=====================================================
implicit none
type(universe)              , intent(inout) :: sys
integer         , optional  , intent(in)    :: N_of_atoms 

! local variables ...
integer :: i , nr , ref , Sol_mol_N_Of_atoms

! define solvent molecule number of atoms ...
If( present(N_of_atoms) ) then
    Sol_mol_N_Of_atoms = N_Of_atoms
else
    If( sys%solvent%N_of_atoms == 0 ) pause ">>> sys%solvent%N_of_atoms not defined in gmx.f : Solvent_residue_groups <<<"
    
    Sol_mol_N_Of_atoms = sys%solvent%N_of_atoms
end If

! define start point to count nresid ...
If( any( sys%atom%fragment /= "S") ) then
    ! solvent residues come always last ...
    ref = maxval( sys%atom%nresid , sys%atom%fragment /= "S") + 1
else
    ref = 1
end If

! start couting nresid ...
nr = 0
do i = 1 , size(sys%atom) 

    if( sys%atom(i)%fragment == 'S' ) then

        sys%atom(i)%nresid = ref + int( nr / Sol_mol_N_Of_atoms )

        nr = nr + 1

    end if

end do

end subroutine Solvent_residue_groups
!
!
!
end module FUNCTION_routines

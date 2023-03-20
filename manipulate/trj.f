module Trajectory_routines

use types_m
use constants_m
use RW_driver           , only : WritingRoutines
use Read_Parms          , only : MMSymbol_2_Symbol , Symbol_2_AtNo , Pack_Residues , Identify_Residues
use RW_routines         , only : view_XYZ , Initialize_System
use Function_routines   , only : ad_hoc_tuning
use diagnosis_m         , only : diagnosis
use EDIT_routines       , only : Eliminate_Fragment , ReGroup , Replicate , Translation
use Statistics_routines , only : Most_Representative_Configuration
use Topology_routines   , only : connect , dump_topol

public :: Read_Trajectories

private

contains
!
!
!========================================
subroutine Read_Trajectories( trj , sys )
!========================================
implicit none
type(universe)  , allocatable   , intent(out)   :: trj(:)
type(universe)                  , intent(out)   :: sys 

! local varibles ...
character(len=1)                :: file_format , YorN , wait
integer                         :: frame , n_frames , i , j , choice , option
real*8                          :: delta_t 
type(universe)  , allocatable   :: frozen_trj(:)

CALL system( "clear" )

write(*,'(/a)') ' Choose the format : '
write(*,'(/a)') ' (p) = PDB '
write(*,'(/a)') ' (v) = VASP '
write(*,'(/a)') ' (d) = DICE '
write(*,'(/a)') ' (x) = XYZ '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') file_format

select case( file_format )

    case( 'p' ) 
        CALL Read_PDB_Trajectories(trj) 

    case( 'v' )
        CALL Read_VASP_Trajectories( trj )

    case( 'x' )
        CALL Read_XYZ_Trajectories( trj )

    case( 'd' )
        CALL Read_DICE_Trajectories( trj )

end select

! select frame OR trajectory ...
do
    write(*,'(/a)') ' (1)  = edit trajectory with AD-HOC'
    write(*,'(/a)') ' (2)  = Translation Operation (only on first frame)'
    write(*,'(/a)') ' (3)  = select frame '      
    write(*,'(/a)') ' (4)  = save PDB trajectory '
    write(*,'(/a)') ' (5)  = re-GROUP molecules ( may need to use AD-HOC & Translation before; first frame must be united )'
    write(*,'(/a)') ' (6)  = DELETE fragment ( uses AD-HOC )'
    write(*,'(/a)') ' (7)  = RMSD of frames'
    write(*,'(/a)') ' (8)  = produce trajectory from single PDB frame'
    write(*,'(/a)') ' (9)  = Interpolate between frames'
    write(*,'(/a)') ' (10) = Reverse time direction in trajectory'
    write(*,'(/a)') ' (11) = Replicate structure'
    write(*,'(/a)') ' (0)  = DONE                '
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(I)') choice 

    select case( choice )

        case( 1 ) 

            do i = 1 , size(trj)
                CALL ad_hoc_tuning( trj(i) , i )
            end do

        case( 2 ) 

           call Translation( trj(1) )               

        case( 3 )

            write(*,'(/a)',advance='no') ' Frame number : '
            read (*,'(i5.5)'           ) frame
            allocate( sys%atom(size(trj(1)%atom)) )

            sys = trj(frame)

            sys%Surface_Characteristics = trj(1)%Surface_Characteristics

            deallocate(trj)

            write(*,'(/a)',advance='no') ' Save Frame (y,n) : [y] '
            read (*,'(a)') YorN

            If( YorN /= "n" ) then 
                call WritingRoutines( sys ) 
                stop
            else 
                return
            End If

        case( 4 )

            CALL Save_PDB_Trajectories( trj )

        case( 5 )

            write(*,'(/a)') ' Choose : '
            write(*,'(/a)') ' (1)  = solution in a box'
            write(*,'(/a)') ' (2)  = adsorbate on surface (no need for using ad-hoc)'  
            write(*,'(/a)',advance='no') '>>>   '
            read (*,'(I)') option

            select case ( option ) 

                case ( 1 )
                   do i = 1 , size(trj)
                        CALL ReGroup( trj(i) )
                   end do

                case ( 2 )
                   CALL ReGroup( trj )

            end select

        case( 6 )

            do i = 1 , size(trj)
                CALL Eliminate_Fragment( trj(i) )
            end do

        case( 7 )

            CALL Most_Representative_Configuration( trj , sys )

            return

        case( 8 ) 

            write(*,'(/a)',advance='no') ' Frame number : '
            read (*,'(i3.3)'           ) frame
            allocate( sys%atom(size(trj(1)%atom)) )

            write(*,'(/a)',advance='no') ' number of frame replicas : '
            read (*,'(i6.6)'           ) n_frames

            write(*,'(/a)',advance='no') ' time interval between frames : '
            read (*,'(f8.5)'           ) delta_t

            sys = trj(frame)

            sys%Surface_Characteristics = trj(1)%Surface_Characteristics

            ! allocate array of identical frames ...
            allocate( frozen_trj(n_frames) )

            do i = 1 , n_frames 
        
                allocate( frozen_trj(i) % atom(sys%N_of_atoms) )

                frozen_trj(i)%Surface_Characteristics = sys%Surface_Characteristics 

                forall(j = 1:3 ) frozen_trj(i) % atom % xyz(j) = sys % atom % xyz(j)

                frozen_trj(i) % atom % MMSymbol = sys % atom % MMSymbol
                frozen_trj(i) % atom % resid    = sys % atom % resid
                frozen_trj(i) % atom % nresid   = sys % atom % nresid
                frozen_trj(i) % atom % Symbol   = sys % atom % Symbol
                frozen_trj(i) % atom % fragment = sys % atom % fragment
                frozen_trj(i) % box             = sys % box
                frozen_trj(i) % N_of_atoms      = sys % N_of_atoms

                frozen_trj(i) % time            = delta_t * (i-1)

            end do

            deallocate(trj)

            CALL Save_PDB_Trajectories( frozen_trj )

       case( 9 ) 

            CALL Interpolate( trj )

       case( 10 ) 

            CALL Save_TimeReversed_PDB_Trajectories( trj )

       case( 11 )

            CALL system( "clear" )

            do i = 1 , size( trj )
                CALL Replicate( trj(i) )
            end do

            CALL connect( trj(1) )

       case default
            exit

    end select

end do    

write(*,'(/a)',advance='no') 'press ENTER '
read (*,'(a)') wait

end subroutine Read_Trajectories
!
!
!=====================================
 subroutine Read_PDB_Trajectories(trj)
!=====================================
type(universe)  , allocatable   , intent(out)   :: trj(:)

! local variables ...
integer                         :: openstatus , inputstatus , i , j , k , model , number_of_atoms , dumb_number
real*8                          :: time_1 , time_2 , delta_t 
character(1)                    :: test
character(4)                    :: keyword
character(5)                    :: MMSymbol_char
character(72)                   :: Surface_Characteristics

open(unit = 31, file = 'frames.pdb', status = 'old', action = 'read', iostat = openstatus)
if (openstatus > 0) stop " *** Cannot open the file *** "

read( unit = 31 , fmt = 31 ) Surface_Characteristics

! find the number of model frames ...
model=0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( inputstatus /= 0 ) exit
    if ( keyword == 'MODE' ) model = model + 1
end do

! return to the top of the file ...
rewind 31

! read number the atoms and time ...
read(unit = 31, fmt = 35, iostat = inputstatus) keyword
do
    if ( keyword == 'MODE' ) then
        exit
    else
        if ( keyword == 'TITL' ) then
            backspace 31
            do 
                read(unit = 31, fmt = 42,advance='no',iostat = inputstatus) test
                if(test == "=") then
                    read(unit = 31, fmt = 41, iostat = inputstatus) time_1 
                    exit    ! <== time_1 read
                end if
            end do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
        else
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
        end if
    end if
end do
number_of_atoms = 0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( keyword == 'ATOM' ) then
        backspace 31
        read(unit = 31, fmt = 32, iostat = inputstatus) dumb_number
        number_of_atoms = number_of_atoms + 1
    else
        if ( keyword == 'TITL' ) then
            backspace 31
            do 
                read(unit = 31, fmt = 42,advance='no',iostat = inputstatus) test
                if(test == "=") then
                    read(unit = 31, fmt = 41, iostat = inputstatus) time_2 
                    exit    ! <== time_2 read
                end if
            end do
            exit    ! <== leave outer do loop
        end if
    end if
end do

delta_t = time_2 - time_1

! return to the top of the file ...
rewind 31

! allocate array of frames ...
allocate( trj(model) )

do j = 1 , model
    if( j == 1 ) then
        trj(j)%time = 0.d0
        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( trj(j)%box(i) , i=1,3 )
            end if
            if( keyword == 'MODE' ) exit
        end do
    
        allocate( trj(j)%atom(number_of_atoms) )
    
        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 33, iostat = inputstatus) MMSymbol_char ,             &
                                                            trj(j)%atom(i)%resid ,      &
                                                            trj(j)%atom(i)%nresid ,     &
                                                            ( trj(j)%atom(i)%xyz(k) , k=1,3 ) 

            trj(j)%atom(i)%MMSymbol = adjustl(MMSymbol_char)

        end do
        CALL MMSymbol_2_Symbol(trj(j)%atom)
        CALL Changing_Fragment(trj(j)%atom%resid,trj(j)%atom)
    else
        trj(j)%time = trj(j-1)%time + delta_t
        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( trj(j)%box(i) , i=1,3 )
            end if
            if ( keyword == 'MODE' ) exit
        end do
    
        allocate( trj(j)%atom(number_of_atoms) )
    
        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 37, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
    end if

    print*, j

end do

trj%Surface_Characteristics = Surface_Characteristics
trj%N_of_Atoms = number_of_atoms

forall(i = 2:model )
    trj(i)%atom%MMSymbol  = trj(1)%atom%MMSymbol
    trj(i)%atom%resid     = trj(1)%atom%resid
    trj(i)%atom%nresid    = trj(1)%atom%nresid
    trj(i)%atom%Symbol    = trj(1)%atom%Symbol
    trj(i)%atom%fragment  = trj(1)%atom%fragment
end forall

! get list of residues in trj ...
CALL Identify_Residues( trj(1) )

! GROUP residues ...
do i = 1 , size(trj)
    CALL Pack_Residues( trj(i)%atom , trj(1)%list_of_residues ) 
end do

! Defining N_of_atoms and N_of_Solvent_molecules ...
trj%N_of_atoms = number_of_atoms
trj%N_of_Solvent_molecules = maxval( trj(1)%atom%nresid , trj(1)%atom%fragment == 'S' ) - minval( trj(1)%atom%nresid , trj(1)%atom%fragment == 'S' ) + 1

close(31)

! Formats ...
31 format(a72)
32 format(5x, i6)
33  format(t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3)                
35 format(a4)
36 format(7x, i7)
37 format(t32, f8.3, t40, f8.3, t48, f8.3)
38 format(10x, a1)
39 format(81x, f7.0)
40 format(6x, 3f9.3)
41 format(f10.5)
42 format(a1)

end subroutine Read_PDB_Trajectories
!
!
!=======================================
subroutine Read_VASP_Trajectories( trj )
!=======================================
implicit none
type(universe)  , allocatable   , intent(out)   :: trj(:)

! local variables ....
integer                         :: openstatus , inputstatus , N_of_atoms , lines , model , i , j , k , N_of_elements
integer         , allocatable   :: atom_No(:)
real*8                          :: x0 , y0 , z0 , factor , box(3) 
character(1)                    :: idx
character(72)                   :: System_Characteristics
character(3)    , allocatable   :: elements(:) 

! read System Characteristics and list of chemical elements ...
OPEN(unit = 13, file = 'VASP.trj', status = 'old', action = 'read', iostat = openstatus)
if( openstatus > 0 ) stop '*** Cannot open the file VASP.trj ***'

read(13 , '(A72)') System_Characteristics

! read multiplication factor for coordinates ...
read(13 , *) factor

! reads the unit cell vectors ...
read(13,*) x0 , y0 , z0
box(1) = dsqrt(x0*x0 + y0*y0 + z0*z0)
read(13,*) x0 , y0 , z0
box(2) = dsqrt(x0*x0 + y0*y0 + z0*z0)
read(13,*) x0 , y0 , z0
box(3) = dsqrt(x0*x0 + y0*y0 + z0*z0)

! from console: reads the number of atoms of each species ...
write(*,'(/a)',advance='no') 'Number of different chemical species (elements) in the system = '
read (*,*) N_of_elements

allocate( atom_No(N_of_elements) )
allocate( elements(N_of_elements) )

read( 13 , * ) (elements(i), i=1,N_of_elements)
read( 13 , * ) (atom_No (i), i=1,N_of_elements)

N_of_atoms = sum(atom_No)

! read the number of models (frames) ...
lines = 0
do
    read(unit = 13, fmt = 21, iostat = inputstatus) idx
    if( inputstatus /= 0 ) exit
    lines = lines + 1
end do

! N_of_frames = N_of_lines / (atoms + 1_blanck-line)
model = lines / (N_of_atoms+1)

rewind 13

! read the system ...
allocate( trj(model) )

do j = 1 , model
    if( j == 1 ) then

        ! skip 7 lines and read dummy character in the 8th line ...
        read(13,'(7/a)') idx

        allocate( trj(j)%atom(N_of_atoms) )
        trj(j)%N_of_atoms = N_of_atoms

        CALL Initialize_System( trj(j) )

        trj(j)%box  = box

        do i = 1 , N_of_atoms
            read(unit = 13, fmt = 22, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )

            do k = 1 , size(elements)
                if( (sum(atom_No(:k-1))+1 <= i) .AND. (i <= sum(atom_No(:k))) ) trj(j)%atom(i)%symbol = elements(k)
            end do

        end do

        CALL Symbol_2_AtNo(trj(j)%atom)

        ! MMSymbol = Symbol ...
        trj(j)%atom%MMSymbol  = trj(j)%atom%Symbol

    else

        read(unit = 13, fmt = 21, iostat = inputstatus) idx

        allocate( trj(j)%atom(N_of_atoms) )

        do i = 1 , N_of_atoms
            read(unit = 13, fmt = 22, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do

    end if
    print*,j
end do

! print elements on screen ...
print*, " " ; print*, elements

! rescale coordinates by the unit cell vectors ...
forall( i=1:model , j=1:3 ) trj(i) % atom % xyz(j) = trj(i) % atom % xyz(j) * box(j)

! fixing atoms to the original unit-cell ...
!do i=1,model 
!    where( trj(i)%atom%xyz(3) > 0.8*box(3) ) trj(i)%atom%xyz(3) = box(3) - trj(i)%atom%xyz(3)
!end do

! finishing the fixing ...
forall(i = 2:model )
    trj(i)%N_of_atoms     = trj(1)%N_of_atoms
    trj(i)%atom%Atno      = trj(1)%atom%Atno
    trj(i)%atom%resid     = trj(1)%atom%resid
    trj(i)%atom%nresid    = trj(1)%atom%nresid
    trj(i)%atom%Symbol    = trj(1)%atom%Symbol
    trj(i)%atom%fragment  = trj(1)%atom%fragment
    trj(i)%box            = box
    ! MMSymbol = Symbol
    trj(i)%atom%MMSymbol  = trj(1)%atom%Symbol
end forall

trj(1)%Surface_Characteristics = System_Characteristics

! formats ...
21 format(a1)
22 format(t2,f12.8,t14,f12.8,t26,f12.8)

end subroutine Read_VASP_Trajectories
!
!
!==================================
subroutine Changing_Fragment(FGT,a)
!==================================
implicit none
character(3)    , intent(in)  :: FGT(:)
type(atomic)    , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(FGT(i))
        case( 'CCC') 
            a(i)%fragment = 'C' 
        case( 'FFF') 
            a(i)%fragment = 'F' 
        case( 'ACN') 
            a(i)%fragment = 'S' 
        case( 'BLK') 
            a(i)%fragment = 'B' 
        case( 'SRF') 
            a(i)%fragment = 'I' 
    end select

 END DO

end subroutine Changing_Fragment
!
!
!
!=====================================
 subroutine Read_XYZ_Trajectories(trj)
!=====================================
implicit none
type(universe)  , allocatable   , intent(out) :: trj(:)

! local variables ....
character(1)                    :: idx , dumb
integer                         :: openstatus , inputstatus , atoms , i , j , k , model  , line
real*8                          :: box(3)
character(72)                   :: System_Characteristics

open(unit = 13, file = 'frames.xyz', status = 'old', action = 'read', iostat = openstatus)
if( openstatus > 0 ) stop '*** Cannot open the file frames.xyz , file not found ***'

! reads information about the system ...
read(13,*) atoms 
read(13,*) System_Characteristics

! reads the unit cell vectors for Direct coordinate mode
read(13,*) box(1)
read(13,*) box(2)
read(13,*) box(3)

! read the number of models ...
model = 0
line  = 0
do
    read(unit = 13, fmt = 21 , iostat = inputstatus) idx
    if( inputstatus /= 0 ) exit
    line = line + 1
    if( mod(line,atoms+2) == 0 ) then
        model = model + 1
    end if
end do

rewind 13

! prepare to read the system ...
allocate( trj(model) )
allocate( trj(1)%atom(atoms) )
CALL Initialize_System( trj(1) )

! read coordinate from the frames ...
do j = 1 , model
print*, j
    if( j == 1 ) then

        read(13,*) atoms 
        read(13,*) System_Characteristics
        read(13,*) trj(1)%box(1)
        read(13,*) trj(1)%box(2)
        read(13,*) trj(1)%box(3)

        trj(1) % Surface_Characteristics = System_Characteristics

        do i = 1 , atoms
            ! it is more robust by not using tabulation        
            ! read(unit = 13, fmt = 22, iostat = inputstatus) trj(j)%atom(i)%Symbol, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
            read(13,*, iostat = inputstatus) trj(j)%atom(i)%Symbol, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
        
        CALL Symbol_2_Atno(trj(j)%atom)

    else

        read(13,*) dumb
        read(13,*) dumb

        allocate( trj(j)%atom(atoms) )
        CALL Initialize_System( trj(j) )

        do i = 1 , atoms
            ! it is more robust by not using tabulation        
            ! read(unit = 13, fmt = 22, iostat = inputstatus) trj(j)%atom(i)%Symbol, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
            read(13,*, iostat = inputstatus)  trj(j)%atom(i)%Symbol, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do

    end if

end do

trj(1)%N_of_atoms = atoms

! Copy information from trj(1) to trj(:) ...
forall(i = 2:model )
    trj(i) % N_of_atoms      = trj(1) % N_of_atoms
    trj(i) % atom % Atno     = trj(1) % atom % Atno
    trj(i) % atom % Symbol   = trj(1) % atom % Symbol
    trj(i)%box               = trj(1) % box
end forall

! formats ...
21 format(a1)
22 format(t2, a2, t7, f11.6, t19, f11.6, t31, f11.6)

end subroutine Read_XYZ_Trajectories
!
!
!
!
!======================================
 subroutine Read_DICE_Trajectories(trj)
!======================================
implicit none
type(universe)  , allocatable   , intent(out) :: trj(:)

! local variables ....
character(1)                    :: dumb
integer                         :: openstatus , inputstatus , atoms , i , j , k , model , lines
character(72)                   :: System_Characteristics

open(unit = 13, file = 'DICE.trj', status = 'old', action = 'read', iostat = openstatus)
if( openstatus > 0 ) stop '*** Cannot open the file DICE.trj ***'

! read the number the atoms ...
read(unit = 13, fmt = 20, iostat = inputstatus) atoms


! number of configurations ...
rewind 13
lines = 0
do
    read(unit = 13, fmt = 21, iostat = inputstatus) dumb
    if( inputstatus /= 0 ) exit
    lines = lines + 1
end do
model = lines / (atoms+2)

! prepare to read the system ...
allocate( trj(model) )
allocate( trj(1)%atom(atoms) )
CALL Initialize_System( trj(1) )

! read structure characteristics ...
write(*,'(\a)',advance='no') ' Enter the System Characteristics : '
read (*,'(\a72)'           )   System_Characteristics
trj(1) % Surface_Characteristics = System_Characteristics

! start reading coordinate from the frames ...
rewind 13
do j = 1 , model

    if( j == 1 ) then
        read(unit = 13, fmt = 21, iostat = inputstatus) dumb

        ! read PBC vectors ... 
        read(13,24) trj(j)%box(1) , trj(j)%box(2) , trj(j)%box(3)

        do i = 1 , atoms
            read(unit = 13, fmt = 22, iostat = inputstatus) trj(j)%atom(i)%Symbol, ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
        CALL Symbol_2_AtNo(trj(j)%atom)
        CALL ad_hoc_tuning(trj(j))
    else
        read(unit = 13, fmt = 21, iostat = inputstatus) dumb

        read(13,24) trj(j)%box(1) , trj(j)%box(2) , trj(j)%box(3)

        allocate( trj(j)%atom(atoms) )
        CALL Initialize_System( trj(j) )

        do i = 1 , atoms
            read(unit = 13, fmt = 23, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
    end if

end do

trj(1)%N_of_atoms = atoms

! no MMSymbol is provided in DICE ...
trj(1)%atom%MMSymbol = trj(1)%atom%Symbol

! Copy information from trj(1) to trj(:) ...
forall(i = 2:model )
    trj(i) % N_of_atoms      = trj(1) % N_of_atoms
    trj(i) % atom % Atno     = trj(1) % atom % Atno
    trj(i) % atom % Symbol   = trj(1) % atom % Symbol
    trj(i) % atom % MMSymbol = trj(1) % atom % MMSymbol
    trj(i) % atom % resid    = trj(1) % atom % resid   
    trj(i) % atom % nresid   = trj(1) % atom % nresid   
    trj(i) % atom % fragment = trj(1) % atom % fragment 
end forall

! formats ...
20 format(I)
21 format(a13)
22 format(t4,a2, t10, f10.5, t25, f10.5, t40, f10.5)
23 format(t10, f10.5, t25, f10.5, t40, f10.5)
24 format(t37,3f9.4)

end subroutine Read_DICE_Trajectories
!
!
!
!=====================================
 subroutine Save_PDB_Trajectories(trj)
!=====================================
type(universe)  , allocatable   , intent(inout)   :: trj(:)

! local variables ...
character(len=1)      :: wait
integer               :: i , j , k , frame_step
real*8  , allocatable :: pm(:)

allocate( pm(trj(1)%N_of_atoms) )

write(*,'(/a)',advance='no') ' Saving with Frame step : '
read (*,'(i3.3)'           ) frame_step

OPEN( unit=4 , file='frames-output.pdb' , status='unknown' , action="write" )

write(4,6) trj(1)%Surface_Characteristics

do j = 1 , size(trj) , frame_step

    write(4,4) 'REMARK' , 'manipulated by edview'
    write(4,5) 'TITLE'  , 'manipulated by edview     t= ',trj(j)%time
    write(4,4) 'REMARK' , 'manipulated by edview'
    write(4,1) 'CRYST1' , trj(j)%box(1) , trj(j)%box(2) , trj(j)%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
    write(4,3) 'MODEL' , j

    ! where MMSymbol is not defined MMSymbol = symbol ...
    where( trj(j)%atom%MMSymbol == "XXX" ) trj(j)%atom%MMSymbol = trj(j)%atom%symbol

    do i = 1 , trj(j)%N_of_atoms

            write(4,2)  'ATOM  '                            ,  &    ! <== non-standard atom
                        i                                   ,  &    ! <== global number
                        trj(j)%atom(i)%MMSymbol             ,  &    ! <== atom type
                        ' '                                 ,  &    ! <== alternate location indicator
                        trj(j)%atom(i)%resid                ,  &    ! <== residue name
                        ' '                                 ,  &    ! <== chain identifier
                        trj(j)%atom(i)%nresid               ,  &    ! <== residue sequence number
                        ' '                                 ,  &    ! <== code for insertion of residues
                        ( trj(j)%atom(i)%xyz(k) , k=1,3 )   ,  &    ! <== xyz coordinates 
                        1.00                                ,  &    ! <== occupancy
                        0.00                                ,  &    ! <== temperature factor
                        ' '                                 ,  &    ! <== segment identifier
                        ' '                                 ,  &    ! <== here only for tabulation purposes
                        trj(j)%atom(i)%symbol               ,  &    ! <== chemical element symbol
                        trj(j)%atom(i)%charge                       ! <== charge on the atom
    end do

    ! check and print topological connections ...
    If( allocated(trj(j)%topol) ) CALL dump_topol(trj(j),4)

    write(4,'(a)') 'TER'
    write(4,'(a)') 'ENDMDL'

end do

close(4)

deallocate( pm )

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)              
3 FORMAT(a5,i8)
4 FORMAT(a6,t15,a21)
5 FORMAT(a5,t15,a35,f12.7)
6 FORMAT(a72)

write(*,'(/a)') ' >>> frames-output.pdb : writing done, press any key <<<'
read (*,'(a)') wait

end subroutine Save_PDB_Trajectories
!
!
!
!==================================================
 subroutine Save_TimeReversed_PDB_Trajectories(trj)
!==================================================
type(universe)  , allocatable   , intent(inout)   :: trj(:)

! local variables ...
character(len=1)    :: wait
integer             :: i , j , k 
real*8              :: delta_t

OPEN( unit=4 , file='frames-output.pdb' , status='unknown' , action="write" )

write(4,6) trj(1)%Surface_Characteristics

!determine delta_t ...
delta_t = trj(2)%time - trj(1)%time

do j = size(trj) , 1 , -1

    write(4,4) 'REMARK' , 'manipulated by edview'
    write(4,5) 'TITLE'  , 'manipulated by edview     t= ', (size(trj) - j)* delta_t
    write(4,4) 'REMARK' , 'manipulated by edview'
    write(4,1) 'CRYST1' , trj(j)%box(1) , trj(j)%box(2) , trj(j)%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
    write(4,3) 'MODEL' , j

    ! where MMSymbol is not defined MMSymbol = symbol ...
    where( trj(j)%atom%MMSymbol == "XXX" ) trj(j)%atom%MMSymbol = trj(j)%atom%symbol

    do i = 1 , trj(j)%N_of_atoms

            write(4,2)  'ATOM  '                            ,  &    ! <== non-standard atom
                        i                                   ,  &    ! <== global number
                        trj(j)%atom(i)%MMSymbol             ,  &    ! <== atom type
                        ' '                                 ,  &    ! <== alternate location indicator
                        trj(j)%atom(i)%resid                ,  &    ! <== residue name
                        ' '                                 ,  &    ! <== chain identifier
                        trj(j)%atom(i)%nresid               ,  &    ! <== residue sequence number
                        ' '                                 ,  &    ! <== code for insertion of residues
                        ( trj(j)%atom(i)%xyz(k) , k=1,3 )   ,  &    ! <== xyz coordinates 
                        1.00                                ,  &    ! <== occupancy
                        0.00                                ,  &    ! <== temperature factor
                        ' '                                 ,  &    ! <== segment identifier
                        ' '                                 ,  &    ! <== here only for tabulation purposes
                        trj(j)%atom(i)%symbol               ,  &    ! <== chemical element symbol
                        trj(j)%atom(i)%charge                       ! <== charge on the atom
    end do

    write(4,'(a)') 'TER'
    write(4,'(a)') 'ENDMDL'

end do

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)              
3 FORMAT(a5,i8)
4 FORMAT(a6,t15,a21)
5 FORMAT(a5,t15,a35,f9.4)
6 FORMAT(a72)

write(*,'(/a)') ' >>> frames-output.pdb : done with writing of Time Reversed trajectory, press any key <<<'
read (*,'(a)') wait

end subroutine Save_TimeReversed_PDB_Trajectories
!
!
!
!===========================
 subroutine Interpolate(trj)
!===========================
implicit none
type(universe)  , allocatable   , intent(inout)   :: trj(:)

!local variables ...
integer                      :: n , i , j , k , new_size, N_of_atoms , indx
real                         :: alpha_box(3) , dt
type(atomic)   , allocatable :: alpha(:)
type(universe) , allocatable :: tmp(:)

write(*,'(/a)',advance='no') 'Number of frames to be interpolated within a slice : '
read (*,*) n

N_of_atoms = trj(1) % N_of_atoms

!size of new trajetory ...
new_size = size(trj)*(n+1) - n
allocate( tmp(new_size) )
tmp % N_of_atoms = N_of_atoms

! broadcast basic information ...
do i = 1 , new_size

    allocate( tmp(i) % atom(N_of_atoms) )

    tmp(i)%Surface_Characteristics = trj(1)%Surface_Characteristics

    tmp(i) % atom % MMSymbol = trj(1) % atom % MMSymbol
    tmp(i) % atom % resid    = trj(1) % atom % resid
    tmp(i) % atom % nresid   = trj(1) % atom % nresid
    tmp(i) % atom % Symbol   = trj(1) % atom % Symbol
    tmp(i) % atom % fragment = trj(1) % atom % fragment
    tmp(i) % atom % charge   = trj(1) % atom % charge  
    tmp(i) % N_of_atoms      = trj(1) % N_of_atoms

end do

! copy information of trj into tmp ...
do i = 1 , size(trj) 

    tmp(i+(i-1)*n) % box     = trj(i) % box
    tmp(i+(i-1)*n) % time    = trj(i) % time

    forall(j = 1:3 ) tmp(i+(i-1)*n) % atom % xyz(j) = trj(i) % atom % xyz(j)

end do    

!start LINEAR interpolation ...
allocate( alpha(N_of_atoms) )

do i = 1 , size(trj)-1

    dt = trj(i+1) % time  - trj(i) % time

    forall(k = 1:3) alpha(:)%xyz(k) = ( trj(i+1)%atom(:)%xyz(k) - trj(i)%atom(:)%xyz(k) ) / dt

    alpha_box  = ( trj(i+1)%box - trj(i)%box ) / dt

    do j = 1 , n

        indx = i+j+(i-1)*n

        forall(k=1:3) tmp(indx)%atom(:)%xyz(k) = trj(i)%atom(:)%xyz(k) + alpha(:)%xyz(k)*dt*dfloat(j)/dfloat(n+1)
        
        tmp(indx) % box(:) = trj(i)%box(:) + alpha_box(:)*dt*dfloat(j)/dfloat(n+1)

        tmp(indx) % time   = trj(i)%time + dt*dfloat(j)/dfloat(n+1)

    end do

end do

CALL move_alloc( from=tmp , to=trj )

end subroutine Interpolate
!
!
!
end module Trajectory_routines

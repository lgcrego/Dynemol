module Elastic_Band

use types_m
use ansi_colors
use Read_Parms
use RW_routines       , only : Sort_Chemical_Elements
use FUNCTION_routines
use GMX_routines      , only : Dump_PDB
use EDT_util_m        , only: translation_mask

public :: Build_Configurations , Pick_Configurations

private

contains
!
!
!
!============================================
subroutine Pick_Configurations( trj , f_type)
!============================================
    implicit none
    type(universe)  , allocatable   , intent(inout) :: trj(:)
    integer                         , intent(in)    :: f_type
    
    ! local varibles ...
    integer                 :: n, N_of_frames
    integer, allocatable    :: frames(:), N_of_elements(:)
    integer                 :: n1, n2, step
    character(21)           :: string
    character(4), parameter :: file_format(2) = ["vasp","pdb"]
    
    CALL system( "clear" )
    
    ! ------------------------------------------------------------
    ! Read frame indices
    ! ------------------------------------------------------------
    write(*,'(/a)') bold//blue// ">>> Extracting frames from the pre-existing trajectory <<<"//reset
    
    write(*,'(/a)', advance='no') green//"start: "//reset
    read (*,'(i3)') n1
    
    write(*,'(/a)', advance='no') green//"end  : "//reset
    read (*,'(i3)') n2
    
    write(*,'(/a)', advance='no') green//"step : "//reset
    read (*,'(i3)') step
    
    ! if N_of_frames > size(trj) select 100 snapshots equally spaced ...
    N_of_frames = (n2 - n1) / max(step,1) + 1
    N_of_frames = min(N_of_frames, size(trj))
    
    allocate(frames(N_of_frames))
    
    ! ------------------------------------------------------------
    ! Select frames
    ! ------------------------------------------------------------
    if (N_of_frames < size(trj)) then
    
        write(*,'(/a)') yellow//"Selecting a subset of frames..."//reset
    
        do n = 1, N_of_frames
            frames(n) = n1 + (n-1) * step
        end do
    
    else
    
        write(*,'(/a)') yellow//"Using all available frames."//reset
        frames = [(n, n=1, size(trj))]
    
    end if
    
    ! Pre-process system & environment ...
    select case( file_format(f_type) )

       case( "vasp" )
           CALL Pre_process( trj(1) , N_of_elements )
           ! save frames in elastic band format ...
           do n = 1 , N_of_frames
               write(string,'(A14,i7.7)')  'frame index = ',frames(n)
               CALL save_Confgs( trj(frames(n)) , N_of_elements , string , n-1 , file_format(f_type) )   
           end do

        case( "pdb" ) 
           print*, bold//red//"to be implemented as required: save pdb trajectoty"//reset
     
    end select

end subroutine Pick_Configurations
!
!
!
!================================================
subroutine Build_Configurations( sys, f_type )
!================================================
implicit none 
type(universe), target, intent(inout) :: sys
integer               , intent(in)    :: f_type

! local variables ...
integer                           :: i, n, n_steps 
real*8                           :: versor(3)
real*8                           :: displacement, step
character(len=96)                :: string
character(4)       , parameter   :: file_format(2)=["vasp","pdb"] 
type(universe)     , allocatable :: trj(:)

CALL system( "clear" )

! ------------------------------------------------------------
! Read frame indices
! ------------------------------------------------------------
write(*,'(/a)') bold//blue// &
    ">>> Creating frames from a single structure <<<"//reset

write(*,'(/a)') green//"Purpose:"//reset
write(*,'(/a)') "  Generate a set of configurations by moving a selected atomic fragment along a prescribed path."

write(*,'(/a)') green//"Required steps:"//reset
write(*,'(/a)') "  1) Define the atomic fragment that will move (a single atom or a group of atoms)."
write(*,'(/a)') "  2) Define the path length (minimum and maximum distance, in consistent units)."
write(*,'(/a)') "  3) Specify the number of configurations to be generated."
write(*,'(/a)') "  4) Define the direction of the motion (typically a bond vector or user-defined axis)."

write(*,'(/a)') green//"Notes:"//reset
write(*,'(/a)') bold//"  - The fragment must be placed at the initial position before starting this routine."//reset
write(*,'(/a)')       "  - All generated frames preserve the original system topology."
write(*,'(/a)')       "  - This procedure is typically used to scan bond stretching,"
write(*,'(/a)')       "    adsorption paths, or constrained reaction coordinates."

! ------------------------------------------------------------ 
! Input data
! ------------------------------------------------------------
write(*,'(/a)') bold//green//"Step 1: define fragment that will move."//reset
call translation_mask(sys)

write(*,'(/a)') bold//cyan//"Steps 2,3,4: define reaction path (origin = current position)."//reset

write(*,'(/a)') green//"> Length of the reaction path "//yellow//"(Angs. units):"//reset
write(*,'(a)',advance='no') yellow//">>> "//reset
read (*,*) displacement
write(*,'(/a)') green//"> Enter the number of steps:"//reset
write(*,'(a)',advance='no') yellow//">>> "//reset
read (*,*) n_steps

step = displacement / float(n_steps)

call define_translation_vector( sys, versor )

! generate configurations ...
allocate( trj(0:n_steps) )

do n = 0, n_steps
    print*, "step ", n

    displacement = n * step  
    trj(n) = sys
    
    write(string,'(A55,F7.4)') "reaction path generated by DynEMol, this bond displacement = ", displacement
 
    trj(n)% System_Characteristics = string

    do concurrent ( i=1:sys%N_of_atoms , sys%atom(i)%translate ) 
            trj(n)%atom(i)%xyz = sys%atom(i)%xyz + versor * displacement
    end do

end do

if( file_format(f_type) == "pdb") then
    call Save_PDB_Trajectories(trj) 
end if

end subroutine Build_Configurations
!
!
!
!===============================================================================
subroutine save_Confgs( structure , N_of_elements , title , frame , file_format)
!===============================================================================
implicit none 
type(universe) , intent(inout) :: structure
integer        , intent(in) :: N_of_elements(:)
character(*)   , intent(in) :: title
integer        , intent(in) :: frame
character(*)   , intent(in) :: file_format

!	local variables ...
integer             :: i , j
character(len=3)    :: string
character(len=11)   :: dir
real*8  , parameter :: bottom_spacer = 3.0 !(Angs)
type(universe)      :: sys

write(string,'(i3.3)') frame
dir = "images/"//string//"/"
CALL system("mkdir "//dir)


select case( file_format )

    case('vasp') !-----------------------------------------------------------

        allocate( sys%atom(structure%N_of_atoms) )
        sys = structure

!       sort by chemical element ...
        CALL Sort_Chemical_Elements( sys )

!       rescale coordinates by the unit cell vectors.
        sys%atom%xyz(1) = ( sys%atom%xyz(1) - minval(sys%atom%xyz(1))                 ) / sys%box(1)
        sys%atom%xyz(2) = ( sys%atom%xyz(2) - minval(sys%atom%xyz(2))                 ) / sys%box(2)
        sys%atom%xyz(3) = ( sys%atom%xyz(3) - minval(sys%atom%xyz(3)) + bottom_spacer ) / sys%box(3)

!       ----------------------------------------------
!             generate input file for VASP
!       ----------------------------------------------

        OPEN(unit=4,file=dir//'POSCAR',status='unknown')

        write(4,'(A96)'   ) title
        write(4,'(A2)'    ) '1.' 
        write(4,'(3F12.5)') sys%box(1)     , 0.d0        , 0.d0
        write(4,'(3F12.5)') 0.d0           , sys%box(2)  , 0.d0 
        write(4,'(3F12.5)') 0.d0           , 0.d0        , sys%box(3)  
        write(4,*         ) ( N_of_elements(i),i=1,size(N_of_elements) )
        write(4,'(A18)'   ) 'Selective Dynamics'
        write(4,'(A6)'    ) 'Direct'

        do i = 1 , size(sys%atom)
            write(4,100) (sys%atom(i)%xyz(j),j=1,3) , 'F' , 'F' , 'F' 
        end do

        deallocate( sys%atom )

    case('pdb') !-----------------------------------------------------------

        OPEN(unit=4,file=dir//'seed.pdb',status='unknown')

        CALL Dump_pdb( structure , title )

end select

CLOSE(4)

100   FORMAT(F10.5,F10.5,F10.5,3A3)

end subroutine save_Confgs
!
!
!
!============================================
subroutine Pre_process( sys , N_of_elements )
!============================================
implicit none 
type(universe)                  , intent(inout) :: sys
integer         , allocatable   , intent(out)   :: N_of_elements(:)

! local variables ...
integer                         :: indx
integer     , allocatable       :: temp(:)
type(universe)                  :: dumb
character(len=3)                :: chemical

CALL system( "clear" )

! system information ...
print*, 'Ti = ', count(sys%atom%symbol=='Ti')
print*, 'O  = ', count(sys%atom%symbol=='O')
print*, 'C  = ', count(sys%atom%symbol=='C')
print*, 'H  = ', count(sys%atom%symbol=='H')
print*, 'N  = ', count(sys%atom%symbol=='N')
print*, 'S  = ', count(sys%atom%symbol=='S')
print*, 'Al = ', count(sys%atom%symbol=='Al')

! need to sort by chemical element ...
allocate( dumb%atom(sys%N_of_atoms) )
dumb = sys
CALL Sort_Chemical_Elements( dumb )

print*, dumb%atom%symbol

! find abundance of chemical elements
allocate( temp(20) )
chemical = dumb%atom(1)%symbol
temp(1) = count(dumb%atom%symbol==chemical)
indx = 1
do 
    if( sum(temp(1:indx)) == dumb%N_of_atoms ) EXIT
    chemical = dumb % atom(sum(temp(1:indx))+1) % symbol
    indx = indx + 1
    temp(indx) = count(dumb%atom%symbol==chemical)
end do
allocate( N_of_elements(indx) , source=temp )
deallocate( dumb%atom , temp )

! preparing the environment ...
CALL system("rm -r -f images")
CALL system("mkdir images")

end subroutine pre_process
!
!
!
!--------------------------
subroutine  sort2(ra,ira)
!--------------------------
real*8  , intent(inout) :: ra(:)
integer , intent(inout) :: ira(:)

! local variables
real    :: rra
integer :: irra, l, n, ir, i, j

!----------------------------------------------------------
!  SORT A(I) , SO THAT THE ELEMENTS IA(I) FOLLOW TOGETHER
!----------------------------------------------------------

      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         rra  = ra(l)
         irra = ira(l)
      else
         rra  = ra(ir)
         irra = ira(ir)
         ra(ir)  = ra(1)
         ira(ir) = ira(1)
         ir = ir - 1
         if(ir .eq. 1) then
             ra(1)  = rra
             ira(1) = irra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(ra(j) .lt. ra(j+1)) j = j + 1
        endif
      if(rra .lt. ra(j)) then
        ra(i)  = ra(j)
        ira(i) = ira(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      ra(i)  = rra
      ira(i) = irra
      goto 10

end subroutine sort2
!
!
!
!=======================================================
subroutine define_translation_vector( system, T_versor )
!=======================================================
    implicit none 
    type(universe), intent(inout) :: system
    real*8        , intent(out)   :: T_versor(3)
    
    ! local variables ...
    integer :: option , at1 , at2 , at3
    real*8  , save :: T_vector(3)
    
    ! define translation vector
    write(*,'(/a)') bold//orange//"define translation vector."//reset
    write(*,'(/a)') yellow//"> Use: (1)-cartesian axis , (2)-ad-hoc vector ,  (3)-normal vector? "//reset
    write(*,'(a)',advance='no') yellow//">>> "//reset
    read(*,*) option
    
    select case (option)
    
       case(1)
            ! use cartesian vectors 
            write(*,'(/a)') '> enter translation vector (T_x,T_y,T_z) as a Real number:'
            write(*,'(a)',advance='no') 'T_x = '
            read (*,*) T_vector(1)
            write(*,'(a)',advance='no') 'T_y = '
            read (*,*) T_vector(2)
            write(*,'(a)',advance='no') 'T_z = '
            read (*,*) T_vector(3)
            T_versor = T_vector / sqrt( dot_product(T_vector,T_vector) )
    
       case(2) 
            ! define translation vector 
            write(*,'(/a)') '> define vector: at1 ======> at2'
            write(*,'(a)',advance='no') 'index of atom 1 = '
            read(*,*) at1
            write(*,'(a)',advance='no') 'index of atom 2 = '
            read(*,*) at2
    
            T_versor = (system% atom(at2)% xyz - system% atom(at1)% xyz) / sqrt(sum( (system% atom(at2)% xyz-system% atom(at1)% xyz)**2) )
    
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
            T_versor = vector_product(system,at1,at2,at3)
    
    end select

    print*,"" 

end subroutine define_translation_vector
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
!================================================
 subroutine Save_PDB_Trajectories(trj, file_name)
!================================================
type(universe), allocatable, intent(inout) :: trj(:)
character(*)  , optional   , intent(in)    :: file_name

! local variables ...
integer      :: i, j, k, last_config 
character(1) :: YorN

If( present(file_name) ) then
    OPEN(unit = 4, file = file_name          , status = 'unknown', action = 'write')
else 
    OPEN(unit = 4, file = 'frames-output.pdb', status = 'unknown', action = 'write')
end if

last_config = size(trj) -1

do j = 0 , last_config

    write(4,5) 'REMARK' , trj(j)% System_Characteristics
    write(4,1) 'CRYST1' , trj(j)%box(1) , trj(j)%box(2) , trj(j)%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
    write(4,3) 'MODEL'  , j

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

    write(4,'(a)') 'MASTER'
    write(4,'(a)') 'END'

end do

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)              
3 FORMAT(a5,i8)
4 FORMAT(a6,t15,a21)
5 FORMAT(a6,t10,a72)

write(*,'(/a)') ' >>> frames-output.pdb : writing done, press any key <<<'
write(*,'(/a)') "That's all ? (y/n)"
read (*,'(a)') YorN
if( YorN /= "n" ) stop

end subroutine Save_PDB_Trajectories
!
!
!
end module Elastic_Band

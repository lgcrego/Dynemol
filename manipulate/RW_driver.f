module RW_driver

use RW_routines  , only : View_XYZ, View_Yaehmop, Save_POSCAR, Save_MD_Urht
use types_m      , only : universe
use GMX_routines , only : Save_GROMACS 
use Read_Parms   , only : atom 
use diagnosis_m    

public :: WritingRoutines

private

contains

!======================================
subroutine WritingRoutines( structure ) 
!======================================
implicit none
type(universe) , intent(inout) :: structure

! local variables
integer            :: AtNo, N_of_atom_type
character(len=1)   :: Writing_Method

!-----------------------------------------------------------
!       Writing the output file

CALL system( "clear" )

CALL diagnosis( structure )

write(*,'(/a)') '>>>    Writing the output file     <<<'

write(*,'(/a)') '1 : XYZ format'
write(*,'(/a)') '2 : Yaehmop format '
write(*,'(/a)') '3 : POSCAR format '
write(*,'(/a)') '4 : PDB format '
write(*,'(/a)') '5 : Urht format '
write(*,'(/a)') '6 : DONE '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') Writing_Method

select case ( Writing_Method )

    case ('1')
        CALL View_XYZ( structure )

    case ('2')
        CALL View_Yaehmop( structure )
        CALL system("tbind seed")      !<== uses Yaehmop software tbind

    case ('3')
        CALL Save_POSCAR( structure )
       
    case ('4')
        CALL Save_GROMACS( structure )

    case ('5')
        CALL Save_MD_Urht( structure )

    case default
        
end select

! total number of atoms of given type ...
If( Writing_Method /= '3' ) then
    do AtNo = 1 , size(atom)

        N_of_atom_type = count(structure % atom % AtNo == AtNo)
          
        If( N_of_atom_type /= 0 ) Print 121 , atom(AtNo) % symbol , N_of_atom_type

    end do
end If

121   FORMAT(1x,A2,' atoms  = ',I5)
122   FORMAT(1x,A34,I6)
123   FORMAT(1x,A12,I4)

end subroutine WritingRoutines


end module RW_driver

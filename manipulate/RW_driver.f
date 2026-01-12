module RW_driver

use ansi_colors
use RW_routines       , only : View_XYZ, View_Yaehmop, Save_POSCAR, Save_MD_Urht
use types_m           , only : universe
use GMX_routines      , only : Save_GROMACS
use Topology_routines , only : get_topology
use Read_Parms        , only : atom 
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
integer            :: AtNo, N_of_atom_type, file_type
character(len=1)   :: Writing_Method

!-----------------------------------------------------------
!       Writing the output file

CALL system( "clear" )

CALL diagnosis( structure )

do

     write(*,'(/a)') bold // cyan // '>>>    Writing the output file     <<<' // reset
     
     write(*,'(a)') green // ' 1 :' // reset // ' XYZ format'
     write(*,'(a)') green // ' 2 :' // reset // ' Yaehmop format'
     write(*,'(a)') green // ' 3 :' // reset // ' POSCAR format'
     write(*,'(a)') green // ' 4 :' // reset // ' PDB format'
     write(*,'(a)') green // ' 5 :' // reset // ' Generate topology ' // &
                          '(requires CONECT in PDB; obabel -ipdb input.pdb -opdb -O output.pdb)'
     
     write(*,'(a)') bold // orange // ' 0 :' // reset // ' DONE'
     
     write(*,'(/a)', advance='no') bold // yellow // '>>> ' // reset
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
             write(*,'(/a)') bold // cyan // 'Choose topology format:' // reset
             write(*,'(a)') green // ' 1 ' // reset // '= OPLS'
             write(*,'(a)') green // ' 2 ' // reset // '= GAFF / NAMD / Amber'
             
             write(*,'(/a)', advance='no') bold // yellow // '>>> ' // reset
             read (*,*) file_type
             
             call get_topology(structure, file_type)
             
             write(*,'(a)') bold // green // '>>> DONE <<< ' // reset
     
         case default
             exit
             
     end select

     call system("clear")

end do

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

program edview

use types_m
use diagnosis_m
use RW_routines
use GMX_routines
use Read_parms
use EDIT_routines     
use SOLVENT_routines
use FUNCTION_routines
use Elastic_Band
use Trajectory_routines
use Multi_Trajectory_routines
use Statistics_routines
use Crystal_routines
use Occupation
use Aminoacids                 

implicit none

! local variables
type(universe)                  :: structure
type(universe)  , allocatable   :: trajectories(:)
type(molecular)                 :: sol_mol
character(len=1)                :: Reading_Method , Writing_Method
character(len=2)                :: Editing_Method 
character(len=3)                :: resid
integer                         :: i ,j , ioerr , N_of_atoms, N_of_C, N_of_M, N_of_S, AtNo, N_of_atom_type
integer                         :: file_type


CALL Read_Atomic_Mass
CALL Read_EHT_params

!-----------------------------------------------------------
!       Reading the input file
do 
    write(*,'(/a)') '>>>    Read the input file     <<<'

    write(*,'(/a)') '1  : Read "poscar.dat" file (coordinates from vasp-CONTCAR)'
    write(*,'(/a)') '2  : Read "input.xyz" file '
    write(*,'(/a)') '3  : Read "solvent.dat" file '
    write(*,'(/a)') '4  : Read "input.gro" file '
    write(*,'(/a)') '5  : Read "input.pdb" file '
    write(*,'(/a)') '6  : Build up crystal'
    write(*,'(/a)') '7  : Read trajectories '
    write(*,'(/a)') '8  : Read Multiple trajectories '
    write(*,'(/a)') '9  : Read "Occupation.bin" file'
    write(*,'(/a)') 'A  : ad_hoc tuning '
    write(*,'(/a)') '0  : DONE '
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(a)') Reading_Method

    select case ( Reading_Method )

        case ('1')
            CALL Read_from_poscar( structure )

        case ('2')
            CALL Read_from_XYZ( structure )
            write(*,'(/a)',advance='no') 'Type of main Residue : '
            read (*,'(a)') resid
            structure % atom % resid = resid

        case ('3')
            CALL Include_Solvent( structure )

        case ('4')
            CALL Read_GROMACS( structure , file_type="gro")
            CALL gro_2_pdb   ( structure )
            pause '>>>  Saving seed.pdb  <<<'

        case ('5')
            CALL Read_GROMACS( structure , file_type="pdb")

        case ('6')
            CALL Buildup_Crystal( structure )

        case ('7')
            CALL Read_Trajectories( trajectories , structure )

        case ('8')
            CALL Work_Multiple_Trajectories

        case ('9')
            CALL Post_proccess_Occupation

        case ('A','a')
            CALL ad_hoc_tuning( structure )

        case default
            exit

    end select
end do

! in case N_of_atoms was not defined ...
structure%N_of_atoms = size(structure%atom)

! total number of atoms of given type ...
do AtNo = 1 , size(atom)

    N_of_atom_type = count(structure%atom%AtNo == AtNo)
          
    If( N_of_atom_type /= 0 ) Print 121 , atom(AtNo)%symbol , N_of_atom_type

end do

! check list ...
if( any(structure%atom%Symbol == "@@") ) stop ">> special atoms not defined ! <<"

!-----------------------------------------------------------
!       Editing the structure

CALL diagnosis( structure )

do
    write(*,'(/a)') '>>>    Modifying the  input file     <<<'

    write(*,'(/a)') '1  : Translation'
    write(*,'(/a)') '2  : Rotation '
    write(*,'(/a)') '3  : Reflection '
    write(*,'(/a)') '4  : Copy '      
    write(*,'(/a)') '5  : Statistics '
    write(*,'(/a)') '6  : Sort fragments together '
    write(*,'(/a)') '7  : Delete fragment '
    write(*,'(/a)') '8  : Include fragment '
    write(*,'(/a)') '9  : Include Solvent '
    write(*,'(/a)') '10 : Replication '
    write(*,'(/a)') '11 : nonbonding topology '
    write(*,'(/a)') '12 : Images '
    write(*,'(/a)') '13 : Aminoacid Stuff '
    write(*,'(/a)') '14 : ad_hoc tuning '
    write(*,'(/a)') '15 : DONE '
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(a2)') Editing_Method

!    CALL system( "clear" )

    select case ( Editing_Method )

        case ('1')
            CALL Translation( structure )

        case ('2')
            CALL Rotation( structure )

        case ('3')
            CALL Reflection( structure )

        case ('4')
            CALL Copy( structure )

        case ('5')
            CALL Statistics( trajectories )

        case ('6')
            CALL Sort_Fragments( structure )

        case ('7')
            CALL Eliminate_Fragment( structure )

        case ('8')
            CALL Include_Fragment( structure )
            write(*,'(/a)') '>>>  Saving seed.pdb  <<<'

        case ('9')
            CALL Include_Solvent( structure )

        case ('10')
            CALL Replicate( structure )

        case ('11')
            CALL Nonbonding_Topology( structure )

        case ('12')
            write(*,'(/a)') "Format of image files :   vasp-1   /   pdb-2"
            read (*,'(i1)') file_type

            if( allocated(trajectories) ) then
                CALL Pick_Configurations( trajectories , file_type ) 
            else
                CALL Build_Configurations( structure , file_type )
            end if

        case ('13')
            CALL AminoacidStuff( structure )

        case ('14')
            CALL ad_hoc_tuning( structure )

        case default
            exit

    end select
end do

structure%N_of_atoms = size(structure%atom)

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
        CALL system("tbind seed")

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

end program edview
!
!
!
!

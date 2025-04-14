program edview

use types_m
use util_m            , only: read_file_name
use diagnosis_m
use RW_routines
use GMX_routines
use Read_parms
use EDIT_routines     , only: Copiar , Translation , Rotation , Reflection , Eliminate_Fragment , Bring_into_PBCBox , ReGroup , Replicate , Nonbonding_Topology , Include_fragment
use SOLVENT_routines
use FUNCTION_routines
use Elastic_Band
use Trajectory_routines
use Multi_Trajectory_routines
use Statistics_routines
use Crystal_routines
use Occupation
use Aminoacids                 
use Amber_routines
use RW_driver       
use EDT_util_m
use Alignment_routines

!use dtw_routines , only : dtw_stuff

implicit none

! local variables
type(universe)                  :: structure
type(universe)   , allocatable  :: trajectories(:)
character(len=1)                :: Reading_Method , yn
character(len=2)                :: Editing_Method 
character(len=3)                :: resid
character(len=30)               :: f_name
integer                         :: AtNo, N_of_atom_type
integer                         :: file_type


CALL get_environment_vars

CALL Read_Atomic_Mass
CALL Read_EHT_params
!call dtw_stuff

CALL system( "clear" )

!-----------------------------------------------------------
!       Reading the input file
do 
    write(*,'(/a)') '>>>    Read the input file     <<<'

    write(*,'(/a)') '1  : Read "poscar.dat" file (coordinates from vasp-CONTCAR)'
    write(*,'(/a)') '2  : Read XYZ file '
    write(*,'(/a)') '3  : Read "solvent.dat" file '
    write(*,'(/a)') '4  : AMBER stuff'
    write(*,'(/a)') '5  : Read PDB file '
    write(*,'(/a)') '6  : Build up crystal'
    write(*,'(/a)') '7  : Read trajectories '
    write(*,'(/a)') '8  : Read Multiple trajectories '
    write(*,'(/a)') '9  : Read "Occupation.bin" file'
    write(*,'(/a)') 'A  : ad_hoc tuning '
    write(*,'(/a)') '0  : DONE / ADVANCE'
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(a)') Reading_Method

    select case ( Reading_Method )

        case ('1')
            CALL Read_from_poscar( structure )

        case ('2')
            CALL read_file_name( f_name , file_type="xyz" )
            CALL Read_from_XYZ( structure , f_name )
            write(*,'(/a)',advance='no') 'Type of main Residue : '
            read (*,'(a)') resid
            structure % atom % resid = resid

        case ('3')
            CALL Include_Solvent( structure )

        case ('4')
            CALL amber_stuff 

        case ('5')
            CALL read_file_name( f_name , file_type="pdb" )
            CALL Read_GROMACS( structure , f_name , file_type="pdb")

        case ('6')
            CALL Buildup_Crystal( structure )

        case ('7')
            CALL Read_Trajectories( trajectories , structure )

        case ('8')
            CALL Work_Multiple_Trajectories

        case ('9')
            CALL Post_proccess_Occupation

        case ('A','a')
            write(*,'(/a)',advance='no') ">>> on the fly tuning? (y/n) "
            read (*,'(a)') yn

            if( yn /= "n" ) &
            then
                CALL on_the_fly_tuning (structure )
            else
                CALL ad_hoc_tuning( structure )
            end if

        case default
            exit

    end select
end do

If( allocated(structure%atom) ) then

         ! in case N_of_atoms was not defined ...
         structure%N_of_atoms = size(structure%atom)
         
         ! total number of atoms of given type ...
         do AtNo = 1 , size(atom)
         
             N_of_atom_type = count(structure%atom%AtNo == AtNo)
                   
             If( N_of_atom_type /= 0 ) Print 121 , atom(AtNo)%symbol , N_of_atom_type
         
         end do
         
         ! check list ...
         if( any(structure%atom%Symbol == "@@") ) stop ">> special atoms not defined ! <<"
         
         CALL diagnosis( structure )
         
end if

!-----------------------------------------------------------
!       Editing the structure
do
    write(*,'(/a)') '>>>    Modifying the  input file     <<<'

    write(*,'(/a)') '1  : Translation'
    write(*,'(/a)') '2  : Rotation '
    write(*,'(/a)') '3  : Reflection '
    write(*,'(/a)') '4  : Align Structures '
    write(*,'(/a)') '5  : Copy '      
    write(*,'(/a)') '6  : Delete '
    write(*,'(/a)') '7  : Sort fragments together '
    write(*,'(/a)') '8  : Statistics '
    write(*,'(/a)') '9  : Include fragment '
    write(*,'(/a)') '10 : Include Solvent '
    write(*,'(/a)') '11 : Replication '
    write(*,'(/a)') '12 : nonbonding topology '
    write(*,'(/a)') '13 : Images '
    write(*,'(/a)') '14 : Bring Solvent into the Bounding Box '
    write(*,'(/a)') '15 : ad_hoc tuning '
    write(*,'(/a)') '16 : DONE '
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(a2)') Editing_Method

    select case ( Editing_Method )

        case ('1')
            CALL Translation( structure )

        case ('2')
            CALL Rotation( structure )

        case ('3')
            CALL Reflection( structure )

        case ('4')
            CALL Alignment( )

        case ('5')
            CALL Copiar( structure )

        case ('6')
            CALL Eliminate_Fragment( structure )

        case ('7')
            CALL Sort_Fragments( structure )

        case ('8')
            CALL Statistics( trajectories )

        case ('9')
            CALL Include_Fragment( structure )
            write(*,'(/a)') '>>>  Saving seed.pdb  <<<'

        case ('10')
            CALL Include_Solvent( structure )

        case ('11')
            CALL Replicate( structure )

        case ('12')
            CALL Nonbonding_Topology( structure )

        case ('13')
            write(*,'(/a)') "Format of image files :   vasp-1   /   pdb-2"
            read (*,'(i1)') file_type

            if( allocated(trajectories) ) then
                CALL Pick_Configurations( trajectories , file_type ) 
            else
                CALL Build_Configurations( structure , file_type )
            end if

        case ('14')
            CALL Bring_into_PBCBox( structure )

        case ('15')
            CALL ad_hoc_tuning( structure )

        case default
            exit

    end select
end do

structure%N_of_atoms = size(structure%atom)

!-----------------------------------------------------------
!       Writing the output file

        call WritingRoutines( structure ) 
!-----------------------------------------------------------

121   FORMAT(1x,A2,' atoms  = ',I5)
122   FORMAT(1x,A34,I6)
123   FORMAT(1x,A12,I4)

end program edview
!
!
!
!

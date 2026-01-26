program edview

use types_m
use ansi_colors
use util_m            , only: read_file_name
use diagnosis_m
use RW_routines
use GMX_routines
use Read_parms
use EDIT_routines     , only: Copiar, Translation, Rotation, Reflection, Eliminate_Fragment, Bring_into_PBCBox, ReGroup, Replicate, Nonbonding_Topology, Include_fragment
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
use EDT_util_m        , only: on_the_fly_tuning 
use Alignment_routines

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

CALL system( "clear" )

!-----------------------------------------------------------
!       Reading the input file
do
    write(*,'(/a)') bold // cyan // '>>>    Read the input file     <<<' // reset

    write(*,'(a)') green // ' 1  :' // reset // ' Read "poscar.dat" file (VASP CONTCAR)'
    write(*,'(a)') green // ' 2  :' // reset // ' Read XYZ file'
    write(*,'(a)') green // ' 3  :' // reset // ' Read "solvent.dat" file'
    write(*,'(a)') green // ' 4  :' // reset // ' AMBER stuff'
    write(*,'(a)') green // ' 5  :' // reset // ' Read PDB file'
    write(*,'(a)') green // ' 6  :' // reset // ' Build up crystal'
    write(*,'(a)') green // ' 7  :' // reset // ' Read trajectories'
    write(*,'(a)') green // ' 8  :' // reset // ' Read multiple trajectories'
    write(*,'(a)') green // ' 9  :' // reset // ' Read "Occupation.bin" file'
    write(*,'(a)') green // ' A  :' // reset // ' Ad hoc tuning'
    write(*,'(a)') bold // orange // ' 0  :' // reset // ' DONE / ADVANCE'

    write(*,'(/a)', advance='no') bold // yellow // '>>> ' // reset
    read (*,'(a)') Reading_Method

    select case (Reading_Method)

        case ('1')
            call Read_from_poscar(structure)

        case ('2')
            call read_file_name(f_name, file_type="xyz")
            call Read_from_XYZ(structure, f_name)
            write(*,'(/a)', advance='no') yellow // 'Type of main residue: ' // reset
            read (*,'(a)') resid
            structure%atom%resid = resid

        case ('3')
            call Include_Solvent(structure)

        case ('4')
            call amber_stuff

        case ('5')
            call read_file_name(f_name, file_type="pdb")
            call Read_GROMACS(structure, f_name, file_type="pdb")

        case ('6')
            call Buildup_Crystal(structure)

        case ('7')
            call Read_Trajectories(trajectories, structure)

        case ('8')
            call Work_Multiple_Trajectories

        case ('9')
            call Post_proccess_Occupation

        case ('A','a')
            write(*,'(/a)', advance='no') yellow // '>>> On-the-fly tuning? (y/n): ' // reset
            read (*,'(a)') yn

            if (yn /= 'n') then
                call on_the_fly_tuning(structure)
            else
                call ad_hoc_tuning(structure)
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
    write(*,'(/a)') bold // cyan // '>>>    Modifying the input file     <<<' // reset
    
    write(*,'(a)') green // ' 1  :' // reset // ' Translation'
    write(*,'(a)') green // ' 2  :' // reset // ' Rotation'
    write(*,'(a)') green // ' 3  :' // reset // ' Reflection'
    write(*,'(a)') green // ' 4  :' // reset // ' Align structures'
    write(*,'(a)') green // ' 5  :' // reset // ' Copy'
    write(*,'(a)') green // ' 6  :' // reset // ' Delete'
    write(*,'(a)') green // ' 7  :' // reset // ' Sort fragments together'
    write(*,'(a)') green // ' 8  :' // reset // ' Statistics'
    write(*,'(a)') green // ' 9  :' // reset // ' Include fragment'
    write(*,'(a)') green // '10  :' // reset // ' Include solvent'
    write(*,'(a)') green // '11  :' // reset // ' Replication'
    write(*,'(a)') green // '12  :' // reset // ' Nonbonding topology'
    write(*,'(a)') green // '13  :' // reset // ' Build Reaction Coordinate'
    write(*,'(a)') green // '14  :' // reset // ' Bring solvent into the bounding box'
    write(*,'(a)') green // '15  :' // reset // ' Ad hoc tuning'
    write(*,'(a)') bold // orange // ' 0  :' // reset // ' DONE'
    
    write(*,'(/a)', advance='no') bold // yellow // '>>> ' // reset
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
            write(*,'(/a)') blue//"Format of image files:"//reset// &
                            "  "//green//"vasp-1"//reset//"  /  "// green//"pdb-2"//reset
            read (*,'(i1)') file_type
            

            if( allocated(trajectories) ) then
                CALL Pick_Configurations( trajectories , file_type ) 
            else
                CALL Build_Configurations( structure , file_type )
            end if

        case ('14')
            CALL Bring_into_PBCBox( structure )

        case ('15')
            CALL on_the_fly_tuning( structure )

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

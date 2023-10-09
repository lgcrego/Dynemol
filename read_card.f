MODULE card_reading

use type_m
use EH_parms_module
use MM_parms_module

private

   public :: ReadInputCard , ReadInputCard_ADHOC , electron_fragment , hole_fragment

   ! module variables ...
   character(len=3) :: electron_fragment , hole_fragment

contains
!
!
!=============================
 subroutine ReadInputCard
!=============================
implicit none

! local variables ...
integer            :: ioerr , i , n
real*8             :: bottom, top
character(len=1)   :: equal_sign , separator
character(len=20)  :: keyword , keyword1 , keyword2
character(len=20)  :: command , command1 , command2 , command3
character(len=120) :: line
logical            :: flag

call default_values()

open(33, file='card.inpt', status='old', iostat=ioerr, err=10)

!   file error msg ...
10 if( ioerr > 0 ) stop '"card.inpt" file not found; terminating execution'

!=====================================================================================
!  reading  the input CARD ...

read_loop: do 

    read(33,'(A)',iostat=ioerr) line
    if ( ioerr /= 0 ) exit read_loop
    read(line,*,iostat=ioerr) keyword
    keyword = to_upper_case(keyword)
    if( index(keyword,"!") /= 0 ) cycle read_loop 

    select case ( keyword )
    
        case( "DRIVER" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            DRIVER = command
    
        case( "NUCLEAR_MATTER" )
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            nuclear_matter = command

        case( "SURVIVAL" ) 
            Survival = get_logical(line)
    
        case( "DP_MOMENT" ) 
            DP_moment = get_logical(line)

        case( "QMMM" ) 
            QMMM = get_logical(line)

        case( "RND_SEED" ) 
            rnd_seed = get_logical(line)

        case( "DK_OF_MIXING" )
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            DK_of_mixing = command
            select case ( DK_of_mixing ) 
                case ("local" , "LOCAL" , "Local" )
                   DK_of_mixing = "local"
                case ("global" , "GLOBAL" , "Global" )
                   DK_of_mixing = "global"
                case default
                   Print*, " >>> Check your card.inpt <<< : DK_of_mixing = " , DK_of_mixing 
                   stop
            end select

        case( "OPT_PARMS" ) 
            OPT_parms = get_logical(line)

        case( "AD_HOC" ) 
            ad_hoc = get_logical(line)

            ! skip for now, entry to be fully read later from subroutine ad_hoc_tuning ...
            if( ad_hoc == .true. ) then
                do 
                   read(33,'(A)',iostat=ioerr) line
                   read(line,*,iostat=ioerr) command
                   command = to_upper_case(command)
                   If( command(1:6) /= "AD_HOC" ) exit
                   !this prevents double reading in the case of blank lines ...
                   line = "XXXXXXXXXXXXXXXXXXXXXX"
                end do
            end if

!--------------------------------------------------------------------                                                                                                   
!           READING FILE FORMAT
!
        case( "FILE_TYPE" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            file_type = command

        case( "FILE_FORMAT" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            file_format = command

!--------------------------------------------------------------------
!           DIAGNOSTIC & DATA-ANALYSIS & VISUALIZATION flags
!
        case( "HFP_FORCES" ) 
                HFP_Forces = get_logical(line)

        case( "SPECTRUM" )  
                SPECTRUM = get_logical(line)

        case( "ALPHA_TENSOR" ) 
                Alpha_Tensor = get_logical(line)

        case( "GAUSSIANCUBE" ) 
                GaussianCube = get_logical(line)

        case( "GAUSSIANCUBE_STEP" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') GaussianCube_step

        case( "NETCHARGE" ) 
                NetCharge = get_logical(line)

        case( "CH_AND_DP_STEP" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') CH_and_DP_step

        case( "DENSITYMATRIX" ) 
                DensityMatrix = get_logical(line)

        case( "AUTOCORRELATION" ) 
                AutoCorrelation = get_logical(line)

        case( "VDOS_" , "VDOS" , "vdos_" , "vdos") 
                VDOS_ = get_logical(line)

!------------------------------------------------------------------------
!               POTENTIALS
!
        case( "ENVFIELD_" , "ENVFIELD" ) 
                EnvField_ = get_logical(line)

        case( "ENVIRON_TYPE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                Environ_Type = command

        case( "ENVIRON_STEP" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') Environ_step

        case( "COULOMB_" ) 
                Coulomb_ = get_logical(line)

        case( "INDUCED_" ) 
                Induced_ = get_logical(line)

!------------------------------------------------------------------------
!               SAMPLING parameters
!
        case( "FRAME_STEP" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') frame_step

!------------------------------------------------------------------------
!               SECURITY COPY
!
        case( "RESTART" ) 
                restart = get_logical(line)

        case( "STEP_SECURITY" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') step_security

!------------------------------------------------------------------------
!               QDynamics parameters
!
        case( "T_I" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) t_i

        case( "T_F" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) t_f

        case( "N_T" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') n_t

        case( "CT_DUMP_STEP" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') CT_dump_step

        case( "N_PART" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') n_part

        case( "HOLE_STATE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                n = len_trim(command)
                
                read(command(5:n),'(i)') hole_state
                hole_fragment = to_upper_case(command(1:3))

        case( "ELECTRON_STATE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                n = len_trim(command)
                
                read(command(5:n),'(i)') electron_state
                electron_fragment = to_upper_case(command(1:3))

!------------------------------------------------------------------------
!               STRUCTURAL  parameters
!
        case( "NNX" ) 
                read(line,*,iostat=ioerr) keyword1 , equal_sign , command1 , separator , &
                                          keyword2 , equal_sign , command2
                read(command1,'(i)') nnx
                read(command2,'(i)') nny

        case( "PBC" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , separator , &
                                          command1 , command2 , command3
                read(command1,'(i)') PBC(1)
                read(command2,'(i)') PBC(2)
                read(command3,'(i)') PBC(3)

!------------------------------------------------------------------------
!               DOS parameters
!
        case( "SIGMA" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) sigma

        case( "DOS_RANGE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , keyword1 , &
                                          command1 , command2 
                read(command1,*) bottom
                read(command2,*) top

                DOS_range = real_interval( bottom , top )

!------------------------------------------------------------------------
!               SPECTRUM  parameters
!
        case( "OCCUPIED" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , keyword1 , &
                                          command1 , command2 
                read(command1,*) bottom
                read(command2,*) top

                occupied = real_interval( bottom , top )

        case( "EMPTY" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , keyword1 , &
                                          command1 , command2 
                read(command1,*) bottom
                read(command2,*) top

                empty = real_interval( bottom , top )

        case( "N_OF_MOLECULES" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') MM % N_of_molecules

        case( "N_OF_SPECIES" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') MM % N_of_species

        case( "SPECIES(1)" )

                CALL allocate_species( MM % N_of_species )

                backspace(33)

                do n = 1 , MM % N_of_species

                   i=1
                   do 
                      i = i + 1
                      read(33,'(A)',iostat=ioerr) line
                      read(line,*,iostat=ioerr) keyword
                      keyword = to_upper_case(keyword)
                      If( verify("SPECIES()",keyword) == 0 ) exit
                      if( i > 10 ) then
                          CALL warning("halting: check N_of_species in card.inpt")
                          stop
                      end if
                   end do
                   read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
                   species(n) % residue = command

                   read(33,'(A)',iostat=ioerr) line
                   read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
                   read(command,'(i)',iostat=ioerr) species(n) % N_of_molecules 

                   read(33,'(A)',iostat=ioerr) line
                   read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
                   read(command,'(i)',iostat=ioerr) species(n) % N_of_atoms 

                   read(33,'(A)',iostat=ioerr) line
                   line = to_upper_case(line)
                   read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
                   flag = merge( .true. , .false. , any( [".TRUE.","TRUE","T","T_"] == command ) )
                   species(n) % flex = flag
                end do
        
        case( "SELECTIVE_DYNAMICS" )
                Selective_Dynamics = get_logical(line)
        
        case( "THERMOSTAT" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                thermostat = command

        case( "TEMPERATURE" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) temperature
        
        case( "PRESSURE" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) pressure
        
        case( "THERMAL_RELAXATION_T" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) thermal_relaxation_time
        
        case( "CUTOFF_RADIUS" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) , cutoff_radius 
        
        case( "DAMPING_WOLF" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) , damping_Wolf 
        
        case( "DRIVER_MM" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                driver_MM = command

        case( "READ_VELOCITIES" )
                read_velocities = get_logical(line)

        case( "MM_INPUT_FORMAT" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                mm_input_format = command

        case( "MM_LOG_STEP" )
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') MM_log_step
        
        case( "MM_FRAME_STEP" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') MM_frame_step
        
        case( "UNITS_MM" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                Units_mm = command

!--------------------------------------------------------------------
!           Genetic_Alg and CG OPTIMIZATION parameters
!

        case( "POP_SIZE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') Pop_size
        
        case( "N_GENERATIONS" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') N_generations

        case( "POP_RANGE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) , Pop_range 

        case( "SELECTION_BY" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                selection_by = command

        case( "MUTATION_RATE" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,*) , Mutation_rate 

        case( "ADAPTIVE_" , "ADAPTIVE" ) 
                Adaptive_ = get_logical(line)

        case( "MUTATE_CROSS" , "MUTATE-CROSS" ) 
                Mutate_Cross = get_logical(line)

        case( "CG_" , "CG" ) 
                CG_ = get_logical(line)

        case( "TOP_SELECTION" , "TOP-SELECTION" ) 
                read(line,*,iostat=ioerr) keyword , equal_sign , command
                read(command,'(i)') Top_Selection

        case( "PROFILING" ) 
                profiling = get_logical(line)

    end select

    !this prevents double reading in the case of blank lines ...
    keyword = "XXXXXXXXXXXXXXXXXXXXXX"
       
end do read_loop

close(33)

call Checklist_and_Warnings

!--------------------------------------------------------------------

include 'formats.h'

end subroutine ReadInputCard
!
!
!
!==================================================
 subroutine ReadInputCard_ADHOC( structure , atom )
!==================================================
implicit none
type(universe)  , optional , intent(inout) :: structure
type(MM_atomic) , optional , intent(inout) :: atom(:)

! local variables ...
integer            :: ioerr , start , finale , int_value
real*8             :: real_value
character(len=3)   :: label
character(len=20)  :: keyword , EH_MM , feature
character(len=120) :: line
logical            :: done

open(33, file='card.inpt', status='old', iostat=ioerr, err=10)

!   file error msg ...
10 if( ioerr > 0 ) stop '"card.inpt" file not found; terminating execution'

!==============================
!  reading  the input CARD ...

read_loop: do 

    read(33,'(A)',iostat=ioerr) line
    if ( ioerr /= 0 ) exit read_loop
    read(line,*,iostat=ioerr) keyword ! <== keyword = first contiguous string from line
    keyword = to_upper_case(keyword)

    if( index(keyword,"!") /= 0 ) cycle read_loop 

    if( keyword == "AD_HOC" ) then

        do 
           read(33,'(A)',iostat=ioerr) line

           CALL get_line_apart( line , done , EH_MM , feature , start , finale , label , int_value , real_value)

           if( done ) exit

           if( EH_MM == "QM" .AND. present(structure) ) then

                   select case(feature)

                          case( "ATOM" )

                          case( "RESIDUE" )
                              structure%atom(start:finale) % residue = label 
                          case( "NR" )
                              structure%atom(start:finale) % nr = int_value 
                          case( "V_SHIFT" )
                              structure%atom(start:finale) % v_shift = real_value 
                          case( "QMMM" )
                              structure%atom(start:finale) % QMMM = label 

                          end select
                          end if

           if( EH_MM == "MM" .AND. present(atom) ) then

                   select case(feature)

                          case( "RESIDUE" )
                              atom(start:finale) % residue = label 
                          case( "NR" )
                              atom(start:finale) % nr = int_value 
                          end select
                          endif

           !this prevents double reading in the case of blank lines ...
           line = "XXXXXXXXXXXXXXXXXXXXXX"
           
        end do

    end if
    !this prevents double reading in the case of blank lines ...
    keyword = "XXXXXXXXXXXXXXXXXXXXXX"
   
end do read_loop

close(33)

include 'formats.h'

end subroutine ReadInputCard_ADHOC
!
!
!                                                                                                                                                                       
!
!================================
 subroutine allocate_species( N )
!================================
implicit none
integer , intent(in) :: N

! local variables ...
integer :: i

allocate( species ( N ) )

do i = 1 , N
    species(i) % my_species     = 0
    species(i) % N_of_atoms     = 0
    species(i) % N_of_molecules = 0
    species(i) % cm(3)          = 0.0d0
    species(i) % mass           = 0.0d0
    species(i) % flex           = .false.
    species(i) % residue        = "XXX"
    species(i) % nr             = 0
    species(i) % Nbonds         = 0
    species(i) % Nangs          = 0
    species(i) % Ndiheds        = 0
    species(i) % Nharm          = 0
    species(i) % Nbonds14       = 0
    species(i) % NintraIJ       = 0
end do

end subroutine allocate_species
!
!
!
!=======================================
 subroutine get_line_apart(  line , done , EH_MM , feature , start , finale , label , int_value , real_value )
!=======================================
implicit none
character(*)            , intent(in)  :: line
logical                 , intent(out) :: done 
character(*)            , intent(out) :: EH_MM
character(*)            , intent(out) :: feature
integer                 , intent(out) :: start
integer                 , intent(out) :: finale
character(*) , optional , intent(out) :: label
integer      , optional , intent(out) :: int_value
real*8       , optional , intent(out) :: real_value

! local variables ...
integer            :: n , n1 , n2 , ioerr
character(len=10)  :: interval , string
character(len=40)  :: command , command1 , command2 

done = .false.

read(line,*,iostat=ioerr) command , command1 , command2
command = to_upper_case(command)
If( command(1:6) /= "AD_HOC" ) then
    done = .true. 
    return
    endif

EH_MM = command(8:9) 

n1 = index(command,"(") 
n2 = index(command,")") 

feature = command(11:n1-1)

interval = command(n1+1:n2-1)
n = index(interval,":") 
if( n /= 0 ) then
    read(interval(:n-1),'(i)') start
    read(interval(n+1:),'(i)') finale
else  
    read(interval,'(i)') start
    finale = start
end if

if( index(command,"=") /= 0 ) then
    if( len_trim(command) > index(command,"=") ) then
        n1 = index(command,"=") 
        n2 = len_trim(command) 
        string = trim(command(n1+1:))
    elseif( len_trim(command) == index(command,"=") ) then
        string = trim(command1)
    endif
endif

if( index(command1,"=") /= 0 ) then
    if( len_trim(command1) > index(command1,"=") ) then
        n1 = index(command1,"=") 
        n2 = len_trim(command1) 
        string = trim(command1(n1+1:))
    elseif( len_trim(command1) == index(command1,"=") ) then
        string = trim(command2)
    endif
endif

select case(feature) 

       case( "RESIDUE" , "ATOM" , "QMMM" )
       label = trim(string)

       case( "V_SHIFT" )
       read(string,'(f9.5)') real_value

       case( "NR" )
       read(string,'(i)') int_value

end select

end subroutine get_line_apart
!
!
!
!=======================================
 function get_logical(line) result(flag)
!=======================================
implicit none
character(*) , intent(in) :: line

!local variables ...
integer           :: ioerr 
character(len=1)  :: equal_sign 
character(len=20) :: keyword , command
logical           :: flag

read(line,*,iostat=ioerr) keyword , equal_sign , command
command = to_upper_case(command)
flag = merge( .true. , .false. , any( [".TRUE.","TRUE","T","T_"] == command ) ) 

end function get_logical
!
!
!
!======================================
 pure FUNCTION TO_UPPER_CASE ( STRING )
!======================================
 implicit none
 CHARACTER ( LEN = * )              , INTENT(IN)    :: STRING
 CHARACTER ( LEN = LEN ( STRING ) )                 :: TO_UPPER_CASE

! Local parameters ...
INTEGER, PARAMETER :: BIG_A = ICHAR ( "A" ), LITTLE_A = ICHAR ( "a" ), LITTLE_Z = ICHAR ( "z" )

! Local scalars ...
INTEGER :: I, ICHR

! Loop over the characters in the string ...
DO I = 1,LEN ( STRING )

!   Get the ASCII order for the character to be converted ...
    ICHR = ICHAR ( STRING(I:I) )

!   Use the order to change the case of the character ...
    IF ( ( ICHR >= LITTLE_A ) .AND. ( ICHR <= LITTLE_Z ) ) THEN
        TO_UPPER_CASE(I:I) = CHAR ( ICHR + BIG_A - LITTLE_A )
    ELSE
        TO_UPPER_CASE(I:I) = STRING(I:I)
    END IF
END DO

END FUNCTION TO_UPPER_CASE

!
!
!
!
!=================================
 subroutine Checklist_and_Warnings
!=================================
implicit none

! local variables ...
integer           :: N_of_Configs
character (len=7) :: argument , aux
logical           :: dynamic 

! local parameter ...
logical, parameter :: T_ = .true. , F_ = .false.

!--------------------------------------------------------------------
!  hereafter only CHECKLIST and  WARNINGS !!!

select case( DRIVER )

    case( "q_dynamics" , "slice_Cheb" , "slice_AO" , "slice_FSSH" , "slice_CSDM" )
        
        dynamic = T_ 

    case( "avrg_confgs" )
 
        dynamic = ( F_ .OR. Survival )
 
        If( Top_Selection > Pop_size ) stop ">> Top_Selection > Pop_size; execution aborted"

    case( "Genetic_Alg" , "diagnostic" )
 
        dynamic  = F_
        survival = F_
 
        If( Top_Selection > Pop_size ) stop ">> Top_Selection > Pop_size; execution aborted"

    case( "MM_Dynamics" )

        QMMM = F_
        dynamic = F_
        
    case default
        Print*, " >>> Check your driver options <<< :" , driver
        stop

end select

static = .not. dynamic

! verbose is T_ only if ...
verbose = (DRIVER /= "Genetic_Alg") .AND. (DRIVER /= "slice_AO") .AND. (DRIVER /= "slice_Cheb") .AND. (DRIVER /= "slice_FSSH") .AND. (DRIVER /= "slice_CSDM")

If ( DRIVER(1:5)=="slice" .AND. nuclear_matter=="extended_sys" .AND. file_type=="structure" ) then
    Print*," >>> halting: " 
    Print*,"     for fixed nuclei use DRIVER = q_dynamics; " 
    Print*,"     for slice use either file_type=trajectory or nuclear_matter=MDynamics <<<" 
    stop
End If    

If ( QMMM == T_ .AND. HFP_Forces == F_ ) then
    HFP_Forces = T_
elseif ( QMMM == F_ .AND. HFP_Forces == T_ .AND. driver /= "diagnostic" ) then
    CALL warning("MUST turn off HFP_Forces; execution halted, check input card")
    stop 
end if

If ( nuclear_matter == "MDynamics" ) NetCharge = T_

!-------------------------------------------------------------
! get command line argument to preview input data, then stop  ...
CALL GET_COMMAND_ARGUMENT( 1 , argument )
preview = F_
if( COMMAND_ARGUMENT_COUNT() /= 0 ) then
    select case ( argument )

        case( "preview" )
            preview = .true.

        case( "resume" )
            resume = .true.

            if( driver /= "MM_Dynamics" ) then
                CALL warning("warning: resume argument only with DRIVER = MM_Dynamics")
            endif

        case( "spawn" )

            if( driver /= "MM_Dynamics" ) then
                CALL warning("warning: configuration spawning only with DRIVER = MM_Dynamics")
            endif

            spawn = .true.
            if( COMMAND_ARGUMENT_COUNT() < 2 ) then
                CALL warning("halting: must also specify <# of configurations> as dynemol argument, after <spawn>")
            endif

            CALL GET_COMMAND_ARGUMENT( 2 , aux ) 
            read( aux , '(i3)' ) N_of_Configs 
            spawn_step = n_t / N_of_Configs

    end select
end if
!--------------------------------------------------------------

include 'formats.h'

end subroutine Checklist_and_Warnings
!
!
!
!
!=========================
 subroutine default_values
!=========================
implicit none

! local parameters ...

nuclear_matter = "extended_sys"
Survival = .true.
DP_moment = .false.
QMMM = .false.
rnd_seed = .false.
DK_of_mixing = "local"
OPT_parms = .false.
ad_hoc = .true.
file_type = "structure"
file_format = "pdb"
HFP_Forces = .false.
SPECTRUM = .false.
Alpha_Tensor = .false.
GaussianCube = .false.
GaussianCube_step = 5000000
NetCharge = .false.
CH_and_DP_step = 5000000
DensityMatrix = .false.
AutoCorrelation = .false.
VDOS_ = .false.
EnvField_ = .false.
Environ_step = 5
Coulomb_ = .false.
Induced_ = .false.
frame_step = 1
restart = .false.
step_security = 1000
t_i = 0.d0
CT_dump_step = 1
n_part = 2
nnx = 0
nny = 0 
PBC = [0 , 0 , 0]
sigma = 0.05d0
DOS_range = real_interval( -15.d0 , 0.d0 )
Selective_Dynamics = .false.
MM % N_of_species = 1
thermostat = "Microcanonical"
temperature = 300.d0
pressure = 1.d0
thermal_relaxation_time = 0.25d0
pressure_relaxation_time = infty
cutoff_radius = 50.d0
damping_Wolf = 0.001d0
driver_MM = "MM_Dynamics"
read_velocities = .true.
MM_log_step = 50
MM_frame_step = 50
Units_mm = "eV"

end subroutine default_values
!
!
!
!
end module card_reading

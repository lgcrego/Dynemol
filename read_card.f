MODULE card_reading

use type_m
use EH_parms_module
use MM_parms_module

private

   public :: ReadInputCard

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
    
        case( "DRIVER" , "driver" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            DRIVER = command
    
        case( "NUCLEAR_MATTER" , "nuclear_matter" )
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            nuclear_matter = command

        case( "SURVIVAL" , "survival" , "Survival" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            Survival = flag
    
        case( "DP_Moment" , "dp_moment" , "DP_moment" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            DP_moment = flag

        case( "QMMM" , "qmmm" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            QMMM = flag

        case( "OPT_PARMS" , "opt_parms" , "OPT_parms" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            OPT_parms = flag

        case( "AD_HOC" , "ad_hoc" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            ad_hoc = flag

!--------------------------------------------------------------------                                                                                                   
!           READING FILE FORMAT
!
        case( "FILE_TYPE" , "file_type" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            file_type = command

        case( "FILE_FORMAT" , "file_format" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            file_format = command

!--------------------------------------------------------------------
!           DIAGNOSTIC & DATA-ANALYSIS & VISUALIZATION flags
!
    case( "HFP_FORCES" , "hfp_forces" , "HFP_Forces" )  !######
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            HFP_Forces = flag

    case( "SPECTRUM" , "spectrum" )  
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            SPECTRUM = flag

    case( "Alpha_TENSOR" , "alpha_tensor" , "Alpha_Tensor") 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            Alpha_Tensor = flag

    case( "GAUSSIANCUBE" , "gaussiancube" , "GaussianCube" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            GaussianCube = flag

    case( "GAUSSIANCUBE_STEP" , "gaussiancube_step" , "GaussianCube_step" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') GaussianCube_step

    case( "NetCHARGE" , "netcharge" , "NetCharge" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            NetCharge = flag

    case( "CH_AND_DP_STEP" , "ch_and_dp_step" , "CH_and_DP_step" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') CH_and_DP_step

    case( "DENSITYMATRIX" , "densitymatrix" , "DensityMatrix" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            DensityMatrix = flag

    case( "AUTOCORRELATION" , "autocorrelation" , "AutoCorrelation" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            AutoCorrelation = flag

    case( "VDOS_" , "vdos_" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            VDOS_ = flag

!--------------------------------------------------------------------
!           POTENTIALS
!
    case( "ENVFIELD_" , "envfield_" , "EnvField_" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            EnvField_ = flag

    case( "ENVIRON_TYPE" , "environ_type" , "Environ_Type" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            Environ_Type = command

    case( "ENVIRON_STEP" , "environ_step" , "Environ_step" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') Environ_step

    case( "COULOMB_" , "coulomb_" , "Coulomb_" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            Coulomb_ = flag

    case( "INDUCED_" , "induced_" , "Induced_") 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            Induced_ = flag

!--------------------------------------------------------------------
!           SAMPLING parameters
!
    case( "FRAME_STEP" , "frame_step" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') frame_step

!--------------------------------------------------------------------
!           SECURITY COPY
!
    case( "RESTART" , "restart" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            restart = flag

    case( "STEP_SECURITY" , "step_security" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') step_security

!--------------------------------------------------------------------
!           QDynamics parameters
!
    case( "T_i" , "t_i" , "t_I" , "T_I" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,*) t_i

    case( "T_f" , "t_f" , "t_F" , "T_F" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,*) t_f

    case( "N_t" , "n_t" , "n_T" , "N_T" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') n_t

    case( "CT_DUMP_STEP" , "ct_dump_step" , "CT_dump_step" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') CT_dump_step

    case( "N_PART" , "n_part" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') n_part

    case( "HOLE_STATE" , "hole_state" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') hole_state

    case( "ELECTRON_STATE" , "electron_state" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,'(i)') electron_state

!--------------------------------------------------------------------
!           STRUCTURAL  parameters
!
    case( "NNX" , "nnx" ) 
            read(line,*,iostat=ioerr) keyword1 , equal_sign , command1 , separator , &
                                      keyword2 , equal_sign , command2
            read(command1,'(i)') nnx
            read(command2,'(i)') nny

    case( "PBC" , "pbc" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , separator , &
                                      command1 , command2 , command3
            read(command1,'(i)') PBC(1)
            read(command2,'(i)') PBC(2)
            read(command3,'(i)') PBC(3)

!--------------------------------------------------------------------
!           DOS parameters
!
    case( "SIGMA" , "sigma" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            read(command,*) sigma

    case( "DOS_range" , "dos_range" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , keyword1 , &
                                      command1 , command2 
            read(command1,*) bottom
            read(command2,*) top

            DOS_range = real_interval( bottom , top )

!--------------------------------------------------------------------
!           SPECTRUM  parameters
!
    case( "OCCUPIED" , "occupied" ) 
            read(line,*,iostat=ioerr) keyword , equal_sign , keyword1 , &
                                      command1 , command2 
            read(command1,*) bottom
            read(command2,*) top

            occupied = real_interval( bottom , top )

    case( "EMPTY" , "empty" ) 
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
                      CALL system("sed '11i >>> halting: check N_of_species in card.inpt <<<' warning.signal |cat")
                      stop
                  end if
               end do
               read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
               species(n) % residue = command

               read(33,'(A)',iostat=ioerr) line
               read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
               read(command,'(i)') species(n) % N_of_molecules 

               read(33,'(A)',iostat=ioerr) line
               read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
               read(command,'(i)') species(n) % N_of_atoms 

               read(33,'(A)',iostat=ioerr) line
               read(line,*,iostat=ioerr) ( keyword , i=1,4) , command
               flag = merge( .true. , .false. , command == "T_") 
               species(n) % flex = flag
            end do
    
    case( "SELECTIVE_DYNAMICS" )
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            Selective_Dynamics = flag
    
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
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
            read_velocities = flag

    case( "MM_INPUT_FORMAT" )
            read(line,*,iostat=ioerr) keyword , equal_sign , command
            flag = merge( .true. , .false. , command == "T_") 
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
    
end select

!this prevents double reading in the case of blank lines ...
keyword = "XXXXXXXXXXXXXXXXXXXXXX"
   
end do read_loop

!--------------------------------------------------------------------

include 'formats.h'

end subroutine ReadInputCard
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
    species(i) % NIntraLJ       = 0
end do

end subroutine allocate_species
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
!=========================
 subroutine default_values
!=========================
implicit none

! local parameters ...

nuclear_matter = "extended_sys"
Survival = .true.
DP_moment = .false.
QMMM = .false.
OPT_parms = .true.
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
sigma = 0.04d0
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

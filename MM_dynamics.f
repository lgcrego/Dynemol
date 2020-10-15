module MM_dynamics_m

    use constants_m
    use parameters_m        , only : driver , restart , step_security , QMMM , VDOS_ , n_t
    use MM_input            , only : MM_log_step , MM_frame_step , Units_MM , thermostat
    use MD_read_m           , only : atom , MM
    use setup_m             , only : setup, move_to_box_CM, Molecular_CM
    use MD_dump_m           , only : output , cleanup , saving_MM_frame
    use f_inter_m           , only : FORCEINTER
    use f_intra_m           , only : FORCEINTRA
    use QMMM_m              , only : QMMM_FORCE
    use for_force           , only : pot_total
    use Babel_m             , only : QMMM_key
    use Structure_Builder   , only : Unit_Cell
    use Backup_MM_m         , only : Security_Copy_MM , Restart_MM
    use MM_types            , only : debug_MM
    use VV_Parent           , only : VV
    use Berendsen_m         , only : Berendsen
    use Nose_Hoover_m       , only : Nose_Hoover
    use NH_Reversible_m     , only : NH_Reversible
    use NVE_m               , only : NVE
    use VDOS_m              , only : VDOS_init , VDOS_Correlation , VDOS_end

    public :: MolecularMechanics , preprocess_MM , Saving_MM_Backup, MoveToBoxCM

    private

    ! local variables ...
    logical             :: done = .false.
    type(Berendsen)     :: VV_Berendsen
    type(Nose_Hoover)   :: VV_Nose_Hoover
    type(NH_Reversible) :: VV_NH_Reversible
    type(NVE)           :: VV_NVE

contains
!
!
!
!============================================================
 subroutine MolecularMechanics( t_rate , frame , Net_Charge )
!============================================================
implicit none
real*8             , intent(in) :: t_rate
integer            , intent(in) :: frame
real*8  , optional , intent(in) :: Net_Charge(:)


If( VDOS_ ) CALL VDOS_Correlation( frame )

! selectc thermostat ...
select case( thermostat )

    case( "Berendsen" )
         CALL VelocityVerlet( VV_Berendsen , t_rate , frame , Net_Charge )

    case( "Nose_Hoover" )
         CALL VelocityVerlet( VV_Nose_Hoover , t_rate , frame , Net_Charge )

    case( "NH_Reversible" )
         CALL VelocityVerlet( VV_NH_Reversible , t_rate , frame , Net_Charge )

    case( "Microcanonical" )
         CALL VelocityVerlet( VV_NVE , t_rate , frame , Net_Charge )

    case default
        Print*, " >>> Check your thermostat options <<< :" , thermostat
        stop

end select

If ( VDOS_ ) then
    If ( DRIVER == "MM_Dynamics" ) then
        if (frame == n_t)  call VDOS_end
    else
        if (frame == n_t - 1) then  ! in this case, the dyn. loop ends at n_t - 1
            call VDOS_Correlation( n_t )
            call VDOS_end
        end if
    end if
end If 

end subroutine MolecularMechanics
!
!
!
!==============================================================
subroutine VelocityVerlet( this , t_rate , frame , Net_Charge )
! nuclear velocities in units of m/s in atom%vel
!==============================================================
implicit none
class(VV)          , intent(inout) :: this
real*8             , intent(in)    :: t_rate
integer            , intent(in)    :: frame
real*8  , optional , intent(in)    :: Net_Charge(:)

! local variables ...
real*8  :: dt , Temperature , pressure , density , Kinetic
integer :: i

! time units are PICOseconds in EHT - seconds in MM ; converts picosecond to second ...
dt = t_rate * pico_2_sec

atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge

! Molecular dynamics ...
CALL this % VV1( dt )

CALL move_to_box_CM

CALL Molecular_CM

CALL ForceInter

CALL ForceIntra

! QMMM coupling ...
if( QMMM ) CALL QMMM_FORCE( Net_Charge )

CALL this% VV2( dt )

if( mod(frame,MM_frame_step) == 0 ) CALL Saving_MM_frame( frame , dt )

Unit_Cell% MD_Kin = this% Kinetic * kJmol_2_eV * MM% N_of_Molecules
Unit_Cell% MD_Pot = Pot_total     * kJmol_2_eV * MM% N_of_Molecules

if( mod(frame,MM_log_step) == 0   ) then 

    Temperature = this % Temperature
    Kinetic     = this % Kinetic 
    Pressure    = this % pressure 
    Density     = this % density

    CALL output( Temperature , frame , dt )

    select case (Units_MM)

        case( "eV" )    
        write(*,10) frame, Temperature, Unit_Cell% MD_Kin, Unit_Cell% MD_Pot, Unit_Cell% MD_Kin + Unit_Cell% MD_Pot 

        case( "kj-mol" )
        write(*,10) frame, Temperature, Unit_Cell% MD_Kin*eV_2_kJmol, Unit_Cell% MD_Pot*eV_2_kJmol, (Unit_Cell% MD_Kin + Unit_Cell% MD_Pot)*eV_2_kJmol

        case default
        write(*,10) frame , Temperature , density , pressure , Kinetic , pot_total , Kinetic + pot_total

    end select

end if

! pass nuclear configuration to QM ...
forall(i=1:size(atom)) Unit_Cell % coord(i,:) = atom( QMMM_key(i) ) % xyz(:)

! saving backup stuff ...
if( driver == "MM_Dynamics" ) CALL Saving_MM_Backup( frame , instance = "from_MM" )

10 format(I7,4F15.5)

end subroutine VelocityVerlet
!
!
!
!==================================================
subroutine preprocess_MM( frame_init , Net_Charge )
!==================================================
implicit none
integer , optional , intent(inout) :: frame_init
real*8  , optional , intent(in)    :: Net_Charge(:)

! local variables ...
integer         :: frame , i 

! select thermostat ...
select case( thermostat )

    case( "Berendsen" )
         VV_Berendsen = Berendsen()

    case( "Nose_Hoover" )
         VV_Nose_Hoover = Nose_Hoover()

    case( "NH_Reversible" )
         VV_NH_Reversible = NH_Reversible()

    case( "Microcanonical" )
         VV_NVE = NVE()

end select

! setting up ...
atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge

CALL Setup
If( .not. done ) CALL move_to_box_CM
CALL Molecular_CM

if( restart ) then

    CALL Restart_MM( MM , atom , frame )

    if( present(frame_init) ) frame_init = frame + 1

    ! pass nuclear configuration for QM ...
    forall(i=1:size(atom)) Unit_Cell % coord(i,:) = atom( QMMM_key(i) ) % xyz(:)

else

    CALL Cleanup

    CALL ForceInter
    CALL ForceIntra

    ! QMMM coupling ...
    if( QMMM ) CALL QMMM_FORCE( Net_Charge )

    ! saving the first frame ==> frame 0 = input ...
    CALL Saving_MM_frame( frame=0 , dt=D_zero )

    if( present(frame_init) ) frame_init = 1

    CALL output( Ttrans=D_zero, frame=0, dt=D_zero )

endif

If( VDOS_ ) CALL VDOS_init

end subroutine preprocess_MM
!
!
!
!========================
subroutine MoveToBoxCM( )
!========================
implicit none

! local variables ...
integer :: i 

! just pass nuclear configuration to QM routines and leave ...

CALL move_to_box_CM

forall(i=1:size(atom)) Unit_Cell % coord(i,:) = atom( QMMM_key(i) ) % xyz(:)

done = .true.

end subroutine MoveToBoxCM
!
!
!
!==============================================
subroutine Saving_MM_Backup( frame , instance )
!==============================================
implicit none
integer          , intent(in) :: frame
character(len=*) , intent(in) :: instance

select case( instance )

        case( "from_MM" )

            if( mod(frame,step_security) == 0 ) CALL Security_Copy_MM( MM , atom , frame )

        case( "from_QM" )

            CALL Security_Copy_MM( MM , atom , frame )

end select

end subroutine Saving_MM_Backup
!
!
!
end module MM_dynamics_m

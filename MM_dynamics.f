module MM_dynamics_m

    use constants_m
    use parameters_m        , only : driver , restart , step_security , QMMM 
    use MM_input            , only : MM_log_step , MM_frame_step , Units_MM
    use MD_read_m           , only : atom , MM
    use setup_m             , only : setup, move_to_box_CM, Molecular_CM
    use MD_dump_m           , only : output , cleanup , saving_MM_frame
    use f_inter_m           , only : FORCEINTER
    use f_intra_m           , only : FORCEINTRA
    use QMMM_m              , only : QMMM_FORCE
    use for_force           , only : pot_total
    use verlet_m            , only : VV1 , VV2 , Summat , PRESS_Boundary
    use Babel_m             , only : QMMM_key
    use Structure_Builder   , only : Unit_Cell
    use Backup_MM_m         , only : Security_Copy_MM , Restart_MM
    use MM_types            , only : debug_MM

    public :: MolecularMechanics , preprocess_MM , Saving_MM_Backup, MoveToBoxCM

    private

    !local variables ...
    logical :: done = .false.

contains
!
!
!
!===========================================================
subroutine MolecularMechanics( t_rate , frame , Net_Charge )
!===========================================================
implicit none
real*8             , intent(in) :: t_rate
integer            , intent(in) :: frame
real*8  , optional , intent(in) :: Net_Charge(:)

! local variables ...
real*8  :: dt , Ttrans , pressure , density
real*8  :: kinetic
integer :: i

! time units are PICOseconds in EHT - seconds in MM ; converts picosencond to second ...
dt = t_rate * pico_2_sec

atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge

! Molecuar dynamic ...
CALL VV1( dt )

CALL Molecular_CM

CALL ForceInter

CALL ForceIntra

! QMMM coupling ...
if( QMMM ) CALL QMMM_FORCE( Net_Charge )

CALL VV2 ( Ttrans , kinetic , dt )

CALL Summat( density ) 

CALL Press_Boundary( pressure , dt )

if( mod(frame,MM_frame_step) == 0 ) CALL Saving_MM_frame( frame , dt )

if( mod(frame,MM_log_step) == 0   ) then 

    CALL output( Ttrans , frame , dt )

    select case (Units_MM)

        case( "eV" )    
        write(*,'(I7,6F15.5)') frame , Ttrans , density , pressure , kinetic*kJmol_2_eV , pot_total*kJmol_2_eV , (kinetic+pot_total)*kJmol_2_eV

        case default
        write(*,'(I7,6F15.5)') frame , Ttrans , density , pressure , kinetic , pot_total , kinetic + pot_total

    end select

end if

! pass nuclear configuration to QM ...
forall(i=1:size(atom)) Unit_Cell % coord(i,:) = atom( QMMM_key(i) ) % xyz(:)

! saving backup stuff ...
if( driver == "MM_Dynamics" ) CALL Saving_MM_Backup( frame , instance = "from_MM" )

end subroutine MolecularMechanics
!
!
!
!==================================================
subroutine preprocess_MM( frame_init , Net_Charge )
!==================================================
implicit none
integer , optional , intent(inout) :: frame_init
real*8  , optional , intent(in)    :: Net_Charge(:)

!local variables ...
integer :: frame , i 

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

endif

end subroutine preprocess_MM
!
!
!
!========================
subroutine MoveToBoxCM( )
!========================
implicit none

!local variables ...
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

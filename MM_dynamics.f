module MM_dynamics_m

    use constants_m
    use parameters_m        , only : restart , QMMM
    use MM_input            , only : MM_log_step , MM_frame_step
    use MD_read_m           , only : atom , Reading , restrt , MM
    use setup_m             , only : setup, move_to_box_CM, Molecular_CM
    use MD_dump_m           , only : output , cleanup , saving_MM_frame
    use f_inter_m           , only : FORCEINTER
    use f_intra_m           , only : FORCEINTRA
    use QMMM_m              , only : QMMM_FORCE
    use for_force           , only : pot
    use verlet_m            , only : VV1 , VV2 , Summat , PRESS_Boundary
    use Babel_m             , only : QMMM_key
    use Structure_Builder   , only : Unit_Cell
    use Data_Output         , only : Net_Charge

    public :: MolecularMechanics , preprocess_MM

    private

contains
!
!
!
!==============================================
subroutine MolecularMechanics( t_rate , frame )
!==============================================
implicit none
real*8  , intent(in)    :: t_rate
integer , intent(in)    :: frame

! local variables ...
real*8  :: dt , Ttrans , pressure , density
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
CALL VV2 ( Ttrans , dt )

CALL Summat( density ) 
CALL Press_Boundary( pressure , dt )

if ( mod(frame,MM_frame_step) == 0 ) CALL Saving_MM_frame( frame , dt )

if ( mod(frame,MM_log_step) == 0 ) write (*,'(I7,4F15.5)') frame , Ttrans , pressure , density

if ( mod(frame,MM_log_step) == 0 ) CALL output( Ttrans , frame , dt )

! pass nuclear configuration to QM ...
forall(i=1:size(atom)) Unit_Cell % coord(i,:) = atom( QMMM_key(i) ) % xyz(:)

end subroutine MolecularMechanics
!
!
!
!=======================
subroutine preprocess_MM
!=======================
implicit none

integer :: i

CALL Reading

atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge

CALL Setup
CALL move_to_box_CM
CALL Molecular_CM

if( restart ) then

    CALL Restrt
    CALL move_to_box_CM

else

    CALL Cleanup
    
    CALL ForceInter
    CALL ForceIntra

    ! QMMM coupling ...
    if( QMMM ) CALL QMMM_FORCE( Net_Charge )

endif

! saving the first frame ==> frame 0 = input ...
CALL Saving_MM_frame( frame=0 , dt=D_zero )

end subroutine preprocess_MM
!
!
!
end module MM_dynamics_m

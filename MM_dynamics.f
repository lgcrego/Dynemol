module MM_dynamics_m

    use constants_m
    use parameters_m        , only : restart , QMMM
    use MD_read_m           , only : atom , Reading , restrt , MM , species , molecule , FF
    use setup_m             , only : setup, cmzero, cmass
    use MD_dump_m           , only : output , cleanup , saving
    use F_inter_m           , only : FORCEINTER
    use F_intra_m           , only : FORCEINTRA
    use for_force           , only : pot
    use MD_dynamics_m       , only : VV1 , VV2 , SUMMAT , PRESS_Boundary
    use MM_types            , only : MM_system , MM_molecular , MM_atomic
    use Babel_m             , only : QMMM_key
    use Structure_Builder   , only : Unit_Cell
    use Data_Output         , only : Net_Charge

    public :: MolecularDynamics

    private

    ! module variables ...
    logical , save  :: done = .false.

contains
!
!
!
!=============================================
subroutine MolecularDynamics( t_rate , frame )
!=============================================
implicit none
real*8  , intent(in)    :: t_rate
integer , intent(in)    :: frame

! local variables ...
real*8  :: Ttrans , pressure , density , dt
integer :: i


! time units are PICOseconds in EHT - seconds in MM ; converts picosecond to second ... 
dt = t_rate * pico_2_sec 

if( .NOT. done ) CALL preprocess_DM

If( QMMM ) atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge + Net_Charge(:)

! Molecuar dynamic ...
CALL VV1( dt )
CALL cmass
CALL ForceInter
CALL ForceIntra
CALL VV2 ( Ttrans , frame , dt )

CALL Summat( density ) 
!CALL Press_Boundary( pressure , dt )

if (mod(frame,1) == 0) CALL Saving( frame , dt )
if (mod(frame,1) == 0) write (*,'(I7,4F15.5)') frame, Ttrans, pot*mol*1.d-6 / dfloat(MM % N_of_molecules), pressure, density

CALL output( Ttrans , dt )

! pass nuclear configuration to QM ...
forall(i=1:size(atom)) Unit_Cell % coord(i,:) = atom( QMMM_key(i) ) % xyz(:)

end subroutine MolecularDynamics
!
!
!
!=======================
subroutine preprocess_DM
!=======================
implicit none

integer :: i

CALL Reading

If( QMMM ) atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge + Net_Charge(:)

CALL Setup
CALL cmzero
CALL cmass

if( restart ) then
    CALL Restrt
    CALL cmzero
else
    CALL Cleanup
    CALL ForceInter
    CALL ForceIntra
endif

done = .true.

end subroutine preprocess_DM
!
!
!
end module MM_dynamics_m

#include "GPU.h"

! Subroutine for computing time evolution adiabatic on the AO
module CSDM_adiabatic_m

    use type_m
    use constants_m
    use blas95
    use parameters_m            , only: t_i , n_t , t_f , n_part ,       &
                                        frame_step , nuclear_matter ,    &
                                        EnvField_ , DP_Moment ,          &
                                        Induced_ , QMMM , restart ,      &
                                        GaussianCube , static ,          &
                                        GaussianCube_step , preview ,    &
                                        hole_state , electron_state ,    &
                                        DensityMatrix, AutoCorrelation,  &
                                        CT_dump_step, Environ_step,      &
                                        step_security
    use Babel_m                 , only: Coords_from_Universe, trj, MD_dt                            
    use Allocation_m            , only: Allocate_UnitCell ,              &
                                        DeAllocate_UnitCell ,            &
                                        DeAllocate_Structures ,          &
                                        Allocate_Brackets                
    use Structure_Builder       , only: Unit_Cell ,                      &
                                        Extended_Cell ,                  &
                                        Generate_Structure ,             &
                                        Basis_Builder
    use FMO_m                   , only: FMO_analysis ,                   &
                                        orbital , eh_tag                 
    use DP_main_m               , only: Dipole_Matrix ,                  &
                                        Dipole_Moment
    use TD_Dipole_m             , only: wavepacket_DP                                        
    use Polarizability_m        , only: Build_Induced_DP
    use Solvated_M              , only: Prepare_Solvated_System 
    use QCModel_Huckel          , only: EigenSystem , S_root_inv 
    use Schroedinger_m          , only: DeAllocate_QDyn
    use Psi_Squared_Cube_Format , only: Gaussian_Cube_Format
    use Data_Output             , only: Populations ,                    &
                                        Net_Charge                       
    use Backup_m                , only: Security_Copy ,                  &
                                        Restart_state ,                  &
                                        Restart_Sys                      
    use MM_dynamics_m           , only: MolecularMechanics ,             &
                                        preprocess_MM , MoveToBoxCM
    use Auto_Correlation_m      , only: MO_Occupation
    use Dielectric_Potential    , only: Environment_SetUp
    use F_intra_m               , only: BcastQMArgs
    use Ehrenfest_CSDM          , only: PST 
    use decoherence_m           , only: apply_decoherence , DecoherenceForce
                                        

    public :: CSDM_adiabatic

    private

    ! module variables ...
    type(STO_basis) , allocatable , dimension(:)   :: ExCell_basis
    Complex*16      , allocatable , dimension(:,:) :: MO_TDSE_bra , MO_TDSE_ket , DUAL_TDSE_bra , DUAL_TDSE_ket ! <== TDSE wvpckt
    Complex*16      , allocatable , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra   ! <== CSDM wvpckt
    Complex*16      , allocatable , dimension(:)   :: phase
    real*8          , allocatable , dimension(:)   :: Net_Charge_MM
    type(R_eigen)   :: UNI , el_FMO , hl_FMO
    real*8          :: t
    integer         :: it , mm , nn

contains
!
!
!
!============================================
 subroutine CSDM_adiabatic( Qdyn , final_it )
!============================================
implicit none
type(f_time)    , intent(out)   :: QDyn
integer         , intent(out)   :: final_it

! local variables ...
integer        :: frame , frame_init , frame_final , frame_restart
real*8         :: t_rate 
type(universe) :: Solvated_System
logical        :: triggered

it = 1
t  = t_i

!--------------------------------------------------------------------------------
! time slicing H(t) : Quantum Dynamics & All that Jazz ...

! time is PICOseconds in EHT & seconds in MM ...
If( nuclear_matter == "MDynamics" ) then
    t_rate      = t_f / float(n_t)
    frame_final = n_t
else
    t_rate      = merge( t_f / float(n_t) , MD_dt * frame_step , MD_dt == epsilon(1.d0) )
    frame_final = size(trj)
end If

If( restart ) then
    CALL Restart_stuff( QDyn , frame_restart )
else
    CALL Preprocess( QDyn )
    triggered = yes
end If

frame_init = merge( frame_restart+1 , frame_step+1 , restart )

do frame = frame_init , frame_final , frame_step

    t = t + t_rate 

    If( (it >= n_t) .OR. (t >= t_f) ) exit    

    it = it + 1
if(it == 40000) stop
    ! calculate for use in MM ...
    If( QMMM ) then
        Net_Charge_MM = Net_Charge
        CALL BcastQMArgs( Extended_Cell, ExCell_basis,  MO_bra, MO_ket, MO_TDSE_bra, MO_TDSE_ket, UNI )
    end If

    ! propagate t -> (t + t_rate) with UNI%erg(t) ...
    CALL U_ad(t_rate) ! <== adiabatic component of the propagation ; 1 of 2 ... 

    ! DUAL representation for efficient calculation of survival probabilities ...
    ! also, intermediary step for NAD-TDSE ...
    CALL DUAL_wvpckts
 
    ! save populations(t + t_rate)  and  update Net_Charge ...
    QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t )

    if( mod(it,CT_dump_step) == 0 ) CALL dump_Qdyn( Qdyn )

    If( GaussianCube .AND. mod(frame,GaussianCube_step) < frame_step ) CALL  Send_to_GaussianCube( frame )

    If( DP_Moment ) CALL DP_stuff( "DP_moment" )

    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    ! build new UNI(t + t_rate) ...
    !============================================================================

    select case ( nuclear_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL DeAllocate_UnitCell ( Unit_Cell )
            CALL Coords_from_Universe( Unit_Cell , Solvated_System )

        case( "extended_sys" )

            CALL DeAllocate_UnitCell ( Unit_Cell )
            CALL Coords_from_Universe( Unit_Cell , trj(frame) )

        case( "MDynamics" )

            ! MM preprocess ...
            if( frame == frame_step+1 ) CALL preprocess_MM( Net_Charge = Net_Charge_MM )   

            ! IF QM_erg < 0 => turn off QMMM ; IF QM_erg > 0 => turn on QMMM ...
            QMMM = (triggered == yes)
print*, QMMM

            ! MM precedes QM ; notice calling with frame -1 ...
            CALL MolecularMechanics( t_rate , frame - 1 , Net_Charge = Net_Charge_MM )   

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
            stop

    end select

    CALL Generate_Structure( frame )
    
    CALL Basis_Builder( Extended_Cell , ExCell_basis )

    If( EnvField_ ) CALL DP_stuff( "EnvField" )

    If( Induced_ )  CALL DP_stuff( "Induced_DP" )

    Deallocate ( UNI%R , UNI%L , UNI%erg )

    CALL EigenSystem( Extended_Cell , ExCell_basis , UNI , it )

    CALL U_nad() ! <== NON-adiabatic component of the propagation ; 2 of 2 ... 

    if( mod(frame,step_security) == 0 ) CALL Security_Copy( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame )

    If( DensityMatrix ) then
        If( n_part == 1 ) CALL MO_Occupation( t, MO_bra, MO_ket, UNI )
        If( n_part == 2 ) CALL MO_Occupation( t, MO_bra, MO_ket, UNI, UNI )
    End If

    CALL QMMM_erg_status( frame , t_rate , triggered )

!    CALL apply_decoherence( MO_bra , MO_ket , UNI%erg , PST , t_rate )

    CALL DecoherenceForce( Extended_Cell , MO_bra , MO_ket , UNI%erg , PST )

    CALL apply_decoherence( MO_bra , MO_ket , UNI%erg , PST , t_rate )
end do

deallocate( MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )
deallocate( MO_TDSE_bra , MO_TDSE_ket , DUAL_TDSE_bra , DUAL_TDSE_ket )

final_it = it

include 'formats.h'

end subroutine CSDM_adiabatic
!
!
!
!=============================
 subroutine Preprocess( QDyn )
!=============================
implicit none
type(f_time)    , intent(out)   :: QDyn

! local variables
integer         :: hole_save , n 
type(universe)  :: Solvated_System

! preprocessing stuff .....................................................

CALL DeAllocate_QDyn( QDyn , flag="alloc" )

select case ( nuclear_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , 1 )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(1) )

    case( "MDynamics" )
    
        CALL MoveToBoxCM
      
    case default

        Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
        stop

end select

CALL Generate_Structure( 1 )

CALL Basis_Builder( Extended_Cell , ExCell_basis )

mm = size(ExCell_basis) ; nn = n_part

CALL Allocate_Brackets( mm , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

If( Induced_ ) CALL Build_Induced_DP( basis = ExCell_basis , instance = "allocate" )

If( EnvField_ ) then

    hole_save  = hole_state
    hole_state = 0
    static     = .true. 

    ! Environ potential in the static GS configuration ...
    CALL Environment_SetUp  ( Extended_Cell )

    hole_state = hole_save
    static     = .false.

    CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

end If

CALL EigenSystem( Extended_Cell , ExCell_basis , UNI , it )

! building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
do n = 1 , n_part                         
    select case( eh_tag(n) )

        case( "el" )

            CALL FMO_analysis( Extended_Cell , ExCell_basis , UNI , el_FMO , instance="E" )

            MO_bra( : , n ) = el_FMO%L( orbital(n) , : )    
            MO_ket( : , n ) = el_FMO%R( : , orbital(n) )   

            Print 591, orbital(n) , el_FMO%erg(orbital(n))
       
        case( "hl" )

            CALL FMO_analysis( Extended_Cell , ExCell_basis , UNI , hl_FMO , instance="H" )

            MO_bra( : , n ) = hl_FMO%L( orbital(n) , : )    
            MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )   

            Print 592, orbital(n) , hl_FMO%erg(orbital(n))
            If( (orbital(n) > hl_FMO%Fermi_State) ) print*,'>>> warning: hole state above the Fermi level <<<'

        end select
end do

! stop here to preview and check input and system info ...
If( preview ) stop

UNI% Fermi_state = Extended_Cell% N_of_Electrons/TWO + mod( Extended_Cell% N_of_Electrons , 2 )

! setup coherent component of the CSDM wavepacket ... 
allocate( MO_TDSE_bra , source=MO_bra )
allocate( MO_TDSE_ket , source=MO_ket )

allocate( Dual_TDSE_bra , mold=Dual_bra )
allocate( Dual_TDSE_ket , mold=Dual_ket )

! DUAL representation for efficient calculation of survival probabilities ...
CALL DUAL_wvpckts
 
! save populations ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t_i )
CALL dump_Qdyn( Qdyn )

If( GaussianCube ) CALL Send_to_GaussianCube( it )

If( DP_Moment    ) CALL DP_stuff( "DP_matrix" )

If( DP_Moment    ) CALL DP_stuff( "DP_moment" )

If( DensityMatrix ) then
    If( n_part == 1 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI )
    If( n_part == 2 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI, UNI )
End If

If( Induced_ ) CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

CALL BcastQMArgs( mm , Extended_Cell%atoms )

Unit_Cell% QM_erg = update_QM_erg()
!..........................................................................

include 'formats.h'

deallocate( el_FMO%L , el_FMO%R , el_FMO%erg )
deallocate( hl_FMO%L , hl_FMO%R , hl_FMO%erg )

end subroutine Preprocess
!
!
!
!
!=========================
 subroutine U_ad( t_rate )
!=========================
implicit none
real*8  , intent(in) :: t_rate 

! local variables ...
integer :: j

phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

! adiabatic component of the propagation ...
forall( j=1:n_part )   
    MO_TDSE_bra(:,j) = merge( conjg(phase(:)) * MO_TDSE_bra(:,j) , C_zero , eh_tag(j) /= "XX" )
    MO_TDSE_ket(:,j) = merge(       phase(:)  * MO_TDSE_ket(:,j) , C_zero , eh_tag(j) /= "XX" )

    MO_bra(:,j)      = merge( conjg(phase(:)) * MO_bra(:,j)      , C_zero , eh_tag(j) /= "XX" )
    MO_ket(:,j)      = merge(       phase(:)  * MO_ket(:,j)      , C_zero , eh_tag(j) /= "XX" )
end forall

end subroutine U_ad
!
!==================
 subroutine U_nad()
!==================
implicit none

! NON-adiabatic component of the propagation ...
! project back to MO_basis with UNI(t + t_rate)
CALL dzgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%R , mm , Dual_TDSE_bra , mm , C_zero , MO_TDSE_bra , mm )
CALL dzgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%L , mm , Dual_TDSE_ket , mm , C_zero , MO_TDSE_ket , mm )

CALL dzgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%R , mm , Dual_bra , mm , C_zero , MO_bra , mm )
CALL dzgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%L , mm , Dual_ket , mm , C_zero , MO_ket , mm )

end subroutine U_nad
!
!
!
!=======================
 subroutine DUAL_wvpckts
!=======================
implicit none

CALL dzgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%R , mm , MO_TDSE_ket , mm , C_zero , DUAL_TDSE_ket , mm )
CALL dzgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_TDSE_bra , mm , C_zero , DUAL_TDSE_bra , mm )

CALL dzgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%R , mm , MO_ket , mm , C_zero , DUAL_ket , mm )
CALL dzgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_bra , mm , C_zero , DUAL_bra , mm )

end subroutine DUAL_wvpckts
!
!
!
!
!========================================
 subroutine Send_to_GaussianCube( frame )
!========================================
implicit none
integer , intent(in) :: frame

! local variables ...
integer :: n

! LOCAL representation for film STO production ...

! coefs of <k(t)| in AO basis 
AO_bra = DUAL_bra

! coefs of |k(t)> in AO basis 
CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

do n = 1 , n_part
    if( eh_tag(n) == "XX" ) cycle
    CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , frame ,t , eh_tag(n) )
end do

!----------------------------------------------------------

end subroutine Send_to_GaussianCube
!
!
!
!
!===============================
 subroutine DP_stuff( instance )
!===============================
implicit none
character(*)  , intent(in)    :: instance

!local variables ...
integer :: i
real*8  :: Total_DP(3)

!----------------------------------------------------------
!       LOCAL representation for DP calculation ...

! coefs of <k(t)| in AO basis 
AO_bra = DUAL_bra

! coefs of |k(t)> in AO basis 
CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

select case( instance )

    case( "DP_matrix" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

    case( "EnvField" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        ! wavepacket component of the dipole vector ...
        ! decide what to do with this ############ 
        !CALL wavepacket_DP( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_ket )

        If( mod(it-1,Environ_step) == 0 ) CALL Environment_SetUp( Extended_Cell )

    case( "DP_moment" )

        CALL Dipole_Moment( Extended_Cell , ExCell_basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

        If( t == t_i ) then
            open( unit = 51 , file = "dyn.trunk/dipole_dyn.dat" , status = "replace" )
        else
            open( unit = 51 , file = "dyn.trunk/dipole_dyn.dat" , status = "unknown", action = "write" , position = "append" )
        end If
        write(51,'(5F10.5)') t , (Total_DP(i) , i=1,3) , sqrt( sum(Total_DP*Total_DP) )
        close(51)

    case( "Induced_DP" ) 

        If( .NOT. EnvField_ ) CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

end select

!----------------------------------------------------------

end subroutine DP_stuff
!
!
!
!============================
 subroutine dump_Qdyn( Qdyn )
!============================
implicit none
type(f_time)    , intent(in) :: QDyn

! local variables ...
integer    :: nf , n
complex*16 :: wp_energy(n_part)

do n = 1 , n_part

    if( eh_tag(n) == "XX" ) cycle

    wp_energy(n) = sum(MO_bra(:,n)*UNI%erg(:)*MO_ket(:,n)) 

    If( it == 1 ) then

        open( unit = 52 , file = "dyn.trunk/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
        write(52,11) "#" ,( nf+1 , nf=0,size(QDyn%fragments)+1 )  ! <== numbered columns for your eyes only ...
        write(52,12) "#" , QDyn%fragments , "total"

        open( unit = 53 , file = "dyn.trunk/"//eh_tag(n)//"_wp_energy.dat" , status = "replace" , action = "write" , position = "append" )

    else

        open( unit = 52 , file = "dyn.trunk/"//eh_tag(n)//"_survival.dat"  , status = "unknown", action = "write" , position = "append" )
        open( unit = 53 , file = "dyn.trunk/"//eh_tag(n)//"_wp_energy.dat" , status = "unknown", action = "write" , position = "append" )

    end If
 
    ! dumps el-&-hl populations ...
    write(52,13) ( QDyn%dyn(it,nf,n) , nf=0,size(QDyn%fragments)+1 ) 

    ! dumps el-&-hl wavepachet energies ...
    write(53,14) QDyn%dyn(it,0,n) , real( wp_energy(n) ) , dimag( wp_energy(n) )

    close(52)
    close(53)

end do

11 FORMAT(A,I9,14I10)
12 FORMAT(/15A10)
13 FORMAT(F11.6,14F10.5)
14 FORMAT(3F12.6)

end subroutine dump_Qdyn
!
!
!
!========================================================
 subroutine QMMM_erg_status( frame , t_rate , triggered )
!========================================================
use MM_input , only : Units_MM , MM_log_step
implicit none
integer, intent(in)    :: frame
real*8 , intent(in)    :: t_rate
logical, intent(inout) :: triggered

! local variables ...

! QM_erg = E_occ - E_empty ; TO BE USED IN "MM_dynamics" FOR ENERGY BALANCE ...
Unit_Cell% QM_erg    =  update_QM_erg( t_rate , triggered )
Unit_Cell% Total_erg =  Unit_Cell% MD_Kin + Unit_Cell% MD_Pot + Unit_Cell% QM_erg

if( mod(frame,MM_log_step) == 0 ) then 

  open( unit = 13 , file = "dyn.trunk/classical_E.dat" , status = "unknown", action = "write" , position = "append" )
  open( unit = 16 , file = "dyn.trunk/quantum_E.dat"   , status = "unknown", action = "write" , position = "append" )

  select case (Units_MM)
    case( "eV" )    
        write(13,'(I7,4F15.5)') frame, Unit_Cell% MD_Kin, Unit_Cell% MD_Pot, Unit_Cell% MD_Kin + Unit_Cell% MD_Pot 
        write(16,'(I7,2F15.5)') frame, Unit_Cell% QM_erg, Unit_Cell% Total_erg

    case( "kj-mol" )
        write(13,'(I7,4F15.5)') frame, Unit_Cell% MD_Kin*eV_2_kJmol, Unit_Cell% MD_Pot*eV_2_kJmol, (Unit_Cell% MD_Kin + Unit_Cell% MD_Pot)*eV_2_kJmol
        write(16,'(I7,2F15.5)') frame, Unit_Cell% QM_erg*eV_2_kJmol, Unit_Cell% Total_erg*eV_2_kJmol

    end select

    close(13)
    close(16)

end if

end subroutine  QMMM_erg_status 
!
!
!
!==============================================
 function update_QM_erg( t_rate , triggered ) &
 result(QM_erg)
!==============================================
implicit none
real*8 , optional , intent(in)    :: t_rate
logical, optional , intent(inout) :: triggered

! local variables ...
integer    :: n
real*8     :: QM_erg , Frontier_gap
complex*16 :: wp_energy(n_part)
logical    :: jump


!            If( it == 1) then
!                QM_erg = UNI%erg(electron_state) - UNI%erg(hole_state)
!                return
!            else
!                QM_erg = UNI%erg(PST(1)) - UNI%erg(PST(2))
!            end If
!
!            If( PST(1) == PST(2) ) then
!            
!               call verify_FSSH_jump( UNI%R , MO_bra , MO_ket , t_rate , jump , method = "Dynemol" ) 
!               
!               Frontier_Gap = UNI% erg( UNI%Fermi_state + 1 ) - UNI% erg( UNI%Fermi_state )
!               
!               if( (Unit_Cell% MD_Kin > Frontier_Gap) .AND. (jump == .true.) ) triggered = yes
!            
!            end if

!===============================================
           do n = 1 , n_part
              wp_energy(n) = sum(MO_bra(:,n)*UNI%erg(:)*MO_ket(:,n)) 
              end do
              QM_erg = real( wp_energy(1) ) - real( wp_energy(2) )

           If( it == 1) then
              return
           elseif( (QM_erg <= d_zero) .OR. (PST(1) == PST(2)) ) then
              PST(:)    = UNI % Fermi_state
              QM_erg    = d_zero
              triggered = NO
           else
              triggered = YES
           end if

!                
!           If( PST(1) == PST(2) ) then

!               Frontier_Gap = UNI% erg( UNI%Fermi_state + 1 ) - UNI% erg( UNI%Fermi_state )
!               
!               if( (Unit_Cell% MD_Kin > Frontier_Gap) .AND. (jump == .true.) ) triggered = yes
!           
!           end if

end function update_QM_erg
!
!
!
!
!===============================================
subroutine Restart_stuff( QDyn , frame_restart )
!===============================================
implicit none
type(f_time)    , intent(out) :: QDyn
integer         , intent(out) :: frame_restart

CALL DeAllocate_QDyn ( QDyn , flag="alloc" )

CALL Restart_State   ( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame_restart )

allocate( phase(size(MO_bra(:,1))) )

CALL Restart_Sys( Extended_Cell , ExCell_basis , Unit_Cell , DUAL_ket , AO_bra , AO_ket , frame_restart , UNI )

mm = size(ExCell_basis)
nn = n_part

allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

If( Induced_ ) then
     CALL Build_Induced_DP( instance = "allocate" )
     CALL DP_stuff ( "Induced_DP" )
end If

end subroutine Restart_stuff
!
!
end module CSDM_adiabatic_m

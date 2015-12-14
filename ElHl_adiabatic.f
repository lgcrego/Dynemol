! Subroutine for computing time evolution adiabatic on the AO
module ElHl_adiabatic_m

    use type_m
    use constants_m
    use blas95
    use parameters_m                , only : t_i , n_t , t_f , n_part ,     &
                                             frame_step , nuclear_matter ,  &
                                             DP_Field_ , DP_Moment ,        &
                                             Induced_ , QMMM ,              &
                                             GaussianCube , static ,        &
                                             GaussianCube_step ,            &
                                             hole_state , initial_state ,   &
                                             restart , Coulomb_ 
    use Babel_m                     , only : Coords_from_Universe ,         &
                                             trj ,                          &
                                             MD_dt
    use Allocation_m                , only : Allocate_UnitCell ,            &
                                             DeAllocate_UnitCell ,          &
                                             DeAllocate_Structures ,        &
                                             Allocate_Brackets 
    use Structure_Builder           , only : Unit_Cell ,                    &
                                             Extended_Cell ,                &
                                             Generate_Structure ,           &
                                             Basis_Builder ,                &
                                             ExCell_basis
    use FMO_m                       , only : FMO_analysis ,                 &
                                             orbital , eh_tag
    use DP_main_m                   , only : Dipole_Matrix ,                &
                                             Dipole_Moment
    use TD_Dipole_m                 , only : wavepacket_DP                                        
    use DP_potential_m              , only : Molecular_DPs                                              
    use Polarizability_m            , only : Build_Induced_DP
    use Solvated_M                  , only : Prepare_Solvated_System 
    use QCModel_Huckel_ElHl         , only : EigenSystem_ElHl                                            
    use Schroedinger_m              , only : DeAllocate_QDyn
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format
    use Data_Output                 , only : Populations ,                  &
                                             Net_Charge
    use Backup_m                    , only : Security_Copy ,                &
                                             Restart_state ,                &
                                             Restart_Sys
    use MD_read_m                   , only : atom                                             
    use MM_dynamics_m               , only : MolecularMechanics ,           &
                                             preprocess_MM , MoveToBoxCM
    use Ehrenfest_Builder           , only : EhrenfestForce


    public :: ElHl_adiabatic

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra , phase
    real*8     , ALLOCATABLE , dimension(:)   :: Net_Charge_MM(:)
    type(R_eigen)                             :: UNI , el_FMO , hl_FMO
    type(R_eigen)                             :: UNI_el , UNI_hl
    integer                                   :: mm , nn
    logical                                   :: static_hole = .false.

contains
!
!
!======================================
 subroutine ElHl_adiabatic( Qdyn , it )
!======================================
implicit none
type(f_time)    , intent(out)   :: QDyn
integer         , intent(out)   :: it

! local variables ...   
real*8                 :: t , t_rate 
integer                :: j , frame , frame_init , frame_final , frame_restart
character(11)          :: argument
type(universe)         :: Solvated_System

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
    CALL Restart_stuff( QDyn , t , it , frame_restart )
else
    CALL Preprocess( QDyn , it )
end If

!-------------------------------------------------------------
! get command line argument to turn off hole dynamics ...
CALL GET_COMMAND_ARGUMENT( 1 , argument ) 

if( COMMAND_ARGUMENT_COUNT() /= 0 ) then
    select case ( argument )

        case( "static_hole" ) 
        static_hole = .true.

        case default 
        write(*,'(a11)', advance="no") argument
        stop " >>> wrong input argument in ElHl_dynamics <<< "

    end select
end if
!-------------------------------------------------------------

frame_init = merge( frame_restart+1 , frame_step+1 , restart )

do frame = frame_init , frame_final , frame_step

    t = t + t_rate 

    If( (it >= n_t) .OR. (t >= t_f) ) exit    

    it = it + 1

    ! propagate t -> (t + t_rate) with UNI%erg(t) ...
    !============================================================================
    phase(:,1) = cdexp(- zi * UNI_el%erg(:) * t_rate / h_bar)
    phase(:,2) = cdexp(- zi * UNI_hl%erg(:) * t_rate / h_bar)

    if( static_hole ) phase(:,2) = (1.d0,0.d0)

    forall( j=1:n_part )   
        MO_bra(:,j) = conjg(phase(:,j)) * MO_bra(:,j) 
        MO_ket(:,j) =       phase(:,j)  * MO_ket(:,j) 
    end forall

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_bra(:,1) , mm , C_zero , DUAL_bra(:,1) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , MO_ket(:,1) , mm , C_zero , DUAL_ket(:,1) , mm )

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_bra(:,2) , mm , C_zero , DUAL_bra(:,2) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , MO_ket(:,2) , mm , C_zero , DUAL_ket(:,2) , mm )

    If( QMMM ) Net_Charge_MM = Net_Charge

    ! save populations(t + t_rate)  and  update Net_Charge ...
    QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t )

    CALL dump_Qdyn( Qdyn , it )

    If( GaussianCube .AND. mod(frame,GaussianCube_step) == 0 ) CALL  Send_to_GaussianCube( frame , t )

    If( DP_Moment ) CALL DP_stuff( t , "DP_moment" )

    If( QMMM ) CALL EhrenfestForce( Extended_Cell , ExCell_basis , MO_bra , MO_ket , UNI_el , UNI_hl )

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

            CALL MolecularMechanics( t_rate , frame - 1 , Net_Charge = Net_Charge_MM )    ! <== MM precedes QM ...

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
            stop

    end select

    CALL Generate_Structure  ( frame )

    CALL Basis_Builder       ( Extended_Cell , ExCell_basis )

    If( DP_field_ )          CALL DP_stuff ( t , "DP_field"   )

    If( Induced_ .OR. QMMM ) CALL DP_stuff ( t , "Induced_DP" )

    ! LOCAL representation for Coulomb calculation ...
    AO_bra = DUAL_bra

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_ket(:,1) , mm , C_zero , AO_ket(:,1) , mm )
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_ket(:,2) , mm , C_zero , AO_ket(:,2) , mm )

    Deallocate              ( UNI_el%R , UNI_el%L , UNI_el%erg )
    Deallocate              ( UNI_hl%R , UNI_hl%L , UNI_hl%erg )

    CALL EigenSystem_ElHl   ( Extended_Cell , ExCell_basis , AO_bra , AO_ket , UNI_el , UNI_hl , flag2=it )

    ! project back to MO_basis with UNI(t + t_rate)
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , Dual_bra(:,1) , mm , C_zero , MO_bra(:,1) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , Dual_ket(:,1) , mm , C_zero , MO_ket(:,1) , mm )

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , Dual_bra(:,2) , mm , C_zero , MO_bra(:,2) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , Dual_ket(:,2) , mm , C_zero , MO_ket(:,2) , mm )

    !============================================================================

    CALL Security_Copy( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame )

    print*, frame 

end do

deallocate( MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

include 'formats.h'

end subroutine ElHl_adiabatic
!
!
!
!==================================
 subroutine Preprocess( QDyn , it )
!==================================
implicit none
type(f_time)    , intent(out)    :: QDyn
integer         , intent(in)     :: it

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

CALL Generate_Structure ( 1 )

CALL Basis_Builder ( Extended_Cell , ExCell_basis )    

If( Induced_ .OR. QMMM ) CALL Build_Induced_DP ( basis = ExCell_basis , instance = "allocate" )

If( DP_field_ ) then
    hole_save  = hole_state
    hole_state = 0
    static     = .true. 

    ! DP potential in the static GS configuration ...
    CALL Molecular_DPs  ( Extended_Cell )

    hole_state = hole_save
    static     = .false.
end If

If( DP_Field_ .OR. Induced_ .OR. QMMM ) CALL Dipole_Matrix ( Extended_Cell , ExCell_basis )

CALL EigenSystem_ElHl   ( Extended_Cell , ExCell_basis , QM_el=UNI_el , QM_hl=UNI_hl , flag2=it )

! build the initial electron-hole wavepacket ...
CALL FMO_analysis( Extended_Cell , ExCell_basis , UNI_el%R , el_FMO , instance="E" )

CALL FMO_analysis( Extended_Cell , ExCell_basis , UNI_hl%R , hl_FMO , instance="H" )

CALL Allocate_Brackets  ( size(ExCell_basis)  ,       & 
                          MO_bra   , MO_ket   ,       &
                          AO_bra   , AO_ket   ,       &
                          DUAL_bra , DUAL_ket ,       &
                          phase )

If( QMMM ) allocate( Net_Charge_MM(Extended_Cell%atoms) , source = D_zero )

mm = size(ExCell_basis)                         

! initial state of the isolated molecule ...
Print 56 , initial_state     

! building up the electron and hole wavepackets with expansion coefficients at t = 0 ...
do n = 1 , n_part                         
    select case( eh_tag(n) )

        case( "el" )

            MO_bra( : , n ) = el_FMO%L( : , orbital(n) )    
            MO_ket( : , n ) = el_FMO%R( : , orbital(n) )   

            Print 591, orbital(n) , el_FMO%erg(orbital(n))
        
        case( "hl" )

            If( (orbital(n) > hl_FMO%Fermi_State) ) print*, '>>> warning: hole state above the Fermi level <<<'

            MO_bra( : , n ) = hl_FMO%L( : , orbital(n) )    
            MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )   

            Print 592, orbital(n) , hl_FMO%erg(orbital(n))

        end select
end do

! DUAL representation for efficient calculation of survival probabilities ...
CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_bra(:,1) , mm , C_zero , DUAL_bra(:,1) , mm )
CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , MO_ket(:,1) , mm , C_zero , DUAL_ket(:,1) , mm )

CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_bra(:,2) , mm , C_zero , DUAL_bra(:,2) , mm )
CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , MO_ket(:,2) , mm , C_zero , DUAL_ket(:,2) , mm )

! save populations ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t_i )

CALL dump_Qdyn( Qdyn , it )

If( GaussianCube       ) CALL Send_to_GaussianCube  ( it , t_i )

If( DP_Moment          ) CALL DP_stuff ( t_i , "DP_matrix" )

If( DP_Moment          ) CALL DP_stuff ( t_i , "DP_moment" )

If( Induced_ .OR. QMMM ) CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

!..........................................................................

include 'formats.h'

end subroutine Preprocess
!
!
!
! 
!=========================================
 subroutine Send_to_GaussianCube( it , t )
!=========================================
implicit none
integer     , intent(in)    :: it
real*8      , intent(in)    :: t

! local variables ...
integer :: n

! LOCAL representation for film STO production ...

! coefs of <k(t)| in AO basis 
AO_bra = DUAL_bra

! coefs of |k(t)> in AO basis 
CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_ket(:,1) , mm , C_zero , AO_ket(:,1) , mm )
CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_ket(:,2) , mm , C_zero , AO_ket(:,2) , mm )

do n = 1 , n_part
    CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
end do

!----------------------------------------------------------

end subroutine Send_to_GaussianCube
!
!
!
!
!===================================
 subroutine DP_stuff( t , instance )
!===================================
implicit none
real*8          , intent(in)    :: t
character(*)    , intent(in)    :: instance

!local variables ...
integer :: i
real*8  :: Total_DP(3)

!----------------------------------------------------------
!       LOCAL representation for DP calculation ...

! coefs of <k(t)| in AO basis 
AO_bra = DUAL_bra

! coefs of |k(t)> in AO basis 
CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_ket(:,1) , mm , C_zero , AO_ket(:,1) , mm )
CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_ket(:,2) , mm , C_zero , AO_ket(:,2) , mm )

select case( instance )

    case( "DP_matrix" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

    case( "DP_field" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        ! wavepacket component of the dipole vector ...
        CALL wavepacket_DP( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_ket )

        CALL Molecular_DPs( Extended_Cell )

    case( "DP_moment" )

!        CALL Dipole_Moment( Extended_Cell , ExCell_basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

        If( t == t_i ) then
            open( unit = 51 , file = "tmp_data/dipole_dyn.dat" , status = "replace" )
        else
            open( unit = 51 , file = "tmp_data/dipole_dyn.dat" , status = "unknown", action = "write" , position = "append" )
        end If
        write(51,'(5F10.5)') t , (Total_DP(i) , i=1,3) , sqrt( sum(Total_DP*Total_DP) )
        close(51)

    case( "Induced_DP" ) 

        If( .NOT. DP_field_ ) CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        CALL Build_Induced_DP( ExCell_basis , Dual_bra , Dual_ket )

end select

!----------------------------------------------------------

end subroutine DP_stuff
!
!
!
!
!=================================
 subroutine dump_Qdyn( Qdyn , it )
!=================================
implicit none
type(f_time)    , intent(in) :: QDyn
integer         , intent(in) :: it 

! local variables ...
integer     :: nf , n
complex*16  :: wp_energy

do n = 1 , n_part

    select case( n_part )
     
        case( 1 ) 
        wp_energy = sum(MO_bra(:,n)*UNI_el%erg(:)*MO_ket(:,n))

        case( 2 ) 
        wp_energy = sum(MO_bra(:,n)*UNI_hl%erg(:)*MO_ket(:,n))

    end select

    If( it == 1 ) then

        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
        write(52,12) "#" , QDyn%fragments , "total"

        open( unit = 53 , file = "tmp_data/"//eh_tag(n)//"_wp_energy.dat" , status = "replace" , action = "write" , position = "append" )
        write(53,14) QDyn%dyn(it,0,n) , real(wp_energy) , dimag(wp_energy)

    else

        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat"  , status = "unknown", action = "write" , position = "append" )
        open( unit = 53 , file = "tmp_data/"//eh_tag(n)//"_wp_energy.dat" , status = "unknown", action = "write" , position = "append" )

    end If

    ! dumps el-&-hl populations ...
    write(52,13) ( QDyn%dyn(it,nf,n) , nf=0,size(QDyn%fragments)+1 ) 

    ! dumps el-&-hl wavepachet energies ...
    write(53,14) QDyn%dyn(it,0,n) , real(wp_energy) , dimag(wp_energy)

    close(52)
    close(53)

end do

12 FORMAT(15A10)
13 FORMAT(15F10.5)
14 FORMAT(3F12.6)

end subroutine dump_Qdyn
!
!
!
!
!========================================================
subroutine Restart_stuff( QDyn , t , it , frame_restart )
!========================================================
implicit none
type(f_time)    , intent(out)   :: QDyn
real*8          , intent(inout) :: t
integer         , intent(inout) :: it
integer         , intent(inout) :: frame_restart

CALL DeAllocate_QDyn ( QDyn , flag="alloc" )

CALL Restart_State   ( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame_restart )

allocate( phase(size(MO_bra(:,1)),size(MO_bra(1,:))) )

CALL Restart_Sys     ( Extended_Cell , ExCell_basis , Unit_Cell , DUAL_ket , AO_bra , AO_ket , frame_restart , it , UNI_el , UNI_hl )

mm = size(ExCell_basis)
nn = n_part

If( QMMM ) then 

    allocate( Net_Charge_MM (Extended_Cell%atoms) , source = D_zero )

    CALL Build_Induced_DP( instance = "allocate" )

    CALL DP_stuff ( t , "Induced_DP" )

end If

end subroutine Restart_stuff
!
!
!
! 
end module ElHl_adiabatic_m

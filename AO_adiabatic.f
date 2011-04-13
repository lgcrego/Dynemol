! Subroutine for computing time evolution adiabatic on the AO
module AO_adiabatic_m

    use type_m
    use constants_m
    use mkl95_blas
    use parameters_m                , only : t_i , n_t , t_f , n_part ,     &
                                             frame_step , state_of_matter , &
                                             DP_Field_ , DP_Moment ,        &
                                             GaussianCube , static ,        &
                                             GaussianCube_step ,            &
                                             hole_state , initial_state ,   &
                                             restart           
    use Babel_m                     , only : Coords_from_Universe ,         &
                                             trj ,                          &
                                             MD_dt
    use Allocation_m                , only : Allocate_UnitCell ,            &
                                             DeAllocate_UnitCell ,          &
                                             DeAllocate_Structures ,        &
                                             Allocate_Brackets ,            &
                                             ReAllocate_Brackets
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
    use Solvated_M                  , only : Prepare_Solvated_System 
    use QCModel_Huckel              , only : EigenSystem                                                 
    use Schroedinger_m              , only : DeAllocate_QDyn
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format
    use Data_Output                 , only : Populations
    use Backup_m                    , only : Security_Copy ,                &
                                             Restart_state ,                &
                                             Restart_Sys

    public :: AO_adiabatic

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Complex*16 , ALLOCATABLE , dimension(:)   :: bra , ket , phase
    type(C_eigen)                             :: UNI , el_FMO , hl_FMO

contains
!
!
!====================================
 subroutine AO_adiabatic( Qdyn , it )
!====================================
implicit none
type(f_time)    , intent(out)   :: QDyn
integer         , intent(out)   :: it

! local variables ...
integer                :: j , frame , frame_init , frame_restart
real*8                 :: t , t_rate 
type(universe)         :: Solvated_System

it = 1
t  = t_i

If( restart ) then
    CALL Restart_stuff( QDyn , t , it , frame_restart )
else
    CALL Preprocess( QDyn , it )
end IF

frame_init = merge( frame_restart+1 , frame_step+1 , restart )

!--------------------------------------------------------------------------------
! time slicing H(t) : Quantum Dynamics & All that Jazz ...

t_rate = MD_dt * frame_step

do frame = frame_init , size(trj) , frame_step

    t = t + t_rate 

    if( (it >= n_t) .OR. (t >= t_f) ) exit    

    it = it + 1

    ! propagate t -> (t + t_rate) with UNI%erg(t) ...
    !============================================================================
    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    forall( j=1:n_part )   
        MO_bra(:,j) = conjg(phase(:)) * MO_bra(:,j) 
        MO_ket(:,j) =       phase(:)  * MO_ket(:,j) 
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    CALL gemm( UNI%L , MO_bra , DUAL_bra , 'T' , 'N' , C_one , C_zero )
    CALL gemm( UNI%R , MO_ket , DUAL_ket , 'N' , 'N' , C_one , C_zero )

    ! save populations(t + t_rate) ...
    QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t )

    CALL dump_Qdyn( Qdyn , it )

    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) CALL  Send_to_GaussianCube( it , t )

    If( DP_Moment ) CALL DP_stuff( t , "DP_moment" )

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    ! build new UNI(t + t_rate) ...
    !============================================================================

    select case ( state_of_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

        case( "extended_sys" )

            CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

        case default

            Print*, " >>> Check your state_of_matter options <<< :" , state_of_matter
            stop

    end select

    CALL Generate_Structure ( frame )

    CALL Basis_Builder      ( Extended_Cell , ExCell_basis )

    If( DP_field_ )         CALL DP_stuff ( t , "DP_field" )

    Deallocate              ( UNI%R , UNI%L , UNI%erg )

    CALL EigenSystem        ( Extended_Cell , ExCell_basis , UNI , flag2=it )

    ! project back to MO_basis with UNI(t + t_rate)
    CALL gemm( UNI%R , DUAL_bra , MO_bra , 'T' , 'N' , C_one , C_zero )
    CALL gemm( UNI%L , DUAL_ket , MO_ket , 'N' , 'N' , C_one , C_zero )

    !============================================================================

    CALL Security_Copy( MO_bra , MO_ket , DUAL_bra , DUAL_ket , bra , ket , t , it , frame )

    print*, frame 

end do

deallocate( MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

include 'formats.h'

end subroutine AO_adiabatic
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
logical         :: el_hl_
type(universe)  :: Solvated_System

! preprocessing stuff .....................................................

CALL DeAllocate_QDyn( QDyn , flag="alloc" )

select case ( state_of_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , 1 )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System , 1 )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(1) , 1 )

    case default

        Print*, " >>> Check your state_of_matter options <<< :" , state_of_matter
        stop

end select

el_hl_ = any( Unit_Cell%fragment == "H")
 
CALL Generate_Structure ( 1 )

CALL Basis_Builder      ( Extended_Cell , ExCell_basis )

If( DP_field_ ) then
    hole_save  = hole_state
    hole_state = 0
    static     = .true. 

    ! DP potential in the static GS configuration ...
    CALL Molecular_DPs  ( Extended_Cell )

    hole_state = hole_save
    static     = .false.
end If

CALL EigenSystem        ( Extended_Cell , ExCell_basis , UNI , flag2=it )

CALL FMO_analysis       ( Extended_Cell , ExCell_basis , UNI%R , el_FMO , fragment="D" )

If( el_hl_ ) CALL FMO_analysis ( Extended_Cell , ExCell_basis , UNI%R , hl_FMO , fragment="H" )

CALL Allocate_Brackets  ( size(ExCell_basis)  ,       & 
                          MO_bra   , MO_ket   ,       &
                          AO_bra   , AO_ket   ,       &
                          DUAL_bra , DUAL_ket ,       &
                          bra      , ket      , phase )

! initial state of the isolated molecule ...
Print 56 , initial_state     

! building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
! assuming non-interacting electrons ...
do n = 1 , n_part                         
    select case( eh_tag(n) )

        case( "el" )

            MO_bra( : , n ) = el_FMO%L( : , orbital(n) )    
            MO_ket( : , n ) = el_FMO%R( : , orbital(n) )   

            Print 591, orbital(n) , el_FMO%erg(orbital(n))
        
        case( "hl" )

            If( (orbital(n) > hl_FMO%Fermi_State) ) pause '>>> quit: hole state above the Fermi level <<<'

            MO_bra( : , n ) = hl_FMO%L( : , orbital(n) )    
            MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )   

            Print 592, orbital(n) , hl_FMO%erg(orbital(n))

        end select
end do

! DUAL representation for efficient calculation of survival probabilities ...
CALL gemm(UNI%L , MO_bra , DUAL_bra , 'T' , 'N' , C_one , C_zero )
CALL gemm(UNI%R , MO_ket , DUAL_ket , 'N' , 'N' , C_one , C_zero )

! save populations ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t_i )

CALL dump_Qdyn( Qdyn , it )

If( GaussianCube ) CALL Send_to_GaussianCube  ( it , t_i )

If( DP_Moment    ) CALL DP_stuff ( t_i , "DP_matrix" )

If( DP_Moment    ) CALL DP_stuff ( t_i , "DP_moment" )

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

!----------------------------------------------------------
!     LOCAL representation for film STO production ...

! coefs of <k(t)| in AO basis 
AO_bra = DUAL_bra

! coefs of |k(t)> in AO basis 
CALL gemm( UNI%L , MO_ket , AO_ket , 'T' , 'N' , C_one , C_zero )

bra(:) = AO_bra(:,1)
ket(:) = AO_ket(:,1)
   
CALL Gaussian_Cube_Format( bra , ket , it ,t )

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
CALL gemm( UNI%L , MO_ket , AO_ket , 'T' , 'N' , C_one , C_zero )

bra(:) = AO_bra(:,1)
ket(:) = AO_ket(:,1)

select case( instance )

    case( "DP_matrix" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

    case( "DP_field" )

        CALL Dipole_Matrix( Extended_Cell , ExCell_basis )

        ! wavepacket component of the dipole vector ...

        CALL wavepacket_DP( Extended_Cell , ExCell_basis , bra , ket , Dual_ket(:,1) )

        CALL Molecular_DPs( Extended_Cell )

    case( "DP_moment" )

        CALL Dipole_Moment( Extended_Cell , ExCell_basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

        If( t == t_i ) then
            open( unit = 51 , file = "tmp_data/dipole_dyn.dat" , status = "replace" )
        else
            open( unit = 51 , file = "tmp_data/dipole_dyn.dat" , status = "unknown", action = "write" , position = "append" )
        end If
        write(51,'(F9.4,4F10.5)') t , (Total_DP(i) , i=1,3) , sqrt( sum(Total_DP*Total_DP) )
        close(51)

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
integer :: nf , n

do n = 1 , n_part

    If( it == 1 ) then
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
        write(52,12) "#" , QDyn%fragments , "total"
    else
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "unknown", action = "write" , position = "append" )
    end If

    write(52,13) ( QDyn%dyn(it,nf,n) , nf=0,size(QDyn%fragments)+1 ) 

    close(52)

end do

12 FORMAT(10A9)
13 FORMAT(10F9.4)

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

CALL DeAllocate_QDyn     ( QDyn , flag="alloc" )

CALL Restart_State       ( MO_bra , MO_ket , DUAL_bra , DUAL_ket , bra , ket , t , it , frame_restart )

CALL ReAllocate_Brackets ( size(MO_bra(:,1)) , AO_bra , AO_ket , phase )

CALL Restart_Sys         ( Extended_Cell , ExCell_basis , Unit_Cell , UNI , DUAL_ket , bra , ket , frame_restart , it )

end subroutine Restart_stuff
!
!
!
!
end module AO_adiabatic_m

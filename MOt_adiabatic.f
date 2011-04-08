! Subroutine for computing time evolution adiabatic on the MO_t
module MOt_adiabatic_m

    use type_m
    use constants_m
    use parameters_m                , only : t_i , t_f , n_part , frame_step ,          &
                                             DP_Moment , DP_Field_ , initial_state
    use mkl95_blas
    use Data_Output                 , only : Populations 
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
    use DP_potential_m              , only : Molecular_DPs                                              
    use QCModel_Huckel              , only : EigenSystem                                                 
    use Schroedinger_m              , only : DeAllocate_QDyn

    public :: MOt_adiabatic

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Complex*16 , ALLOCATABLE , dimension(:)   :: bra , ket , phase
    type(C_eigen)                             :: el_FMO , hl_FMO , UNI

contains
!
!
!====================================
 subroutine MOt_adiabatic( QDyn , it )
!====================================
implicit none
type(f_time)    ,   intent(out) :: QDyn
integer         ,   intent(out) :: it

! local variables ...
integer  :: j , frame 
real*8   :: t , t_rate 

it = 1
t  = t_i
CALL Preprocess( QDyn , it )

!--------------------------------------------------------------------------------
! time slicing H(t) : Quantum Dynamics & All that Jazz ...

t_rate = MD_dt * frame_step

do frame = (1 + frame_step) , size(trj) , frame_step

    t = t + t_rate 

    if( t >= t_f ) exit

    it = it + 1

    ! propagate t -> (t + t_rate) with UNI%erg(t) ...
    !============================================================================

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    forall( j=1:n_part )   
        MO_bra(:,j) = conjg(phase(:)) * MO_bra(:,j) 
        MO_ket(:,j) =       phase(:)  * MO_ket(:,j) 
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)
    CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

    ! save populations(t + t_rate) ...
    QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t )

    If( DP_Moment ) Print*, ">>> DP_Moment not implemented for this routine <<<"

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    ! build new UNI(t + t_rate) ...
    !============================================================================

    CALL Coords_from_Universe   ( Unit_Cell , trj(frame) , frame )

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) &
    CALL Molecular_DPs          ( Extended_Cell )

    CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it )

    !============================================================================

    print*, frame 

end do

deallocate( MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

include 'formats.h'

end subroutine MOt_adiabatic
!
!
!
!========================================
 subroutine Preprocess( QDyn , it )
!========================================
implicit none
type(f_time)                    ,   intent(out) :: QDyn
integer                         ,   intent(in)  :: it

! local variables ...
integer :: n
logical :: el_hl_

! preprocessing stuff .....................................................

CALL DeAllocate_QDyn        ( QDyn , flag="alloc" )

CALL Coords_from_Universe   ( Unit_Cell , trj(1) , 1 )

el_hl_ = any( Unit_Cell%fragment == "H")
 
CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) &
CALL Molecular_DPs          ( Extended_Cell )

CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it )

CALL FMO_analysis           ( Extended_Cell , ExCell_basis , UNI%R , el_FMO , fragment="D" )

If( el_hl_ ) CALL FMO_analysis ( Extended_Cell , ExCell_basis , UNI%R , hl_FMO , fragment="H" )

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(UNI%L(1,:))     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

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
CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)
CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

! save populations ...
QDyn%dyn(it,:,:) = Populations( QDyn%fragments , ExCell_basis , DUAL_bra , DUAL_ket , t_i )

If( DP_Moment ) Print*, ">>> DP_Moment not implemented for this routine <<<"

!..........................................................................

include 'formats.h'

end subroutine Preprocess
!
!
!
end module MOt_adiabatic_m

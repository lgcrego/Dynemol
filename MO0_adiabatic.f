! Subroutine for computing time evolution adiabatic on the MO_0
module MO0_adiabatic_m

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
                                             orbital
    use DP_potential_m              , only : Molecular_DPs                                              
    use QCModel_Huckel              , only : EigenSystem                                                 
    use Schroedinger_m              , only : DeAllocate_QDyn

    public :: MO0_adiabatic

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Complex*16 , ALLOCATABLE , dimension(:)   :: bra , ket , phase
    type(C_eigen)                             :: FMO , UNI_0

contains
!
!
!====================================
 subroutine MO0_adiabatic( QDyn , it )
!====================================
implicit none
type(f_time)    ,   intent(out) :: QDyn
integer         ,   intent(out) :: it

! local variables ...
integer                :: j , frame 
real*8                 :: t , t_rate 
real*8  , allocatable  :: erg(:) , QDyn_temp(:,:) 

it = 1
t  = t_i
CALL Preprocess( erg , QDyn , it )

!--------------------------------------------------------------------------------
! time slicing H(t) : Quantum Dynamics & All that Jazz ...

t_rate = MD_dt * frame_step

do frame = (1 + frame_step) , size(trj) , frame_step

    t = t + t_rate 

    if( t >= t_f ) exit

    it = it + 1

    ! propagate t -> (t + t_rate) with UNI%erg(t) ...
    !============================================================================

    phase(:) = cdexp(- zi * erg(:) * t_rate / h_bar)

    forall( j=1:n_part )   
        MO_bra(:,j) = conjg(phase(:)) * MO_bra(:,j) 
        MO_ket(:,j) =       phase(:)  * MO_ket(:,j) 
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    ! coefs of <k(t)| in DUAL basis ...
    CALL gemm(UNI_0%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

    ! coefs of |k(t)> in DUAL basis ...
    CALL gemm(UNI_0%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

    ! save populations(t + t_rate) ...
    bra(:) = DUAL_bra(:,1)
    ket(:) = DUAL_ket(:,1)

    QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t )

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

    CALL EigenSystem            ( Extended_Cell , ExCell_basis , erg )

    !============================================================================

    print*, frame 

end do

deallocate( MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase , erg )

include 'formats.h'

end subroutine MO0_adiabatic
!
!
!
!========================================
 subroutine Preprocess( erg , QDyn , it )
!========================================
implicit none
real*8        ,  allocatable    ,   intent(out) :: erg(:)
type(f_time)                    ,   intent(out) :: QDyn
integer                         ,   intent(in)  :: it

! local variables ...
integer :: i

! preprocessing stuff .....................................................

CALL DeAllocate_QDyn        ( QDyn , flag="alloc" )

CALL Coords_from_Universe   ( Unit_Cell , trj(1) , 1 )

CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) &
CALL Molecular_DPs          ( Extended_Cell )

CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI_0 , flag2=0 )

CALL FMO_analysis           ( Extended_Cell , ExCell_basis , UNI_0%R , FMO )

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(UNI_0%L(1,:))   ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

MO_bra = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0
MO_ket = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0

! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI_0%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI_0%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

! save populations ...
bra(:) = DUAL_bra(:,1)
ket(:) = DUAL_ket(:,1)

QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t_i )

If( DP_Moment ) Print*, ">>> DP_Moment not implemented for this routine <<<"

!..........................................................................

allocate( erg(size(ExCell_basis)) )
forall( i=1:size(ExCell_basis) ) erg(i) = UNI_0%erg(i)

include 'formats.h'

end subroutine Preprocess
!
!
!
!
end module MO0_adiabatic_m

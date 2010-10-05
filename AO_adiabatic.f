! Subroutine for computing time evolution adiabatic on the AO
module AO_adiabatic_m

    use type_m
    use constants_m
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
    use Solvated_M                  , only : Prepare_Solvated_System 
    use Dipole_potential_m          , only : Molecular_DPs                                              
    use QCModel_Huckel              , only : EigenSystem                                                 
    use Schroedinger_m              , only : DeAllocate_QDyn
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format

    public :: AO_adiabatic

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:) :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Complex*16 , ALLOCATABLE , dimension(:)   :: bra , ket , phase
    type(C_eigen)                             :: UNI , FMO 

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
integer                :: j , frame 
real*8                 :: t , t_rate 
real*8  , allocatable  :: QDyn_temp(:,:)
type(universe)         :: Solvated_System

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
    ! coefs of <k(t)| in DUAL basis ...
    CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

    ! coefs of |k(t)> in DUAL basis ...
    CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

    ! save populations(t + t_rate) ...
    bra(:) = DUAL_bra(:,1)
    ket(:) = DUAL_ket(:,1)

    QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t )

    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) CALL  Send_to_GaussianCube( UNI%L , it , t )

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )
    Deallocate                  ( UNI%R , UNI%L , UNI%erg )

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

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) &
    CALL Molecular_DPs          ( Extended_Cell )

    CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it )

    ! project back to MO_basis with UNI(t + t_rate)
    CALL gemm(UNI%R,DUAL_bra,MO_bra,'T','N',C_one,C_zero)
    CALL gemm(UNI%L,DUAL_ket,MO_ket,'N','N',C_one,C_zero)

    !============================================================================

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
integer         :: frame
type(universe)  :: Solvated_System

! preprocessing stuff .....................................................

CALL DeAllocate_QDyn        ( QDyn , flag="alloc" )

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

CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) &
CALL Molecular_DPs          ( Extended_Cell )

CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it )

CALL FMO_analysis           ( Extended_Cell , ExCell_basis , UNI%R , FMO )

Print 56 , initial_state     ! <== initial state of the isolated molecule 

CALL Allocate_Brackets( size(UNI%L(1,:))     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

MO_bra = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0
MO_ket = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0

! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

! save populations ...
bra(:) = DUAL_bra(:,1)
ket(:) = DUAL_ket(:,1)

QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t_i )

If( GaussianCube ) CALL Send_to_GaussianCube( UNI%L , it , t_i )

!..........................................................................

include 'formats.h'

end subroutine Preprocess
!
!
!
! 
!================================================
subroutine Send_to_GaussianCube( UNI_L , it , t )
!================================================
implicit none
complex*16  , intent(in)    :: UNI_L(:,:)
integer     , intent(in)    :: it
real*8      , intent(in)    :: t

!----------------------------------------------------------
! LOCAL representation for film STO production ...

! coefs of <k(t)| in AO basis 
CALL gemm(UNI_L,MO_bra,AO_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in AO basis 
CALL gemm(UNI_L,MO_ket,AO_ket,'T','N',C_one,C_zero)

bra(:) = AO_bra(:,1)
ket(:) = AO_ket(:,1)
   
CALL Gaussian_Cube_Format(bra,ket,it,t)

!----------------------------------------------------------

end subroutine Send_to_GaussianCube
!
!
!
!
end module AO_adiabatic_m

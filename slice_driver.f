! Subroutine for computing time evolution through time slices
module EigenSlice_m

    use type_m
    use constants_m
    use mkl95_blas
    use Data_Output                 , only : Populations , Dump_stuff 
    use Babel_m                     , only : Coords_from_Universe ,         &
                                             trj
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
    use Dipole_potential_m          , only : Solvent_Molecule_DP                                              
    use QCModel_Huckel              , only : EigenSystem                                                 
    use Schroedinger_m              , only : DeAllocate_QDyn

    public :: Eigen_Slice

    private

contains
!
!
!=======================
 subroutine Eigen_Slice
!=======================
implicit none

! local variables ...
integer                     :: j , it , frame 
real*8                      :: t 
real*8                      :: t_rate = MD_dt * frame_step
real*8       , allocatable  :: QDyn_temp(:,:)
type(f_time)                :: QDyn
type(C_eigen)               :: UNI , FMO

! more local variables ...
complex*16 , ALLOCATABLE :: zG_L(:,:)     , MO_bra(:,:) 
complex*16 , ALLOCATABLE :: zG_R(:,:)     , MO_ket(:,:)
complex*16 , ALLOCATABLE :: AO_bra(:,:)   , AO_ket(:,:) 
complex*16 , ALLOCATABLE :: DUAL_ket(:,:) , DUAL_bra(:,:) 
complex*16 , ALLOCATABLE :: phase(:)      , bra(:)        , ket(:)

! preprocessing stuff .....................................................
it = 1

CALL DeAllocate_QDyn        ( QDyn , flag="alloc" )

CALL Coords_from_Universe   ( Unit_Cell , trj(1) , 1 )

CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) CALL Solvent_Molecule_DP( Extended_Cell )

CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it )

CALL FMO_analysis           ( Extended_Cell , ExCell_basis , UNI%R , FMO )

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(UNI%L(1,:))     ,      &
                         zG_L     , zG_R     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

zG_L = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0
zG_R = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0

! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,zG_L,DUAL_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,zG_R,DUAL_ket,'N','N',C_one,C_zero)

! save populations(time) ...
bra(:) = DUAL_bra(:,1)
ket(:) = DUAL_ket(:,1)

QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t_i )

!..........................................................................
! time slicing H(t) : Quantum Dynamics & All that Jazz ...
it = 2
t  = t_i

do frame = (1 + frame_step) , size(trj) , frame_step

    ! propagate until from t -> (t + t_rate) with UNI%erg(t) ...
    t  =  t + t_rate 

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    forall(j=1:n_part)   
        MO_bra(:,j) = conjg(phase(:)) * zG_L(:,j) 
        MO_ket(:,j) =       phase(:)  * zG_R(:,j) 
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    ! coefs of <k(t)| in DUAL basis ...
    CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

    ! coefs of |k(t)> in DUAL basis ...
    CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

    ! save populations(time) ...
    bra(:) = DUAL_bra(:,1)
    ket(:) = DUAL_ket(:,1)

    QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t )

    zG_L = MO_bra       ! <== updating expansion coefficients at t 
    zG_R = MO_ket       ! <== updating expansion coefficients at t

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    if( t >= t_f ) exit

    CALL Coords_from_Universe   ( Unit_Cell , trj(frame) , frame )

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) &
    CALL Solvent_Molecule_DP    ( Extended_Cell )

    CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it )

    ! coefs of <k(t)| in DUAL basis ...
    CALL gemm(UNI%R,DUAL_bra,zG_L,'T','N',C_one,C_zero)

    ! coefs of |k(t)> in DUAL basis ...
    CALL gemm(UNI%L,DUAL_ket,zG_R,'N','N',C_one,C_zero)

    it = it + 1

    print*, frame 
end do

! prepare data for survival probability ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn )

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

deallocate( zG_L , zG_R , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

include 'formats.h'

end subroutine Eigen_Slice
!
!
!
end module EigenSlice_m

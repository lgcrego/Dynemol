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
    use Semi_Empirical_Parms        , only : Read_Atomic_Mass ,             &
                                             Atomic_Mass

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
complex*16 , ALLOCATABLE :: AO_bra(:,:)   , AO_ket(:,:) , Psi(:) , Psi_A(:) , Psi_temp(:)
complex*16 , ALLOCATABLE :: DUAL_ket(:,:) , DUAL_bra(:,:) 
complex*16 , ALLOCATABLE :: phase(:)      , bra(:)        , ket(:)
complex*16 , ALLOCATABLE :: A_r(:)
real*8     , allocatable :: H(:,:) , S(:,:)
integer    , allocatable :: ind(:)
real*8                   :: var
integer :: i , counter

! preprocessing stuff .....................................................
it = 1

CALL DeAllocate_QDyn        ( QDyn , flag="alloc" )

CALL Coords_from_Universe   ( Unit_Cell , trj(1) , 1 )

CALL Generate_Structure     ( 1 )

CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

If( DP_field_ ) CALL Solvent_Molecule_DP( Extended_Cell )

CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it , flag3=H , flag4=S)

CALL FMO_analysis           ( Extended_Cell , ExCell_basis , UNI%R , FMO )

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(UNI%L(1,:))     ,      &
                         zG_L     , zG_R     ,      &
                         MO_bra   , MO_ket   ,      &
                         AO_bra   , AO_ket   ,      &
                         DUAL_bra , DUAL_ket ,      &
                         bra      , ket      , phase)

zG_L = FMO%L( : , orbital(0:n_part-1) )    ! <== expansion coefficients at t = 0
zG_R = FMO%R( : , orbital(0:n_part-1) )    ! <== expansion coefficients at t = 0

CALL Read_Atomic_Mass

! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,zG_L,DUAL_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,zG_R,DUAL_ket,'N','N',C_one,C_zero)

! save populations(time) ...
bra(:) = DUAL_bra(:,1)
ket(:) = DUAL_ket(:,1)

allocate( A_r      ( size(bra) ) )
allocate( Psi      ( size(bra) ) )
allocate( Psi_A    ( size(bra) ) )
allocate( Psi_temp ( size(bra) ) )

Psi_A  = bra
Psi(:) = zG_L(:,1)

A_r = matmul(S,Psi_A)

write(70,*), 0.0d0 , real( sum( conjg(Psi_A(561:596)) * A_r(561:596) ) )
write(71,*), 0.0d0 , real( sum( conjg(Psi_A(1157:1192)) * A_r(1157:1192) ) )
write(72,*), 0.0d0 , real( sum( conjg(Psi_A(1753:1788)) * A_r(1753:1788) ) )
write(73,*), 0.0d0 , real( sum( conjg(Psi_A(1:560)) * A_r(1:560) ) ) + real( sum( conjg(Psi_A(597:1156)) * A_r(597:1156) ) ) + &
                   + real( sum( conjg(Psi_A(1193:1752)) * A_r(1193:1752) ) )
!..........................................................................
! time slicing H(t) : Quantum Dynamics & All that Jazz ...
it = 2
t  = t_i

!frame = 1

!do
do frame = (1 + frame_step) , size(trj) , frame_step

!    write(35,*), t , real( sum( conjg(Psi_A) * matmul( H , Psi_A ) ) )

    write(1200,*), t , UNI%erg(831) - UNI%erg(830)
    write(1201,*), t , UNI%erg(830) - UNI%erg(829)

    ! propagate until from t -> (t + t_rate) with UNI%erg(t) ...
    t  =  t + t_rate 

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    Psi = phase * Psi

    forall(j=1:n_part)   
        MO_bra(:,j) = conjg(phase(:)) * zG_L(:,j) 
        MO_ket(:,j) =       phase(:)  * zG_R(:,j) 
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    ! coefs of <k(t)| in DUAL basis ...
    CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',C_one,C_zero)

    ! coefs of |k(t)> in DUAL basis ...
    CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

    ! Psi_A(t) ...
    CALL gemv(UNI%L,Psi,Psi_A,C_one,C_zero,'T')

!   ===================== Populations =============================================
    A_r   = matmul(S,Psi_A)
    write(70,*), t , real( sum( conjg(Psi_A(561:596)) * A_r(561:596) ) )
    write(71,*), t , real( sum( conjg(Psi_A(1157:1192)) * A_r(1157:1192) ) )
    write(72,*), t , real( sum( conjg(Psi_A(1753:1788)) * A_r(1753:1788) ) )
    write(73,*), t , real( sum( conjg(Psi_A(1:560)) * A_r(1:560) ) ) + real( sum( conjg(Psi_A(597:1156)) * A_r(597:1156) ) ) + &
                   + real( sum( conjg(Psi_A(1193:1752)) * A_r(1193:1752) ) )

    ! save populations(time) ...
    bra(:) = DUAL_bra(:,1)
    ket(:) = DUAL_ket(:,1)

    QDyn%dyn(it,:) = Populations( QDyn%fragments , ExCell_basis , bra , ket , t )

    zG_L = MO_bra       ! <== updating expansion coefficients at t 
    zG_R = MO_ket       ! <== updating expansion coefficients at t

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )
    Deallocate( S , H )

    if( t >= t_f ) exit

    CALL Coords_from_Universe   ( Unit_Cell , trj(frame) , frame )

    CALL Generate_Structure     ( frame )

    CALL Basis_Builder          ( Extended_Cell , ExCell_basis )

    If( DP_field_ ) CALL Solvent_Molecule_DP ( Extended_Cell )

    CALL EigenSystem            ( Extended_Cell , ExCell_basis , UNI , flag2=it , flag3=H , flag4=S )

    ! Psi_A(t) ...
    CALL gemv(UNI%R,Psi_A,Psi,C_one,C_zero,'T')

    ! coefs of <k(t)| in DUAL basis ...
    CALL gemm(UNI%R,DUAL_bra,zG_L,'T','N',C_one,C_zero)

    ! coefs of |k(t)> in DUAL basis ...
    CALL gemm(UNI%L,DUAL_ket,zG_R,'N','N',C_one,C_zero)

    it = it + 1

    print*, frame 
end do

!CALL DeAllocate_UnitCell    ( Unit_Cell     )
!CALL DeAllocate_Structures  ( Extended_Cell )
!DeAllocate                  ( ExCell_basis  )
!deallocate( S , H )

! prepare data for survival probability ...
allocate ( QDyn_temp( it , 0:size(QDyn%fragments)+1 ) , source=QDyn%dyn( 1:it , 0:size(QDyn%fragments)+1 ) )
CALL move_alloc( from=QDyn_temp , to=QDyn%dyn )

CALL Dump_stuff( QDyn=QDyn )

! final procedures ...
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

deallocate( zG_L , zG_R , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase , A_r )

include 'formats.h'

end subroutine Eigen_Slice
!
!
!
end module EigenSlice_m

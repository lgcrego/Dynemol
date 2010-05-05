! Subroutine for computing time evolution through time slices
module TimeSlice_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Babel_m                     , only : Coords_from_Universe ,         &
                                             trj
    use Allocation_m                , only : Allocate_UnitCell ,            &
                                             DeAllocate_UnitCell ,          &
                                             DeAllocate_Structures ,        &
                                             Allocate_Brackets
    use QCModel_Huckel              , only : EigenSystem
    use FMO_m                       , only : FMO_analysis ,                 &
                                             orbital
    use Structure_Builder           , only : Unit_Cell ,                    &
                                             Extended_Cell ,                &
                                             Generate_Structure ,           &
                                             Basis_Builder ,                &
                                             ExCell_basis
    use Schroedinger_m              , only : DeAllocate_QDyn
    use Data_Output                 , only : Dump_stuff ,                   &
                                             get_Populations
    use dipole_potential_m          , only : Solvent_Molecule_DP
    use Psi_Squared_Cube_Format     , only : Gaussian_Cube_Format
    use Chebyshev_m                 , only : Chebyshev


    public :: Time_Slice

    private

    ! module variables ...
    integer :: it
    real*8  :: t , t_rate
    logical :: done = .false.

    ! Huckel_Slice_Dynamics variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:) :: zG_L , zG_R , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Complex*16 , ALLOCATABLE , dimension(:)   :: phase , bra , ket

    interface preprocess
        module procedure preprocess_Huckel_Slice
        module procedure preprocess_Chebyshev
    end interface

contains
!
!
!=======================
 subroutine Time_Slice
!=======================
implicit none

! local variables ...
integer                     :: frame , P1 , P2 
real*8       , allocatable  :: QDyn(:,:) , QDyn_temp(:,:)
complex*16   , allocatable  :: MO(:)
character(1) , allocatable  :: QDyn_fragments(:)

CALL DeAllocate_QDyn( QDyn , QDyn_fragments , flag="alloc" )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! time slicing H(t) : Quantum Dynamics & All that Jazz ...

do frame = 1 , size(trj) , frame_step

    CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

    CALL Generate_Structure( frame )
   
    CALL Basis_Builder( Extended_Cell , ExCell_basis )

    If( DP_field_ ) CALL Solvent_Molecule_DP( Extended_Cell )

    select case ( DRIVER )

        case( "chebyshev" )

            if( .NOT. done ) CALL preprocess_Chebyshev( Extended_Cell , ExCell_basis , MO , P1 , P2 )
    
            CALL Chebyshev( Extended_Cell , ExCell_basis , MO , it , P1 , P2 )

        case( "eigen_slice" )

            If( .NOT. done ) CALL preprocess( Extended_Cell , ExCell_basis , QDyn )

            CALL Huckel_Slice_Dynamics( Extended_Cell , ExCell_basis , QDyn )

    end select

    CALL DeAllocate_UnitCell   ( Unit_Cell     )
    CALL DeAllocate_Structures ( Extended_Cell )
    DeAllocate                 ( ExCell_basis  )

    if( t >= t_f ) exit

print*, frame
end do

! prepare data for survival probability ...
allocate( QDyn_temp(it-1,size(QDyn_fragments)+1) )
QDyn_temp(:,:) = QDyn(1:it-1,:)
CALL move_alloc( from=QDyn_temp , to=QDyn )

CALL Dump_stuff( QDyn=QDyn , list_of_fragments=QDyn_fragments , t_rate=t_rate)

! turn off the ligths ...

CALL DeAllocate_QDyn( QDyn , QDyn_fragments , flag="dealloc" )

deallocate(zG_L , zG_R , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra , phase , bra , ket )

end subroutine Time_Slice
!
!
!
!=========================================================
 subroutine Huckel_Slice_Dynamics( system , basis , QDyn )
!=========================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(inout) :: QDyn(:,:)

! local variables ...
integer       :: j
type(C_eigen) :: UNI 

! local parameters ...
complex*16 , parameter :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

!------------------------------------------------------------------------
! Q-DYNAMICS ... 

CALL EigenSystem( system , basis, UNI )

phase(:) = exp(- zi * UNI%erg(:) * t_rate / h_bar)

IF( t == t_i ) phase = one

forall( j=1:n_part )   
    MO_bra(:,j) = conjg(phase(:)) * zG_L(:,j) 
    MO_ket(:,j) =       phase(:)  * zG_R(:,j) 
end forall

! DUAL representation for efficient calculation of survival probabilities ...
! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,MO_bra,DUAL_bra,'T','N',one,zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',one,zero)

QDyn(it,:) = get_Populations(system,basis,DUAL_bra,DUAL_ket)

zG_L = MO_bra       ! <== updating expansion coefficients at t 
zG_R = MO_ket       ! <== updating expansion coefficients at t

t  = t  + t_rate
it = it + 1

!------------------------------------------------------------------------
! . LOCAL representation for film STO production ...

IF( GaussianCube ) then

    ! coefs of <k(t)| in AO basis 
    CALL gemm(UNI%L,MO_bra,AO_bra,'T','N',one,zero)

    ! coefs of |k(t)> in AO basis 
    CALL gemm(UNI%L,MO_ket,AO_ket,'T','N',one,zero)

    bra(:) = AO_bra(:,1)
    ket(:) = AO_ket(:,1)
   
    CALL Gaussian_Cube_Format(bra,ket,it,t)

end IF    
!------------------------------------------------------------------------

include 'formats.h'

end subroutine Huckel_Slice_Dynamics
!
!
!
!================================================================
 subroutine preprocess_Chebyshev( system , basis , MO , P1 , P2 )
!================================================================
implicit none
type(structure)                 , intent(in)    :: system
type(STO_basis)                 , intent(in)    :: basis(:)
complex*16      , allocatable   , intent(inout) :: MO(:)
integer                         , intent(inout) :: P1
integer                         , intent(inout) :: P2

!local variables ...
integer                         :: i , li , N
real*8          , allocatable   :: wv_FMO(:)
type(C_eigen)                   :: UNI , FMO

! local parameter ...
complex*16      , parameter     :: z0 = ( 0.0d0 , 0.0d0 )

CALL EigenSystem( system , basis, UNI )

CALL FMO_analysis( system , basis, UNI%R, FMO , wv_FMO )

li = minloc( basis%indx , DIM = 1 , MASK = basis%fragment == "D" )
N  = size(wv_FMO)
P1 = li
P2 = li + N

allocate( MO(size(basis)) , source=z0 )
MO(li:li+N-1) = cmplx( wv_FMO(:) )

deallocate( wv_FMO )

it =  1

done = .true. 

end subroutine preprocess_Chebyshev
!
!
!
!===========================================================
 subroutine preprocess_Huckel_Slice( system , basis , Qdyn )
!===========================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(inout) :: QDyn(:,:)

!local variables ...
type(C_eigen) :: UNI , FMO
complex*16 , parameter :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

CALL EigenSystem( system , basis, UNI )

CALL FMO_analysis( system , basis, UNI%R, FMO )

CALL Allocate_Brackets( size(basis) , zG_L , zG_R , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

zG_L = FMO%L( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 
zG_R = FMO%R( : , orbital(1:n_part) )    ! <== expansion coefficients at t = 0 

it     =  1
t      =  t_i              
t_rate =  MD_dt * frame_step

Print 56 , initial_state     ! <== initial state of the isolated molecule 

! DUAL representation for efficient calculation of survival probabilities ...
! coefs of <k(t)| in DUAL basis ...
CALL gemm(UNI%L,ZG_L,DUAL_bra,'T','N',one,zero)

! coefs of |k(t)> in DUAL basis ...
CALL gemm(UNI%R,ZG_R,DUAL_ket,'N','N',one,zero)

QDyn(it,:) = get_Populations(system,basis,DUAL_bra,DUAL_ket)

t  = t  + t_rate
it = it + 1

done = .true. 

include 'formats.h'

end subroutine preprocess_Huckel_Slice
!
!
!
!
end module TimeSlice_m

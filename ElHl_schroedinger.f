 module Schroedinger_m

 use type_m
 use constants_m
 use f95_precision
 use blas95
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube ,          &
                                         GaussianCube_step ,  DP_Moment , electron_state ,  &
                                         Coulomb_ , restart , DensityMatrix , CT_dump_step
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use FMO_m                      , only : FMO_analysis , orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations , Net_Charge
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format
 use Backup_m                   , only : Security_Copy , Restart_state
 use Auto_Correlation_m         , only : MO_Occupation


    public :: Simple_dynamics , DeAllocate_QDyn

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:)   :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Real*8     , ALLOCATABLE , dimension(:,:,:) :: Pops(:,:,:)
    type(R_eigen)                               :: UNI , el_FMO , hl_FMO

 contains
!
!
!=====================================================
 subroutine Simple_dynamics(system, basis, UNI, QDyn )
!=====================================================
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(inout) :: basis(:)
 type(R_eigen)   , intent(in)    :: UNI
 type(f_time)    , intent(inout) :: QDyn

! local variables ...
integer                          :: j , nn , mm
integer                          :: it , n , it_init
real*8                           :: t , t_rate
real*8                           :: Total_DP(3)
complex*16      , ALLOCATABLE    :: phase(:)
character(11)                    :: argument

! ------------------ preprocess stuff --------------------

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) )

mm = size(basis) ; nn = n_part

CALL Allocate_Brackets( mm , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

t  = t_i
it = 1

!   building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
do n = 1 , n_part
    select case( eh_tag(n) )

        case( "el" )

            CALL FMO_analysis ( system , basis , UNI%R , el_FMO , instance="E" )

            MO_bra( : , n ) = el_FMO%L( : , orbital(n) )
            MO_ket( : , n ) = el_FMO%R( : , orbital(n) )

            Print 591, orbital(n) , el_FMO%erg(orbital(n))

        case( "hl" )

            CALL FMO_analysis ( system , basis , UNI%R , hl_FMO , instance="H" )

            MO_bra( : , n ) = hl_FMO%L( : , orbital(n) )
            MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )

            Print 592, orbital(n) , hl_FMO%erg(orbital(n))
            If( (orbital(n) > hl_FMO%Fermi_State) ) write(*,"(/a)") '>>> warning: hole state above the Fermi level <<<'

        end select
end do

! deallocate after use ...
deallocate( el_FMO%L , el_FMO%R , el_FMO%erg , hl_FMO%L , hl_FMO%R , hl_FMO%erg )

! DUAL representation for efficient calculation of survival probabilities ...
CALL DZgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%R , mm , MO_ket , mm , C_zero , DUAL_ket , mm )
CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_bra , mm , C_zero , DUAL_bra , mm )

! save populations ...
Pops(1,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

QDyn%dyn(1,:,:) = Pops(1,:,:)
CALL dump_QDyn( QDyn , 1 )

If( DensityMatrix ) then
    If( n_part == 1 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI )
    If( n_part == 2 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI, UNI )
End If

!   save the initial GaussianCube file ...
If( GaussianCube ) then

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

    do n = 1 , n_part
        if( eh_tag(n) == "XX" ) cycle
        CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
    end do

end If
!-------------------------------------------------------------
!                       Q-DYNAMICS

it_init = it + 1

t_rate = (t_f - t_i) / float(n_t)

DO it = it_init , n_t

    t = t + t_rate

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    forall( j=1:n_part)
        MO_bra(:,j) = merge( conjg(phase(:)) * MO_bra(:,j) , C_zero , eh_tag(j) /= "XX" )
        MO_ket(:,j) = merge(       phase(:)  * MO_ket(:,j) , C_zero , eh_tag(j) /= "XX" )
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    CALL dzgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%R , mm , MO_ket , mm , C_zero , DUAL_ket , mm )
    CALL dzgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_bra , mm , C_zero , DUAL_bra , mm )

    ! save populations ...
    Pops(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

    QDyn%dyn(it,:,:) = Pops(it,:,:)
    if( mod(it,CT_dump_step) == 0 ) CALL dump_QDyn( QDyn , it )

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then

        do n = 1 , n_part
            if( eh_tag(n) == "XX" ) cycle
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
        end do

    end If

    if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

END DO

! sum population dynamics over frames ...
QDyn%dyn = D_zero
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

include 'formats.h'

end subroutine Simple_dynamics
!
!
!
!
!========================================
 subroutine dump_Qdyn( Qdyn , it )
!========================================
implicit none
type(f_time)    , intent(in) :: QDyn
integer         , intent(in) :: it

! local variables ...
integer     :: nf , n
complex*16  :: wp_energy

do n = 1 , n_part

    if( eh_tag(n) == "XX" ) cycle

    wp_energy = sum(MO_bra(:,n) * UNI%erg(:) * MO_ket(:,n))

    If( it == 1 ) then

        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
        write(52,15) "#" ,( nf+1 , nf=0,size(QDyn%fragments)+1 )  ! <== numbered columns for your eyes only ...
        write(52,12) "#" , QDyn%fragments , "total"

        open( unit = 53 , file = "tmp_data/"//eh_tag(n)//"_wp_energy.dat" , status = "replace" , action = "write" , position = "append" )

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

12 FORMAT(/15A10)
13 FORMAT(F11.6,14F10.5)
14 FORMAT(3F12.6)
15 FORMAT(A,I9,14I10)

end subroutine dump_Qdyn
!
!
!
!
!=========================================
 subroutine DeAllocate_QDyn( QDyn , flag )
!=========================================
implicit none
type(f_time)  , intent(inout) :: QDyn
character(*)  , intent(in)    :: flag

! local variable ...
integer      :: N_of_fragments
character(1) :: first_in_line
logical      :: A_flag

select case( flag )

    case( "alloc" )

        if( allocated(trj) ) then

            CALL Coords_from_Universe( Unit_Cell, trj(2) )          ! <== use number 2 to avoid verbose
            CALL Generate_Structure( 2 )
            N_of_fragments = size( Extended_Cell%list_of_fragments )

        else

            CALL Generate_Structure( 2 )                            ! <== use number 2 to avoid verbose
            N_of_fragments = size( Extended_Cell%list_of_fragments )

        end if

        ! for the sake of having the DONOR or ACCEPTOR survival probability in the first column at output ...
        A_flag = any(Extended_Cell%list_of_fragments == "A")
        first_in_line = Extended_Cell%list_of_fragments(1)
        If( A_flag ) then
            where( Extended_Cell%list_of_fragments == "A" ) Extended_Cell%list_of_fragments = first_in_line
        else
            where( Extended_Cell%list_of_fragments == "D" ) Extended_Cell%list_of_fragments = first_in_line
        end If
        Extended_Cell%list_of_fragments(1) = merge( "A" , "D" , A_flag )

        ! QDyn%dyn = ( time ; fragments ; all fragments ) ...
        allocate( QDyn%fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )
        allocate( QDyn%dyn      ( n_t , 0:N_of_fragments+1 , n_part       ) , source = 0.d0                              )

        ! allocatating Net_Charte for future use ...
        allocate( Net_Charge(Extended_Cell%atoms) , source = D_zero )

        ! cleaning the mess ...
        CALL DeAllocate_Structures( Extended_Cell )

    case( "dealloc" )

        deallocate( QDyn%dyn , QDyn%fragments )

    case( "update" )  ! <== used for saving populations of atomic orbitals ...

        ! start re-building ...
        deallocate( Qdyn%fragments , Qdyn%dyn )

        N_of_fragments = size( Extended_Cell%list_of_fragments )

        ! for the sake of having the DONOR or ACCEPTOR survival probability in the first column at output ...
        A_flag = any(Extended_Cell%list_of_fragments == "A")
        first_in_line = Extended_Cell%list_of_fragments(1)
        If( A_flag ) then
            where( Extended_Cell%list_of_fragments == "A" ) Extended_Cell%list_of_fragments = first_in_line
        else
            where( Extended_Cell%list_of_fragments == "D" ) Extended_Cell%list_of_fragments = first_in_line
        end If
        Extended_Cell%list_of_fragments(1) = merge( "A" , "D" , A_flag )

        ! QDyn%dyn = ( time ; fragments ; all fragments ) ...
        allocate( QDyn%fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )
        allocate( QDyn%dyn      ( n_t , 0:N_of_fragments+1 , n_part       ) , source = 0.d0                              )

end select

end subroutine DeAllocate_QDyn
!
!
!
end module Schroedinger_m

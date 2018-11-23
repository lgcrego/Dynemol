 module Schroedinger_m

 use type_m
 use constants_m
 use f95_precision
 use blas95
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube ,          &
                                         GaussianCube_step ,  DP_Moment , initial_state ,   &
                                         Coulomb_ , restart , DensityMatrix
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use QCModel_Huckel_ElHl        , only : EigenSystem_ElHl
 use FMO_m                      , only : orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations , Net_Charge
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format
 use PDOS_tool_m                , only : Partial_DOS
 use Backup_m                   , only : Security_Copy , Restart_state
 use Auto_Correlation_m         , only : MO_Occupation


    public :: ElHl_dynamics , Huckel_dynamics , DeAllocate_QDyn

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:)   :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Real*8     , ALLOCATABLE , dimension(:,:,:) :: Pops(:,:,:)
    type(R_eigen)                               :: UNI_el , UNI_hl
    logical                                     :: static_hole = .false.

 contains
!
!
!=====================================================================
 subroutine ElHl_dynamics(system, basis, UNI, el_FMO , hl_FMO , QDyn )
!=====================================================================
 implicit none
 type(structure) , intent(inout)    :: system
 type(STO_basis) , intent(in)       :: basis(:)
 type(R_eigen)   , intent(in)       :: UNI
 type(R_eigen)   , intent(inout)    :: el_FMO
 type(R_eigen)   , intent(inout)    :: hl_FMO
 type(f_time)    , intent(inout)    :: QDyn

! local variables ...
integer                             :: mm
integer                             :: it , n , it_init
real*8                              :: t , t_rate
real*8                              :: Total_DP(3)
complex*16      , ALLOCATABLE       :: phase(:,:)
character(11)                       :: argument

! ------------------ preprocess stuff --------------------

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) )

Print 56 , initial_state     ! <== initial state of the isolated molecule

CALL Allocate_Brackets( size(basis) , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

mm = size(basis)

! at t=t_i UNI = UNI_el = UNI_hl ...
allocate( UNI_el%L(mm,mm) , UNI_el%R(mm,mm) , UNI_el%erg(mm) )
allocate( UNI_hl%L(mm,mm) , UNI_hl%R(mm,mm) , UNI_hl%erg(mm) )

UNI_el = UNI
UNI_hl = UNI

! Restart or build up wavepackets ...
if( restart .AND. Coulomb_ ) then

    deallocate( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket )

    CALL Restart_State( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it )

    Deallocate ( UNI_el%R , UNI_el%L , UNI_el%erg )
    Deallocate ( UNI_hl%R , UNI_hl%L , UNI_hl%erg )

    CALL EigenSystem_ElHl( system , basis , AO_bra , AO_ket , UNI_el , UNI_hl )

    it_init = it

else

    t  = t_i
    it = 1

!   building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
    do n = 1 , n_part
        select case( eh_tag(n) )

            case( "el" )

                MO_bra( : , n ) = el_FMO%L( : , orbital(n) )
                MO_ket( : , n ) = el_FMO%R( : , orbital(n) )

                Print 591, orbital(n) , el_FMO%erg(orbital(n))

            case( "hl" )

                    If( (orbital(n) > hl_FMO%Fermi_State) ) write(*,"(/a)") '>>> warning: hole state above the Fermi level <<<'

                MO_bra( : , n ) = hl_FMO%L( : , orbital(n) )
                MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )

                Print 592, orbital(n) , hl_FMO%erg(orbital(n))

            end select
    end do

!   deallocate after use ...
    deallocate( el_FMO%L , el_FMO%R , el_FMO%erg , hl_FMO%L , hl_FMO%R , hl_FMO%erg )

!   DUAL representation for efficient calculation of survival probabilities ...
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_bra(:,1) , mm , C_zero , DUAL_bra(:,1) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , MO_ket(:,1) , mm , C_zero , DUAL_ket(:,1) , mm )

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_bra(:,2) , mm , C_zero , DUAL_bra(:,2) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , MO_ket(:,2) , mm , C_zero , DUAL_ket(:,2) , mm )

!   save populations ...
    Pops(1,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

    QDyn%dyn(1,:,:) = Pops(1,:,:)
    CALL dump_QDyn( QDyn , 1 )

    If( DensityMatrix ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI_el, UNI_hl )

!   save the initial GaussianCube file ...
    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then

        ! LOCAL representation for film STO production ...
        AO_bra = DUAL_bra

        CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_ket(:,1) , mm , C_zero , AO_ket(:,1) , mm )
        CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_ket(:,2) , mm , C_zero , AO_ket(:,2) , mm )

        do n = 1 , n_part
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
        end do

    end If

    it_init = it + 1

end if

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
!                       Q-DYNAMICS

t_rate = (t_f - t_i) / float(n_t)

DO it = it_init , n_t

    t = t + t_rate

    phase(:,1) = cdexp(- zi * UNI_el%erg(:) * t_rate / h_bar)
    phase(:,2) = cdexp(- zi * UNI_hl%erg(:) * t_rate / h_bar)

    if( static_hole ) phase(:,2) = (1.d0,0.d0)

    forall( n=1:n_part)
        MO_bra(:,n) = conjg(phase(:,n)) * MO_bra(:,n)
        MO_ket(:,n) =       phase(:,n)  * MO_ket(:,n)
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_bra(:,1) , mm , C_zero , DUAL_bra(:,1) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , MO_ket(:,1) , mm , C_zero , DUAL_ket(:,1) , mm )

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_bra(:,2) , mm , C_zero , DUAL_bra(:,2) , mm )
    CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , MO_ket(:,2) , mm , C_zero , DUAL_ket(:,2) , mm )

    ! save populations ...
    Pops(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

    QDyn%dyn(it,:,:) = Pops(it,:,:)
    CALL dump_QDyn( QDyn , it )

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra

    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_ket(:,1) , mm , C_zero , AO_ket(:,1) , mm )
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_ket(:,2) , mm , C_zero , AO_ket(:,2) , mm )

    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then

        do n = 1 , n_part
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
        end do

    end If

    if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

    If( Coulomb_ ) then

        ! saving a little memory space ...
        Deallocate ( UNI_el%R , UNI_el%L , UNI_el%erg )
        Deallocate ( UNI_hl%R , UNI_hl%L , UNI_hl%erg )

        CALL EigenSystem_ElHl( system , basis , AO_bra , AO_ket , UNI_el , UNI_hl )

        ! project back to MO_basis with UNI(t + t_rate)
        CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , Dual_bra(:,1) , mm , C_zero , MO_bra(:,1) , mm )
        CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , Dual_ket(:,1) , mm , C_zero , MO_ket(:,1) , mm )

        CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , Dual_bra(:,2) , mm , C_zero , MO_bra(:,2) , mm )
        CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , Dual_ket(:,2) , mm , C_zero , MO_ket(:,2) , mm )

        CALL Security_Copy( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it )

        If( DensityMatrix ) CALL MO_Occupation( t, MO_bra, MO_ket, UNI_el, UNI_hl )

    end If

END DO

! sum population dynamics over frames ...
QDyn%dyn = D_zero
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

include 'formats.h'

end subroutine ElHl_dynamics
!
!
!
!==============================================================
 subroutine Huckel_dynamics(system, basis, UNI, el_FMO , QDyn )
!==============================================================
 implicit none
 type(structure) , intent(inout)    :: system
 type(STO_basis) , intent(in)       :: basis(:)
 type(R_eigen)   , intent(in)       :: UNI
 type(R_eigen)   , intent(inout)    :: el_FMO
 type(f_time)    , intent(inout)    :: QDyn

! local variables ...
integer                             :: mm
integer                             :: it , n
real*8                              :: t , t_rate
real*8                              :: Total_DP(3)
complex*16      , ALLOCATABLE       :: phase(:)

! ------------------ preprocess stuff --------------------

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) )

Print 56 , initial_state     ! <== initial state of the isolated molecule

CALL Allocate_Brackets( size(basis) , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

mm = size(basis)

MO_bra( : , 1 ) = el_FMO%L( : , orbital(1) )
MO_ket( : , 1 ) = el_FMO%R( : , orbital(1) )

Print 591, orbital(1) , el_FMO%erg(orbital(1))

! deallocate after use ...
deallocate( el_FMO%L , el_FMO%R , el_FMO%erg )

!---------------------------------------------------------
!                    Q-DYNAMICS

t = t_i

t_rate = (t_f - t_i) / float(n_t)

DO it = 1 , n_t

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    If( t == t_i ) phase = C_one

    forall( n=1:n_part)
        MO_bra(:,n) = conjg(phase(:)) * MO_bra(:,n)
        MO_ket(:,n) =       phase(:)  * MO_ket(:,n)
    end forall

!--------------------------------------------------------------------------
! LOCAL representation for film STO production ...
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI%L , mm , MO_bra , mm , C_zero , AO_bra , mm )
    CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

    if( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then
        do n = 1 , n_part
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it , t , eh_tag(n) )
        end do
    end if
!--------------------------------------------------------------------------
! DUAL representation for efficient calculation of survival probabilities ...
   DUAL_bra = AO_bra
   CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI%R , mm , MO_ket , mm , C_zero , DUAL_ket , mm )

   Pops(it,:,1) = Populations( QDyn%fragments , basis , DUAL_bra(:,1) , DUAL_ket(:,1) , t )

   if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

   t = t + t_rate

END DO

! sum population dynamics over frames ...
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

include 'formats.h'

end subroutine Huckel_dynamics
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
13 FORMAT(F11.6,14F10.5)
14 FORMAT(3F12.6)

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
logical      :: E_flag

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

        ! for the sake of having the DONOR or EXCITON survival probability in the first column at output ...
        E_flag = any(Extended_Cell%list_of_fragments == "E")
        first_in_line = Extended_Cell%list_of_fragments(1)
        If( E_flag ) then
            where( Extended_Cell%list_of_fragments == "E" ) Extended_Cell%list_of_fragments = first_in_line
        else
            where( Extended_Cell%list_of_fragments == "D" ) Extended_Cell%list_of_fragments = first_in_line
        end If
        Extended_Cell%list_of_fragments(1) = merge( "E" , "D" , E_flag )

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

        ! for the sake of having the DONOR or EXCITON survival probability in the first column at output ...
        E_flag = any(Extended_Cell%list_of_fragments == "E")
        first_in_line = Extended_Cell%list_of_fragments(1)
        If( E_flag ) then
            where( Extended_Cell%list_of_fragments == "E" ) Extended_Cell%list_of_fragments = first_in_line
        else
            where( Extended_Cell%list_of_fragments == "D" ) Extended_Cell%list_of_fragments = first_in_line
        end If
        Extended_Cell%list_of_fragments(1) = merge( "E" , "D" , E_flag )

        ! QDyn%dyn = ( time ; fragments ; all fragments ) ...
        allocate( QDyn%fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )
        allocate( QDyn%dyn      ( n_t , 0:N_of_fragments+1 , n_part       ) , source = 0.d0                              )

end select

end subroutine DeAllocate_QDyn
!
!
!
end module Schroedinger_m

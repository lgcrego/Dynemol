 module Schroedinger_m

 use type_m
 use constants_m
 use mkl95_precision
 use mkl95_blas
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube ,          &
                                         GaussianCube_step ,  DP_Moment , initial_state ,   &
                                         Coulomb_
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use QCModel_Huckel_ElHl        , only : EigenSystem_ElHl 
 use FMO_m                      , only : orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format
 use PDOS_tool_m                , only : Partial_DOS


    public :: ElHl_dynamics , Huckel_dynamics , DeAllocate_QDyn 

    private

    ! module variables ...
    Complex*16 , ALLOCATABLE , dimension(:,:)   :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra 
    Complex*16 , ALLOCATABLE , dimension(:)     :: bra , ket
    Real*8     , ALLOCATABLE , dimension(:,:,:) :: Pops(:,:,:)

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
 type(R_eigen)   , intent(in)       :: el_FMO 
 type(R_eigen)   , intent(in)       :: hl_FMO 
 type(f_time)    , intent(inout)    :: QDyn

! local variables ...
integer                             :: mm 
integer                             :: it , i , n 
real*8                              :: t , t_rate 
real*8                              :: Total_DP(3)
complex*16      , ALLOCATABLE       :: phase(:,:)
type(R_eigen)                       :: UNI_el , UNI_hl

real*8 :: E_el , E_hl

! ------------------ preprocess stuff --------------------

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) ) 

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(basis) , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

mm = size(basis)                          

! at t=t_i UNI = UNI_el = UNI_hl ...                        
allocate( UNI_el%L(mm,mm) , UNI_el%R(mm,mm) , UNI_el%erg(mm) )
allocate( UNI_hl%L(mm,mm) , UNI_hl%R(mm,mm) , UNI_hl%erg(mm) )

UNI_el = UNI
UNI_hl = UNI

! building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
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

! deallocate after use ...
deallocate( el_FMO%L , el_FMO%R , el_FMO%erg , hl_FMO%L , hl_FMO%R , hl_FMO%erg )

! DUAL representation for efficient calculation of survival probabilities ...
CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_el%L , mm , MO_bra(:,1) , mm , C_zero , DUAL_bra(:,1) , mm )
CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_el%R , mm , MO_ket(:,1) , mm , C_zero , DUAL_ket(:,1) , mm )

CALL DZgemm( 'T' , 'N' , mm , 1 , mm , C_one , UNI_hl%L , mm , MO_bra(:,2) , mm , C_zero , DUAL_bra(:,2) , mm )
CALL DZgemm( 'N' , 'N' , mm , 1 , mm , C_one , UNI_hl%R , mm , MO_ket(:,2) , mm , C_zero , DUAL_ket(:,2) , mm )

! save populations ...
Pops(1,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

CALL dump_QDyn( QDyn , 1 )

!-------------------------------------------------------------
!                       Q-DYNAMICS  

t = t_i              

t_rate = (t_f - t_i) / float(n_t)
    
DO it = 2 , n_t    

    t = t + t_rate

    phase(:,1) = cdexp(- zi * UNI_el%erg(:) * t_rate / h_bar)
    phase(:,2) = cdexp(- zi * UNI_hl%erg(:) * t_rate / h_bar)

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

!    if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

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

    end If

END DO

! sum population dynamics over frames ...
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

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
 type(R_eigen)   , intent(in)       :: el_FMO 
 type(f_time)    , intent(inout)    :: QDyn

! local variables ...
integer                             :: mm 
integer                             :: it , i , n 
real*8                              :: t , t_rate 
real*8                              :: Total_DP(3)
complex*16      , ALLOCATABLE       :: phase(:)

! ------------------ preprocess stuff --------------------

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) ) 

Print 56 , initial_state     ! <== initial state of the isolated molecule 
 
CALL Allocate_Brackets( size(basis) , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase ) 

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

   Pops(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

   if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

   t = t + t_rate

END DO

! sum population dynamics over frames ...
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , bra , ket , phase )

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
integer :: nf , n

do n = 1 , n_part

    If( it == 1 ) then
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
        write(52,12) "#" , QDyn%fragments , "total"
    else
        open( unit = 52 , file = "tmp_data/"//eh_tag(n)//"_survival.dat" , status = "unknown", action = "write" , position = "append" )
    end If

    write(52,13) ( Pops(it,nf,n) , nf=0,size(QDyn%fragments)+1 ) 

    close(52)

end do

12 FORMAT(10A9)
13 FORMAT(10F10.5)

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
integer      :: i , N_of_fragments
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

        ! cleaning the mess ...
        CALL DeAllocate_Structures( Extended_Cell )

    case( "dealloc" )

        deallocate( QDyn%dyn , QDyn%fragments )

end select

end subroutine DeAllocate_QDyn
!
!
!
! 
end module Schroedinger_m

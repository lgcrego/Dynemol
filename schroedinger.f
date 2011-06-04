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
 use FMO_m                      , only : orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format
 use Coulomb_m                  , only : wormhole_to_Coulomb , Build_Coulomb_potential

 public :: Huckel_dynamics , DeAllocate_QDyn

 private

 contains
!
! 
!=======================================================================
 subroutine Huckel_dynamics(system, basis, UNI, el_FMO , hl_FMO , QDyn )
!=======================================================================
 implicit none
 type(structure)            , intent(inout)    :: system
 type(STO_basis)            , intent(in)       :: basis(:)
 type(C_eigen)              , intent(in)       :: UNI 
 type(C_eigen)              , intent(in)       :: el_FMO 
 type(C_eigen)   , optional , intent(in)       :: hl_FMO 
 type(f_time)               , intent(inout)    :: QDyn

! local variables ...
integer                             :: it , i , n 
real*8                              :: t , t_rate 
real*8                              :: Total_DP(3)
real*8          , ALLOCATABLE       :: Pops(:,:,:)
complex*16      , ALLOCATABLE       :: MO_bra(:,:)   , MO_ket(:,:)
complex*16      , ALLOCATABLE       :: AO_bra(:,:)   , AO_ket(:,:) 
complex*16      , ALLOCATABLE       :: DUAL_ket(:,:) , DUAL_bra(:,:) 
complex*16      , ALLOCATABLE       :: phase(:)      , bra(:)        , ket(:)

! preprocessing stuff ..........................................................

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) ) 

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

!...............................................................................

!=========================================================
!                       Q-DYNAMICS  
!=========================================================
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

! coefs of <k(t)| in AO basis 
    CALL gemm(UNI%L,MO_bra,AO_bra,'T','N',C_one,C_zero)

! coefs of |k(t)> in AO basis 
    CALL gemm(UNI%L,MO_ket,AO_ket,'T','N',C_one,C_zero)

    if( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then
        do n = 1 , n_part
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it , t , eh_tag(n) )
        end do
    end if 

    CALL Coulomb_stuff( AO_bra , AO_ket , basis )

!--------------------------------------------------------------------------
! DUAL representation for efficient calculation of survival probabilities ...

! coefs of <k(t)| in DUAL basis ...
   DUAL_bra = AO_bra

! coefs of |k(t)> in DUAL basis ...
   CALL gemm(UNI%R,MO_ket,DUAL_ket,'N','N',C_one,C_zero)

   Pops(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )

   if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

   t = t + t_rate

END DO

! sum population dynamics over frames ...
QDyn%dyn = QDyn%dyn + Pops

deallocate( Pops )

include 'formats.h'

end subroutine Huckel_dynamics
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
!===================================================
 subroutine Coulomb_stuff( AO_bra , AO_ket , basis )
!===================================================
implicit none
complex*16      , intent(in) :: AO_bra(:,:) , AO_ket(:,:) 
type(STO_basis) , intent(in) :: basis(:)

! save coefficients in modulo Coulomb_m ...
CALL wormhole_to_Coulomb( AO_bra , AO_ket )

end subroutine Coulomb_stuff
!
!
!
end module Schroedinger_m

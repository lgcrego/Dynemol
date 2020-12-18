! Subroutine for computing average properties over nuclear configurations ...

module Sampling_m

    use type_m
    use constants_m
    use parameters_m         , only : frame_step , spectrum ,        &
                                      survival , nuclear_matter ,    &
                                      DP_Moment , EnvField_ ,        &
                                      file_type , sigma , n_part ,   &
                                      restart , SOC
    use Babel_m              , only : System_Characteristics ,       &
                                      Coords_from_Universe ,         &
                                      trj
    use Allocation_m         , only : Allocate_UnitCell ,            &
                                      DeAllocate_UnitCell ,          &
                                      DeAllocate_Structures 
    use Solvated_M           , only : Prepare_Solvated_System ,      &
                                      DeAllocate_TDOS ,              &
                                      DeAllocate_PDOS ,              &
                                      DeAllocate_SPEC 
    use QCModel_Huckel       , only : EigenSystem
    use FMO_m                , only : FMO_analysis , eh_tag
    use Structure_Builder    , only : Unit_Cell ,                    &
                                      Extended_Cell ,                &
                                      Generate_Structure ,           &
                                      Basis_Builder ,                &
                                      ExCell_basis
    use DOS_m
    use Oscillator_m         , only : Optical_Transitions
    use DP_main_m            , only : Dipole_Matrix 
    use Dielectric_Potential , only : Environment_SetUp
    use Schroedinger_m       , only : Simple_dynamics ,              &
                                      DeAllocate_QDyn ,              &
                                      RunningStat 
    use Data_Output          , only : Dump_stuff , FileName


    public :: Avrg_Confgs 

    private

contains
!
!
!=======================
 subroutine Avrg_Confgs
!=======================
 implicit none

! local variables ...
integer                         :: frame , frame_init , nr , N_of_residues
real*8                          :: internal_sigma
logical                         :: FMO_ , DIPOLE_
type(R_eigen)                   :: UNI 
type(f_grid)                    :: TDOS , SPEC
type(f_grid)    , allocatable   :: PDOS(:) 
type(f_time)                    :: QDyn
type(universe)                  :: Solvated_System

! preprocessing stuff .....................................................

FMO_    = ( spectrum .OR. survival  )
DIPOLE_ = ( spectrum .OR. DP_Moment )

internal_sigma = sigma / ( float(size(trj))/float(frame_step) )      

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!..........................................................................

! Average over samples : Quantum Dynamics & All that Jazz ...

if( restart ) then
    CALL restart_AVRG( frame , TDOS , PDOS , SPEC )
    frame_init = frame + frame_step
else
    frame_init = 1
end if

do frame = frame_init , size(trj) , frame_step

    select case ( nuclear_matter )

        case( "solvated_sys" )

            CALL Prepare_Solvated_System( Solvated_System , frame )

            CALL Coords_from_Universe( Unit_Cell , Solvated_System )

        case( "extended_sys" )

            CALL Coords_from_Universe( Unit_Cell , trj(frame) )

        case default

            Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
            stop

    end select

    CALL Generate_Structure( frame )

    CALL Basis_Builder( Extended_Cell , ExCell_basis )

    If( EnvField_ ) CALL Environment_SetUp( Extended_Cell )

    CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

    ! TDOS += avrg
    CALL Total_DOS( UNI , ExCell_basis , TDOS , internal_sigma )                                             

    do nr = 1 , N_of_residues
        ! PDOS += avrg
        CALL Partial_DOS( Extended_Cell , ExCell_basis , UNI , PDOS , nr , internal_sigma )            
    end do

    If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )

    ! SPEC += avrg
    If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC , internal_sigma )

    If( survival ) CALL Simple_dynamics( Extended_Cell, ExCell_basis, UNI, QDyn )

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    ! temporary data storage [frame in trj] ...
    CALL dump_stuff( TDOS , PDOS , SPEC )

    print*, frame

end do

! finish calculation of STD ...
CALL RunningStat( QDyn , instance="getSTD" ) 

! final data storage, replacing temporary data ...
! histograms:{TDOS , PDOS , SPEC} , average survival:(QDyn} ...
CALL Dump_stuff( TDOS , PDOS , SPEC , QDyn )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

end subroutine Avrg_Confgs
!
!
!
!===========================================================
subroutine restart_AVRG( frame , TDOS , PDOS , SPEC , QDyn )
!===========================================================
implicit none
integer                     , intent(out)   :: frame
type(f_grid)    , optional  , intent(out)   :: TDOS
type(f_grid)    , optional  , intent(out)   :: PDOS(:)
type(f_grid)    , optional  , intent(out)   :: SPEC
type(f_time)    , optional  , intent(out)   :: QDyn

! local variables ...
integer          :: i , nr , nf , np , N_of_residues , N_of_fragments , n_spin , spin
character(26)    :: string
character(len=:) ,  allocatable :: f_name

! save TDOS ...
If( present(TDOS) ) then
    OPEN( unit = 33 , file = "dyn.trunk/TDOS.dat" , status = "unknown" , action = "read" )
    read(33,14) frame
    do i = 1 , size(TDOS%func)
        read(33,10) TDOS%grid(i) , TDOS%average(i) , TDOS%peaks(i) ,  TDOS%occupation(i)
    end do
    CLOSE(33)
end if

! save PDOS ...
If( present(PDOS) ) then
    N_of_residues = size( PDOS )
    do nr = 1 , N_of_residues
        string = "dyn.trunk/PDOS-"//PDOS(nr)%residue//".dat"
        OPEN( unit = 33 , file=string , status='unknown' )
        read(33,14) frame
            do i = 1 , size(PDOS(nr)%func)
                read(33,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i) , PDOS(nr)%peaks(i) , PDOS(nr)%occupation(i)
            end do
        CLOSE(33)
    end do
end if

! save peak and broadened specs ...
If( present(SPEC) ) then
    OPEN( unit = 33 , file='dyn.trunk/spectrum.dat' , status='unknown' )
        read(33,14) frame
        do i = 1 , size(SPEC%func)
            read(33,11) SPEC%grid(i) , SPEC%average(i) , SPEC%peaks(i)
        end do
    CLOSE(33)
end if

! save time-dependent electron or hole populations ...
If( present(QDyn) ) then
   N_of_fragments = size( QDyn%fragments )
   n_spin         = merge(2,1,SOC)
   
   do spin = 1 , n_spin
   
      do np = 1 , n_part
   
           If( eh_tag(np) == "XX" ) cycle
   
           call FileName( f_name , np , spin , instance="dens" )
           open( unit = 33 , file = f_name , status = "unknown" )
               read(33,14) frame           
               read(33,12) QDyn%fragments
               DO i = 1 , size( QDyn%dyn(:,1,1,1) )
                   read(33,13) ( QDyn%dyn(i,nf,np,spin) , nf=0,N_of_fragments+1 )
               end do
           close(33)

           end do
           end do

           end if

10   FORMAT(4F12.5)
11   FORMAT(3F13.9)
12   FORMAT(10A10)
13   FORMAT(F11.6,9F10.5)
14   FORMAT(I4)

end subroutine restart_AVRG
!
!
!
end module Sampling_m

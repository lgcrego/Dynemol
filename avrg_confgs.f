! Subroutine for computing average properties over nuclear configurations ...

module Sampling_m

    use type_m
    use constants_m
    use parameters_m        , only : frame_step , spectrum ,        &
                                     survival , nuclear_matter ,    &
                                     DP_Moment , DP_Field_ ,        &
                                     file_type , sigma , n_part ,   &
                                     restart
    use Babel_m             , only : System_Characteristics ,       &
                                     Coords_from_Universe ,         &
                                     trj
    use Allocation_m        , only : Allocate_UnitCell ,            &
                                     DeAllocate_UnitCell ,          &
                                     DeAllocate_Structures 
    use Solvated_M          , only : Prepare_Solvated_System ,      &
                                     DeAllocate_TDOS ,              &
                                     DeAllocate_PDOS ,              &
                                     DeAllocate_SPEC 
    use QCModel_Huckel      , only : EigenSystem
    use FMO_m               , only : FMO_analysis , eh_tag
    use Structure_Builder   , only : Unit_Cell ,                    &
                                     Extended_Cell ,                &
                                     Generate_Structure ,           &
                                     Basis_Builder ,                &
                                     ExCell_basis
    use DOS_m
    use Oscillator_m        , only : Optical_Transitions
    use DP_main_m           , only : Dipole_Matrix 
    use DP_potential_m      , only : Molecular_DPs
    use Schroedinger_m      , only : Huckel_dynamics ,              &
                                     DeAllocate_QDyn
    use Data_Output         , only : Dump_stuff


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
type(R_eigen)                   :: UNI , FMO
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

    If( DP_field_ ) CALL Molecular_DPs( Extended_Cell )

    CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

    CALL Total_DOS( UNI%erg , TDOS , internal_sigma )                                             

    do nr = 1 , N_of_residues
        CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr , internal_sigma )            
    end do

    If( FMO_     ) CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, FMO , instance="E" )

    If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )

    If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC , internal_sigma )

    If( survival ) CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, FMO , QDyn=QDyn )

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    CALL dump_QDyn( frame , TDOS , PDOS , SPEC )

    print*, frame

end do

! average over configurations ...
If( file_type == "trajectory" ) QDyn%dyn = QDyn%dyn / ( float(size(trj))/float(frame_step) )

CALL Dump_stuff( TDOS , PDOS , SPEC , QDyn )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )
CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

end subroutine Avrg_Confgs
!
!
!
!========================================================
subroutine dump_Qdyn( frame , TDOS , PDOS , SPEC , QDyn )
!========================================================
implicit none
integer                     , intent(in) :: frame
type(f_grid)    , optional  , intent(in) :: TDOS
type(f_grid)    , optional  , intent(in) :: PDOS(:)
type(f_grid)    , optional  , intent(in) :: SPEC
type(f_time)    , optional  , intent(in) :: QDyn

! local variables ...
integer         :: i , n , nr , nf , np , N_of_residues , N_of_fragments
character(26)   :: string

! save TDOS ...
If( present(TDOS) ) then
    OPEN( unit = 33 , file = "tmp_data/TDOS.dat" , status = "unknown" , action = "write" )
    write(33,14) frame
    do i = 1 , size(TDOS%func)
        write(33,10) TDOS%grid(i) , TDOS%average(i) , TDOS%peaks(i) ,  TDOS%occupation(i)
    end do
    CLOSE(33)
end if

! save PDOS ...
If( present(PDOS) ) then
    N_of_residues = size( PDOS )
    do nr = 1 , N_of_residues
        string = "tmp_data/PDOS-"//PDOS(nr)%residue//".dat"
        OPEN( unit = 33 , file=string , status='unknown' )
            write(33,14) frame
            do i = 1 , size(PDOS(nr)%func)
                write(33,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i) , PDOS(nr)%peaks(i) , PDOS(nr)%occupation(i)
            end do
        CLOSE(33)
    end do
end if

! save peak and broadened specs ...
If( present(SPEC) ) then
    OPEN( unit = 33 , file='tmp_data/spectrum.dat' , status='unknown' )
        write(33,14) frame
        do i = 1 , size(SPEC%func)
            write(33,11) SPEC%grid(i) , SPEC%average(i) , SPEC%peaks(i)
        end do
    CLOSE(33)
end if

! save time-dependent electron or hole populations ...
If( present(QDyn) ) then

    N_of_fragments = size( QDyn%fragments )

        do np = 1 , n_part

            OPEN( unit = 33 , file="tmp_data/"//eh_tag(n)//"_survival.dat" , status="unknown" )

            write(33,14) frame
            write(33,12) "#" , QDyn%fragments , "total"
            do i = 1 , size( QDyn%dyn(:,1,1) )
                write(33,13) ( QDyn%dyn(i,nf,np) , nf=0,N_of_fragments+1 )
            end do

            CLOSE(33)

        end do
end if

10   FORMAT(4F12.5)
11   FORMAT(3F13.9)
12   FORMAT(10A10)
13   FORMAT(F11.6,9F10.5)
14   FORMAT(I4)

end subroutine dump_Qdyn
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
integer         :: i , n , nr , nf , np , N_of_residues , N_of_fragments
character(26)   :: string

! save TDOS ...
If( present(TDOS) ) then
    OPEN( unit = 33 , file = "tmp_data/TDOS.dat" , status = "unknown" , action = "read" )
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
        string = "tmp_data/PDOS-"//PDOS(nr)%residue//".dat"
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
    OPEN( unit = 33 , file='tmp_data/spectrum.dat' , status='unknown' )
        read(33,14) frame
        do i = 1 , size(SPEC%func)
            read(33,11) SPEC%grid(i) , SPEC%average(i) , SPEC%peaks(i)
        end do
    CLOSE(33)
end if

! save time-dependent electron or hole populations ...
If( present(QDyn) ) then

    N_of_fragments = size( QDyn%fragments )

        do np = 1 , n_part

            OPEN( unit = 33 , file="tmp_data/"//eh_tag(n)//"_survival.dat" , status="unknown" )

            read(33,14) frame           
            read(33,12) QDyn%fragments
            do i = 1 , size( QDyn%dyn(:,1,1) )
                read(33,13) ( QDyn%dyn(i,nf,np) , nf=0,N_of_fragments+1 )
            end do

            CLOSE(33)

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

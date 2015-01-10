module Backup_m

    use type_m
    use mkl95_blas
    use parameters_m        , only : driver                     , &
                                     QMMM                       , &
                                     nuclear_matter             , &
                                     DP_field_                  , &
                                     DP_Moment                  , &
                                     Coulomb_                   , &
                                     restart
    use Solvated_M          , only : Prepare_Solvated_System
    use Babel_m             , only : Coords_from_Universe       , &
                                     trj
    use Structure_Builder   , only : Generate_Structure         , &
                                     Basis_Builder
    use FMO_m               , only : orbital                    , &
                                     eh_tag    
    use QCModel_Huckel      , only : EigenSystem
    use QCModel_Huckel_ElHl , only : EigenSystem_ElHl    
    use DP_potential_m      , only : Molecular_DPs
    use TD_Dipole_m         , only : wavepacket_DP
    use DP_main_m           , only : Dipole_Matrix   
    use MM_dynamics_m       , only : preprocess_MM              , &
                                     Saving_MM_Backup
    use Data_Output         , only : Net_Charge
                                     

    public  :: Security_Copy , Restart_State , Restart_Sys

contains
!
!
!
!==============================================================================================================================
 subroutine Restart_Sys( Extended_Cell , ExCell_basis , Unit_Cell , DUAL_ket , AO_bra , AO_ket , frame , it , UNI_el , UNI_hl )
!==============================================================================================================================
implicit none
type(structure)                 , intent(out)   :: Extended_Cell
type(STO_basis) , allocatable   , intent(out)   :: ExCell_basis(:)
type(structure)                 , intent(inout)   :: Unit_Cell
complex*16                      , intent(in)    :: DUAL_ket (:,:)
complex*16                      , intent(in)    :: AO_bra   (:,:)
complex*16                      , intent(in)    :: AO_ket   (:,:)
integer                         , intent(in)    :: frame
integer                         , intent(in)    :: it
type(R_eigen)                   , intent(out)   :: UNI_el
type(R_eigen)   , optional      , intent(out)   :: UNI_hl

! local variables ...
type(universe) :: Solvated_System

select case ( nuclear_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , frame )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

    case( "MDynamics" )

        ! MM preprocess ...
        CALL preprocess_MM( Net_Charge = Net_Charge )   

    case default

        Print*, " >>> Check your nuclear_matter options <<< :" , nuclear_matter
        stop

end select

CALL Generate_Structure ( frame )

CALL Basis_Builder      ( Extended_Cell , ExCell_basis )

if( DP_Moment ) &
CALL Dipole_Matrix      ( Extended_Cell , ExCell_basis )

if( DP_field_ ) then

    CALL Dipole_Matrix  ( Extended_Cell , ExCell_basis )

    CALL wavepacket_DP  ( Extended_Cell , ExCell_basis , AO_bra , AO_ket , Dual_ket )

    CALL Molecular_DPs  ( Extended_Cell )

end If

if( driver == "slice_ElHl") then

    CALL EigenSystem_ElHl( Extended_Cell , ExCell_basis , AO_bra , AO_ket , UNI_el , UNI_hl , flag2=it )

else

    CALL EigenSystem( Extended_Cell , ExCell_basis , UNI_el , flag2=it )

end if

end subroutine Restart_Sys
!
!                                                                    
!
!===================================================================================================
subroutine Security_Copy( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame )
!===================================================================================================
implicit none
complex*16              , intent(in)    :: MO_bra   (:,:)
complex*16              , intent(in)    :: MO_ket   (:,:)
complex*16              , intent(in)    :: DUAL_bra (:,:)
complex*16              , intent(in)    :: DUAL_ket (:,:)
complex*16              , intent(in)    :: AO_bra   (:,:)
complex*16              , intent(in)    :: AO_ket   (:,:)
real*8                  , intent(in)    :: t
integer                 , intent(in)    :: it
integer     , optional  , intent(in)    :: frame

! local variables ...
integer         :: i , j , basis_size , file_err , n_part
logical , save  :: first_time = .true. 

! check whether restart is properly set ...
if( (.NOT. restart) .AND. first_time ) then
    open(unit=3, file="Restart_copy.dat", status="new", iostat=file_err , err=11 )
    11 if( file_err > 0 ) stop " <Restart_copy.dat> exists; check restart parameter "
    close(3)
    CALL system( "rm Restart_copy.dat" )
    first_time = .false.
end if

if( QMMM ) CALL Saving_MM_Backup( frame , instance = "from_QM" )

open(unit=33, file="Security_copy.dat", status="unknown", form="unformatted", action="write")

if( present(frame) ) write(33) frame
write(33) it
write(33) t
write(33) size(MO_bra(:,1))
write(33) size(MO_bra(1,:))
write(33) size(eh_tag)

basis_size = size(MO_bra(:,1))
n_part     = size(MO_bra(1,:))

write(33) ( orbital(i) , eh_tag(i) , i=1,n_part )

do j = 1 , n_part

    write(33) ( MO_bra(i,j)   , MO_ket   (i,j) , i=1,basis_size )

    write(33) ( DUAL_bra(i,j) , DUAL_ket (i,j) , i=1,basis_size )

    write(33) ( AO_bra(i,j)   , AO_ket   (i,j) , i=1,basis_size )

end do

write(33) ( Net_Charge , i=1,size(Net_Charge) )

close( 33 )

! erase restart_copy file ...
CALL system( "rm Restart_copy.dat 2> qdynamo.err" )

end subroutine Security_Copy
!
!
!
!===================================================================================================
subroutine Restart_State( MO_bra , MO_ket , DUAL_bra , DUAL_ket , AO_bra , AO_ket , t , it , frame )
!===================================================================================================
implicit none
complex*16  , allocatable   , intent(out) :: MO_bra     (:,:)
complex*16  , allocatable   , intent(out) :: MO_ket     (:,:)
complex*16  , allocatable   , intent(out) :: DUAL_bra   (:,:)
complex*16  , allocatable   , intent(out) :: DUAL_ket   (:,:)
complex*16  , allocatable   , intent(out) :: AO_bra     (:,:)
complex*16  , allocatable   , intent(out) :: AO_ket     (:,:)
real*8                      , intent(out) :: t
integer                     , intent(out) :: it
integer     , optional      , intent(out) :: frame

! local variables ...
integer :: i , j , size_r , size_c , file_err , size_eh_tag

open(unit=33, file="Restart_copy.dat", form="unformatted", status="old", action="read" , iostat=file_err , err=11 )

if( present(frame) ) read(33) frame
read(33) it
read(33) t
read(33) size_r
read(33) size_c
read(33) size_eh_tag

allocate( MO_bra   ( size_r , size_c ) )
allocate( MO_ket   ( size_r , size_c ) )
allocate( DUAL_bra ( size_r , size_c ) )
allocate( DUAL_ket ( size_r , size_c ) )
allocate( AO_bra   ( size_r , size_c ) )
allocate( AO_ket   ( size_r , size_c ) )

if( .NOT. allocated( orbital) ) allocate( orbital(size_eh_tag) )
if( .NOT. allocated( eh_tag ) ) allocate( eh_tag(size_eh_tag) )

read(33) ( orbital(i) , eh_tag(i) , i=1,size_eh_tag )

do j = 1 , size_c

    read(33) ( MO_bra(i,j)   , MO_ket   (i,j) , i=1,size_r )

    read(33) ( DUAL_bra(i,j) , DUAL_ket (i,j) , i=1,size_r )

    read(33) ( AO_bra(i,j)   , AO_ket   (i,j) , i=1,size_r )

end do

read(33) ( Net_Charge , i=1,size(Net_Charge) )

close( 33 )

11 if( file_err > 0 ) stop " <Restart_copy.dat> file not found; terminating execution"

end subroutine Restart_State   
!
!
!
end module Backup_m

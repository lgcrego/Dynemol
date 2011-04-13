module Backup_m

    use type_m
    use mkl95_blas
    use parameters_m        , only : state_of_matter            , &
                                     DP_field_                  , &
                                     DP_Moment                  , &
                                     restart
    use Solvated_M          , only : Prepare_Solvated_System
    use Babel_m             , only : Coords_from_Universe       , &
                                     trj
    use Structure_Builder   , only : Generate_Structure         , &
                                     Basis_Builder
    use QCModel_Huckel      , only : EigenSystem
    use DP_potential_m      , only : Molecular_DPs
    use TD_Dipole_m         , only : wavepacket_DP
    use DP_main_m           , only : Dipole_Matrix   
                                     

    public  :: Security_Copy , Restart_State , Restart_Sys

contains
!
!
!
!============================================================================================================
 subroutine Restart_Sys( Extended_Cell , ExCell_basis , Unit_Cell , UNI , DUAL_ket , bra , ket , frame , it )
!============================================================================================================
implicit none
type(structure)                 , intent(out)   :: Extended_Cell
type(STO_basis) , allocatable   , intent(out)   :: ExCell_basis(:)
type(structure)                 , intent(out)   :: Unit_Cell
type(C_eigen)                   , intent(out)   :: UNI
complex*16                      , intent(in)    :: DUAL_ket(:,:)
complex*16                      , intent(in)    :: bra(:)
complex*16                      , intent(in)    :: ket(:)
integer                         , intent(in)    :: frame
integer                         , intent(in)    :: it

! local variables ...
type(universe) :: Solvated_System

select case ( state_of_matter )

    case( "solvated_sys" )

        CALL Prepare_Solvated_System( Solvated_System , frame )

        CALL Coords_from_Universe( Unit_Cell , Solvated_System , frame )

    case( "extended_sys" )

        CALL Coords_from_Universe( Unit_Cell , trj(frame) , frame )

    case default

        Print*, " >>> Check your state_of_matter options <<< :" , state_of_matter
        stop

end select

CALL Generate_Structure ( frame )

CALL Basis_Builder      ( Extended_Cell , ExCell_basis )

if( DP_Moment ) &
CALL Dipole_Matrix      ( Extended_Cell , ExCell_basis )

if( DP_field_ ) then

    CALL Dipole_Matrix  ( Extended_Cell , ExCell_basis )

    CALL wavepacket_DP  ( Extended_Cell , ExCell_basis , bra , ket , Dual_ket(:,1) )

    CALL Molecular_DPs  ( Extended_Cell )

end If

CALL EigenSystem( Extended_Cell , ExCell_basis , UNI , flag2=it )

end subroutine Restart_Sys
!
!                                                                    
!
!=============================================================================================
subroutine Security_Copy( MO_bra , MO_ket , DUAL_bra , DUAL_ket , bra , ket , t , it , frame )
!=============================================================================================
implicit none
complex*16  , intent(in)    :: MO_bra(:,:)
complex*16  , intent(in)    :: MO_ket(:,:)
complex*16  , intent(in)    :: DUAL_bra(:,:)
complex*16  , intent(in)    :: DUAL_ket(:,:)
complex*16  , intent(in)    :: bra(:)
complex*16  , intent(in)    :: ket(:)
real*8      , intent(in)    :: t
integer     , intent(in)    :: it
integer     , intent(in)    :: frame

! local variables ...
integer         :: i , file_err
logical , save  :: first_time = .true. 

! check whether restart is properly set ...
if( (.NOT. restart) .AND. first_time ) then
    open(unit=3, file="Restart_copy.dat", status="new", iostat=file_err , err=11 )
    11 if( file_err > 0 ) stop " <Restart_copy.dat> exists; check restart parameter "
    close(3)
    CALL system( "rm Restart_copy.dat" )
    first_time = .false.
end if

open(unit=33, file="Security_copy.dat", status="unknown", form="unformatted", action="write")

write(33) frame
write(33) it
write(33) t
write(33) size(MO_bra(:,1))
write(33) size(MO_bra(1,:))

do i = 1 , size(MO_bra(:,1))

    write(33) MO_bra(i,1) , MO_ket(i,1)

end do

do i = 1 , size(DUAL_bra(:,1))

    write(33) DUAL_bra(i,1) , DUAL_ket(i,1)

end do

do i = 1 , size(bra)

    write(33) bra(i) , ket(i)

end do

close( 33 )

! erase restart_copy file ...
CALL system( "rm Restart_copy.dat 2> qdynamo.err" )

end subroutine Security_Copy
!
!
!
!=============================================================================================
subroutine Restart_State( MO_bra , MO_ket , DUAL_bra , DUAL_ket , bra , ket , t , it , frame )
!=============================================================================================
implicit none
complex*16  , allocatable   , intent(out) :: MO_bra(:,:)
complex*16  , allocatable   , intent(out) :: MO_ket(:,:)
complex*16  , allocatable   , intent(out) :: DUAL_bra(:,:)
complex*16  , allocatable   , intent(out) :: DUAL_ket(:,:)
complex*16  , allocatable   , intent(out) :: bra(:)
complex*16  , allocatable   , intent(out) :: ket(:)
real*8                      , intent(out) :: t
integer                     , intent(out) :: it
integer                     , intent(out) :: frame

! local variables ...
integer :: i , size_r , size_c , file_err

open(unit=33, file="Restart_copy.dat", form="unformatted", status="old", action="read" , iostat=file_err , err=11 )

read(33) frame
read(33) it
read(33) t
read(33) size_r
read(33) size_c

allocate( MO_bra ( size_r , size_c ) )
allocate( MO_ket ( size_r , size_c ) )

do i = 1 , size_r

    read(33) MO_bra(i,1) , MO_ket(i,1)

end do

allocate( DUAL_bra ( size_r , size_c ) )
allocate( DUAL_ket ( size_r , size_c ) )

do i = 1 , size_r

    read(33) DUAL_bra(i,1) , DUAL_ket(i,1)
    
end do

allocate( bra ( size_r ) )
allocate( ket ( size_r ) )

do i = 1 , size_r

    read(33) bra(i) , ket(i)
    
end do

close( 33 )

11 if( file_err > 0 ) stop " <Restart_copy.dat> file not found; terminating execution"

end subroutine Restart_State   
!
!
!
end module Backup_m

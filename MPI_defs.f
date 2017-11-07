module MPI_definitions_m

    use MPI

    public :: launch_MPI, myid , world , master , slave , np , EigenComm , EigenCrew , ForceCrew

    private

    ! module variables ...
    integer    :: world , EigenComm , myid , np 
    logical    :: master = .false. , slave = .true. 
    logical    :: EigenCrew = .false. , ForceCrew = .false.
 
contains
!
!
!
!=====================
 subroutine launch_MPI
!=====================
 implicit none

! local variables ...
 integer :: err , my_color
 real*8 :: test 

! define MPI_COMM_WORLD ...
 world = MPI_COMM_WORLD

 CALL MPI_INIT(err)                   ! <== initiates MPI environment
 CALL MPI_COMM_RANK(world,myid,err)   ! <== sets the rank of processes
 CALL MPI_COMM_SIZE(world,np,err)     ! <== gets the total number of processes

 if( myid == 0 ) master = .true.
 if( myid == 0 ) slave  = .false.

! define sub_groups and new communicators ...
 if( myid <= 1 ) then
        my_color = 0
        EigenCrew = .true.
 else 
        my_color = MPI_undefined
        ForceCrew = .true.
 end if

 CALL MPI_Comm_split( world , my_color , myid , EigenComm  , err )

 end subroutine launch_MPI
!
!
end module MPI_definitions_m

module MPI_definitions_m

    use MPI

    public :: launch_MPI
    public :: world      , myid     , master , slave , np 
    public :: EigenComm  , EigenCrew  , myEigen 
    public :: KernelComm , KernelCrew , myKernel 
    public :: ForceComm  , ForceCrew  , myForce  , npForce

    private

    ! module variables ...
    integer    :: world , myid , np 
    integer    :: KernelComm , myKernel 
    integer    :: EigenComm , myEigen
    integer    :: ForceComm , myForce , npForce
    logical    :: master = .false. , slave = .true. 
    logical    :: EigenCrew  = .false. 
    logical    :: ForceCrew  = .false. 
    logical    :: KernelCrew = .false.
 
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
 logical :: drafted = .false.

! define MPI_COMM_WORLD ...
 world = MPI_COMM_WORLD

 CALL MPI_INIT(err)                   ! <== initiates MPI environment
 CALL MPI_COMM_RANK(world,myid,err)   ! <== sets the rank of processes
 CALL MPI_COMM_SIZE(world,np,err)     ! <== gets the total number of processes

 if( myid == 0 ) master = .true.
 if( myid == 0 ) slave  = .false.

! define sub_groups and new communicators ...
!------------------------------------------------------------------------
! KernelComm group = (0,1,2) ...
 select case ( myid )
    case (0)
        my_color   =  0
    case (1:2)
        my_color   =  0
        KernelCrew =  .true.
    case (3:)
        my_color   =  MPI_undefined
 end select
 CALL MPI_Comm_split( world , my_color , myid , KernelComm  , err )
 If( KernelCrew ) CALL MPI_COMM_RANK ( KernelComm , myKernel , err )   ! <== sets the rank of processes in KernelComm
 
!------------------------------------------------------------------------
! EigenComm group = (0,3)  ...
 select case ( myid )
    case (0,3)
        my_color = 0
        EigenCrew = .true.
    case (1,2,4:) 
        my_color = MPI_undefined
 end select
 CALL MPI_Comm_split( world , my_color , myid , EigenComm  , err )
 If( EigenCrew ) CALL MPI_COMM_RANK ( EigenComm , myEigen , err )      ! <== sets the rank of processes in EigenComm

!------------------------------------------------------------------------
! ForceComm group = KernelCrew + ForceCrew  ...
 select case ( myid )
    case (0:2)
        my_color = 0
        drafted  = .true.
    case (4:)
        my_color   = 0
        drafted    = .true.
        ForceCrew  = .true.
    case (3) 
        my_color = MPI_undefined
 end select
 CALL MPI_Comm_split( world , my_color , myid , ForceComm  , err )
 If( drafted ) then
    CALL MPI_COMM_RANK (ForceComm , myForce , err )   ! <== sets the rank of processes
    CALL MPI_COMM_SIZE (ForceComm , npForce , err )   ! <== gets the total number of processes
 end If

 end subroutine launch_MPI
!
!
end module MPI_definitions_m

module MPI_definitions_m

    use MPI
    use parameters_m     , only : driver , EnvField_

    public :: world , myid , master , slave , np 
    public :: ChebyComm  , ChebyCrew  , myCheby 
    public :: KernelComm , KernelCrew , myKernel 
    public :: ForceComm  , ForceCrew  , myForce  , npForce
    public :: EnvComm    , EnvCrew    , myEnvId  , npEnv
    public :: ChebyKernelComm , myChebyKernel

    public :: launch_MPI , setup_MPI_labor_force

    private

    ! module variables ...
    integer :: world , myid , np 
    integer :: KernelComm , myKernel 
    integer :: ChebyComm , myCheby
    integer :: ChebyKernelComm , myChebyKernel
    integer :: ForceComm , myForce , npForce
    integer :: EnvComm , myEnvId , npEnv
    logical :: master = .false. , slave = .true. 
    logical :: ChebyCrew  = .false. 
    logical :: ForceCrew  = .false. 
    logical :: KernelCrew = .false.
    logical :: EnvCrew    = .false.

    ! module parameters ...
    integer :: EnvProcs = 6  ! <== MPI procs dedicated to EnvFields; default = 6 ...   
 
contains
!
!
!
!=====================
 subroutine launch_MPI
!=====================
 implicit none

! local variables ...
 integer :: err 

! define MPI_COMM_WORLD ...
 world = MPI_COMM_WORLD

 CALL MPI_INIT(err)                   ! <== initiates MPI environment
 CALL MPI_COMM_RANK(world,myid,err)   ! <== sets the rank of processes
 CALL MPI_COMM_SIZE(world,np,err)     ! <== gets the total number of processes

 if( myid == 0 ) master = .true.
 if( myid == 0 ) slave  = .false.

 end subroutine launch_MPI
!
!
!
!================================
 subroutine setup_MPI_labor_force
!================================
implicit none

select case ( driver )

       case ( "slice_AO" , "slice_ElHl" )
           CALL setup_Adiabatic_Crew

       case ( "slice_Cheb" )
           CALL setup_Chebyshev_Crew

       case default

end select

end subroutine setup_MPI_labor_force
!
!
!
!===============================
 subroutine setup_Adiabatic_Crew
!===============================
 implicit none

! local variables ...
 integer :: err , my_color , ForceCrewLimit
 logical :: drafted = .false.

! define sub_groups and new communicators ...
!------------------------------------------------------------------------
! KernelComm group = (0,1,2) ; to work in "EhrenfestForce" ...
 select case ( myid )
    case (0)
        my_color   = 0
    case (1:2)
        my_color   = 0
        KernelCrew = .true.
    case (3:)
        my_color   = MPI_undefined
 end select

 CALL MPI_Comm_split( world , my_color , myid , KernelComm  , err )
 If( KernelCrew ) CALL MPI_COMM_RANK ( KernelComm , myKernel , err )   ! <== sets the rank of processes in KernelComm
 
!------------------------------------------------------------------------
! ForceComm group = KernelCrew + ForceCrew ; to work in "EhrenfestForce" ...

 ForceCrewLimit = (np-1) - merge( EnvProcs , 0 , EnvField_ )

 If( myid <= 2) then                                        ! <== case(0:2)
        my_color = 0
        drafted  = .true.
 ElseIf( (3 <= myid) .AND. (myid <= ForceCrewLimit) ) then  ! <== case(3:ForceCrewLimit)
        my_color  = 0
        drafted   = .true.
        ForceCrew = .true.
 Else                                                       ! <== case(ForceCrewLimit+1:)
        my_color = MPI_undefined
 EndIf

 CALL MPI_Comm_split( world , my_color , myid , ForceComm  , err )

 If( drafted ) then
    CALL MPI_COMM_RANK (ForceComm , myForce , err )   ! <== sets the rank of processes
    CALL MPI_COMM_SIZE (ForceComm , npForce , err )   ! <== gets the total number of processes
 end If 

! processes released for next drafting ...
 drafted = .false.
!------------------------------------------------------------------------
! EnvComm group = (0,[EnvProcs]) ; to work in "even_more_extended_huckel" ...

 If( EnvField_ ) then
 
     IF( myid == 0 ) then                    ! <== case(0)
            my_color = 0
            drafted  = .true.
     ElseIf( myid > ForceCrewLimit ) then    ! <== case(ForceCrewLimit+1:)
            my_color = 0
            drafted  = .true.
            EnvCrew  = .true.
     Else
            my_color = MPI_undefined
     EndIf
 
     CALL MPI_Comm_split( world , my_color , myid , EnvComm  , err )
 
     If( drafted ) then
        CALL MPI_COMM_RANK (EnvComm , myEnvId , err )   ! <== sets the rank of processes
        CALL MPI_COMM_SIZE (EnvComm , npEnv   , err )   ! <== gets the total number of processes
     end If
 
 EndIf
!------------------------------------------------------------------------

 end subroutine setup_Adiabatic_Crew
!
!
!
!===============================
 subroutine setup_Chebyshev_Crew
!===============================
 implicit none

! local variables ...
 integer :: err , my_color , ForceCrewLimit
 logical :: drafted = .false.

! define sub_groups and new communicators ...
!------------------------------------------------------------------------
! ChebyComm group = (0,1)  ...
 select case ( myid )
    case (0,1)
        my_color = 0
        ChebyCrew = .true.
    case (2:) 
        my_color = MPI_undefined
 end select
 CALL MPI_Comm_split( world , my_color , myid , ChebyComm  , err )
 If( ChebyCrew ) CALL MPI_COMM_RANK ( ChebyComm , myCheby , err )      ! <== sets the rank of processes in ChebyComm

! processes released for next drafting ...
 drafted = .false.
!------------------------------------------------------------------------
! KernelComm group = (0,2) ...
 select case ( myid )
    case (0,2)
        my_color   =  0
        drafted    = .true.
    case (1,3:)
        my_color   =  MPI_undefined
 end select
 CALL MPI_Comm_split( world , my_color , myid , KernelComm  , err )
 If( drafted ) CALL MPI_COMM_RANK ( KernelComm , myKernel , err )   ! <== sets the rank of processes in KernelComm

! processes released for next drafting ...
 drafted = .false.
!------------------------------------------------------------------------
! ChebyKernelComm group = (0,1,2) = ChebyComm + KernelComm ...
 select case ( myid )
    case (0:2)
        my_color   =  0
        drafted    = .true.
    case (3:)
        my_color   =  MPI_undefined
 end select
 CALL MPI_Comm_split( world , my_color , myid , ChebyKernelComm  , err )
 If( drafted ) CALL MPI_COMM_RANK ( ChebyKernelComm , myChebyKernel , err )   ! <== sets the rank of processes in KernelComm

! processes released for next drafting ...
 drafted = .false.
!------------------------------------------------------------------------
! ForceComm group = KernelCrew + ForceCrew ...

 ForceCrewLimit = (np-1) - merge( EnvProcs , 0 , EnvField_ )

 If( myid == 0) then                                        ! <== case(0)
        my_color = 0
        drafted  = .true.
 ElseIf( (2 <= myid) .AND. (myid <= ForceCrewLimit) ) then  ! <== case(2:ForceCrewLimit)
        my_color  = 0
        drafted   = .true.
        ForceCrew = .true.
 Else                                                       ! <== case(1,ForceCrewLimit+1:)
        my_color = MPI_undefined
 EndIf

 CALL MPI_Comm_split( world , my_color , myid , ForceComm  , err )
 If( drafted ) then
    CALL MPI_COMM_RANK (ForceComm , myForce , err )   ! <== sets the rank of processes
    CALL MPI_COMM_SIZE (ForceComm , npForce , err )   ! <== gets the total number of processes
 end If

! processes released for next drafting ...
 drafted = .false.
!------------------------------------------------------------------------
! EnvComm group = (0,[EnvProcs]) ; to work in "even_more_extended_huckel" ...

 If( EnvField_ ) then
 
     IF( myid == 0 ) then                    ! <== case(0)
            my_color = 0
            drafted  = .true.
     ElseIf( myid > ForceCrewLimit ) then    ! <== case(ForceCrewLimit+1:)
            my_color = 0
            drafted  = .true.
            EnvCrew  = .true.
     Else
            my_color = MPI_undefined
     EndIf
 
     CALL MPI_Comm_split( world , my_color , myid , EnvComm  , err )
 
     If( drafted ) then
        CALL MPI_COMM_RANK (EnvComm , myEnvId , err )   ! <== sets the rank of processes
        CALL MPI_COMM_SIZE (EnvComm , npEnv   , err )   ! <== gets the total number of processes
     end If
 
 EndIf
!------------------------------------------------------------------------

 end subroutine setup_Chebyshev_Crew
!
!
!
end module MPI_definitions_m

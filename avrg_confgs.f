! Subroutine for computing average properties over nuclear configurations ...

module Sampling_m

    use mpi
    use type_m
    use constants_m
    use MPI_definitions_m    , only : master , world , myid , np , slave

    use parameters_m         , only : n_t , frame_step , spectrum ,  &
                                      survival , nuclear_matter ,    &
                                      DP_Moment , EnvField_ ,        &
                                      file_type , sigma , n_part ,   &
                                      restart
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
    use FMO_m                , only : eh_tag
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
    use Data_Output          , only : Dump_stuff


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
integer                         :: k , frame , frame_init , nr , N_of_residues , err
real*8                          :: internal_sigma 
logical                         :: DIPOLE_
type(R_eigen)                   :: UNI 
type(f_grid)                    :: TDOS , SPEC
type(f_grid)    , allocatable   :: PDOS(:) 
type(f_time)                    :: QDyn
type(universe)                  :: Solvated_System

! preprocessing stuff .....................................................

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

k = 0

do frame = frame_init + (myid*frame_step) , size(trj) , (frame_step*np)

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

    ! TDOS += this TDOS
    CALL Total_DOS( UNI%erg , TDOS , internal_sigma )                                             

    do nr = 1 , N_of_residues
        ! PDOS += this PDOS
        CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr , internal_sigma )            
    end do

    If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )

    ! SPEC += this SPEC
    If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC , internal_sigma )

    If( survival ) CALL Simple_dynamics( Extended_Cell, ExCell_basis, UNI, QDyn )

    CALL DeAllocate_UnitCell    ( Unit_Cell     )
    CALL DeAllocate_Structures  ( Extended_Cell )
    DeAllocate                  ( ExCell_basis  )

    k = k + 1

end do

! finish calculation of STD ...
CALL Overall_RunningStat( k , TDOS , PDOS , SPEC , QDyn )

If( slave ) return

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
!
!==============================================================
subroutine Overall_RunningStat( k , TDOS , PDOS , SPEC , QDyn )
!==============================================================
implicit none
integer                    , intent(in)    :: k
type(f_grid)               , intent(inout) :: TDOS
type(f_grid) , allocatable , intent(inout) :: PDOS(:)
type(f_grid)               , intent(inout) :: SPEC
type(f_time)               , intent(inout) :: QDyn

!local variables ...
integer               :: i , err , TDOS_size , PDOS_size , SPEC_size , Qdyn_size , N_of_residues , n_f , all_k
integer               :: mpi_D_R = mpi_double_precision
real*8                :: dummy , factor
real*8  , allocatable :: recv_TDOS(:) , recv_SPEC(:) , time(:)
real*8  , allocatable :: recv_PDOS(:,:) 
real*8  , allocatable :: aux_Qdyn(:,:,:) , recv_std(:,:,:)

N_of_residues = size( PDOS            )
TDOS_size     = size( TDOS%average    )
PDOS_size     = size( PDOS(1)%average )
SPEC_size     = size( SPEC%average    )

!..........................................................................
! master gathers TDOS/PDOS/SPEC data ...
If( master ) then

     ! TDOS
     allocate( recv_TDOS(TDOS_size) , source = 0.d0 )
     CALL MPI_reduce( TDOS%average , recv_TDOS , TDOS_Size , MPI_D_R , mpi_SUM , 0 , world , err )
     TDOS%average = recv_TDOS
  
     ! PDSO(:)
     allocate( recv_PDOS(PDOS_size,N_of_residues) , source = 0.d0 )
     do i = 1 , N_of_residues
        CALL MPI_reduce( PDOS(i)%average , recv_PDOS(:,i) , PDOS_Size , MPI_D_R , mpi_SUM , 0 , world , err )
        PDOS(i)%average(:) = recv_PDOS(:,i)
     end do
   
     ! SPEC
     allocate( recv_SPEC(SPEC_size) , source = 0.d0 )
     CALL MPI_reduce( SPEC%average , recv_SPEC , SPEC_Size , MPI_D_R , mpi_SUM , 0 , world , err )
     SPEC%average = recv_SPEC
  
     deallocate( recv_TDOS , recv_PDOS, recv_SPEC ) 

else ! slaves send data ...

     ! TDOS 
     CALL MPI_reduce( TDOS%average , dummy , TDOS_Size , MPI_D_R , mpi_SUM , 0 , world , err )
  
     ! PDOS(:)
     do i = 1 , N_of_residues
        CALL MPI_reduce( PDOS(i)%average , dummy , PDOS_Size , MPI_D_R , mpi_SUM , 0 , world , err )
     end do
  
     !SPEC
     CALL MPI_reduce( SPEC%average , dummy , SPEC_Size , MPI_D_R , mpi_SUM , 0 , world , err )
  
     CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
     CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
     CALL DeAllocate_SPEC( SPEC , flag="dealloc" )
  
end If
!..........................................................................


!..........................................................................
! master gathers Qdyn data ...

!number of fragments ...
n_f       = size( Qdyn%fragments )
Qdyn_size = n_t * (n_f+2) * n_part

allocate( aux_Qdyn( n_t , 0:n_f+1 , n_part ) , source = 0.d0 )

If( master ) then

     allocate( time(n_t) , source = Qdyn% dyn(:,0,1) )

     CALL MPI_AllReduce( k , all_k , 1 , mpi_Integer , mpi_SUM , world , err )
     factor = dfloat(k)/dfloat(all_k)
  
     aux_Qdyn = Qdyn%dyn
  
     ! QDyn% dyns
     CALL MPI_reduce( factor * aux_Qdyn , Qdyn%dyn , Qdyn_size , MPI_D_R , mpi_SUM , 0 , world , err )
  
     CALL MPI_BCAST( Qdyn% dyn, Qdyn_size , mpi_D_R , 0 , world , err )
  
     ! QDyn% std
     ! for master: aux_Qdyn is the LOCAL average
     Qdyn% std = ( (k-1)*Qdyn% std + k*(aux_Qdyn - Qdyn%dyn)**2 ) / dfloat( all_k - 1 )
  
     CALL MPI_reduce( Qdyn% std , aux_Qdyn , Qdyn_size , MPI_D_R , mpi_SUM , 0 , world , err )

     ! finally, overall STD
     Qdyn% std(:,1:n_f+1,:) = sqrt( aux_Qdyn(:,1:n_f+1,:) )

     forall( i=1:n_part) Qdyn% dyn(:,0,i) = time
     forall( i=1:n_part) Qdyn% std(:,0,i) = time
  
     deallocate( aux_Qdyn ) 

else ! slaves process and send data ...

     CALL MPI_AllReduce( k , all_k , 1 , mpi_Integer , mpi_SUM , world , err )
     factor = dfloat(k)/dfloat(all_k)
  
     ! QDyn% dyn
     CALL MPI_reduce( factor * Qdyn%dyn , dummy , Qdyn_size , MPI_D_R , mpi_SUM , 0 , world , err )
  
     CALL MPI_BCAST( aux_Qdyn, Qdyn_size , mpi_D_R , 0 , world , err )
  
     ! QDyn% std
     ! for slaves: aux_Qdyn is the OVERALL average
     Qdyn% std = ( (k-1)*Qdyn% std + k*(Qdyn%dyn - aux_Qdyn)**2 ) / dfloat( all_k - 1 )
  
     CALL MPI_reduce( Qdyn% std , aux_Qdyn , Qdyn_size , MPI_D_R , mpi_SUM , 0 , world , err )
  
     deallocate( aux_Qdyn ) 
     CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

end If

! mission accomplished; slaves go home

end subroutine Overall_RunningStat
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
    OPEN( unit = 33 , file = "dyn.trunk/TDOS.dat" , status = "unknown" , action = "read" )
    read(33,14) frame
    do i = 1 , size(TDOS%func)
        read(33,11) TDOS%grid(i) , TDOS%average(i) 
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
                read(33,11) PDOS(nr)%grid(i) , PDOS(nr)%average(i) 
            end do
        CLOSE(33)
    end do
end if

! save peak and broadened specs ...
If( present(SPEC) ) then
    OPEN( unit = 33 , file='dyn.trunk/spectrum.dat' , status='unknown' )
        read(33,14) frame
        do i = 1 , size(SPEC%func)
            read(33,11) SPEC%grid(i) , SPEC%average(i) 
        end do
    CLOSE(33)
end if

! save time-dependent electron or hole populations ...
If( present(QDyn) ) then

    N_of_fragments = size( QDyn%fragments )

        do np = 1 , n_part

            if( eh_tag(n) == "XX" ) cycle

            OPEN( unit = 33 , file="dyn.trunk/"//eh_tag(n)//"_survival.dat" , status="unknown" )

            read(33,14) frame           
            read(33,12) QDyn%fragments
            do i = 1 , size( QDyn%dyn(:,1,1) )
                read(33,13) ( QDyn%dyn(i,nf,np) , nf=0,N_of_fragments+1 )
            end do

            CLOSE(33)

        end do
end if

11   FORMAT(2F13.9)
12   FORMAT(10A10)
13   FORMAT(F11.6,9F10.5)
14   FORMAT(I4)

end subroutine restart_AVRG
!
!
!
end module Sampling_m

! Program for computing Ehrenfest forces from Huckel Hamiltonian with Coherent-Switch-Decay-of-Mixing
module CSDM_Master

    use MPI
    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use MD_read_m         , only: atom
    use Structure_Builder , only: Unit_Cell
    use MPI_definitions_m , only: KernelComm

    public :: Ehrenfest_Master , PST , dNA_El , dNA_Hl , NewPointerState

    private

    !module variables ...
    integer                                            :: dim_N , dim_E , PST(2)
    real*8            , allocatable , dimension(:,:)   :: tmp_El, tmp_Hl
    type(d_NA_vector) , allocatable , dimension(:,:)   :: dNA_El, dNA_Hl              
    real*8            , allocatable , dimension(:,:,:) :: tmp_El_xyz , tmp_Hl_xyz

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

contains
!
!
!
!==================================================
 subroutine Ehrenfest_Master( system , size_basis )
!==================================================
implicit none
type(structure) , intent(inout) :: system
integer         , intent(in)    :: size_basis

! local variables ... 
real*8  , allocatable :: Force(:) , Force_xyz(:,:)
integer :: i , j , xyz , err
integer :: mpi_D_R = mpi_double_precision
integer :: displs(4) , recv_counts(4)

dim_E = size_basis

CALL preprocess( system )

displs = system%atoms * [ 0 , 0 , 1 , 2 ]
recv_counts = system%atoms * [ 0 , 1 , 1 , 1 ]

allocate( Force(0) )
allocate( Force_xyz(system%atoms,3) , source = D_zero )

CALL MPI_GatherV( Force , 0 , mpi_D_R , Force_xyz , recv_counts , displs , mpi_D_R , 0 , KernelComm , err ) 

do concurrent ( xyz=1:3 , i=1:system% atoms )
   atom(i)% Ehrenfest(xyz) = Force_xyz(i,xyz)
   end do
   deallocate( Force_xyz )

displs = dim_N*dim_E * [0 , 0 , 1 , 2 ]
recv_counts = dim_N*dim_E * [0 , 1 , 1 , 1 ]

CALL MPI_GatherV( tmp_El , 0 , mpi_D_R , tmp_El_xyz , recv_counts , displs , mpi_D_R , 0 , KernelComm , err ) 
CALL MPI_GatherV( tmp_Hl , 0 , mpi_D_R , tmp_Hl_xyz , recv_counts , displs , mpi_D_R , 0 , KernelComm , err ) 

do xyz=1,3 
do j=1,dim_E
do i=1,dim_N 
   dNA_El(i,j)%vec(xyz) = tmp_El_xyz(i,j,xyz)
   dNA_Hl(i,j)%vec(xyz) = tmp_Hl_xyz(i,j,xyz)
   end do
   end do
   end do
   deallocate( tmp_El_xyz , tmp_Hl_xyz )

deallocate( Force , tmp_El , tmp_Hl )

include 'formats.h'

end subroutine Ehrenfest_Master
!
!
!
!============================
 subroutine preprocess( sys )
!============================
implicit none
type(structure) , intent(in) :: sys

! local variables ...
integer :: i , j 
logical , save :: first_time = .true.                                                                                                                           

dim_N = count( sys%QMMM == "QM" .AND. sys%flex == T_ )

allocate( tmp_El (0,0) )
allocate( tmp_Hl (0,0) )

allocate( tmp_El_xyz (dim_N,dim_E,3) , source = d_zero ) 
allocate( tmp_Hl_xyz (dim_N,dim_E,3) , source = d_zero ) 

if( first_time ) then

    allocate( dNA_El (dim_N,dim_E) )
    allocate( dNA_Hl (dim_N,dim_E) )
    do concurrent( i=1:dim_N , j=1:dim_E )
       dNA_El (i,j)% vec(:) = d_zero
       dNA_Hl (i,j)% vec(:) = d_zero
       enddo
    
    CALL init_random_seed()
    first_time = F_

end if

end subroutine Preprocess
!
!
!
!
!==============================================================
 subroutine NewPointerState( system , bra , ket , QM , t_rate )
!==============================================================
implicit none
! args
type(structure), intent(in):: system
complex*16     , intent(in):: bra(:,:)
complex*16     , intent(in):: ket(:,:)
type(R_eigen)  , intent(in):: QM
real*8         , intent(in):: t_rate

! local variables ...
integer             :: i , j , Fermi , newPST(2)
real*8              :: rn , EH_jump
real*8, allocatable :: rho(:,:) , base(:,:) , P_switch(:,:) , B_kl(:,:) , Omega(:,:)
character(len=7)    :: method

Fermi = QM%Fermi_state

allocate( rho(dim_E, 2) )

! this loop: Symm. Re(rho_ij)/rho_ii, j=1(el), 2(hl) ...
do j = 1 , 2
   rho(:,j) = real( ket(:,j)*bra(PST(j),j) + ket(PST(j),j)*bra(:,j) ) / TWO
   rho(:,j) = rho(:,j) / rho( PST(j) , j )
   enddo

!============================================
! choose method = "Dynemol" or "Tully" ...
  method = "Dynemol"
!============================================
! both methods are equivalent ...
allocate(P_switch(dim_E,2))
if ( method == "Dynemol" ) then

   call Dynemol_way(QM,Omega)
   P_switch(:,:) = two * rho * Omega
   deallocate(Omega)

else

   call Tully_way( system , rho , B_kl )
   forall( j=1:2 ) P_switch(:,j) = t_rate * B_kl(:,j) 
   deallocate(B_kl)

end if
!============================================

call random_number(rn)

newPST = PST
allocate( base(0:dim_E,2) , source=d_zero )
base(0,:) = d_zero
do j = 1 , 2 
   do i = 1 , dim_E
      base(i,j) = base(i-1,j) + max(d_Zero,P_switch(i,j)) 
      if( rn > base(i-1,j) .AND. rn <= base(i,j) ) then
          newPST(j) = i     
          exit
          end if
      end do
      end do

EH_jump = (QM%erg(newPST(1)) - QM%erg(PST(1))) - (QM%erg(newPST(2)) - QM%erg(PST(2)))

if( EH_jump > Unit_Cell% MD_Kin ) then
   !transitions are not allowed ; energy forbidden 
   return
   endif

if( newPST(1) > Fermi .AND. newPST(2) <= Fermi ) then
         ! do nothing, transitions are allowed
   elseif( newPST(1) == newPST(2) ) then
         ! electron/hole annihilation
         ! system returns to GS
         newPST(1:2) = Fermi 
   elseif( (newPST(1) == PST(2)) .AND. (newPST(2) == PST(1)) ) then
         ! electron/hole exchange transition
         ! system returns to GS
         newPST(1:2) = Fermi 
   else
         ! transitions not allowed
         newPST = PST  
   endif

If( newPST(1) < newPST(2) ) then
    CALL warning("ATTENTION: electron below hole state")
    stop 
    end If

deallocate( rho , base , P_switch ) 

PST = newPST

end subroutine NewPointerState
!
!
!
!====================================
 subroutine Dynemol_way( QM , Omega ) 
!====================================
implicit none
type(R_eigen)        , intent(in)  :: QM
real*8 , allocatable , intent(out) :: Omega(:,:)

! local variables ...
integer               :: i , j
real*8  , allocatable :: newQ(:,:) 
logical               :: flip 

!local saved variables ...
logical               , save :: done = F_
real*8  , allocatable , save :: pastQ(:,:)

allocate( newQ  (dim_E,2) )
allocate( Omega (dim_E,2) )

if( .not. done ) then
    ! setup environment ...
    allocate( pastQ (dim_E,dim_E) )
    pastQ = QM%R
    done = T_
else
    ! used to calculate P_switch via Scattering Matrix (Omega): DynEMol method ...
    do i = 1 , dim_E 
       flip = dot_product( QM%L(i,:) , pastQ(:,i) ) < 0
       if(flip) pastQ(:,i) = -pastQ(:,i)
       end do

    forall(j=1:2) newQ(:,j) = QM%L(PST(j),:)

    call gemm( pastQ , newQ , Omega , 'T' )    

    !change sign for hole wvpckt ...
    Omega(:,2) = -Omega(:,2)

    do j=1,2
       Omega(PST(j),j) = d_zero
       end do

    pastQ = QM%R
end if

deallocate( newQ ) 

end subroutine Dynemol_way
!
!
!
!===========================================
 subroutine Tully_way( system , rho , B_kl ) 
!===========================================
implicit none
type(structure)               , intent(in)  :: system
real*8          , allocatable , intent(in)  :: rho(:,:)
real*8          , allocatable , intent(out) :: B_kl(:,:)

! local parameters ...
real*8, parameter :: V_factor  = 1.d-2   ! <== converts nuclear velocity: m/s (MM) to Ang/ps (QM)

!local variables ...
integer             :: i , j , n
real*8, allocatable :: v_x_dNA(:,:)

!local saved variables ...
logical               , save :: done = F_
real*8  , allocatable , save :: past_rho(:,:) 

allocate( v_x_dNA (dim_E,2) , source = d_zero )
allocate( B_kl    (dim_E,2) )

if( .not. done ) then
    ! setup environment ...
    allocate( past_rho (dim_E,2) )
    past_rho = rho
    done = T_
else
    do i = 1 , dim_E
         do n = 1 , system%atoms
            If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
            v_x_dNA(i,1) = v_x_dNA(i,1) + dot_product( atom(n)% vel(:) , dNA_El(n,i)% vec(:) )
            v_x_dNA(i,2) = v_x_dNA(i,2) - dot_product( atom(n)% vel(:) , dNA_Hl(n,i)% vec(:) )
            end do
            enddo
    v_x_dNA = v_x_dNA * V_factor

    forall( j=1:2 ) B_kl(:,j) = - two * past_rho(:,j) * v_x_dNA(:,j)

    do j=1,2
       B_kl(PST(j),j) = d_zero
       end do

    past_rho = rho
end if

deallocate( v_x_dNA ) 

end subroutine Tully_way
!
!
!
!
!==============================
 subroutine init_random_seed ()
!==============================
implicit none

!local variables ...
integer :: seed(5)

seed = [10051965,27092004,2092002,22021967,-76571]

call random_seed(put=seed(1:5))
    
end subroutine init_random_seed
!
!
!
!
end module CSDM_Master

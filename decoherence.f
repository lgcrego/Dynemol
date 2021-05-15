module decoherence_m

    use type_m
    use constants_m
    use parameters_m      , only: n_part
    use Structure_Builder , only: Unit_Cell

    public :: apply_decoherence , DecoherenceRate , DecoherenceForce

    private

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

    !module variables ...
    real*8 , allocatable :: tau_inv(:,:) , veloc(:,:) , coord(:,:)

    interface apply_decoherence
        module procedure FirstOrderDecoherence
    end interface apply_decoherence

contains
!
!
!
!========================================================================
 subroutine FirstOrderDecoherence( MO_bra , MO_ket , erg , PES , t_rate )
!========================================================================
implicit none
complex*16 , intent(inout) :: MO_bra(:,:)
complex*16 , intent(inout) :: MO_ket(:,:)
real*8     , intent(in)    :: erg(:)
integer    , intent(in)    :: PES(:)
real*8     , intent(in)    :: t_rate

! local variables ...
integer :: n , i 
real*8  :: dt , coeff , summ(2)

! J. Chem. Phys. 126, 134114 (2007)
CALL  DecoherenceRate( erg , PES )
! because wavefunction tau(wvpckt) = 2.0*tau(rho) ...
dt = t_rate * HALF

summ = d_zero
do n = 1 , n_part 
do i = 1 , size(erg) 
     if( i == PES(n) ) cycle
     MO_bra(i,n) = MO_bra(i,n) * exp(-dt*tau_inv(i,n))
     MO_ket(i,n) = MO_ket(i,n) * exp(-dt*tau_inv(i,n))
     summ(n) = summ(n) + MO_bra(i,n)*MO_ket(i,n)
     end do
     end do

do n = 1 , n_part
     coeff = MO_bra(PES(n),n) * MO_ket(PES(n),n)
     coeff = (d_one - summ(n)) / coeff
     coeff = sqrt(coeff)
     MO_bra(PES(n),n) = MO_bra(PES(n),n) * coeff
     MO_ket(PES(n),n) = MO_ket(PES(n),n) * coeff
     end do

end subroutine FirstOrderDecoherence
!
!
!
!=======================================
 subroutine DecoherenceRate( erg , PES )
!=======================================
implicit none
real*8  , intent(in) :: erg(:)
integer , intent(in) :: PES(:)

!local parameters ...
real*8 , parameter :: C = 0.5 * Hartree_2_eV   ! <== eV units

!local variables ...   
integer :: i , j
real*8  :: Const , tau , dE

! kinetic energy in eV units ...
Const = d_one + C/Unit_Cell%MD_Kin  !   E_kin()

allocate( tau_inv(size(erg),2) , source = d_zero )

do j = 1 , n_part 
do i = 1 , size(erg)
     if( i == PES(j) ) cycle 
        dE = abs(erg(i) - erg(PES(j)))
        tau_inv(i,j) = h_bar * Const / dE
        tau_inv(i,j) = d_one/tau_inv(i,j)
        end do
        end do

end subroutine DecoherenceRate
!
!
!
!
!===================================================================
 subroutine DecoherenceForce( system , MO_bra , MO_ket , erg , PST )
!===================================================================
 use MD_read_m   , only: atom
 use MM_input    , only: read_velocities
 implicit none
 type(structure) , intent(in) :: system
 complex*16      , intent(in) :: MO_bra(:,:)
 complex*16      , intent(in) :: MO_ket(:,:)
 real*8          , intent(in) :: erg(:)
 integer         , intent(in) :: PST(:)

! local parameters ...
real*8  , parameter :: eVAngs_2_Newton = 1.602176565d-9 
real*8  , parameter :: V_factor  = 1.d-2   ! <== converts nuclear velocity: m/s (MM) to Ang/ps (QM)

! local variables ...
integer  :: i, j, n, N_atoms, space
real*8   :: aux
real*8  , allocatable , dimension(:,:)     :: rho, v_x_s, tmp
real*8  , allocatable , dimension(:,:,:)   :: ForceN
real*8  , allocatable , dimension(:,:,:,:) :: s_N_ik

CALL preprocess( system )

N_atoms = system%atoms
space   = size(erg)

! set atomic forces to zero beforehand ...
forall( i=1:N_atoms ) atom(i)% f_CSDM(:) = d_zero

if( Unit_Cell% MD_Kin < mid_prec ) return
     
s_N_ik = get_S_versor( system , PST , space )

allocate( v_x_s(space,n_part) )
allocate( tmp(space,n_part) )
do j = 1 , n_part
     !!$OMP parallel do private(n) firstprivate(j) default(shared) reduction(+:v_x_s)
     do i = 1 , space
          v_x_s(i,j) = d_zero
          tmp(i,j) = d_zero
          do n = 1 , N_atoms
             If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
             v_x_s(i,j) = v_x_s(i,j) + sum(veloc(:,n)*s_N_ik(:,n,i,j)) 
             tmp(i,j) = tmp(i,j) + sum(s_N_ik(:,n,i,j)*s_N_ik(:,n,i,j)) 
             end do  
     end do
     !!$OMP end parallel do
end do
v_x_s = v_x_s * V_factor  

print*, tmp(18,1), tmp(17,1)
do concurrent( i=1:space , j=1:n_part , v_x_s(i,j)/=d_zero ) 
     v_x_s(i,j) = d_one/v_x_s(i,j)
     end do

allocate( rho(space,n_part) )
forall(j=1:2) rho(:,j) = MO_ket(:,j)*MO_bra(:,j) 

allocate( ForceN(3,N_atoms,n_part) , source=d_zero )

do i = 1 , space 
     !===================================================================
     ! electron ...
     If( i == PST(1) ) cycle     
     aux = rho(i,1)*tau_inv(i,1)*(erg(i)-erg(PST(1)))*v_x_s(i,1)
     do n = 1 , N_atoms
          If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
          ForceN(:,n,1) = ForceN(:,n,1) + aux*s_N_ik(:,n,i,1)
          end do
     !===================================================================
end do

print*, rho(16,1), rho(17,1)
print*, v_x_s(16,1), v_x_s(17,1)
 



do i = 1 , space 
     !===================================================================
     ! hole ...
     If( i == PST(2) ) cycle     
     aux = rho(i,2)*tau_inv(i,2)*(erg(i)-erg(PST(2)))*v_x_s(i,2)
     do n = 1 , N_atoms
          If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
          ForceN(:,n,2) = ForceN(:,n,2) + aux*s_N_ik(:,n,i,2)
          end do
     !===================================================================
end do

do n = 1 , N_atoms
     If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
     atom(n)% f_CSDM(:) = ( ForceN(:,n,1) - ForceN(:,n,2) ) * eVAngs_2_Newton 
     end do 
print*, atom(2)%f_CSDM(1)
deallocate( rho , tau_inv , v_x_s , s_N_ik , ForceN )  

end subroutine DecoherenceForce
!
!
!
!===============================================
 function get_S_versor( system , PST , space ) &
 result(s_N_ik)
!===============================================
use Ehrenfest_CSDM, only: d_NA_El , d_NA_Hl

implicit none
type(structure) , intent(in) :: system
integer         , intent(in) :: PST(:)
integer         , intent(in) :: space

! local parameters ...
real*8  , parameter :: V_factor  = 1.d-2   ! <== converts nuclear velocity: m/s (MM) to Ang/ps (QM)

! local variables ...
integer :: i , n , N_atoms
real*8  :: norm_S , R2 , v_x_R , aux , D_k(2)
real*8 , allocatable , dimension(:,:)     :: V_vib , v_x_dNA , norm_d
real*8 , allocatable , dimension(:,:,:,:) :: s_N_ik

N_atoms = system%atoms

! V_vib = (\vec{V}\cdot\hat{R})\hat{R} ... units=Ang/ps
allocate( V_vib(3,N_atoms) , source=d_zero )
R2 = d_zero
v_X_R = d_zero
do n = 1 , N_atoms
     R2 = R2 + dot_product(coord(:,n),coord(:,n))
     v_X_R = v_X_R + dot_product(veloc(:,n),coord(:,n))*V_factor
     enddo
do n = 1 , N_atoms
     V_vib(:,n) = (v_X_R/R2)*coord(:,n)
     enddo

! MIND: d_NA_EL and d_NA_HL vectors are zero for "fixed" or "MM" atoms ...
allocate( v_x_dNA (space,n_part) , source=d_zero )
allocate( norm_d  (space,n_part) , source=d_zero )
!!$OMP parallel do private(n,xyz) firstprivate(i) default(shared) reduction(+:v_x_dNA,norm_d)
do i = 1 , space
     do n = 1 , N_atoms
        If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
        ! electrons ... 
        v_x_dNA(i,1) = v_x_dNA(i,1) + sum(veloc(:,n)*d_NA_El(:,n,i)) 
         norm_d(i,1) =  norm_d(i,1) + sum(d_NA_El(:,n,i)*d_NA_El(:,n,i)) 
        ! holes ...
        v_x_dNA(i,2) = v_x_dNA(i,2) + sum(veloc(:,n)*d_NA_Hl(:,n,i)) 
         norm_d(i,2) =  norm_d(i,2) + sum(d_NA_Hl(:,n,i)*d_NA_Hl(:,n,i)) 
        enddo  
        enddo
!!$OMP end parallel do
v_x_dNA = a_Bohr * v_x_dNA * V_factor  ! <== units = Ang/ps ...
norm_d  = sqrt(d_one/norm_d)


! building decoherent force versor s_N_ik ...
allocate( s_N_ik(3,N_atoms,space,n_part) , source = d_zero )

do i = 1 , space
     If( i == PST(1) ) cycle     
     !=============================================================================================================
     ! electron ...
     do n = 1 , N_atoms
          If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
          ! the vector ...
          s_N_ik(:,n,i,1) = v_x_dNA(i,1)*d_NA_El(:,n,i)*norm_d(i,1) + V_vib(:,n)
          enddo
     norm_s = d_zero
     do n = 1 , n_atoms
          norm_s = norm_s + sum( s_n_ik(:,n,i,1) * s_n_ik(:,n,i,1) )
          enddo
          norm_s = sqrt(d_one/norm_s)
     do n = 1 , n_atoms
          s_n_ik(:,n,i,1) = s_n_ik(:,n,i,1) * norm_s
          enddo
     !========================================================================================================
     enddo

do i = 1 , space
     If( i == PST(2) ) cycle     
     !=============================================================================================================
     ! hole ...
     do n = 1 , N_atoms
          If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
          ! the vector ...
          s_N_ik(:,n,i,2) = v_x_dNA(i,2)*d_NA_Hl(:,n,i)*norm_d(i,2) + V_vib(:,n)
          enddo
     norm_s = d_zero
     do n = 1 , n_atoms
          norm_s = norm_s + sum( s_n_ik(:,n,i,2) * s_n_ik(:,n,i,2) )
          enddo
          norm_s = sqrt(d_one/norm_s)
     do n = 1 , n_atoms
          s_n_ik(:,n,i,2) = s_n_ik(:,n,i,2) * norm_s
          enddo
     !=========================================================================================================
     enddo

deallocate( V_vib , v_x_dNA , norm_d )

end function get_S_versor
!
!
!
!=================================
 subroutine preprocess( system )
!=================================
use MD_read_m   , only: atom
implicit none
type(structure) , intent(in) :: system

! local variables ...
integer :: n , xyz

If(.NOT. allocated(coord)) then
   allocate( coord (3,system%atoms) )
   allocate( veloc (3,system%atoms) )
   end if

do concurrent( n=1:system%atoms , xyz=1:3 )
   coord(xyz,n) = system% coord(n,xyz)
   veloc(xyz,n) = atom(n)% vel(xyz)
   end do

end subroutine preprocess
!
!
!
!================================
 function E_kin() result(kinetic)
!================================
use MD_read_m, only: MM, atom
implicit none

! local variables ...
integer :: i
real*8  :: kinetic

! calculation of the kinetic energy ...                                                                                                                                 
kinetic = d_zero
do i = 1 , MM % N_of_atoms
    kinetic = kinetic + atom(i)% mass * sum( atom(i)% vel(:) * atom(i)% vel(:) ) * half   ! <== J/kmol
end do
kinetic = kinetic * micro * kJmol_2_eV   ! <== eV

end function E_kin
!
!
!
end module decoherence_m

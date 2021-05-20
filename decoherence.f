module decoherence_m

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m      , only: n_part
    use Structure_Builder , only: Unit_Cell

    public :: apply_decoherence , DecoherenceRate , DecoherenceForce

    private

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

    !module variables ...
    integer              :: dim_3N
    real*8 , allocatable :: tau_inv(:,:) , veloc(:) , coord(:)

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
real*8 , parameter :: C = 0.1 * Hartree_2_eV   ! <== eV units

!local variables ...   
integer :: i , j
real*8  :: Const , dE

! kinetic energy in eV units ...
Const = d_one + C/Unit_Cell%MD_Kin  !   E_kin()

if( .not. allocated(tau_inv) ) then
    allocate( tau_inv(size(erg),2) , source = d_zero )
    endif

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
use MD_read_m  , only: atom
implicit none
type(structure), intent(in):: system
complex*16     , intent(in):: MO_bra(:,:)
complex*16     , intent(in):: MO_ket(:,:)
real*8         , intent(in):: erg(:)
integer        , intent(in):: PST(:)

! local parameters ...
real*8, parameter:: eVAngs_2_Newton = 1.602176565d-9 

! local variables ...
integer  :: i, j, h, n, N_atoms, space, xyz
real*8   :: f_ik 
real*8  , allocatable , dimension(:,:)   :: rho, v_x_s, s_El_ik, s_Hl_ik
real*8  , allocatable , dimension(:,:) :: Force3N 

CALL preprocess( system )
if( .not. allocated(tau_inv) ) then
    allocate( tau_inv(size(erg),2) , source = d_zero )
    endif

N_atoms = system%atoms
space   = size(erg)

! reset decoherence force to zero ...
forall( i=1:N_atoms ) atom(i)% f_CSDM(:) = d_zero

if( Unit_Cell% MD_Kin < mid_prec ) return
     
allocate( rho(space,n_part) )
forall(j=1:2) rho(:,j) = MO_ket(:,j)*MO_bra(:,j) 

CALL get_S_versor( s_El_ik , s_Hl_ik , system , PST , space )

allocate( v_x_s(space,n_part) )
CALL gemv( s_El_ik , veloc , v_x_s(:,1) , 1.d0 , 0.d0 , 'T' )
CALL gemv( s_Hl_ik , veloc , v_x_s(:,2) , 1.d0 , 0.d0 , 'T' )

do concurrent( i=1:space , j=1:n_part , v_x_s(i,j)/=d_zero ) 
     v_x_s(i,j) = d_one/v_x_s(i,j)
     enddo

allocate( Force3N(dim_3N,n_part) , source=d_zero )

do i = 1 , space 
     !===================================================================
     ! electron = 1
     If( i == PST(1) ) cycle     
     f_ik = rho(i,1)*tau_inv(i,1)*(erg(i)-erg(PST(1)))*v_x_s(i,1)
     Force3N(:,1) = Force3N(:,1) + f_ik*s_El_ik(:,i)
     !===================================================================
end do

do i = 1 , space 
     !===================================================================
     ! hole = 2
     If( i == PST(2) ) cycle     
     f_ik = rho(i,2)*tau_inv(i,2)*(erg(i)-erg(PST(2)))*v_x_s(i,2)
     Force3N(:,2) = Force3N(:,2) + f_ik*s_Hl_ik(:,i)
     !===================================================================
end do

h = 0
do n = 1 , N_atoms 
do xyz = 1 , 3 
   If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
   h = h + 1
   atom(n)% f_CSDM(xyz) = ( Force3N(h,1) - Force3N(h,2) ) * eVAngs_2_Newton 
   enddo 
   enddo

deallocate( rho , tau_inv , v_x_s , s_El_ik , s_Hl_ik , Force3N )  

end subroutine DecoherenceForce
!
!
!
!===================================================================
 subroutine get_S_versor( s_El_ik , s_Hl_ik , system , PST , space ) 
!===================================================================
use Ehrenfest_CSDM, only: d_NA_El , d_NA_Hl
implicit none
type(structure)     , intent(in) :: system
integer             , intent(in) :: PST(:)
integer             , intent(in) :: space
real*8 , allocatable, intent(out):: s_El_ik(:,:)
real*8 , allocatable, intent(out):: s_Hl_ik(:,:)

! local variables ...
integer :: i , N_atoms
real*8  :: norm , R2 , v_x_R 
real*8 , allocatable , dimension(:)   :: V_vib
real*8 , allocatable , dimension(:,:) :: v_x_dNA , versor_d_EL , versor_d_Hl 

N_atoms = system%atoms

! V_vib, units=Ang/ps
allocate( V_vib(dim_3N) )
R2    = dot_product( coord , coord )
v_X_R = dot_product( veloc , coord ) 
V_vib = v_X_R / R2 * coord

! MIND: d_NA_EL and d_NA_HL vectors are NOT defined for "fixed" or "MM" atoms ...
allocate( v_x_dNA     (space,n_part) , source=d_zero )
allocate( versor_d_El (dim_3N,space) , source=d_zero )
allocate( versor_d_Hl (dim_3N,space) , source=d_zero )

do i = 1 , space
     If( i == PST(1) ) cycle     
     !========================================================
     ! electron = 1 
     v_x_dNA(i,1)     = dot_product( veloc(:)     , d_NA_El(:,i) ) 
     norm             = dot_product( d_NA_El(:,i) , d_NA_El(:,i) ) 
     norm             = sqrt(d_one/norm)
     versor_d_El(:,i) = d_NA_El(:,i) * norm
     !========================================================
     enddo
do i = 1 , space
     If( i == PST(2) ) cycle     
     !========================================================
     ! hole  = 2 
     v_x_dNA(i,2)     = dot_product( veloc(:)     , d_NA_Hl(:,i) ) 
     norm             = dot_product( d_NA_Hl(:,i) , d_NA_Hl(:,i) ) 
     norm             = sqrt(d_one/norm)
     versor_d_Hl(:,i) = d_NA_Hl(:,i) * norm
     !========================================================
     enddo

v_x_dNA = a_Bohr * v_x_dNA ! <== units = Ang/ps ...

! building decoherent force versor s_ik ...

allocate( s_El_ik(dim_3N,space) , source = d_zero )
do i = 1 , space
     If( i == PST(1) ) cycle     
     !========================================================
     ! electron = 1 
     s_El_ik(:,i) = v_x_dNA(i,1)*versor_d_El(:,i) + V_vib(:)
     norm = dot_product( s_El_ik(:,i) , s_El_ik(:,i) )
     norm = sqrt(d_one/norm)
     s_El_ik(:,i) = s_El_ik(:,i) * norm
     !========================================================
     enddo

allocate( s_Hl_ik(dim_3N,space) , source = d_zero )
do i = 1 , space
     If( i == PST(2) ) cycle     
     !========================================================
     ! hole = 2
     s_Hl_ik(:,i) = v_x_dNA(i,2)*versor_d_Hl(:,i) + V_vib(:)
     norm = dot_product( s_Hl_ik(:,i) , s_Hl_ik(:,i) )
     norm = sqrt(d_one/norm)
     s_Hl_ik(:,i) = s_Hl_ik(:,i) * norm
     !========================================================
     enddo

deallocate( V_vib , v_x_dNA , versor_d_El , versor_d_Hl )

end subroutine get_S_versor
!
!
!
!=================================
 subroutine preprocess( system )
!=================================
use MD_read_m   , only: atom
implicit none
type(structure) , intent(in) :: system

! local parameters ...
real*8, parameter:: V_factor  = 1.d-2   ! <== converts nuclear velocity: m/s (MM) to Ang/ps (QM)

! local variables ...
integer :: k , n , xyz

If(.NOT. allocated(coord)) then
   dim_3N = 3*count( system%QMMM == "QM" .AND. system%flex == T_ )
   allocate( coord (dim_3N) )
   allocate( veloc (dim_3N) )
   endif

k = 0
do n = 1 , system%atoms 
   If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
   do xyz = 1 , 3 
      k = k + 1
      coord(k) = system% coord(n,xyz)
      veloc(k) = atom(n)% vel(xyz) * V_factor
      enddo
      enddo

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

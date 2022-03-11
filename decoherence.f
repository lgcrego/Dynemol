module decoherence_m

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m      , only: n_part
    use Structure_Builder , only: Unit_Cell

    public :: apply_decoherence , DecoherenceRate , DecoherenceForce , AdjustNuclearVeloc

    private

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

    !module variables ...
    integer                       :: dim_N
    real*8          , allocatable :: tau_inv(:,:)
    type(R3_vector) , allocatable :: nucleus(:)

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

! using kinetic energy in eV units ...
Const = d_one + C/Unit_Cell%MD_Kin  

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
integer:: i, j, h, n, N_atoms, dim_E, xyz
real*8 :: f_ik , aux , tmp
real*8           , allocatable, dimension(:,:):: rho, v_x_s
type(d_NA_vector), allocatable, dimension(:,:):: s_El_ik, s_Hl_ik, Force 

CALL preprocess( system )
if( .not. allocated(tau_inv) ) then
    allocate( tau_inv(size(erg),2) , source = d_zero )
    endif

N_atoms = system%atoms
dim_E   = size(erg)

if( Unit_Cell% MD_Kin < mid_prec ) return
     
allocate( rho(dim_E,n_part) )
forall(j=1:2) rho(:,j) = MO_ket(:,j)*MO_bra(:,j) 

CALL get_S_versor( s_El_ik , s_Hl_ik , system , PST , dim_E )

write(32,'(5F9.6)') rho(51,2), rho(50,2), rho(49,2), rho(48,2), rho(47,2)

allocate( v_x_s(dim_E,n_part) , source = d_zero )
do i = 1 , dim_E
     do n = 1 , dim_N
          aux = dot_product( nucleus(n)% v(:) , s_EL_ik(n,i)% vec(:) )
          v_x_s(i,1) = v_x_s(i,1) + aux
     
          aux = dot_product( nucleus(n)% v(:) , s_HL_ik(n,i)% vec(:) )
          v_x_s(i,2) = v_x_s(i,2) + aux
          end do
end do

do concurrent( i=1:dim_E , j=1:n_part , v_x_s(i,j)/=d_zero ) 
     v_x_s(i,j) = d_one/v_x_s(i,j)
     enddo

allocate( Force(dim_N,n_part) )
do n = 1 , dim_N

     Force(n,1)%vec = d_zero
     do i = 1 , dim_E 
          !===================================================================
          ! electron = 1
          If( i == PST(1) ) cycle     
          f_ik = rho(i,1)*tau_inv(i,1)*(erg(i)-erg(PST(1)))*v_x_s(i,1)
          Force(n,1)%vec(:) = Force(n,1)%vec(:) + f_ik * s_El_ik(n,i)%vec(:)
          !===================================================================
     end do
     
     Force(n,2)%vec = d_zero
     do i = 1 , dim_E 
          !===================================================================
          ! hole = 2
          If( i == PST(2) ) cycle     
          f_ik = rho(i,2)*tau_inv(i,2)*(erg(i)-erg(PST(2)))*v_x_s(i,2)
          Force(n,2)%vec(:) = Force(n,2)%vec(:) + f_ik * s_Hl_ik(n,i)%vec(:)
          !===================================================================
     end do

end do

h = 0
do n = 1 , N_atoms 

     ! reset decoherence force to zero ...
     atom(n)% f_CSDM(:) = d_zero

     If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle

     h = h + 1
     atom(n)% f_CSDM(:) = ( Force(h,1)%vec(:) - Force(h,2)%vec(:) ) * eVAngs_2_Newton 

     enddo

if( PST(1) /= 52 .or. pst(2) /= 50 )  print*, "before:", PST

aux = d_zero
tmp = d_zero
do n = 1 , N_atoms
     aux = aux + sum( abs(Force(n,2)%vec(:)) )
     tmp = tmp + sum( abs(Force(n,1)%vec(:)) )
end do
write(17,*)  tmp , aux


deallocate( rho , tau_inv , v_x_s , s_El_ik , s_Hl_ik , Force )  

end subroutine DecoherenceForce
!
!
!
!===================================================================
 subroutine get_S_versor( s_El_ik , s_Hl_ik , system , PST , dim_E ) 
!===================================================================
use Ehrenfest_CSDM, only: dNA_El , dNA_Hl
implicit none
type(structure)               , intent(in) :: system
integer                       , intent(in) :: PST(:)
integer                       , intent(in) :: dim_E
type(d_NA_vector), allocatable, intent(out):: s_El_ik(:,:)
type(d_NA_vector), allocatable, intent(out):: s_Hl_ik(:,:)

! local variables ...
integer :: i , n , N_atoms
real*8  :: norm , R2 , v_x_R , v_x_dNA 

N_atoms = system%atoms

! V_vib, units=Ang/ps
do n = 1 , dim_N
   R2    = dot_product( nucleus(n)%r , nucleus(n)%r )
   v_X_R = dot_product( nucleus(n)%v , nucleus(n)%r ) 
   nucleus(n)% V_vib = v_X_R / R2 * nucleus(n)%r
   end do

! MIND: dNA_El and dNA_Hl vectors are NOT defined for "fixed" or "MM" atoms ...
allocate( s_El_ik (dim_N,dim_E) )
allocate( s_Hl_ik (dim_N,dim_E) )
do concurrent( n=1:dim_N , i=1:dim_E )
   s_El_ik(n,i)% vec(:) = d_zero
   s_Hl_ik(n,i)% vec(:) = d_zero
   enddo

do n = 1 , dim_N
     do i = 1 , dim_E
          If( i == PST(1) ) cycle     
          !========================================================
          ! electron = 1 
          v_x_dNA      = dot_product( nucleus(n)% v(:)    , dNA_El(n,i)% vec(:) ) 
          norm         = dot_product( dNA_El(n,i)% vec(:) , dNA_El(n,i)% vec(:) ) 
          v_x_dNA      = v_x_dNA / sqrt(norm)

          s_El_ik(n,i)% vec = a_Bohr * v_x_dNA * dNA_El(n,i)%vec 

          s_El_ik(n,i)% vec = s_El_ik(n,i)% vec + nucleus(n)% V_vib     ! <== units = Ang/ps ...

          norm = dot_product( s_El_ik(n,i)% vec , s_El_ik(n,i)% vec )

          ! building decoherence force versor s_ik ...
          s_El_ik(n,i)% vec = s_El_ik(n,i)% vec  / sqrt(norm)
          !========================================================
          enddo

     do i = 1 , dim_E
          If( i == PST(2) ) cycle     
          !========================================================
          ! hole  = 2 

          v_x_dNA      = dot_product( nucleus(n)% v(:)    , dNA_Hl(n,i)% vec(:) ) 
          norm         = dot_product( dNA_Hl(n,i)% vec(:) , dNA_Hl(n,i)% vec(:) ) 
          v_x_dNA      = v_x_dNA / sqrt(norm)

          s_Hl_ik(n,i)% vec = a_Bohr * v_x_dNA * dNA_Hl(n,i)%vec 

          s_Hl_ik(n,i)% vec = s_Hl_ik(n,i)% vec + nucleus(n)% V_vib     ! <== units = Ang/ps ...

          norm = dot_product( s_Hl_ik(n,i)% vec , s_Hl_ik(n,i)% vec )
          norm = sqrt(d_one/norm)

          ! building decoherence force versor s_ik ...
          s_Hl_ik(n,i)% vec = s_Hl_ik(n,i)% vec / sqrt(norm)
          !========================================================
          enddo
end do

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

If(.NOT. allocated(nucleus)) then
   dim_N = count( system%QMMM == "QM" .AND. system%flex == T_ )
   allocate( nucleus (dim_N) )
   endif

k = 0
do n = 1 , system%atoms 
   If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
   do xyz = 1 , 3 
      nucleus(n)% r(xyz) = system% coord(n,xyz)
      nucleus(n)% v(xyz) = atom(n)% vel(xyz) * V_factor
      enddo
      enddo

end subroutine preprocess
!
!
!
!===============================================
 subroutine AdjustNuclearVeloc( system , QM_erg)
!===============================================
use MD_read_m   , only: atom
implicit none
type(structure), intent(in):: system
real*8         , intent(in):: QM_erg

! local variables ...
integer:: n , Nactive
real*8 :: erg_per_part , V_adjustment 

! update atomic kinetic energy ...
do n = 1 , system%atoms 
   atom(n)%kinetic = atom(n)%mass * sum(atom(n)%vel(:)*atom(n)%vel(:)) * half   ! <== J/kmol
   enddo
   atom%kinetic = atom%kinetic * kJmol_2_eV * micro    ! <== eV

! return negative QM_erg to the nuclei ...
Nactive = count( system%QMMM == "QM" .AND. system%flex == T_ )
erg_per_part = QM_erg/float(Nactive)

! reset nuclear velocities for GS ...
do n = 1 , system%atoms 
   If( system%QMMM(n) == "MM" .OR. system%flex(n) == F_ ) cycle
   V_adjustment = dsqrt(d_one + erg_per_part/atom(n)%kinetic)
   atom(n)%vel  = atom(n)%vel * V_adjustment
   enddo

! reset kinetic energy and forces for GS ...
do n = 1 , system%atoms 
   atom(n)%kinetic = atom(n)%mass * sum(atom(n)%vel(:)*atom(n)%vel(:)) * half   ! <== J/kmol
   atom(n)%kinetic = atom(n)%kinetic * kJmol_2_eV * micro                       ! <== eV
   atom(n)%ftotal  = atom(n)%f_MM
   enddo
Unit_Cell% MD_kin = sum(atom%kinetic)

end subroutine AdjustNuclearVeloc
!
!
!
end module decoherence_m

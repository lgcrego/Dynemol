module decoherence_m

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m      , only: n_part
    use Structure_Builder , only: Unit_Cell

    public :: apply_decoherence , DecoherenceRate , DecoherenceForce , AdjustNuclearVeloc , Bcast_H_matrix , Bcast_EigenVecs

    private

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

    !module variables ...
    integer                       :: dim_N
    real*8          , allocatable :: tau_inv(:,:) , H_ij(:,:) , S_ij(:,:) , Q_ij(:,:)
    type(R3_vector) , allocatable :: nucleus(:)

    interface apply_decoherence
        module procedure FirstOrderDecoherence
    end interface apply_decoherence

contains
!
!
!
!=======================================================================================
 subroutine FirstOrderDecoherence( basis , bra , ket , erg , PES , t_rate , atenuation )
!=======================================================================================
implicit none
type(STO_basis)      , intent(in)    :: basis(:)
complex*16           , intent(inout) :: bra(:,:)
complex*16           , intent(inout) :: ket(:,:)
real*8               , intent(in)    :: erg(:)
integer              , intent(in)    :: PES(:)
real*8               , intent(in)    :: t_rate
real      , optional , intent(in)    :: atenuation

! local variables ...
integer :: n , i 
real*8  :: dt , coeff , gauge , summ(2)

! J. Chem. Phys. 126, 134114 (2007)
!CALL  DecoherenceRate( erg , PES )

!default is local decoherence ...
CALL  LocalDecoherenceRate( basis , erg , PES )

! because wavefunction tau(wvpckt) = 2.0*tau(rho) ...
dt = t_rate * HALF

if(.not. present(atenuation)) then
    summ = d_zero
    do n = 1 , n_part 
    do i = 1 , size(erg) 
         if( i == PES(n) ) cycle
         bra(i,n) = bra(i,n) * exp(-dt*tau_inv(i,n))
         ket(i,n) = ket(i,n) * exp(-dt*tau_inv(i,n))
         summ(n) = summ(n) + bra(i,n)*ket(i,n)
         end do
         end do
else

    gauge = 1.d0 / atenuation

    summ = d_zero
    do n = 1 , n_part 
    do i = 1 , size(erg) 
         if( i == PES(n) ) cycle
         bra(i,n) = bra(i,n) * exp(-dt*tau_inv(i,n)*gauge)
         ket(i,n) = ket(i,n) * exp(-dt*tau_inv(i,n)*gauge)
         summ(n) = summ(n) + bra(i,n)*ket(i,n)
         end do
         end do
endif

do n = 1 , n_part
     coeff = bra(PES(n),n) * ket(PES(n),n)
     coeff = (d_one - summ(n)) / coeff
     coeff = sqrt(coeff)
     bra(PES(n),n) = bra(PES(n),n) * coeff
     ket(PES(n),n) = ket(PES(n),n) * coeff
     end do

end subroutine FirstOrderDecoherence
!
!
!
!====================================================
 subroutine LocalDecoherenceRate( basis , erg , PES )
!====================================================
use Structure_Builder , only: sys => Extended_Cell
implicit none
type(STO_basis) , intent(in) :: basis(:)
real*8          , intent(in) :: erg(:)
integer         , intent(in) :: PES(:)

!local parameters ...
real*8 , parameter :: C = 0.1 * Hartree_2_eV   ! <== eV units

!local variables ...   
integer              :: b , j , k , N , f , N_f
real*8               :: Const , tmp
character(len=1)     :: fragment
real*8 , allocatable :: aux(:) , tau_mtx(:,:) , tau_frag(:,:) , Qij_frag(:,:)

N = size(erg)

! using kinetic energy in eV units ...
Const = d_one + C/Unit_Cell%MD_Kin  

if( .not. allocated(tau_inv) ) then
    allocate( tau_inv(N,2) , source = d_zero )
    endif

allocate( tau_mtx(N,N) )

do k = 1 , n_part 

     tau_mtx(:,:) = (H_ij(:,:) - erg(PES(k))*S_ij(:,:)) / (h_bar * Const) 

     do f = 1 , size(sys%list_of_fragments)

          fragment = sys%list_of_fragments(f)

          N_f = count(basis(:)%fragment == fragment )
          allocate( tau_frag(N_f,N_f) , Qij_frag(N_f,N) , aux(N_f) )

          b=0
          do j=1,N
               if(basis(j)%fragment /= fragment ) cycle
               b=b+1
               tau_frag(:,b) = pack( tau_mtx(:,j) , basis(:)%fragment == fragment )
               enddo

          do concurrent(j=1:N)
               Qij_frag(:,j) = pack( Q_ij(:,j) , basis(:)%fragment == fragment )
               enddo

          do concurrent(j=1:N) local( aux , tmp ) shared( tau_frag , Qij_frag , tau_inv)
               do b=1,N_f 
                    aux(b) = sum( tau_frag(:,b)*Qij_frag(:,j) )
                    enddo ! <== b_loop
               tmp = sum( aux(:)*Qij_frag(:,j) )
               tau_inv(j,k) = tau_inv(j,k) + abs(tmp)
               enddo !<== j_loop

          deallocate( tau_frag , Qij_frag , aux )

     enddo !<== f_loop
enddo !<== k_loop

deallocate( tau_mtx )

end subroutine LocalDecoherenceRate
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
        tau_inv(i,j) = dE / (h_bar * Const)
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
integer:: i, j, h, n, N_atoms, dim_E
real*8 :: f_ik , aux
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

          s_El_ik(n,i)% vec = a_Bohr * v_x_dNA * dNA_El(n,i)% vec 

          s_El_ik(n,i)% vec = s_El_ik(n,i)% vec + nucleus(n)% V_vib     ! <== units = Ang/ps ...

          norm = dot_product( s_El_ik(n,i)% vec , s_El_ik(n,i)% vec )

          ! building decoherence force versor s_ik ...
          s_El_ik(n,i)% vec = s_El_ik(n,i)% vec / sqrt(norm)
         !========================================================
         enddo

     do i = 1 , dim_E
          If( i == PST(2) ) cycle     
          !========================================================
          ! hole  = 2 

          v_x_dNA      = dot_product( nucleus(n)% v(:)    , dNA_Hl(n,i)% vec(:) ) 
          norm         = dot_product( dNA_Hl(n,i)% vec(:) , dNA_Hl(n,i)% vec(:) ) 
          v_x_dNA      = v_x_dNA / sqrt(norm)

          s_Hl_ik(n,i)% vec = a_Bohr * v_x_dNA * dNA_Hl(n,i)% vec 

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
   k = k + 1
   do xyz = 1 , 3 
      nucleus(k)% r(xyz) = system% coord(n,xyz)
      nucleus(k)% v(xyz) = atom(n)% vel(xyz) * V_factor
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
!====================================================
 subroutine Bcast_H_matrix( H_matrix , S_matrix , N )
!====================================================
 implicit none
 real*8  , intent(in) :: H_matrix(:,:)
 real*8  , intent(in) :: S_matrix(:,:)
 integer , intent(in) :: N

! local variables ... 
integer :: i , j

if( .not. allocated(H_ij)) allocate( H_ij(N,N) , S_ij(N,N) )

H_ij = H_matrix 
do j = 1 , N
  do i = j+1 , N
        H_ij(j,i) = H_ij(i,j)
    end do
end do

S_ij = S_matrix

end subroutine Bcast_H_matrix
!
!
!
!==============================================
 subroutine Bcast_EigenVecs( Eigen_matrix , N )
!==============================================
 implicit none
 real*8  , intent(in) :: Eigen_matrix(:,:)
 integer , intent(in) :: N

! local variables ... 

if( .not. allocated(Q_ij)) allocate( Q_ij(N,N) )
Q_ij = Eigen_matrix

end subroutine Bcast_EigenVecs
!
!
!
end module decoherence_m

module decoherence_m

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m      , only: n_part
    use MD_read_m         , only: atom
    use Structure_Builder , only: Unit_Cell

    public :: apply_decoherence , DecoherenceForce , AdjustNuclearVeloc , Bcast_Matrices

    private

    !module parameters ...
    logical , parameter :: T_ = .true. , F_ = .false.

    !module variables ...
    integer                       :: dim_N
    real*8          , allocatable :: d_rho_ii_dt(:,:) , S_ij(:,:) , QR_ij(:,:) , QL_ij(:,:)
    type(R3_vector) , allocatable :: nucleus(:)

    interface apply_decoherence
        module procedure Local_CSDM
        module procedure Global_CSDM
    end interface apply_decoherence

contains
!
!
!
!=======================================================================================
 subroutine Local_CSDM( basis , dual_bra , PST , t_rate , MO_bra , MO_ket , slow_Decoh )
!=======================================================================================
use Structure_Builder    , only: sys => Extended_Cell
use Semi_empirical_parms , only: ChemAtom => atom
implicit none
type(STO_basis)        , intent(in)    :: basis(:)
complex*16             , intent(inout) :: dual_bra(:,:)
integer                , intent(in)    :: PST(:)
real*8                 , intent(in)    :: t_rate
complex*16             , intent(out)   :: MO_bra(:,:)
complex*16             , intent(out)   :: MO_ket(:,:)
logical     , optional , intent(in)    :: slow_Decoh

! local variables ...
integer    :: n , i , ia , a , L , dim_E
real*8     :: dt , coeff , summ(2)
integer    , allocatable , save :: list(:)
real*8     , allocatable :: decay(:,:)
complex*16 , allocatable :: AO_bra(:,:) , AO_ket(:,:) , dual_ket(:,:) , aux(:,:)
complex*16 , allocatable :: d_AL_dt(:,:) , d_AR_dt(:,:) , d_CL_dt(:,:) , d_CR_dt(:,:)

dim_E = size(basis)

CALL  Local_CSDM_Rate( sys , PST , decay )

If( present(slow_Decoh) .AND. slow_Decoh == T_ ) then
    decay = decay*HALF
    endif

! list of atoms subject to Ehrenfest force ...
if( .not. allocated(list) ) then
    allocate( list , &
    source = pack( [( L , L=1,sys%atoms )] , sys%QMMM(:) == "QM" .AND. sys%flex(:) == T_ ) &
    )
end if
 
!####################################################
! get AO_brackets ...
allocate( AO_bra(dim_E,2) , AO_ket(dim_E,2) )
AO_bra = dual_bra
do concurrent (a=1:dim_E)
       AO_ket(a,1) = sum( QL_ij(:,a)*MO_ket(:,1) )
       AO_ket(a,2) = sum( QL_ij(:,a)*MO_ket(:,2) )
       enddo

!####################################################
! decoherence of AO_brackets ... 

dt = t_rate 
! because wavefunction tau(wvpckt) = 2.0*tau) ...

allocate( d_AL_dt(dim_E,2) , source = C_zero )
allocate( d_AR_dt(dim_E,2) , source = C_zero )

do L = 1 , size(list)
     n = list(L)
     do ia = 1 , ChemAtom( sys%AtNo(n) )% DOS
           a = sys% BasisPointer(n) + ia

           AO_bra(a,:)  = AO_bra(a,:) * exp(-dt*decay(L,:))
           AO_ket(a,:)  = AO_ket(a,:) * exp(-dt*decay(L,:))
     
           d_AL_dt(a,:) =  - ( decay(L,:) ) * AO_bra(a,:)
           d_AR_dt(a,:) =  - ( decay(L,:) ) * AO_ket(a,:)
     enddo
enddo

!####################################################
! recover dual_brackets after CSDM decoherence ...

dual_bra = AO_bra

allocate( dual_ket(dim_E,2) )
do concurrent (a=1:dim_E)
       dual_ket(a,1) = sum( S_ij(:,a)*AO_ket(:,1) )
       dual_ket(a,2) = sum( S_ij(:,a)*AO_ket(:,2) )
       enddo

allocate( aux(dim_E,2) , source = d_AR_dt )
do concurrent (a=1:dim_E)
       d_AR_dt(a,1) = sum( S_ij(:,a)*aux(:,1) )
       d_AR_dt(a,2) = sum( S_ij(:,a)*aux(:,2) )
       enddo

deallocate( aux , AO_bra , AO_ket )

!####################################################
! calculating MO_brackets with CSDM decoherence ...

do concurrent (i=1:dim_E)
       MO_bra(i,1) = sum( QR_ij(:,i)*dual_bra(:,1) )
       MO_bra(i,2) = sum( QR_ij(:,i)*dual_bra(:,2) )

       MO_ket(i,1) = sum( QL_ij(i,:)*dual_ket(:,1) )
       MO_ket(i,2) = sum( QL_ij(i,:)*dual_ket(:,2) )
enddo
deallocate( dual_ket )

summ = d_zero
do n = 1 , n_part 
do i = 1 , dim_E 
     if( i == PST(n) ) cycle
     summ(n) = summ(n) + MO_bra(i,n)*MO_ket(i,n)
     end do
     end do

do n = 1 , n_part
     coeff = MO_bra(PST(n),n) * MO_ket(PST(n),n)
     coeff = (d_one - summ(n)) / coeff
     coeff = sqrt(coeff)
     MO_bra(PST(n),n) = MO_bra(PST(n),n) * coeff
     MO_ket(PST(n),n) = MO_ket(PST(n),n) * coeff
     end do

deallocate( decay )

!####################################################
! calculating d_rho_dt ...

if( .not. present(slow_Decoh) ) then
     allocate( d_CL_dt(dim_E,2) , d_CR_dt(dim_E,2) )
     
     do concurrent (i=1:dim_E)
            d_CL_dt(i,1) = sum( QR_ij(:,i)*d_AL_dt(:,1) )
            d_CL_dt(i,2) = sum( QR_ij(:,i)*d_AL_dt(:,2) )
     
            d_CR_dt(i,1) = sum( QL_ij(i,:)*d_AR_dt(:,1) )
            d_CR_dt(i,2) = sum( QL_ij(i,:)*d_AR_dt(:,2) )
            enddo
            deallocate( d_AL_dt , d_AR_dt )
     
     allocate( d_rho_ii_dt(dim_E,2) )
     do n = 1 , 2 
          d_rho_ii_dt(:,n) = real( d_CL_dt(:,n)*MO_ket(:,n) + MO_bra(:,n)*d_CR_dt(:,n) )
     enddo
     d_rho_ii_dt( PST(1) , 1 ) = d_zero
     d_rho_ii_dt( PST(2) , 2 ) = d_zero
     
     deallocate( d_CL_dt , d_CR_dt )
endif

end subroutine Local_CSDM
!
!
!
!===============================================
 subroutine Local_CSDM_Rate( sys , PST , decay )
!===============================================
use CSDM_master , only: dNA_El , dNA_Hl
implicit none
type(structure)       , intent(in)  :: sys
integer               , intent(in)  :: PST(:)
real*8  , allocatable , intent(out) :: decay(:,:)

!local parameters ...
real*8 , parameter :: V_factor = 1.d-2   ! <== converts nuclear velocity: m/s (MM) to Ang/ps (QM)

!local variables ...   
integer :: dim_E , i , k , n , xyz
real*8  :: aux

dim_N = size(dNA_El(:,1))
dim_E = size(dNA_El(1,:))

allocate( decay( dim_N , 2 ) , source = d_zero )

CALL preprocess( sys )

k = 0
do n = 1 , sys%atoms 
   If( sys%QMMM(n) == "MM" .OR. sys%flex(n) == F_ ) cycle
   k = k + 1
   do xyz = 1 , 3 
      nucleus(k)% v(xyz) = atom(n)% vel(xyz) * V_factor
      enddo
enddo

do n = 1 , dim_N

     aux = 0.d0
     do i = 1 , dim_E
        If( i == PST(1) ) cycle     
        ! electron = 1 
        decay(n,1) = aux + abs( dot_product(dNA_El(n,i)% vec(:) , nucleus(n)% v(:)) )
        enddo

     aux = 0.d0
     do i = 1 , dim_E
        If( i == PST(2) ) cycle     
        ! hole  = 2 
        decay(n,2) = aux + abs( dot_product(dNA_Hl(n,i)% vec(:) , nucleus(n)% v(:)) )
        enddo
end do

end subroutine Local_CSDM_Rate
!
!
!
!
!===================================================================
 subroutine DecoherenceForce( system , MO_bra , MO_ket , erg , PST )
!===================================================================
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
real*8           , allocatable, dimension(:,:):: v_x_s
type(d_NA_vector), allocatable, dimension(:,:):: s_El_ik, s_Hl_ik, Force 

CALL preprocess( system )

N_atoms = system%atoms
dim_E   = size(erg)

if( Unit_Cell% MD_Kin < mid_prec ) return

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
          f_ik = - d_rho_ii_dt(i,1)*(erg(i)-erg(PST(1)))*v_x_s(i,1)
          Force(n,1)%vec(:) = Force(n,1)%vec(:) + f_ik * s_El_ik(n,i)%vec(:)
          !===================================================================
     end do
     
     Force(n,2)%vec = d_zero
     do i = 1 , dim_E 
          !===================================================================
          ! hole = 2
          If( i == PST(2) ) cycle     
          f_ik = - d_rho_ii_dt(i,2)*(erg(i)-erg(PST(2)))*v_x_s(i,2)
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

deallocate( d_rho_ii_dt , v_x_s , s_El_ik , s_Hl_ik , Force )  

end subroutine DecoherenceForce
!
!
!
!===================================================================
 subroutine get_S_versor( s_El_ik , s_Hl_ik , system , PST , dim_E ) 
!===================================================================
use CSDM_master , only: dNA_El , dNA_Hl
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
!==============================================
 subroutine Bcast_Matrices( A , B , C , N )
!==============================================
 implicit none
 real*8  , intent(in) :: A(:,:)
 real*8  , intent(in) :: B(:,:)
 real*8  , intent(in) :: C(:,:)
 integer , intent(in) :: N

! local variables ... 

if( .not. allocated(QR_ij)) allocate( QR_ij(N,N) )
QR_ij = A
if( .not. allocated(QL_ij)) allocate( QL_ij(N,N) )
QL_ij = B
if( .not. allocated(S_ij) ) allocate( S_ij(N,N)  )
S_ij  = C

end subroutine Bcast_Matrices
!
!
!
!
!
!
!
!=====================================================================
 subroutine Global_CSDM( bra , ket , erg , PST , t_rate , slow_Decoh )
!=====================================================================
implicit none
complex*16           , intent(inout) :: bra(:,:)
complex*16           , intent(inout) :: ket(:,:)
real*8               , intent(in)    :: erg(:)
integer              , intent(in)    :: PST(:)
real*8               , intent(in)    :: t_rate
logical   , optional , intent(in)    :: slow_Decoh

! local variables ...
integer :: n , i 
real*8  :: dt , coeff , gauge , summ(2) 
real*8, allocatable :: decay(:,:)

! J. Chem. Phys. 126, 134114 (2007)
CALL  Global_CSDM_Rate( erg , PST , decay )
! because wavefunction tau(wvpckt) = 2.0*tau(rho) ...
dt = t_rate

If( present(slow_Decoh) .AND. slow_Decoh == T_ ) then
     decay = decay*HALF
     endif

summ = d_zero
do n = 1 , n_part 
do i = 1 , size(erg) 
     if( i == PST(n) ) cycle
     bra(i,n) = bra(i,n) * exp(-dt*decay(i,n) * HALF)
     ket(i,n) = ket(i,n) * exp(-dt*decay(i,n) * HALF)
     summ(n) = summ(n) + bra(i,n)*ket(i,n)
     end do
     end do

do n = 1 , n_part
     coeff = bra(PST(n),n) * ket(PST(n),n)
     coeff = (d_one - summ(n)) / coeff
     coeff = sqrt(coeff)
     bra(PST(n),n) = bra(PST(n),n) * coeff
     ket(PST(n),n) = ket(PST(n),n) * coeff
     end do

!####################################################
! calculating d_rho_dt ...

if( .not. present(slow_Decoh) ) then

     allocate( d_rho_ii_dt(size(erg),2) )
     
     forall(n=1:2) d_rho_ii_dt(:,n) = -decay(:,n) * bra(:,n)*ket(:,n)
     
     d_rho_ii_dt( PST(1) , 1 ) = d_zero
     d_rho_ii_dt( PST(2) , 2 ) = d_zero

endif

deallocate( decay )

end subroutine Global_CSDM
!
!
!
!================================================
 subroutine Global_CSDM_Rate( erg , PST , decay )
!================================================
implicit none
real*8  , intent(in)                :: erg(:)
integer , intent(in)                :: PST(:)
real*8  , allocatable , intent(out) :: decay(:,:)

!local parameters ...
real*8 , parameter :: C = 0.1 * Hartree_2_eV   ! <== eV units

!local variables ...   
integer :: i , j
real*8  :: Const , dE

! using kinetic energy in eV units ...
Const = d_one + C/Unit_Cell%MD_Kin  

allocate( decay( size(erg) , 2 ) , source = d_zero )

do j = 1 , n_part 
do i = 1 , size(erg)
     if( i == PST(j) ) cycle 
        dE = abs(erg(i) - erg(PST(j)))
        decay(i,j) = dE / (h_bar * Const)
        end do
        end do

end subroutine Global_CSDM_Rate
!
!
!
end module decoherence_m

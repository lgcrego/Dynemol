module DWFF

    use constants_m
    use omp_lib
    use parameters_m       , only : PBC
    use MM_parms_module    , only : DWFF_type
    use Berendsen_Barostat , only : virial_tensor
    use syst               , only : using_barostat
    use md_read_m          , only : atom, MM, molecule, special_pair_mtx
    use for_force          , only : rcut, rcut2, vscut, fscut, KAPPA, DWFF_erg
    use Build_DWFF         , only : HOH => HOH_diss_parms
    use DWFF_QMMM          , only : qd_qd, mix_q_qd

    public :: f_DWFF

    private

    ! module variables ...
    real*8  :: bond_erg
    integer :: nOX, nHX
    integer, allocatable :: H_of_O(:,:), nH_of_O(:), O_ptr(:), H_ptr(:)
    real*8 , allocatable :: f_bond_aux(:,:,:), f_ang_aux(:,:,:), ang_erg(:)

            !-----------------------------------------------------------!
            ! Legacy Conversion procedure for Electrostatic Interaction ! 
            !                                                           ! 
            !        e^2                                                !
            !   --------------- =  2.3071 * 10^(-28)  [N.m^2]           !
            !   4.pi.epsilon_0                                          !
            !                                                           ! 
            !   Therefore:                                              !
            !        e^2          1              10^(-28)               !
            !   -------------- * ---- = 2.3071 * --------  [N.m^2]      !
            !   4.pi.epsilon_0   Angs              Angs                 !
            !                                                           !
            !        e^2          1              10^(-28)  [N.m^2]      !
            !   -------------- * ---- = 2.3071 * --------  -------      !
            !   4.pi.epsilon_0   Angs            10^(-10)    [m]        !
            !                                                           !
            !        e^2          1                                     !
            !   -------------- * ---- = 2.3071 * 10^(-18)  [N.m]        !
            !   4.pi.epsilon_0   Angs                                   !
            !                                                           !
            !        e^2          1                                     !
            !   -------------- * ---- = 230.71 * 10^(-20)  [J]          !
            !   4.pi.epsilon_0   Angs                                   !
            !                                                           !
            !        e^2          1                                     !
            !   -------------- * ---- = 230.71 * factor3  [J]           !
            !   4.pi.epsilon_0   Angs                                   !
            !                                                           !
            !        e^2          1                                     !
            !   -------------- * ---- = coulomb * factor3 [J]           !
            !   4.pi.epsilon_0   Angs     |          |                  !
            !                             |          |                  !
            !                            \|/         |                  !
            !         mantissa significant figures   |                  !
            !                                       \|/                 !
            !                          applied after force calculation  !
            !                                                           ! 
            !   See parameter definitions in modulo header              ! 
            !-----------------------------------------------------------!    
contains
!
!
!===================
 subroutine f_DWFF()
!===================
    implicit none
   
    ! local variables
    integer :: i, j
    integer :: nOX, nHX
    logical, save :: done = .false.

    if( .not. done ) then
        call preprocess
        done = .true.
    end if
    call HOH_bond_topology

    do i = 1 , MM % N_of_atoms
        atom(i)% f_DWFF(:) = D_zero  
    end do

    call calculate_DWFF

    ! force units = J/Angs ...
    ! manual reduction (+: f_bond , f_ang) ...
    do i = 1, MM % N_of_atoms
        do j = 1,3
            atom(i) % f_DWFF(j) = sum(f_bond_aux(i,j,:)) + sum(f_ang_aux(i,j,:))
        end do
    end do
   
    ! energy 
    DWFF_erg = ( bond_erg + sum(ang_erg) )*factor3 
 
    deallocate( f_bond_aux , f_ang_aux , ang_erg )

end subroutine f_DWFF
!
!
!
!=========================
 subroutine calculate_DWFF
!=========================
    implicit none
    
    !local variables ...
    real*8  :: rkl(3)
    real*8  :: rkl2 , force , erg
    real*8  :: virial_private(3,3)
    integer :: i, j, k, l, O_idx, pair_of_kind
    integer :: ithr, numthr
    logical :: DWFF_special_pair
    character(len=2) :: type1, type2
    
    numthr = OMP_get_max_threads()  !safe if runtime ≤ max_threads, which is usually true

    allocate( f_bond_aux ( MM % N_of_atoms , 3 , numthr) , source = D_zero )
    allocate( f_ang_aux  ( MM % N_of_atoms , 3 , numthr) , source = D_zero )

    bond_erg = D_zero
    allocate( ang_erg(numthr) , source = D_zero )
    
    !##############################################################################
    ! INTER-MOLECULAR DWFF calculations ...

!$OMP parallel default (shared) &
!$OMP private (i, j, k, l, O_idx, rkl, rkl2, force, erg, DWFF_special_pair, type1, type2, pair_of_kind, ithr, virial_private)  &
!$OMP reduction (+: bond_erg)
                           
    ! initialize thread-local variables
    ithr = OMP_get_thread_num() + 1
    virial_private = D_zero

    !$OMP do schedule(dynamic,4)
    do k = 1 , MM % N_of_atoms - 1
       do l = k+1 , MM % N_of_atoms
       
            ! only for DWFF special pairs ...
            DWFF_special_pair = (special_pair_mtx(k,l) == 3)
            if ( .not. DWFF_special_pair ) cycle
       
            rkl(:) = atom(k) % xyz(:) - atom(l) % xyz(:)
            rkl(:) = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)
       
            rkl2 = sum( rkl(:)**2 )
       
            ! only inside cutoff radius ... 
            if( rkl2 > rcut2 ) cycle
       
            type1 = atom(k)% MMSymbol
            type2 = atom(l)% MMSymbol
   
            select case (trim(type1)//'-'//trim(type2))
            case ('HX-HX')
                pair_of_kind = 3

            case ('OX-OX')
                pair_of_kind = 2
                ! 3body does not apply 
       
            case ('HX-OX' , 'OX-HX')
                pair_of_kind = 1

            end select

            ! evaluate 2-body interaction (force and energy)
            call evaluate_2body_DWFF ( k , l , pair_of_kind , rkl2 , force , erg )
       
            f_bond_aux(k,1:3,ithr) = f_bond_aux(k,1:3,ithr) + force * rkl(1:3)
            f_bond_aux(l,1:3,ithr) = f_bond_aux(l,1:3,ithr) - force * rkl(1:3)
            
            bond_erg = bond_erg + erg
            
            !-------------------------------------------------------------------------------
            if( using_barostat% anyone ) then
                do i=1,3 ; do j=i,3
                   virial_private(i,j) = virial_private(i,j) + rkl(i) * force * rkl(j)
                end do; end do
            end if
            !---------------------------------------------------------------------------------
       end do
    end do
    !$OMP end do

    !$OMP do schedule(dynamic,4)
    do i = 1, size(O_ptr)
        O_idx = O_ptr(i)
        ! every pair of hydrogens bonded to this oxygen forms one angle
        do k = 1, nH_of_O(i) - 1
            do l = k+1, nH_of_O(i)
                call DWFF_3body ( H_of_O(i,k) , O_idx , H_of_O(i,l) , ithr , virial_private )
            end do
        end do
    end do
    !$OMP end do

    ! reduce thread-local virial into the shared virial_tensor safely
    !$OMP critical
       virial_tensor = virial_tensor + virial_private
    !$OMP end critical    

!$OMP end parallel 
    !##############################################################################

end subroutine calculate_DWFF
!
!
!
!========================================================
 subroutine DWFF_3body( atj , ati , atk , ithr , virial )
!========================================================
    implicit none
    integer , intent(in)    :: atj , ati , atk , ithr
    real*8  , intent(inout) :: virial(3,3)
    
    ! local_variables ...
    real*8 , dimension(3) :: rij, rik, f_atj, f_atk 
    real*8  :: rij_norm, rik_norm
    real*8  :: r0, cos_theta, U3, U03, exp_arg, exponential
    real*8  :: a1, a2, a3, f_ij, f_ik, inv_delta_0ij, inv_delta_0ik
    integer :: i , j
    
    !================================
    !          Angle potential ...
    !
    !             J     K
    !              \   /
    !               \ / 
    !                I 
    !
    ! MIND: the atomic sequence is JIK
    !================================

     r0 = HOH%Angle(1,2)
    
     ! rij = r_j - r_i
     rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
     rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
     rij_norm  = norm2(rij)
     if ( (rij_norm+milli) > r0 ) return

     ! rik = r_k - r_i 
     rik(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
     rik(:) = rik(:) - MM % box(:)*DNINT( rik(:) * MM % ibox(:) ) * PBC(:)
     rik_norm  = norm2(rik)
     if ( (rik_norm+milli) > r0) return

     cos_theta = dot_product(rij,rik) / ( rij_norm * rik_norm )
    
     inv_delta_0ij = 1.d0/(r0-rij_norm)
     inv_delta_0ik = 1.d0/(r0-rik_norm)
    
     exp_arg = HOH%Angle(1,3)*( inv_delta_0ij + inv_delta_0ik )
     exponential = exp(-exp_arg)
    
     U03 = HOH%Angle(1,1) * (cos_theta - HOH%Angle(1,4)) * exponential
     U3  = U03 * (cos_theta - HOH%Angle(1,4))

     ! energy of the triplet
     ang_erg(ithr) = ang_erg(ithr) + U3
    
     a1 = U3*HOH%Angle(1,3)    
     a2 = two*U03
     a3 = a2*cos_theta
     
     ! forces on each atom of the triplet 
     f_ij  = a1*inv_delta_0ij**2 + a3/rij_norm
     f_ik  = a2/rij_norm
     f_atj = f_ij*rij(:)/rij_norm - f_ik*rik(:)/rik_norm
     f_ang_aux(atj,:,ithr) = f_ang_aux(atj,:,ithr) + f_atj
    
     f_ik  = a1*inv_delta_0ik**2 + a3/rik_norm
     f_ij  = a2/rik_norm
     f_atk = f_ik*rik(:)/rik_norm - f_ij*rij(:)/rij_norm
     f_ang_aux(atk,:,ithr) = f_ang_aux(atk,:,ithr) + f_atk
    
     f_ang_aux(ati,:,ithr) = f_ang_aux(ati,:,ithr) - (f_atj + f_atk)
    
     ! inside DWFF_3body, after computing f_atj, f_atk:
     if( using_barostat% anyone ) then
         do i = 1,3 ; do j = i,3
             virial(i,j) = virial(i,j) + rij(i)*f_atj(j) + rik(i)*f_atk(j)
         end do ; end do
     end if

end subroutine DWFF_3body
!
!
!======================================================================
 subroutine evaluate_2body_DWFF( k , l , m , rkl2 , force , erg )
!======================================================================
    implicit none
    integer , intent(in)  :: k , l , m
    real*8  , intent(in)  :: rkl2 
    real*8  , intent(out) :: force
    real*8  , intent(out) :: erg
    
    ! local parameters ...
    real*8, parameter :: a1 = 1.1283791671d0  ! <== 2/sqrt(PI)
    
    ! local variables ...
    integer :: atk , atl
    real*8 :: irkl , ir2 , ir6 , ir8 , rkl
    real*8 :: zeta , erfc_zeta , arg , exp_arg2
    real*8 :: arg_Wolf , decay_Wolf , exp_Wolf
    real*8 :: Ecoul , Fcoul , f_sr , E_sr
    real*8 :: A, B, C, a2, a3, a4, d1, d2, U0
    
    rkl  = SQRT(rkl2)
    irkl = D_one / rkl
    ir2  = D_one / rkl2
    
    !----------------------------
    ! SR (short-range) only for:
    ! O-H ==> m = 1
    ! O-O ==> m = 2
    !----------------------------
    if ( any( m == [1,2]) ) then
        A = HOH% SR(m,1)
        B = HOH% SR(m,2)
        C = HOH% SR(m,3)
        
        zeta = rkl * B
        erfc_zeta = erfc(zeta) / zeta
        
        ir6 = ir2 * ir2 * ir2
        ir8 = ir6 * ir2
        
        ! SR Energy
        E_sr = A*erfc_zeta - C*ir6
        ! SR Force
        f_sr = A*( erfc_zeta + a1*exp(-zeta**2) )*ir2 - SIX*C*ir8
    else
        E_sr = 0.d0
        f_sr = 0.d0 
    end if

    !-----------------------------
    ! Coulomb electrostatic
    !-----------------------------
    a2  = irsqPI * (two * HOH% Coul(m,4))   
    arg = rkl * HOH%Coul(m,4)
    exp_arg2 = EXP(-arg**2)

    arg_Wolf   = KAPPA * rkl
    decay_Wolf = erfc(arg_Wolf)
    exp_Wolf   = EXP(-arg_Wolf**2)
    
    select case( DWFF_type )
        case("DIFFUSE")
            d1 = HOH%Coul(m,2)*erf(arg) + HOH%Coul(m,3)*erf(arg*sqrt2)
            d2 = HOH%Coul(m,2) + sqrt2*HOH%Coul(m,3)*exp_arg2 
        case("SPC_LIKE")
            d1 = D_zero
            d2 = D_zero
        case("QMMM")
            d1 = qd_qd(k,l)*erf(arg) + mix_q_qd(k,l)*erf(arg*sqrt2)
            d2 = qd_qd(k,l) + sqrt2*mix_q_qd(k,l)*exp_arg2
    end select

    ! Energy
    U0 = HOH%Coul(m,1) + d1
    Ecoul = coulomb * U0 * decay_Wolf * irkl
    
    ! Force
    ! Fcoul (damped)
    a3 = a2 * d2
    a4 = decay_Wolf + TWO*irsqPI*KAPPA*rkl*exp_Wolf
    Fcoul = coulomb * (U0*a4*ir2*irkl - a3*exp_arg2*decay_Wolf*ir2)
    
    !----------------------------------------------------
    ! total: intent(out)
    !----------------------------------------------------
    atk = atom(k)% my_intra_species_id
    atl = atom(l)% my_intra_species_id

    erg   = E_sr + Ecoul - vscut(atk,atl) + fscut(atk,atl)*( rkl - rcut )
    force = f_sr + Fcoul - fscut(atk,atl)*irkl

end subroutine evaluate_2body_DWFF
!
!
!
!====================
subroutine preprocess
!====================
    implicit none

    ! Local variables
    integer :: i, j, k, n_atoms

    !------------------------------------
    ! Basic system sizes
    !------------------------------------
    n_atoms = size(atom)
    nOX     = count(atom%MMsymbol == "OX")
    nHX     = count(atom%MMsymbol == "HX")

    !------------------------------------
    ! Build pointer lists for O and H
    !------------------------------------
    allocate(O_ptr(nOX))
    allocate(H_ptr(nHX))

    i = 0
    j = 0
    do k = 1, n_atoms
        select case (atom(k)%MMsymbol)
        case ("OX")
            i = i + 1
            O_ptr(i) = k
        case ("HX")
            j = j + 1
            H_ptr(j) = k
        end select
    end do

end subroutine preprocess
!
!
!===========================
subroutine HOH_bond_topology
!===========================
    implicit none

    ! Local parameters
    integer, parameter :: max_coord = 6

    ! Local variables
    integer :: i, j
    real(8) :: rij2, rijlen, OH_bond_cut
    real(8) :: rij(3), Txyz(3)
    real(8), allocatable :: OH_distance_table(:,:)

    Txyz = MM % box(:)

    allocate( OH_distance_table(nOX, nHX) , source = 0.d0 )

    ! --------------------------------------------------------
    ! Compute O–H and O–O distances (minimum image convention)
    ! --------------------------------------------------------
    associate (OX => atom(O_ptr), HX => atom(H_ptr))
        ! O–H distances
        do i = 1, nOX
            do j = 1, nHX
                rij = OX(i)%xyz - HX(j)%xyz
                rij = rij - Txyz * dnint(rij / Txyz) * PBC(:)

                rij2   = dot_product(rij, rij)
                rijlen = sqrt(rij2)
                  
                ! not a symmetric matrix
                OH_distance_table(i, j) = rijlen
            end do
        end do
    end associate

    ! --------------------------------------------------------
    !          build oxygen-centered adjacency
    ! --------------------------------------------------------
    if( .not. allocated(H_of_O) ) then
        allocate( H_of_O (nOX , max_coord) )
        allocate( nH_of_O(nOX) )
    end if
    nH_of_O = 0

    OH_bond_cut = HOH%Angle(1,2)
    do i = 1, nOX
        do j = 1, nHX
            if ( OH_distance_table(i,j) < OH_bond_cut ) then
                if ( nH_of_O(i) < max_coord ) then
                     nH_of_O(i) = nH_of_O(i) + 1
                     H_of_O(i, nH_of_O(i)) = H_ptr(j)
                     end if
            end if
        end do
    end do

    deallocate( OH_distance_table )

end subroutine HOH_bond_topology
!
!
!
end module DWFF

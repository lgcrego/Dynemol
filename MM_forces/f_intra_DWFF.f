module F_intra_DWFF

    use constants_m
    use omp_lib
    use parameters_m , only: PBC 
    use MD_read_m    , only: atom , molecule , MM 
    use for_force    , only: rcut, rcut2, vscut, fscut, KAPPA, DWFF_intra
    use Build_DWFF   , only: HOH => HOH_diss_parms
                                

    public :: DW_f_intra 

    private

    ! module variables
    logical              :: done = .false.
    real*8 , save        :: A, B, C
    real*8               :: bond_erg, ang_erg
    real*8 , allocatable :: f_bond(:,:), f_ang(:,:)

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
!========================
 subroutine DW_f_intra ()
!========================
    implicit none
    
    !local variables
    integer :: i, j
    
    if( .not. done ) call set_local_parameters
    
    do j = 1 , MM % N_of_atoms
        atom(j)% f_intra_DWFF(:) = D_zero  
    end do
    
    call intra_2body_DWFF
    
    call intra_3body_DWFF
    
    ! local force (units = J/Angs) ...
    do i = 1, MM % N_of_atoms
       atom(i)% f_intra_DWFF(:) = f_bond(i,:) + f_ang(i,:)
    end do

    ! energy (Joule units)
    DWFF_intra = (bond_erg + ang_erg)*factor3 
    
    deallocate( f_bond , f_ang )
    
end subroutine DW_f_intra
!
!
!
!==============================
 subroutine intra_2body_DWFF ()
!==============================
    implicit none
    
    ! local variables
    real*8  :: rkl(3)
    real*8  :: rkl2, force, erg
    integer :: i, j, k, l, pair_of_kind
    character(len=2) :: type1, type2
    
    allocate( f_bond (MM%N_of_atoms,3) , source = D_zero )
    
    bond_erg = D_zero
    
    do i = 1 , MM % N_of_molecules
    
         if ( .not. molecule(i)% DWFF ) cycle
    
         do j = 1 , molecule(i) % NintraIJ
    
              k = molecule(i) % IntraIJ(j,1) 
              l = molecule(i) % IntraIJ(j,2) 

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
                  pair_of_kind = 2  !<== not for intra_DWFF
                  cycle
              case ('HX-OX' , 'OX-HX')
                  pair_of_kind = 1
              end select
    
              call evaluate_2body_DWFF( k , l , pair_of_kind , rkl2 , force , erg )
    
              f_bond(k,1:3) = f_bond(k,1:3) + force * rkl(:)
              f_bond(l,1:3) = f_bond(l,1:3) - force * rkl(:)
    
              bond_erg = bond_erg + erg
    
         end do
    end do

end subroutine intra_2body_DWFF
!
!
!
!======================================================================
 subroutine evaluate_2body_DWFF( k , l , m , rkl2 , DW_force , DW_erg )
!======================================================================
    implicit none
    integer , intent(in)  :: k, l, m
    real*8  , intent(in)  :: rkl2 
    real*8  , intent(out) :: DW_force
    real*8  , intent(out) :: DW_erg
    
    ! local parameters ...
    real*8, parameter :: a1 = 1.1283791671d0  ! <== 2/sqrt(PI)
    
    ! local variables ...
    integer :: atk , atl
    real*8 :: irkl , ir2 , ir6 , ir8 , rkl
    real*8 :: zeta , erfc_zeta , arg , exp_arg2
    real*8 :: arg_Wolf , decay_Wolf , exp_Wolf
    real*8 :: a2, a3, a4, d1, d2, U0
    real*8 :: Ecoul , Fcoul , f_sr , E_sr
    
    rkl  = SQRT(rkl2)
    irkl = D_one / rkl
    ir2  = D_one / rkl2
    
    !----------------------------------------------------
    ! SR (short-range) applies only to O-H pair ...
    !----------------------------------------------------
    if (m==1) then 
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

    !----------------------------------------------------
    ! Coulomb applies to O-H and H-H pairs
    ! O-H ==> m = 1
    ! H-H ==> m = 3
    !----------------------------------------------------
    a2  = irsqPI * ( two * HOH% Coul(m,4) )   
    arg = rkl * HOH%Coul(m,4)
    exp_arg2 = EXP(-arg**2)

    arg_Wolf   = KAPPA * rkl
    decay_Wolf = erfc(arg_Wolf)
    exp_Wolf   = EXP(-arg_Wolf**2)   

    if( HOH%contain_diffuse ) then
        d1 = HOH%Coul(m,2)*erf(arg) + HOH%Coul(m,3)* erf(arg*sqrt2)
        d2 = HOH%Coul(m,2) + sqrt2*HOH%Coul(m,3)* exp_arg2 
    else
        d1 = D_zero
        d2 = D_zero
    end if
    
    ! Energy
    U0 = HOH%Coul(m,1) + d1
    Ecoul = coulomb * U0 * decay_Wolf * irkl
    
    ! Force
    ! Fcoul (not damped)
    a3 = a2 * d2
    a4 = decay_Wolf + TWO*irsqPI*KAPPA*rkl*exp_Wolf
    Fcoul = coulomb * (U0*a4*ir2*irkl - a3*exp_arg2*decay_Wolf*ir2)
    
    !----------------------------------------------------
    ! total 
    !----------------------------------------------------
    
    atk = atom(k)% my_intra_species_id
    atl = atom(l)% my_intra_species_id

    DW_erg   = E_sr + Ecoul - vscut(atk,atl) + fscut(atk,atl)*( rkl - rcut )
    DW_force = f_sr + Fcoul - fscut(atk,atl)*irkl

end subroutine evaluate_2body_DWFF
!
!
!
!==============================
 subroutine intra_3body_DWFF( )
!==============================
    implicit none
    
    ! local parameter ...
    real*8 :: r0 = 1.6
    
    ! local_variables ...
    real*8 , dimension (3):: rij, rik 
    real*8  :: rij_norm, rik_norm, fxyz
    real*8  :: cos_theta, U3, U03, exp_arg, exponential
    real*8  :: a1, a2, a3, f_ij, f_ik, inv_delta_0ij, inv_delta_0ik
    integer :: i, j, k, l, n, ati, atj, atk 
    
    !================================
    ! Angle - bending potential ...
    !             J     K
    !              \   /
    !               \ / 
    !                I 
    !================================
    
    allocate( f_ang (MM% N_of_atoms,3) , source = D_zero )
    
    ang_erg = D_zero
    
    do i = 1 , MM % N_of_molecules
    
         if ( .not. molecule(i)% DWFF ) cycle
    
            ! MIND: the atomic sequence is JIK, with ATOM I AT THE VERTEX
            call get_bond_triplet(i, ati, atj, atk)
    
            if(.not. any([atom(atj)%flex, atom(ati)%flex, atom(atk)%flex])) cycle
    
                ! MIND: the atomic sequence is JIK, with ATOM I AT THE VERTEX
                ! rij = r_j - r_i
                rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
                rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
                rij_norm = norm2(rij)
                if ( (rij_norm+milli) > r0 ) cycle

                ! rik = r_k - r_i 
                rik(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
                rik(:) = rik(:) - MM % box(:)*DNINT( rik(:) * MM % ibox(:) ) * PBC(:)
                rik_norm = norm2(rik)
                if ( (rik_norm+milli) > r0) cycle
    
                cos_theta = dot_product(rij,rik) / ( rij_norm * rik_norm )
    
                inv_delta_0ij = 1.d0/(r0-rij_norm)
                inv_delta_0ik = 1.d0/(r0-rik_norm)
    
                exp_arg = HOH%Angle(1,3)*( inv_delta_0ij + inv_delta_0ik )
                exponential = exp(-exp_arg)
    
                U03 = HOH%Angle(1,1) * (cos_theta - HOH%Angle(1,4)) * exponential
                U3  = U03 * (cos_theta - HOH%Angle(1,4))
    
                ! energy of the triplet
                ang_erg = ang_erg + U3
    
                a1 = U3*HOH%Angle(1,3)    
                a2 = two*U03
                a3 = a2*cos_theta
                
                ! forces on each atom of the triplet 
                f_ij = a1*inv_delta_0ij**2 + a3/rij_norm
                f_ik = a2/rij_norm
                f_ang(atj,:) = f_ij*rij(:)/rij_norm - f_ik*rik(:)/rik_norm
    
                f_ik = a1*inv_delta_0ik**2 + a3/rik_norm
                f_ij = a2/rik_norm
                f_ang(atk,:) = f_ik*rik(:)/rik_norm - f_ij*rij(:)/rij_norm
    
                f_ang(ati,:) = -f_ang(atj,:) -f_ang(atk,:)
    end do
    
end subroutine intra_3body_DWFF
!
!
!
!======================================================
 subroutine get_bond_triplet(i, ati, atj, atk)
! Given intraIJ array, determine the OX/HX atom indices
! and return them in ati, atj, atk.
!======================================================
    implicit none
    integer, intent(in)  :: i
    integer, intent(out) :: ati, atj, atk ! selected indices
    
    ! local variables
    integer :: idx1, idx2, idx3, idx4
    
    associate( intraIJ => molecule(i)% intraIJ )
          ! Unpack once
          idx1 = intraIJ(1,1)
          idx2 = intraIJ(1,2)
          idx3 = intraIJ(2,1)
          idx4 = intraIJ(2,2)
    
          if (atom(idx1)%MMSymbol == "OX") then
              ati = idx1
              atj = idx2
              atk = merge(idx3, idx4, atom(idx3)%MMSymbol == "HX")
    
          elseif (atom(idx2)%MMSymbol == "OX") then
              ati = idx2
              atj = idx1
              atk = merge(idx3, idx4, atom(idx3)%MMSymbol == "HX")
    
          else
              atj = idx1
              atk = idx2
              ati = merge(idx3, idx4, atom(idx3)%MMSymbol == "OX")
          end if
    end associate
    
end subroutine get_bond_triplet
!
!
!
!===================
 function ERFC ( X )
!===================
    ! ERFC = (1 - ERF), complementary error function ...
    implicit none
    real*8 :: ERFC
    real*8 :: A1, A2, A3, A4, A5, P, T, X, XSQ, TP
    parameter ( A1 = 0.254829592d0, A2 = -0.284496736d0 ) 
    parameter ( A3 = 1.421413741d0, A4 = -1.453122027d0 ) 
    parameter ( A5 = 1.061405429d0, P  =  0.3275911d0   ) 
    
    T    = d_one / ( d_one + P * X )
    XSQ  = X * X
    TP   = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
    ERFC = TP * EXP ( -XSQ )

end function ERFC
!
!
!
!================================
subroutine set_local_parameters()
!================================
    implicit none
    
    !local variables
    real*8 :: chrg_decay
    
    A = HOH% SR(1,1)
    B = HOH% SR(1,2)
    C = HOH% SR(1,3)
    
    done = .true.
    
end subroutine set_local_parameters
!
!
!
end module F_intra_DWFF

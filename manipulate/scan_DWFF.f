module DWFF_stuff

    use types_m
    use constants_m
    use Build_DWFF  , only : HOH => HOH_diss_parms

    private

    public :: configuration_energy

    !module variables
    integer, allocatable:: H2_list(:, :), O_list(:)

    !module parameters
    real(8), parameter:: milli   = 1.0d-3
    real(8), parameter:: J_2_eV  = 6.242d18
    real(8), parameter:: factor3 = 1.0d-20
    real(8), parameter:: coulomb = 230.7113d0   ! e^2/(4.pi.epsilon_0)*(1/Angs) [Joule], mantissa part

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
!====================================================
 function configuration_energy(atom) result(DWFF_erg)
!====================================================
    implicit none
    type(atomic), intent(in):: atom(:)
   
    ! local variables
    real(8) :: two_body_erg, three_body_erg, DWFF_erg
    real(8) :: bond, ang

    call build_H2_list(atom) 
    
    call DWFF( atom, two_body_erg, three_body_erg )

    ! energy (Joule units)
    bond = two_body_erg * factor3 
    ang  = three_body_erg * factor3 

    !energy (eV units)
    bond = bond * J_2_eV 
    ang  = ang  * J_2_eV 

    !total energy (eV units)
    DWFF_erg = bond + ang 
 
    deallocate( H2_list, O_list )

end function configuration_energy
!
!
!
!=====================================================
 subroutine DWFF( atom, two_body_erg, three_body_erg )
!=====================================================
    implicit none
    type(atomic), intent(in)  :: atom(:)
    real*8      , intent(out) :: two_body_erg
    real*8      , intent(out) :: three_body_erg

    !local variables ...
    real*8  :: rkl(3), rkl2 , E_kl
    integer :: k, l, pair_of_kind, n_atoms, n_O, n_H2
    logical :: DWFF_pair
    character(len=2) :: type1, type2
    
    n_atoms = size(atom)

    two_body_erg = D_zero

    do k = 1, n_atoms - 1
       do l = k+1, n_atoms
       
            ! only for DWFF atoms ...
            DWFF_pair = (atom(k)%DWFF .and. atom(l)%DWFF)
            if ( .not. DWFF_pair ) cycle
       
            rkl(:) = atom(k) % xyz(:) - atom(l) % xyz(:)
            rkl2   = dot_product(rkl,rkl)
       
            type1 = atom(k)% MMSymbol
            type2 = atom(l)% MMSymbol
   
            select case (trim(type1)//'-'//trim(type2))
            case ('HX-HX')
                pair_of_kind = 3
                ! 3body does not apply 
       
            case ('OX-OX')
                pair_of_kind = 2
                ! 3body does not apply 
       
            case ('HX-OX' , 'OX-HX')
                pair_of_kind = 1
                ! not usgin 3body for the moment

            case default
                cycle

            end select

            ! evaluate 2-body interaction 
            call evaluate_2body_DWFF (pair_of_kind, rkl2, E_kl )
       
            two_body_erg = two_body_erg + E_kl
            
       end do
    end do

    !---------------------------------------------------------
    ! 3body calculations: 3body_DWFF( H, O, H )
    !---------------------------------------------------------
    n_O  = size(O_list)
    n_H2 = size(H2_list, dim=1)

    three_body_erg = D_zero

    do k = 1, n_O
       do l = 1, n_H2
          
          if( .not. all( [atom(H2_list(l,1))%DWFF, &
                          atom( O_list(k)  )%DWFF, &
                          atom(H2_list(l,2))%DWFF] ) ) cycle

          three_body_erg = three_body_erg + &
                         + evaluate_3body_DWFF ( atom, H2_list(l,1), O_list(k), H2_list(l,2) )
       end do
    end do

end subroutine DWFF
!
!
!
!=================================================
 subroutine evaluate_2body_DWFF( m , rkl2 , E_kl )
!=================================================
    implicit none
    integer , intent(in)  :: m
    real*8  , intent(in)  :: rkl2 
    real*8  , intent(out) :: E_kl
    
    ! local variables ...
    real*8 :: irkl, ir2, ir6, rkl
    real*8 :: zeta, erfc_zeta, arg
    real*8 :: Ecoul, E_sr, A, B, C, U0
    
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
        
        ! SR Energy
        E_sr = A*erfc_zeta - C*ir6
    else
        E_sr = 0.d0
    end if

    !-----------------------------
    ! Coulomb contribution
    !-----------------------------
    arg = rkl * HOH%Coul(m,4)

    ! Energy
    U0 = HOH%Coul(m,1)                     & 
       + HOH%Coul(m,2) * erf(arg)          & 
       + HOH%Coul(m,3) * erf(arg * sqrt2)

    Ecoul = coulomb * U0 * irkl
    
    E_kl = E_sr + Ecoul 

end subroutine evaluate_2body_DWFF
!
!
!
!======================================================================
 function evaluate_3body_DWFF( atom, atj , ati , atk )  result(ang_erg)
!======================================================================
    implicit none
    type(atomic), intent(in):: atom(:)
    integer     , intent(in):: atj , ati , atk
    real(8):: ang_erg
    
    ! local parameter ...
    real*8 :: r0 = 1.6d0
    
    ! local_variables ...
    real*8 , dimension (3):: rij, rik 
    real*8  :: rij_norm, rik_norm
    real*8  :: cos_theta, U3, U03, exp_arg, exponential
    real*8  :: inv_delta_0ij, inv_delta_0ik

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

     !initializing result
     ang_erg = 0.d0 
    
     ! rij = r_j - r_i
     rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
     rij_norm  = norm2(rij)
     if ( (rij_norm+milli) > r0 ) return

     ! rik = r_k - r_i 
     rik(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
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
     ang_erg = U3
    
end function evaluate_3body_DWFF
!
!
!
!==============================
 subroutine build_H2_list(atom)
!==============================
    implicit none
    type(atomic), intent(in) :: atom(:)

    !local variables
    integer :: i, j, k, n_atoms
    character(len=2), allocatable :: MMSymbol(:)

    n_atoms = size(atom)

    allocate( MMSymbol(n_atoms) )
    forall(i=1:n_atoms) MMSymbol(i) = trim(adjustL(atom(i)% MMSymbol))

    !========================================
    !               H2_list 
    !========================================
    ! ---------- First pass: count ----------
    k = 0
    do i = 1, n_atoms - 1
        if (MMSymbol(i) /= "HX") cycle
        do j = i + 1, n_atoms
            if (MMSymbol(j) /= "HX") cycle
            k = k + 1
        end do
    end do

    allocate(H2_list(k, 2))

    ! ---------- Second pass: fill ----------
    k = 0
    do i = 1, n_atoms - 1
        if (MMSymbol(i) /= "HX") cycle
        do j = i + 1, n_atoms
            if (MMSymbol(j) /= "HX") cycle
            k = k + 1
            H2_list(k, :) = [i, j]
        end do
    end do

    !========================================
    !               O_list 
    !========================================
    k = count( MMSymbol == "OX" )
    
    allocate( O_list(k) )

    k=0
    do i = 1, n_atoms
        if (MMSymbol(i) == "OX") then
           k = k + 1
           O_list(k) = i
        end if
    end do

end subroutine build_H2_list
!
!
!
end module DWFF_stuff

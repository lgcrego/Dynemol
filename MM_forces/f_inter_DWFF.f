module F_inter_DWFF

    use constants_m
    use omp_lib
    use syst         , only : using_barostat
    use parameters_m , only : PBC
    use md_read_m    , only : atom, MM, molecule, special_pair_mtx
    use Data_Output  , only : Net_Charge
    use for_force    , only : rcut, vrecut, frecut, rcutsq, vscut, fscut, KAPPA, DWFF_inter
    use Build_DWFF   , only : HOH => HOH_diss_parms

    public :: f_DWFF_inter , virial_tensor

    private

    ! module variables ...
    real*8 , allocatable :: f_bond(:,:), f_ang(:,:)
    real*8 , save        :: virial_tensor(3,3)
    real*8               :: bond_erg, ang_erg, f_angle
    
contains
!
!
!=========================
 subroutine f_DWFF_inter()
!=========================
    implicit none
   
    ! local variables
    integer :: i, j
    
    do j = 1 , MM % N_of_atoms
        atom(j)% f_inter_DWFF(:) = D_zero  
    end do
   
    call inter_DWFF
    
    ! local force units = J/Angs ...
    do i = 1, MM % N_of_atoms
       atom(i)% f_inter_DWFF(:) = f_bond(i,:) + f_ang(i,:)
    end do
   
    ! energy 
    DWFF_inter = (bond_erg + ang_erg)*factor3 
   
    deallocate( f_bond , f_ang )

end subroutine f_DWFF_inter
!
!
!
!=====================
 subroutine inter_DWFF
!=====================
    implicit none
    
    !local variables ...
    real*8  :: rkl(3) , cm_kl(3)
    real*8  :: rkl2 , force , erg
    integer :: i , j , k , l , pair_of_kind
    integer :: nresidl , nresidk
    logical :: DWFF_special_pair
    character(len=2) :: type1, type2
    

    bond_erg = D_zero
    ang_erg  = D_zero
    allocate( f_bond (MM% N_of_atoms,3) , source = D_zero )
    allocate( f_ang  (MM% N_of_atoms,3) , source = D_zero )
    
    !##############################################################################
    ! INTER-MOLECULAR DWFF calculations ...
                            
         do k = 1 , MM % N_of_atoms - 1
         do l = k+1 , MM % N_of_atoms
    
              ! only for DWFF special pairs ...
              DWFF_special_pair = (special_pair_mtx(k,l) == 3)
              if ( .not. DWFF_special_pair ) cycle
         
              ! only for different molecules ...
              if ( atom(k)% nr == atom(l)% nr ) cycle
    
              rkl(:) = atom(k) % xyz(:) - atom(l) % xyz(:)
              rkl(:) = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)
         
              rkl2 = sum( rkl(:)**2 )
        
              ! only inside cutoff radius ... 
              if( rkl2 > rcutsq ) cycle
    
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
                  !------------------------------------------------------
                  ! 3body calculations
                  if ( atom(k)% MMSymbol == 'HX' ) then
                       call inter_3body_DWFF ( k , l , HOH%H_ptr(1) )
                       call inter_3body_DWFF ( k , l , HOH%H_ptr(2) )
                  else
                       call inter_3body_DWFF ( l , k , HOH%H_ptr(1) )
                       call inter_3body_DWFF ( l , k , HOH%H_ptr(2) )
                  end if
                  !------------------------------------------------------
     
              end select
    
              call evaluate_2body_inter_DWFF ( pair_of_kind , rkl2 , force , erg )
         
              f_bond(k,1:3) = f_bond(k,:) + force * rkl(:)
              f_bond(l,1:3) = f_bond(l,:) - force * rkl(:)
              
              bond_erg = bond_erg + erg
              
              !-------------------------------------------------------------------------------
              if( using_barostat% inter ) &
              then
                    nresidk = atom(k)% nr
                    nresidl = atom(l)% nr
                    cm_kl(:) = molecule(nresidk) % cm(:) - molecule(nresidl) % cm(:)
                    cm_kl(:) = cm_kl(:) - MM % box * DNINT( cm_kl(:) * MM % ibox(:) ) * PBC(:)
                    do i=1,3 ; do j=i,3
                       virial_tensor(i,j) = virial_tensor(i,j) + cm_kl(i) * force * rkl(j)
                    end do; end do
              end if
              !---------------------------------------------------------------------------------

         end do
         end do

end subroutine inter_DWFF
!
!
!
!==============================================
 subroutine inter_3body_DWFF( atj , ati , atk )
!==============================================
    implicit none
    integer , intent(in) :: atj , ati , atk
    
    ! local parameter ...
    real*8 :: r0 = 1.6
    
    ! local_variables ...
    real*8 , dimension (3):: rij, rik 
    real*8  :: rijq, rijsq, rikq, riksq
    real*8  :: cos_theta, U3, U03, exp_arg, exponential
    real*8  :: a1, a2, a3, f_ij, f_ik, inv_delta_0ij, inv_delta_0ik
    
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
    
     ! rij = r_j - r_i
     rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
     rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
     rijq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
     rijsq  = SQRT(rijq)
     if ( (rijsq+milli) > r0 ) return
    
     ! rik = r_k - r_i 
     rik(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
     rik(:) = rik(:) - MM % box(:)*DNINT( rik(:) * MM % ibox(:) ) * PBC(:)
     rikq   = rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3)
     riksq  = SQRT(rikq)
     if ( (riksq+milli) > r0) return
    
     cos_theta = (rij(1)*rik(1) + rij(2)*rik(2) + rij(3)*rik(3)) / ( rijsq * riksq )
    
     inv_delta_0ij = 1.d0/(r0-rijsq)
     inv_delta_0ik = 1.d0/(r0-riksq)
    
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
     f_ij = a1*inv_delta_0ij**2 + a3/rijsq
     f_ik = a2/rijsq
     f_ang(atj,:) = f_ang(atj,:) + f_ij*rij(:)/rijsq - f_ik*rik(:)/riksq
    
     f_ik = a1*inv_delta_0ik**2 + a3/riksq
     f_ij = a2/riksq
     f_ang(atk,:) = f_ang(atk,:) + f_ik*rik(:)/riksq - f_ij*rij(:)/rijsq
    
     f_ang(ati,:) = f_ang(ati,:) - (f_ang(atj,:) + f_ang(atk,:))
    
end subroutine inter_3body_DWFF
!
!
!==============================================================
 subroutine evaluate_2body_inter_DWFF( k , rkl2 , force , erg )
!==============================================================
    implicit none
    integer , intent(in)  :: k
    real*8  , intent(in)  :: rkl2 
    real*8  , intent(out) :: force
    real*8  , intent(out) :: erg
    
    ! local parameters ...
    real*8, parameter :: a1 = 1.1283791671d0  ! <== 2/sqrt(PI)
    
    ! local variables ...
    real*8 :: irkl , ir2 , ir6 , ir8 , rkl
    real*8 :: zeta , erfc_zeta , arg , exp_arg2
    real*8 :: a2, a3, U0 , Ecoul , Fcoul , f_sr , E_sr
    real*8 :: A, B, C
    
    rkl  = SQRT(rkl2)
    irkl = D_one / rkl
    ir2  = D_one / rkl2
    
    !----------------------------
    ! SR (short-range) only for:
    ! O-H ==> k = 1
    ! O-O ==> k = 2
    !----------------------------
    if ( any( k == [1,2]) ) then
        A = HOH% SR(k,1)
        B = HOH% SR(k,2)
        C = HOH% SR(k,3)
        
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
    a2  = irsqPI * (two * HOH% Coul(k,4))   
    arg = rkl * HOH%Coul(k,4)
    exp_arg2 = EXP(-arg**2)
    
    ! Energy
    U0 = HOH%Coul(k,1) + HOH%Coul(k,2)*erf(arg) + HOH%Coul(k,3)* erf(arg*sqrt2)
    Ecoul = coulomb * U0 * irkl
    
    ! Force
    ! Fcoul (not damped)
    a3 = a2 * ( HOH%Coul(k,2) + sqrt2*HOH%Coul(k,3) * exp_arg2 )
    Fcoul = coulomb * (U0*ir2*irkl - a3*exp_arg2*ir2)
    
    !----------------------------------------------------
    ! total: intent(out)
    !----------------------------------------------------
    erg   = E_sr + Ecoul
    force = f_sr + Fcoul

end subroutine evaluate_2body_inter_DWFF
!
!
!
!
end module F_inter_DWFF

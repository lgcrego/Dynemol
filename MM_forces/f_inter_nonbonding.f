module F_inter_nonbond

    use constants_m
    use omp_lib
    use syst               , only : using_barostat
    use type_m             , only : warning
    use parameters_m       , only : PBC, QMMM
    use md_read_m          , only : atom, MM, molecule, special_pair_mtx
    use Data_Output        , only : Net_Charge
    use gmx2mdflex         , only : SpecialPairs
    use Berendsen_Barostat , only : virial_tensor
    use for_force          , only : rcut, vrecut, frecut, rcutsq, pot_INTER, Coul_inter, evdw, vscut, fscut, KAPPA

    public :: f_inter_nonbonding

    private

contains
!
!
!
!==============================
 subroutine f_inter_nonbonding
!==============================
implicit none

!local variables ...
real*8 , allocatable :: tmp_f_inter(:,:)
real*8               :: rkl(3) , cm_kl(3)
real*8               :: Ecoul_damped , Fcoul , fs , vsr  , rkl2 
integer              :: i , j , k , l , atk , atl, nresidk, nresidl

allocate( tmp_f_inter (MM % N_of_atoms, 3), source = D_zero )
evdw = D_zero
Coul_inter = D_zero

!##############################################################################
! INTER-MOLECULAR vdW and Coulomb forces ...
!$OMP parallel DO &
!$OMP default (shared) &
!$OMP private (k, l, atk, atl, rkl, rkl2, fs, vsr, Ecoul_damped, Fcoul, i, j, nresidk, nresidl, cm_kl)  &
!$OMP reduction (+: tmp_f_inter, Coul_inter, evdw, virial_tensor) 
do k = 1 , MM % N_of_atoms - 1
do l = k+1 , MM % N_of_atoms

     ! DWFF special pairs are treated elsewhere ...
     if ( special_pair_mtx(k,l) == 3 ) cycle

     ! for different molecules only ...
     if ( atom(k)% nr == atom(l)% nr ) cycle

     rkl(:) = atom(k) % xyz(:) - atom(l) % xyz(:)
     rkl(:) = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)

     rkl2 = sum( rkl(:)**2 )

     if( rkl2 > rcutsq ) cycle

     !-------------------------------------------------------------
     atk = atom(k)% my_intra_species_id
     atl = atom(l)% my_intra_species_id

     select case ( special_pair_mtx(k,l) )

            case(0) ! <== this is a general pair (most common)
                    ! selecting case
                    if( atom(k)%LJ .AND. atom(l)%LJ ) then

                        call Lennard_Jones( k , l , atk , atl , rkl2 , fs , vsr )

                    elseif( atom(k)%Buck .AND. atom(l)%Buck ) then

                        call Buckingham( k , l , atk , atl , rkl2 , fs , vsr )

                    endif

            case(1) !  <== this is a LJ special-pair
                    call Lennard_Jones( k , l , atk , atl , rkl2 , fs , vsr )

            case(2) !  <== this is a BUCK special-pair
                    call Buckingham( k , l , atk , atl , rkl2 , fs , vsr )

            case(3) !  <== this is a DWFF special-pair
                    CALL warning("f_inter.f: DWFF molecule reached a forbidden kernel")
                    stop 
                    
     end select

     ! Coulomb Interaction
     call Electrostatic( k , l , rkl2 , Fcoul , Ecoul_damped )

     ! attention: calculation subject to non-associative floating-point arithemtic
     tmp_f_inter(k,:) = tmp_f_inter(k,:) + (fs + Fcoul) * rkl(:)
     tmp_f_inter(l,:) = tmp_f_inter(l,:) - (fs + Fcoul) * rkl(:)

     !Energy
     evdw       = evdw + vsr
     Coul_inter = Coul_inter + Ecoul_damped
     
     !-------------------------------------------------------------------------------                                                                    
     if( using_barostat% inter ) &
     then
           nresidk = atom(k)% nr
           nresidl = atom(l)% nr
           cm_kl(:) = molecule(nresidk)% cm(:) - molecule(nresidl) % cm(:)
           cm_kl(:) = cm_kl(:) - MM % box * DNINT( cm_kl(:) * MM % ibox(:) ) * PBC(:)
           do i=1,3 ; do j=i,3
              virial_tensor(i,j) = virial_tensor(i,j) + cm_kl(i)*(fs + Fcoul)*rkl(j)
           end do; end do
     end if
     !---------------------------------------------------------------------------------

end do
end do
!$OMP end parallel do

!##############################################################################
evdw       = evdw * factor3
Coul_inter = Coul_inter * factor3
pot_INTER  = evdw + Coul_inter

! force units = J/mts = Newtons ...
do i = 1, MM % N_of_atoms
    atom(i) % f_inter_nonbond(:) = tmp_f_inter(i,:)
end do

deallocate ( tmp_f_inter )

end subroutine f_inter_nonbonding
!
!
!
!===============================================================
 subroutine Lennard_Jones( k , l , atk , atl , rkl2 , fs , vsr )
!===============================================================
implicit none
integer , intent(in)  :: k 
integer , intent(in)  :: l 
integer , intent(in)  :: atk 
integer , intent(in)  :: atl 
real*8  , intent(in)  :: rkl2 
real*8  , intent(out) :: fs
real*8  , intent(out) :: vsr

! local variables ...
integer :: n 
real*8  :: sr2 , sr6 , sr12 , eps , r_kl
logical :: flag1 , flag2 

! Lennard Jones ...

if( special_pair_mtx(k,l) == 0 ) then 

      select case ( MM % CombinationRule )
          case (2) 
              ! AMBER FF :: GMX COMB-RULE 2
              sr2 = (atom(k)%sig + atom(l)%sig)**2 / rkl2
      
          case (3)
              ! OPLS  FF :: GMX COMB-RULE 3
              sr2 = (atom(k)%sig * atom(l)%sig )**2 / rkl2
      end select
      eps = atom(k)%eps * atom(l)%eps

else
      
      read_loop: do  n = 1, size(SpecialPairs)
      
         flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                 ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(l) % MMSymbol ) )
         flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                 ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(l) % MMSymbol ) )
      
         if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
            sr2 = SpecialPairs(n)%Parms(1)**2 / rkl2
            eps = SpecialPairs(n)%Parms(2) 
            exit read_loop
         end if
      
      end do read_loop
endif

r_kl = SQRT(rkl2)
sr6  = sr2 * sr2 * sr2
sr12 = sr6 * sr6

! Forces
fs = 24.d0 * eps * ( two*sr12 - sr6 )
fs = fs/rkl2 - fscut(atk,atl)/r_kl

! LJ energy
vsr = 4.d0 * eps * ( sr12 - sr6 )
vsr = vsr - vscut(atk,atl) + fscut(atk,atl)*( r_kl - rcut )

end subroutine Lennard_Jones
!
!
!
!============================================================
 subroutine Buckingham( k , l , atk , atl , rkl2 , fs , vsr )
!============================================================
implicit none
integer , intent(in)  :: k 
integer , intent(in)  :: l 
integer , intent(in)  :: atk 
integer , intent(in)  :: atl 
real*8  , intent(in)  :: rkl2 
real*8  , intent(out) :: fs
real*8  , intent(out) :: vsr

! local variables ...
integer :: n 
real*8  :: ir2 , ir6 , ir8 , r_kl
real*8  :: A , B , C 
logical :: flag1 , flag2 

! Bukingham Potential and Forces ...

if( special_pair_mtx(k,l) == 0 ) then 

      ! Combination Rules
      A = atom(k)% BuckA * atom(l)% BuckA
      B = atom(k)% BuckB + atom(l)% BuckB
      C = atom(k)% BuckC * atom(l)% BuckC
      
else
      
      read_loop: do  n = 1, size(SpecialPairs)
      
         flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                 ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(l) % MMSymbol ) )
         flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                 ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(l) % MMSymbol ) )
      
         if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
            A = SpecialPairs(n)% Parms(1) 
            B = SpecialPairs(n)% Parms(2)
            C = SpecialPairs(n)% Parms(3)
            exit read_loop
         end if
      
      end do read_loop
end if

r_kl = SQRT(rkl2)

ir2 = 1.d0 / rkl2
ir6 = ir2 * ir2 * ir2
ir8 = ir2 * ir2 * ir2 * ir2

!Force
fs = A*B*exp(-B*r_kl)/r_kl - SIX*C*ir8
fs = fs - fscut(atk,atl)/r_kl

!Buckingham Energy
vsr = A*exp(-B*r_kl) - C*ir6
vsr = vsr - vscut(atk,atl) + fscut(atk,atl)*( r_kl - rcut )

end subroutine Buckingham
!
!
            !-----------------------------------------------------------!
            !                 Legacy Conversion procedure               ! 
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
            !                  significant figures   |                  !
            !                                       \|/                 !
            !                          applied after force calculation  !
            !-----------------------------------------------------------!
!
!===============================================================
 subroutine Electrostatic( k , l , rkl2 , Fcoul , Ecoul_damped )
!===============================================================
implicit none
integer , intent(in)  :: k 
integer , intent(in)  :: l 
real*8  , intent(in)  :: rkl2 
real*8  , intent(out) :: Fcoul
real*8  , intent(out) :: Ecoul_damped

! local variables ...
real*8 :: chrgk , chrgl , ir2 , KRIJ , expar , r_kl

! Electrostatic Coulomb Interaction ...

chrgk = atom(k)% charge
chrgl = atom(l)% charge

if( QMMM ) &
then
     chrgk = atom(k)% charge + Net_Charge(k)
     chrgl = atom(l)% charge + Net_Charge(l)
endif

r_kl  = SQRT(rkl2)

ir2   = 1.d0 / rkl2
KRIJ  = KAPPA * r_kl
expar = EXP(-KRIJ**2)

!Force
Fcoul = coulomb * chrgk * chrgl * (ir2/r_kl)
! damped Fcoul
Fcoul = Fcoul * ( ERFC(KRIJ) + TWO*irsqPI*KAPPA*r_kl*expar )
! shifted force: F_sf(R) = F(R) - F(Rc) ...
Fcoul = Fcoul - (frecut * chrgk * chrgl / r_kl)

!Energy
Ecoul_damped = (coulomb*chrgk*chrgl/r_kl) * ERFC(KRIJ)
! Coulomb shifted potential: V_sf(R) = V(R) - V(Rc) - (dV/dR)[Rc]x(R-Rc) ...
Ecoul_damped = Ecoul_damped - (vrecut*chrgk*chrgl) + (frecut*chrgk*chrgl*( r_kl-rcut ))

end subroutine Electrostatic
!
!
end module F_inter_nonbond

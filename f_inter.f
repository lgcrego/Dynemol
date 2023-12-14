module F_inter_m

    use constants_m
    use omp_lib
    use type_m       , only : warning
    use parameters_m , only : PBC , QMMM
    use for_force    , only : rcut , vrecut , frecut , rcutsq , pot_INTER , Coul_inter , &
                              Vself , evdw , vscut , fscut , KAPPA
    use md_read_m    , only : atom , MM , molecule , special_pair_mtx
    use MM_types     , only : MM_system , MM_molecular , MM_atomic , debug_MM
    use setup_m      , only : offset
    use Data_Output  , only : Net_Charge
    use gmx2mdflex   , only : SpecialPairs

    public :: FORCEINTER
    
    ! module variables ...
    real*8, public, save :: stresSR(3,3), stresRE(3,3)

contains
!
!
!
!=====================
 subroutine FORCEINTER
!=====================
implicit none

!local variables ...
real*8  , allocatable :: tmp_fsr(:,:,:) , tmp_fch(:,:,:) 
integer , allocatable :: species_offset(:)
real*8                :: rkl(3) , cm_kl(3)
real*8                :: total_chrg , Ecoul_damped , Fcoul , fs , vsr  , rkl2 
real*8                :: stresSR11 , stresSR22 , stresSR33 , stresSR12 , stresSR13 , stresSR23
real*8                :: stresRE11 , stresRE22 , stresRE33 , stresRE12 , stresRE13 , stresRE23 
integer               :: i , k , l , atk , atl
integer               :: OMP_get_thread_num , ithr , numthr , nresidl , nresidk

CALL offset( species_offset )

stresSR(:,:) = D_zero
stresRE(:,:) = D_zero

!upper triangle
stresSR11 = D_zero; stresSR12 = D_zero; stresSR13 = D_zero
stresSR22 = D_zero; stresSR23 = D_zero; stresSR33 = D_zero

!upper triangle
stresRE11 = D_zero; stresRE12 = D_zero; stresRE13 = D_zero
stresRE22 = D_zero; stresRE23 = D_zero; stresRE33 = D_zero

numthr = OMP_get_max_threads()

allocate( tmp_fsr ( MM % N_of_atoms , 3 , numthr ) , source = D_zero )
allocate( tmp_fch ( MM % N_of_atoms , 3 , numthr ) , source = D_zero )

do i = 1, MM % N_of_atoms 
    atom(i) % fsr(:) = D_zero
    atom(i) % fch(:) = D_zero
end do

evdw = D_zero
Coul_inter = D_zero

! vself part of the Coulomb calculation
total_chrg = sum(atom(:)% charge)
vself = (HALF*vrecut + rsqPI*KAPPA*coulomb) * total_chrg**2
vself = vself*factor3

!##############################################################################
! vdW and Coulomb calculations ...

!$OMP parallel DO &
!$OMP default (shared) &
!$OMP private (k, l, atk, atl, rkl2, cm_kl, rkl, fs, vsr, Ecoul_damped, Fcoul, nresidk, nresidl , ithr)  &
!$OMP reduction (+: Coul_inter, evdw, stresSR11, stresSR22, stresSR33, stresSR12, stresSR13, stresSR23,  &
!$OMP                                 stresRE11, stresRE22, stresRE33, stresRE12, stresRE13, stresRE23)
                   
do k = 1 , MM % N_of_atoms - 1
    do l = k+1 , MM % N_of_atoms

       ! for different molecules ...
        if ( atom(k)% nr /= atom(l)% nr ) then
        
            ithr = OMP_get_thread_num() + 1

            nresidk = atom(k)% nr
            nresidl = atom(l)% nr

            cm_kl(:) = molecule(nresidk) % cm(:) - molecule(nresidl) % cm(:)
            cm_kl(:) = cm_kl(:) - MM % box * DNINT( cm_kl(:) * MM % ibox(:) ) * PBC(:)

            rkl(:) = atom(k) % xyz(:) - atom(l) % xyz(:)
            rkl(:) = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)

            rkl2 = sum( rkl(:)**2 )

            if( rkl2 < rcutsq ) then

                    atk = atom(k)% my_intra_id + species_offset(atom(k)% my_species)
                    atl = atom(l)% my_intra_id + species_offset(atom(l)% my_species)

                    select case ( special_pair_mtx(k,l) )

                           case(0) ! <== not a SpecialPair
                                   if( atom(k)%LJ .AND. atom(l)%LJ ) then

                                       call Lennard_Jones( k , l , atk , atl , rkl2 , fs , vsr )

                                   elseif( atom(k)%Buck .AND. atom(l)%Buck ) then

                                       call Buckingham( k , l , atk , atl , rkl2 , fs , vsr )

                                   endif
                           case(1) 
                                   call Lennard_Jones( k , l , atk , atl , rkl2 , fs , vsr )
                           case(2)
                                   call Buckingham( k , l , atk , atl , rkl2 , fs , vsr )
                    end select

                    tmp_fsr(k,1:3,ithr) = tmp_fsr(k,1:3,ithr) + fs * rkl(1:3)
                    tmp_fsr(l,1:3,ithr) = tmp_fsr(l,1:3,ithr) - fs * rkl(1:3)
                    
                    stresSR11 = stresSR11 + cm_kl(1) * fs * rkl(1) 
                    stresSR22 = stresSR22 + cm_kl(2) * fs * rkl(2)
                    stresSR33 = stresSR33 + cm_kl(3) * fs * rkl(3)
                    stresSR12 = stresSR12 + cm_kl(1) * fs * rkl(2)
                    stresSR13 = stresSR13 + cm_kl(1) * fs * rkl(3)
                    stresSR23 = stresSR23 + cm_kl(2) * fs * rkl(3)
                    
                    ! Coulomb Interaction
                    call Electrostatic( k , l , rkl2 , Fcoul , Ecoul_damped )

                    tmp_fch(k,1:3,ithr) = tmp_fch(k,1:3,ithr) + Fcoul * rkl(1:3)
                    tmp_fch(l,1:3,ithr) = tmp_fch(l,1:3,ithr) - Fcoul * rkl(1:3)

                    stresRE11 = stresRE11 + cm_kl(1) * Fcoul * rkl(1)
                    stresRE22 = stresRE22 + cm_kl(2) * Fcoul * rkl(2)
                    stresRE33 = stresRE33 + cm_kl(3) * Fcoul * rkl(3)
                    stresRE12 = stresRE12 + cm_kl(1) * Fcoul * rkl(2)
                    stresRE13 = stresRE13 + cm_kl(1) * Fcoul * rkl(3)
                    stresRE23 = stresRE23 + cm_kl(2) * Fcoul * rkl(3)

                    !Energy
                    evdw       = evdw + vsr
                    Coul_inter = Coul_inter + Ecoul_damped
         
            end if
        end if
    end do
end do
!$OMP end parallel do

evdw       = evdw * factor3
Coul_inter = Coul_inter * factor3
pot_INTER  = evdw + Coul_inter

! ################################################################################3
!--------------------------------------------------------
stresSR11 = stresSR11 * factor3
stresSR22 = stresSR22 * factor3
stresSR33 = stresSR33 * factor3
stresSR12 = stresSR12 * factor3
stresSR13 = stresSR13 * factor3
stresSR23 = stresSR23 * factor3

stresSR(1,1) = stresSR11; stresSR(2,2) = stresSR22
stresSR(3,3) = stresSR33; stresSR(1,2) = stresSR12
stresSR(1,3) = stresSR13; stresSR(2,3) = stresSR23

stresSR(2,1) = stresSR(1,2); stresSR(3,1) = stresSR(1,3)
stresSR(3,2) = stresSR(2,3) 
!--------------------------------------------------------

!--------------------------------------------------------
stresRE11 = stresRE11 * factor3
stresRE22 = stresRE22 * factor3
stresRE33 = stresRE33 * factor3
stresRE12 = stresRE12 * factor3
stresRE13 = stresRE13 * factor3
stresRE23 = stresRE23 * factor3

stresRE(1,1) = stresRE11; stresRE(2,2) = stresRE22
stresRE(3,3) = stresRE33; stresRE(1,2) = stresRE12
stresRE(1,3) = stresRE13; stresRE(2,3) = stresRE23

stresRE(2,1) = stresRE(1,2); stresRE(3,1) = stresRE(1,3)
stresRE(3,2) = stresRE(2,3)
!--------------------------------------------------------

! force units = J/mts = Newtons ...
do i = 1, MM % N_of_atoms
    do k = 1 , numthr
        atom(i) % fsr(1:3) = atom(i) % fsr(1:3) + tmp_fsr(i,1:3,k)
        atom(i) % fch(1:3) = atom(i) % fch(1:3) + tmp_fch(i,1:3,k)
    end do
    atom(i) % f_MM(1:3) = ( atom(i) % fsr(1:3) + atom(i) % fch(1:3) ) * Angs_2_mts
end do

deallocate ( tmp_fsr , tmp_fch )

end subroutine FORCEINTER
!
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
!
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

!if( QMMM ) &
!then
!     chrgk = atom(k)% charge + Net_Charge(k)
!     chrgl = atom(l)% charge + Net_Charge(l)
!endif

r_kl  = SQRT(rkl2)

ir2   = 1.d0 / rkl2
KRIJ  = KAPPA * r_kl
expar = EXP(-KRIJ**2)

!Force
Fcoul = coulomb * chrgk * chrgl * (ir2/r_kl)
! damped Fcoul
Fcoul = Fcoul * ( ERFC(KRIJ) + TWO*rsqPI*KAPPA*r_kl*expar )
! shifted force: F_sf(R) = F(R) - F(Rc) ...
Fcoul = Fcoul - (frecut * chrgk * chrgl / r_kl)

!Energy
Ecoul_damped = (coulomb*chrgk*chrgl/r_kl) * ERFC(KRIJ)
! Coulomb shifted potential: V_sf(R) = V(R) - V(Rc) - (dV/dR)[Rc]x(R-Rc) ...
Ecoul_damped = Ecoul_damped - (vrecut*chrgk*chrgl) + (frecut*chrgk*chrgl*( r_kl-rcut ))

end subroutine Electrostatic
!
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
end module F_inter_m

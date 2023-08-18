module F_inter_m

    use constants_m
    use omp_lib
    use type_m       , only : warning
    use parameters_m , only : PBC
    use for_force    , only : forcefield , rcut , vrecut , frecut , rcutsq , pot_INTER , Coul_inter , &
                              Vself , evdw , vscut , fscut , KAPPA
    use md_read_m    , only : atom , MM , molecule , special_pair_mtx
    use MM_types     , only : MM_system , MM_molecular , MM_atomic , debug_MM
    use setup_m      , only : offset
    use gmx2mdflex   , only : SpecialPairs

    public :: FORCEINTER
    
    ! module variables ...
    real*8  , public  , save :: stresSR(3,3), stresRE(3,3)

    real*8  :: rkl2 , dkl , atk , atl
    logical :: there_are_NB_SpecialPairs = .false.

contains
!
!
!
!=====================
 subroutine FORCEINTER
!=====================
implicit none

!local variables ...
real*8  , allocatable   :: tmp_fsr(:,:,:) , tmp_fch(:,:,:) 
real*8  , allocatable   :: erfunc(:,:)
integer , allocatable   :: species_offset(:)
real*8                  :: rij(3) , rjk(3) , rkl(3) , cm_kl(3)
real*8                  :: rij2 , dij , tmp , pikap , erfkrq , chrgk , chrgl , eps , total_chrg
real*8                  :: Ecoul_damped , Fcoul , sr2 , fs , KRIJ , expar , vsr , pot
real*8                  :: stresSR11 , stresSR22 , stresSR33 , stresSR12 , stresSR13 , stresSR23
real*8                  :: stresRE11 , stresRE22 , stresRE33 , stresRE12 , stresRE13 , stresRE23 
integer                 :: i , j , k , l , j1 , j2
integer                 :: OMP_get_thread_num , ithr , numthr , nresid , nresidl , nresidk
logical                 :: flag1 , flag2 

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

allocate( erfunc ( MM % N_of_atoms , MM % N_of_atoms ) , source = D_zero )

do i = 1, MM % N_of_atoms 
    atom(i) % fsr(:) = D_zero
    atom(i) % fch(:) = D_zero
end do

pot  = D_zero
evdw = D_zero
Coul_inter = D_zero

If( allocated(SpecialPairs) ) there_are_NB_SpecialPairs = .true. 

! ##################################################################
! vself part of the Coulomb calculation

!!$OMP parallel do private(i,nresid,j1,j2,j,rij,rij2,dij,tmp) default(shared)
!do i = 1 , MM % N_of_atoms 
!
!    nresid = atom(i) % nr
!    
!    if ( molecule(nresid) % N_of_atoms > 1 ) then
!
!        j1 = sum(molecule(1:nresid-1)% N_of_atoms) + 1
!        j2 = sum(molecule(1:nresid)% N_of_atoms)
!
!        do j =  j1 , j2
!            if ( i /= j ) then
!
!                rij(:) = atom(i) % xyz(:) - atom(j) % xyz(:)
!                rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
!                rij2 = sum( rij(:)**2 )
!                dij = sqrt(rij2)
!                tmp = KAPPA * dij
!                erfunc(i,j) = Coulomb*(D_one - ERFC(tmp))/dij
!
!            end if
!        end do
!    end if
!end do
!!$OMP end parallel do
!
!pikap = (HALF*vrecut) + (rsqPI*KAPPA*coulomb)
!
!!$OMP parallel do private(i,nresid,j1,j2,j,erfkrq) default(shared) reduction( + : vself )
!do i = 1 , MM % N_of_atoms
!
!    vself  = vself + pikap*atom(i)%charge**2
!    nresid = atom(i) % nr
!
!    if ( molecule(nresid) % N_of_atoms > 1 ) then
!
!        j1 = sum(molecule(1:nresid-1)% N_of_atoms) + 1
!        j2 = sum(molecule(1:nresid)% N_of_atoms)       
!
!       do j =  j1 , j2
!          if ( i /= j ) then
!
!             erfkrq = HALF * ( erfunc(i,j) + vrecut )
!             vself  = vself + atom(i)%charge * atom(j)%charge * erfkrq
!
!          endif
!       end do
!    endif
!end do
!!$OMP end parallel do







total_chrg = sum(atom(:)% charge)
vself = (HALF*vrecut + rsqPI*KAPPA*coulomb) * total_chrg**2
vself = vself*factor3

!##############################################################################
!!$OMP parallel DO &
!!$OMP private (k, l, atk, atl, rkl2, dkl, chrgk, chrgl, sr2, KRIJ, cm_kl, rkl, fs, vsr, Ecoul_damped, &
!!$OMP          expar, Fcoul, nresidk, nresidl , ithr )                                         &
!!$OMP reduction (+ : pot, Coul_inter, evdw, stresSR11, stresSR22, stresSR33, stresSR12, stresSR13, stresSR23,  &
!!$OMP                                  stresRE11, stresRE22, stresRE33, stresRE12, stresRE13, stresRE23)
                   
! vdW and Coulomb calculations ...

do k = 1 , MM % N_of_atoms - 1
    do l = k+1 , MM % N_of_atoms

       ! for different molecules ...
        if ( atom(k) % nr /= atom(l) % nr ) then
        
            ithr    = OMP_get_thread_num() + 1

            nresidk = atom(k) % nr
            nresidl = atom(l) % nr

            cm_kl(:)  = molecule(nresidk) % cm(:) - molecule(nresidl) % cm(:)
            cm_kl(:)  = cm_kl(:) - MM % box * DNINT( cm_kl(:) * MM % ibox(:) ) * PBC(:)

            rkl(:)  = atom(k) % xyz(:) - atom(l) % xyz(:)
            rkl(:)  = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)

            rkl2    = sum( rkl(:)**2 )

            if( rkl2 < rcutsq ) then

                    atk = atom(k) % my_intra_id + species_offset(atom(k) % my_species)
                    atl = atom(l) % my_intra_id + species_offset(atom(l) % my_species)

                    dkl = SQRT(rkl2)

                    select case ( special_pair_mtx(k,l) )

                           case(0) ! <== not a SpecialPair
                                   if( atom(k)%LJ .AND. atom(l)%LJ ) then
                                       call Lennard_Jones( k , l , fs , vsr )
                                   elseif( atom(k)%Buck .AND. atom(l)%Buck ) then
                                       call Buckingham( k , l , fs , vsr )
                                   endif
                           case(1) 
                                   call Lennard_Jones( k , l , fs , vsr )
                           case(2)
                                   call Buckingham( k , l , fs , vsr )
                           case default
                                   CALL warning("unknown non-bonding special pair code in special_pair_mtx")
                                   STOP
                    end select

                    tmp_fsr(k,1:3,ithr) = tmp_fsr(k,1:3,ithr) + fs * rkl(1:3)
                    tmp_fsr(l,1:3,ithr) = tmp_fsr(l,1:3,ithr) - fs * rkl(1:3)
                    
                    stresSR11 = stresSR11 + cm_kl(1) * fs * rkl(1) 
                    stresSR22 = stresSR22 + cm_kl(2) * fs * rkl(2)
                    stresSR33 = stresSR33 + cm_kl(3) * fs * rkl(3)
                    stresSR12 = stresSR12 + cm_kl(1) * fs * rkl(2)
                    stresSR13 = stresSR13 + cm_kl(1) * fs * rkl(3)
                    stresSR23 = stresSR23 + cm_kl(2) * fs * rkl(3)
                    
                    vsr  = vsr - vscut(atk,atl) + fscut(atk,atl)*( dkl - rcut )
                    evdw = evdw + vsr*factor3

                    ! Coulomb Interaction
                    chrgk = atom(k) % charge
                    chrgl = atom(l) % charge

                    sr2   = 1.d0 / rkl2
                    KRIJ  = KAPPA * dkl
                    expar = EXP(-KRIJ**2)

                    !Force
                    Fcoul = coulomb * chrgk * chrgl * (sr2/dkl)
                    ! damped Fcoul
                    Fcoul = Fcoul * ( ERFC(KRIJ) + TWO*rsqPI*KAPPA*dkl*expar )
                    ! shifted force: F_sf(R) = F(R) - F(Rc) ...
                    Fcoul = Fcoul - (frecut * chrgk * chrgl / dkl)
                    tmp_fch(k,1:3,ithr) = tmp_fch(k,1:3,ithr) + Fcoul * rkl(1:3)
                    tmp_fch(l,1:3,ithr) = tmp_fch(l,1:3,ithr) - Fcoul * rkl(1:3)

                    stresRE11 = stresRE11 + cm_kl(1) * Fcoul * rkl(1)
                    stresRE22 = stresRE22 + cm_kl(2) * Fcoul * rkl(2)
                    stresRE33 = stresRE33 + cm_kl(3) * Fcoul * rkl(3)
                    stresRE12 = stresRE12 + cm_kl(1) * Fcoul * rkl(2)
                    stresRE13 = stresRE13 + cm_kl(1) * Fcoul * rkl(3)
                    stresRE23 = stresRE23 + cm_kl(2) * Fcoul * rkl(3)

                    !Energy
                    Ecoul_damped = (coulomb*chrgk*chrgl/dkl) * ERFC(KRIJ)
                    ! Coulomb shifted potential: V_sf(R) = V(R) - V(Rc) - (dV/dR)[Rc]x(R-Rc) ...
                    Ecoul_damped = Ecoul_damped - (vrecut*chrgk*chrgl) + (frecut*chrgk*chrgl*( dkl-rcut ))
                    Coul_inter   = Coul_inter + Ecoul_damped*factor3
         
                    pot = pot + (vsr + Ecoul_damped) *factor3
                    
            end if
        end if
    end do
end do
!!$OMP end parallel do

pot_INTER = pot

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

deallocate ( tmp_fsr , tmp_fch , erfunc )

end subroutine FORCEINTER
!
!
!
!
!============================================
 subroutine Lennard_Jones( k , l , fs , vsr )
!============================================
implicit none
integer , intent(in)  :: k 
integer , intent(in)  :: l 
real*8  , intent(out) :: fs
real*8  , intent(out) :: vsr

! local variables ...
integer :: n 
real*8  :: sr2 , sr6 , sr12 , eps
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

sr6  = sr2 * sr2 * sr2
sr12 = sr6 * sr6

! Forces
fs = 24.d0 * eps * ( 2.d0 * sr12 - sr6 )
fs = fs/rkl2 - fscut(atk,atl)/dkl

! LJ energy
vsr = 4.d0 * eps * ( sr12 - sr6 )

end subroutine Lennard_Jones
!
!
!
!
!=========================================
 subroutine Buckingham( k , l , fs , vsr )
!=========================================
implicit none
integer , intent(in)  :: k 
integer , intent(in)  :: l 
real*8  , intent(out) :: fs
real*8  , intent(out) :: vsr

! local variables ...
integer :: n 
real*8  :: sr2 , sr6 , sr12 , eps
real*8  :: A , B , C , sr8
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

sr2 = 1.d0 / rkl2
sr6 = sr2 * sr2 * sr2
sr8 = sr2 * sr2 * sr2 * sr2

!Force
fs = A*B*exp(-B*dkl)/dkl - SIX*C*sr8
fs = fs - fscut(atk,atl)/dkl

!Buckingham Energy
vsr = A*exp(-B*dkl) - C*sr6

end subroutine Buckingham
!
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

module FF_intra_nonbond

    use type_m   
    use omp_lib
    use constants_m
    use parameters_m       , only : PBC
    use syst               , only : using_barostat
    use for_force          , only : rcut, vrecut, frecut, vscut, fscut, KAPPA, LJ_14, LJ_intra, Coul_14, Coul_intra, rcutsq
    use MD_read_m          , only : atom, molecule, MM 
    use gmx2mdflex         , only : SpecialPairs, SpecialPairs14
    use Berendsen_Barostat , only : virial_tensor

    private

    public :: f_intra_nonbonding

    ! module variables ...
    logical :: there_are_NB_SpecialPairs   = .false.
    logical :: there_are_NB_SpecialPairs14 = .false.

contains
!
!
!=============================
subroutine f_intra_nonbonding
!=============================
implicit none

! local_variables ...
real*8 , allocatable :: tmp_vdw(:,:,:), tmp_ele(:,:,:)
real*8 , dimension (3):: rij 
real*8  :: rij2, dij, fs, Fcoul, E_vdw, E_coul 
integer :: i, j, k, l, n, ithr, numthr, ati, atj 
logical :: flag1 , flag2


do j = 1 , MM % N_of_atoms
   atom(j) % fnonbd14(:) = D_zero     ! Non-bonded 1-4
   atom(j) % fnonbd(:)   = D_zero     ! Non-bonded Intramolecular
   atom(j) % fnonch14(:) = D_zero     ! Non-bonded coulomb 1-4
   atom(j) % fnonch(:)   = D_zero     ! Non-bonded coulomb Intramolecular
end do
LJ_14      = D_zero
LJ_intra   = D_zero
Coul_14    = D_zero
Coul_intra = D_zero

numthr = OMP_get_max_threads()

!====================================================================
! Nonbonding Interactions
! parameters:
! rcut = cutoff radius for vdW and Coulomb interactions (defined in card.intp)
! rcutsq = rcut*rcut
! FScut(i,j) = LJ/Bck Force at a spherical Surface of radius rcut for each MM atom pair
! VScut(i,j) = LJ/Bck energy at a spherical Surface of radius rcut for each MM atom pair
! rij(:)/dij is the direction versor of the atomic pair ij
!====================================================================
! Non-bonded 1,4 intramolecular interactions ...

If( allocated(SpecialPairs14) ) there_are_NB_SpecialPairs14 = .true.

do i = 1 , MM % N_of_molecules
    do j = 1 , molecule(i) % Nbonds14

        ati = molecule(i) % bonds14(j,1)
        atj = molecule(i) % bonds14(j,2)

        if ( atom(atj) % flex .OR. atom(ati) % flex ) then

            rij(:) = atom(ati) % xyz(:) - atom(atj) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)

            rij2 = sum( rij(:)**2 )

            call Lennard_Jones_14 (rij2 , ati , atj , fs)
            atom(ati) % fnonbd14(1:3) = atom(ati) % fnonbd14(1:3) + fs * rij(1:3)
            atom(atj) % fnonbd14(1:3) = atom(atj) % fnonbd14(1:3) - fs * rij(1:3)


            call Coulomb14 (rij2 , ati , atj , Fcoul)
            atom(ati) % fnonch14(1:3) = atom(ati) % fnonch14(1:3) + Fcoul * rij(1:3)
            atom(atj) % fnonch14(1:3) = atom(atj) % fnonch14(1:3) - Fcoul * rij(1:3)

        end if
    end do
end do
 
!====================================================================
! Non-bonding intramolecular interactions ...
! Van der Walls , Buckingham , Electrostatic

if( using_barostat% intra ) call InitializeStressMatrix

allocate( tmp_vdw (MM%N_of_atoms,3,numthr) , source = D_zero )
allocate( tmp_ele (MM%N_of_atoms,3,numthr) , source = D_zero )

If( allocated(SpecialPairs) ) there_are_NB_SpecialPairs = .true.

do i = 1 , MM % N_of_molecules
     if ( molecule(i)% NintraIJ==0 .or. molecule(i)% DWFF ) cycle
     !$OMP parallel DO &
     !$OMP default (shared) &
     !$OMP private (j , ithr , ati , atj , rij , rij2 , fs , Fcoul , E_vdw , E_coul)  &
     !$OMP reduction (+: LJ_intra, Coul_intra, virial_tensor)
     do j = 1 , molecule(i) % NintraIJ

         ithr = OMP_get_thread_num() + 1

         ati  = molecule(i) % IntraIJ(j,1) 
         atj  = molecule(i) % IntraIJ(j,2) 

         if ( atom(atj) % flex .OR. atom(ati) % flex ) then

             rij(:) = atom(ati) % xyz(:) - atom(atj) % xyz(:)
             rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)

             rij2 = sum( rij(:)**2 )

             if ( rij2 < rcutsq ) then

                  select case( molecule(i)% intraIJ(j,3) )

                         case(1)
                         call Lennard_Jones( rij2 , ati , atj , fs , E_vdw )

                         case(2) 
                         call Buckingham( rij2 , ati , atj , fs , E_vdw )

                         case(0)  ! <== electrostatic
                         fs    = d_zero
                         E_vdw = d_zero

                  end select

                  tmp_vdw(ati,1:3,ithr) = tmp_vdw(ati,1:3,ithr) + fs * rij(:)
                  tmp_vdw(atj,1:3,ithr) = tmp_vdw(atj,1:3,ithr) - fs * rij(:)

                  LJ_intra = LJ_intra + E_vdw*factor3 ! <== LJ and/or Buck energy

                  ! Coulomb Interaction
                  call Coulomb_DSF( rij2 , ati , atj , Fcoul , E_coul )

                  tmp_ele(ati,1:3,ithr) = tmp_ele(ati,1:3,ithr) + Fcoul * rij(:)
                  tmp_ele(atj,1:3,ithr) = tmp_ele(atj,1:3,ithr) - Fcoul * rij(:)

                  Coul_intra = Coul_intra + E_coul*factor3 ! <== Coulomb energy

                  !-------------------------------------------------------------------------------
                  if( using_barostat% intra ) &
                  then
                       do k=1,3 ; do l=k,3
                          virial_tensor(k,l) = virial_tensor(k,l) + rij(k) * (fs+Fcoul) * rij(l)
                       end do; end do
                  end if
                  !---------------------------------------------------------------------------------

             end if
         end if
     end do
     !$OMP end parallel do
end do

! force units = J/mts = Newtons ...    
! manual reduction (+: fnonbd , fnonch) ...
do i = 1, MM % N_of_atoms
    do k=1,numthr              
        atom(i)% fnonbd(1:3) = atom(i)% fnonbd(1:3) + tmp_vdw(i,1:3,k)
        atom(i)% fnonch(1:3) = atom(i)% fnonch(1:3) + tmp_ele(i,1:3,k)
    enddo
end do    

!====================================================================
! dissociative forces

if( using_barostat% intra ) call ConcludeStressMatrix

deallocate( tmp_vdw , tmp_ele )

end subroutine f_intra_nonbonding
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
!==========================================================
 subroutine Coulomb_DSF (rij2 , ati , atj , fcoul , E_coul)
!==========================================================
implicit none
real*8  , intent(in)  :: rij2
integer , intent(in)  :: ati
integer , intent(in)  :: atj
real*8  , intent(out) :: fcoul
real*8  , intent(out) :: E_coul

! local variables ...
real*8 :: chrgi , chrgj , ir2 , KRIJ , expar , dij

! Coulomb damped shifted field, V_dsf ...
! the damped and cutoff-neutralized Coulomb potential ...
! incorporates both image charges at Rc surface and damping of electrostatic interaction ...
! for details: JPC 124, 234104 (2006)

chrgi = atom(ati)% charge
chrgj = atom(atj)% charge

dij   = sqrt(rij2)
ir2   = 1.d0 / rij2
KRIJ  = KAPPA * dij
expar = EXP(-KRIJ**2)

!Force
Fcoul = coulomb * chrgi * chrgj * (ir2/dij)
! damped Fcoul
Fcoul = Fcoul * ( ERFC(KRIJ) + TWO*irsqPI*KAPPA*dij*expar )
! shifted force: F_sf(R) = F(R) - F(Rc) ...
Fcoul = Fcoul - (frecut * chrgi * chrgj / dij)   

!Energy
E_coul = (coulomb*chrgi*chrgj/dij) * ERFC(KRIJ)
! Coulomb shifted potential: V_sf(R) = V(R) - V(Rc) - (dV/dR)[Rc]x(R-Rc) ...
E_coul = E_coul - (vrecut*chrgi*chrgj) + (frecut*chrgi*chrgj*( dij-rcut ))

end subroutine Coulomb_DSF
!
!
!
!
!===============================================
 subroutine Coulomb14 (rij2 , ati , atj , Fcoul)
!===============================================
implicit none
real*8  , intent(in)  :: rij2
integer , intent(in)  :: ati
integer , intent(in)  :: atj
real*8  , intent(out) :: Fcoul

! local variables ...
real*8 :: chrgi , chrgj , ir2 , KRIJ , expar , dij , E_coul

! Coulomb damped field for 1-4 pairs ...

chrgi = atom(ati)% charge
chrgj = atom(atj)% charge

dij = sqrt(rij2)
ir2   = 1.d0 / rij2
KRIJ  = KAPPA * dij
expar = EXP(-KRIJ**2)

! for preserving 1-4 Coulomb interactions mind the parameters below
!erfc(kappa*dij) ~ 0.98 for kappa = 0.005 and dij = 2.5

!Force
Fcoul = coulomb * chrgi * chrgj * ( ir2/dij )
! damped Fcoul
Fcoul = Fcoul * ( ERFC(KRIJ) + TWO*irsqPI*KAPPA*dij*expar ) * MM%fudgeQQ

!Energy 
E_coul = (coulomb*chrgi*chrgj/dij) * ERFC(KRIJ) * MM%fudgeQQ

Coul_14 = Coul_14 + E_coul*factor3

end subroutine Coulomb14
!
!
!
!
!========================================================
 subroutine Lennard_Jones (rij2 , ati , atj , fs , E_vdw) 
!========================================================
implicit none
real*8  , intent(in)  :: rij2
integer , intent(in)  :: ati
integer , intent(in)  :: atj
real*8  , intent(out) :: fs
real*8  , intent(out) :: E_vdw

! local variables ...  
integer :: n , ati1 , atj1
real*8  :: sr2 , sr6 , sr12 , dij , eps
logical :: flag1 , flag2

! Lennard Jones ...
select case ( MM % CombinationRule )

    case (2)
        ! AMBER FF :: GMX COMB-RULE 2
        sr2 = (atom(ati)%sig + atom(atj)%sig)**2 / rij2

    case (3)
        ! OPLS  FF :: GMX COMB-RULE 3
        sr2 = (atom(ati)%sig * atom(atj)%sig)**2 / rij2

end select
eps = atom(ati)%eps * atom(atj)%eps

If( there_are_NB_SpecialPairs ) then    ! <== check whether (I,J) is a SpecialPair ... 

   read_loop: do  n = 1, size(SpecialPairs)

      flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
              ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(atj) % MMSymbol ) )
      flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
              ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(atj) % MMSymbol ) )

      if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
         sr2 = SpecialPairs(n)%Parms(1)**2 / rij2 
         eps = SpecialPairs(n)%Parms(2) 
         exit read_loop
      end if

   end do read_loop

end if

dij  = sqrt(rij2)
sr6  = sr2 * sr2 * sr2
sr12 = sr6 * sr6

!Forces
fs = 24.d0 * eps * ( TWO * sr12 - sr6 )

! shifted force: F_sf(R) = F(R) - F(Rc) ...
ati1 = atom(ati) % my_intra_species_id
atj1 = atom(atj) % my_intra_species_id 

fs = fs/rij2 - fscut(ati1,atj1)/dij     

!Energy
E_vdw = 4.d0 * eps * ( sr12 - sr6 )
! LJ shifted potential: V_sf(R) = V(R) - V(Rc) - (dV/dR)[Rc]x(R-Rc) ...
E_vdw = E_vdw - vscut(ati1,atj1) + fscut(ati1,atj1)*( dij - rcut ) 

end subroutine Lennard_Jones
!
!
!
!
!===================================================
 subroutine Lennard_Jones_14 (rij2 , ati , atj , fs)
!===================================================
implicit none
real*8  , intent(in)  :: rij2
integer , intent(in)  :: ati
integer , intent(in)  :: atj
real*8  , intent(out) :: fs

! local variables ...
integer :: n 
real*8  :: sr2 , sr6 , sr12 , eps , E_vdw
logical :: flag1 , flag2

! Lennard Jones for 1-4 pairs ...
select case ( MM % CombinationRule )

    case (2)
        ! AMBER FF :: GMX COMB-RULE 2
        sr2 = (atom(ati)%sig14 + atom(atj)%sig14)**2 / rij2

    case (3)
        ! OPLS  FF :: GMX COMB-RULE 3
        sr2 = (atom(ati)%sig14 * atom(atj)%sig14)**2 / rij2

end select
eps = atom(ati)%eps14 * atom(atj)%eps14 

If( there_are_NB_SpecialPairs14 ) then    ! <== check whether (I,J) is a SpecialPair ... 

   read_loop1: do  n = 1, size(SpecialPairs14)

       flag1 = ( adjustl( SpecialPairs14(n) % MMSymbols(1) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
               ( adjustl( SpecialPairs14(n) % MMSymbols(2) ) == adjustl( atom(atj) % MMSymbol ) )
       flag2 = ( adjustl( SpecialPairs14(n) % MMSymbols(2) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
               ( adjustl( SpecialPairs14(n) % MMSymbols(1) ) == adjustl( atom(atj) % MMSymbol ) )

       if ( flag1 .OR. flag2 ) then       ! <== apply SpecialPair parms ... 
           sr2 = SpecialPairs14(n)%Parms(1)**2 / rij2
           eps = SpecialPairs14(n)%Parms(2) 
           exit read_loop1
       end if
    
   end do read_loop1

end if

sr6  = sr2 * sr2 * sr2
sr12 = sr6 * sr6

!Forces
fs = 24.d0 * eps * ( TWO * sr12 - sr6 )
fs = (fs / rij2) * MM % fudgeLJ

!Energy
E_vdw = 4.d0 * eps * ( sr12 - sr6 ) * MM % fudgeLJ

LJ_14 = LJ_14 + E_vdw*factor3  ! <== LJ energy

end subroutine Lennard_Jones_14
!
!
!
!
!=====================================================
 subroutine Buckingham (rij2 , ati , atj , fs , E_vdw)
!=====================================================
implicit none
real*8  , intent(in)  :: rij2
integer , intent(in)  :: ati
integer , intent(in)  :: atj
real*8  , intent(out) :: fs
real*8  , intent(out) :: E_vdw

! local variables ...
real*8  :: Aij , Bij , Cij
integer :: n , ati1 , atj1
real*8  :: ir2 , ir6 , ir8 , dij
logical :: flag1 , flag2
                                
! Bukingham Potential and Forces ...

! Combination Rules

Aij = atom(ati)% BuckA * atom(atj)% BuckA

Bij = atom(ati)% BuckB + atom(atj)% BuckB  

Cij = atom(ati)% BuckC * atom(atj)% BuckC

If( there_are_NB_SpecialPairs ) then    ! <== check whether (I,J) is a SpecialPair ... 

   read_loop: do  n = 1, size(SpecialPairs)

      flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
              ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(atj) % MMSymbol ) )
      flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(ati) % MMSymbol ) ) .AND. &
              ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(atj) % MMSymbol ) )

      if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
         Aij = SpecialPairs(n)% Parms(1) 
         Bij = SpecialPairs(n)% Parms(2)
         Cij = SpecialPairs(n)% Parms(3)
         exit read_loop
      end if

   end do read_loop

end if

dij = sqrt(rij2)
ir2 = 1.d0 / rij2
ir6 = ir2 * ir2 * ir2 
ir8 = ir2 * ir2 * ir2 * ir2

!Force
fs = Aij*Bij*exp(-Bij*dij) / dij - SIX*Cij*ir8

! shifted force: F_sf(R) = F(R) - F(Rc) ...
ati1 = atom(ati) % my_intra_species_id
atj1 = atom(atj) % my_intra_species_id 

fs = fs - fscut(ati1,atj1)/dij     

!Energy
E_vdw = Aij*exp(-Bij*dij) - Cij*ir6
! Buckingham shifted potential: V_sf(R) = V(R) - V(Rc) - (dV/dR)[Rc]x(R-Rc) ...
E_vdw = E_vdw - vscut(ati1,atj1) + fscut(ati1,atj1)*( dij - rcut ) 

end subroutine Buckingham
!
!
!
!                         
!                            
!=================================
 subroutine InitializeStressMatrix
!=================================
  implicit none
  !-----------------------------------------------------
  ! initializing variables for this integration step ...
  !-----------------------------------------------------
  virial_tensor(:,:)   = D_zero
  !-----------------------------------------------------
end subroutine InitializeStressMatrix  
!
!===============================
 subroutine ConcludeStressMatrix
!===============================
  implicit none

  ! local variables
  integer :: i,j
  !-------------------------------------------------
  ! symmetrizing the tensors ...
  !-------------------------------------------------
  virial_tensor   = virial_tensor * factor3

  do concurrent (i = 1:2, j = 1:3, j>i)
    virial_tensor(j,i)   = virial_tensor(i,j) 
  end do
  !-------------------------------------------------
end subroutine ConcludeStressMatrix
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
end module FF_intra_nonbond

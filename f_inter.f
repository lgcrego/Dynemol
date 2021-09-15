module F_inter_m

    use constants_m
    use omp_lib
    use parameters_m , only : PBC
    use for_force    , only : forcefield , rcut , vrecut , frecut , rcutsq , pot_INTER , ecoul , &
                              eintra , evdw , vscut , fscut , KAPPA
    use MD_read_m    , only : atom , MM , molecule
    use MM_types     , only : MM_system , MM_molecular , MM_atomic , debug_MM
    use setup_m      , only : offset
    use gmx2mdflex   , only : SpecialPairs

    public :: FORCEINTER
    
    ! module variables ...
    real*8  , public  , save :: stressr(3,3), stresre(3,3)
    logical                  :: there_are_NB_SpecialPairs = .false.

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
real*8  , allocatable   :: erfkr(:,:)
integer , allocatable   :: species_offset(:)
real*8                  :: rij(3) , rjk(3) , rkl(3)
real*8                  :: rjkq , rklq , rjksq , rklsq , tmp , pikap , erfkrq , chrgk , chrgl , eps 
real*8                  :: vreal , freal , sr2 , sr6 , sr12 , fs , KRIJ , expar , vsr , vself , pot
real*8                  :: stressr11 , stressr22 , stressr33 , stressr12 , stressr13 , stressr23
real*8                  :: stresre11 , stresre22 , stresre33 , stresre12 , stresre13 , stresre23 
integer                 :: i , j , k , l , n , atk , atl , j1 , j2
integer                 :: OMP_get_thread_num , ithr , numthr , nresid , nresidl , nresidk
logical                 :: flag1 , flag2 

CALL offset( species_offset )

stressr(:,:) = D_zero
stresre(:,:) = D_zero

stressr11 = D_zero; stressr22 = D_zero; stressr33 = D_zero
stressr12 = D_zero; stressr13 = D_zero; stressr23 = D_zero
stresre11 = D_zero; stresre22 = D_zero; stresre33 = D_zero
stresre12 = D_zero; stresre13 = D_zero; stresre23 = D_zero

numthr = OMP_get_max_threads()

allocate( tmp_fsr ( MM % N_of_atoms , 3 , numthr ) , source = D_zero )
allocate( tmp_fch ( MM % N_of_atoms , 3 , numthr ) , source = D_zero )

allocate( erfkr ( MM % N_of_atoms , MM % N_of_atoms ) , source = D_zero )

do i = 1, MM % N_of_atoms 
    atom(i) % fsr(:) = D_zero
    atom(i) % fch(:) = D_zero
end do

pot    = D_zero
vself  = D_zero 
ecoul  = D_zero
evdw   = D_zero
eintra = D_zero

If( allocated(SpecialPairs) ) there_are_NB_SpecialPairs = .true. 

! ##################################################################
! vself part of the Coulomb calculation
!$OMP parallel do private(i,nresid,j1,j2,j,rjk,rjkq,rjksq,tmp) default(shared)
do i = 1 , MM % N_of_atoms 

    nresid = atom(i) % nr
    
    if ( molecule(nresid) % N_of_atoms > 1 ) then

        j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
        j2 = sum(molecule(1:nresid) % N_of_atoms)

        do j =  j1 , j2
            if ( i /= j ) then

                rjk(:)     = atom(i) % xyz(:) - atom(j) % xyz(:)
                rjk(:)     = rjk(:) - MM % box(:) * DNINT( rjk(:) * MM % ibox(:) ) * PBC(:)
                rjkq       = sum( rjk(:) * rjk(:) )
                rjksq      = sqrt(rjkq)
                tmp        = KAPPA * rjksq
                erfkr(i,j) = ( D_one - ERFC(tmp) ) / rjksq
                erfkr(i,j) = erfkr(i,j) * coulomb * factor3

            end if
        end do
    end if
end do
!$OMP end parallel do

pikap = HALF * vrecut + rsqpi * KAPPA * coulomb * factor3

!$OMP parallel do private(i,nresid,j1,j2,j,erfkrq) default(shared) reduction( + : vself )
do i = 1 , MM % N_of_atoms

    vself  = vself + pikap * atom(i) % charge * atom(i) % charge
    nresid = atom(i) % nr

    if ( molecule(nresid) % N_of_atoms > 1 ) then

        j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
        j2 = sum(molecule(1:nresid) % N_of_atoms)       

       do j =  j1 , j2
          if ( i /= j ) then

             erfkrq = HALF * ( erfkr(i,j) + vrecut )
             vself  = vself + atom(i) % charge * atom(j) % charge * erfkrq

          endif
       end do
    endif
end do
!$OMP end parallel do

eintra = eintra + vself

!##############################################################################
!$OMP parallel DO &
!$OMP private (k, l, atk, atl, rklq, rklsq, chrgk, chrgl, sr2, sr6, sr12, KRIJ, rij, rkl, fs, vsr, vreal, &
!$OMP          expar, freal, nresidk, nresidl , ithr , eps , n , flag1 , flag2 )                          &
!$OMP reduction (+ : pot, ecoul, evdw, stressr11, stressr22, stressr33, stressr12, stressr13, stressr23,  &
!$OMP                                  stresre11, stresre22, stresre33, stresre12, stresre13, stresre23)
                   
! LJ and Coulomb calculation

do k = 1 , MM % N_of_atoms - 1
    do l = k , MM % N_of_atoms

       ! do it for different molecules ...
        if ( atom(k) % nr /= atom(l) % nr ) then
        
            ithr    = OMP_get_thread_num() + 1

            nresidk = atom(k) % nr
            nresidl = atom(l) % nr

            rij(:)  = molecule(nresidk) % cm(:) - molecule(nresidl) % cm(:)
            rij(:)  = rij(:) - MM % box * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)

            chrgk   = atom(k) % charge
            chrgl   = atom(l) % charge

            rkl(:)  = atom(k) % xyz(:) - atom(l) % xyz(:)
            rkl(:)  = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)

            rklq    = sum( rkl(:) * rkl(:) )

            if( rklq < rcutsq ) then

                select case(forcefield)

                    case( 1 )
                    ! Born-Mayer ; short range ...

                    case( 2 )
                        ! Lennard Jones ; short range ...

                        atk = atom(k) % my_intra_id + species_offset(atom(k) % my_species)
                        atl = atom(l) % my_intra_id + species_offset(atom(l) % my_species)
                       
                    select case ( MM % CombinationRule )

                        case (2) 
                            ! AMBER FF :: GMX COMB-RULE 2

                            sr2   = ( ( atom(k) % sig + atom(l) % sig ) * ( atom(k) % sig + atom(l) % sig ) ) / rklq

                        case (3)
                            ! OPLS  FF :: GMX COMB-RULE 3

                            sr2   = ( ( atom(k) % sig * atom(l) % sig ) * ( atom(k) % sig * atom(l) % sig ) ) / rklq

                    end select
                    eps = atom(k) % eps * atom(l) % eps

                    If( there_are_NB_SpecialPairs ) then    ! <== check whether (K,L) is a SpecialPair ... 

                       read_loop: do  n = 1, size(SpecialPairs)

                          flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                                  ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(l) % MMSymbol ) )
                          flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( atom(k) % MMSymbol ) ) .AND. &
                                  ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( atom(l) % MMSymbol ) )

                          if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
                             sr2 = ( SpecialPairs(n)%Parms(1) * SpecialPairs(n)%Parms(1) ) / rklq
                             eps = SpecialPairs(n) % Parms(2) 
                             exit read_loop
                          end if

                       end do read_loop

                    end if
                   
                    sr6  = sr2 * sr2 * sr2
                    sr12 = sr6 * sr6
                    
                    rklsq   = SQRT(rklq)

                    ! factor3 is used here in the force calculation because fscut was multiplied by it in md_setup ...
                    fs   = 24.d0 * ( eps * factor3 ) * ( 2.d0 * sr12 - sr6 )
                    fs   = fs / rklq - fscut(atk,atl) / rklsq

                    vsr  = 4.d0 * ( eps * factor3 ) * ( sr12 - sr6 )
                    vsr  = vsr - vscut(atk,atl) + fscut(atk,atl) * ( rklsq - rcut )

                    pot  = pot + vsr
                    evdw = evdw + vsr

                    stressr11 = stressr11 + rij(1) * fs * rkl(1)
                    stressr22 = stressr22 + rij(2) * fs * rkl(2)
                    stressr33 = stressr33 + rij(3) * fs * rkl(3)
                    stressr12 = stressr12 + rij(1) * fs * rkl(2)
                    stressr13 = stressr13 + rij(1) * fs * rkl(3)
                    stressr23 = stressr23 + rij(2) * fs * rkl(3)
                   
                    ! compensating the factor3 with 1.0d20 ... 
                    fs = fs * 1.0d20
                    tmp_fsr(k,1:3,ithr) = tmp_fsr(k,1:3,ithr) + fs * rkl(1:3)
                    tmp_fsr(l,1:3,ithr) = tmp_fsr(l,1:3,ithr) - fs * rkl(1:3)

                    ! Coulomb Part
                    sr2   = 1.d0 / rklq
                    KRIJ  = KAPPA * rklsq
                    expar = EXP( -(KRIJ*KRIJ) )

                    freal = coulomb * chrgk * chrgl * (sr2 / rklsq)
                    freal = freal * ( ERFC(KRIJ) + 2.d0 * rsqpi * KAPPA * rklsq * expar )
                    freal = freal - frecut / rklsq * chrgk * chrgl

                    vreal = coulomb * factor3 * chrgk * chrgl * ERFC(KRIJ)/rklsq
                    vreal = vreal - vrecut * chrgk * chrgl + frecut * chrgk * chrgl * ( rklsq-rcut ) * factor3

                    pot   = pot + vreal
                    ecoul = ecoul + vreal
         
                    stresre11 = stresre11 + rij(1) * freal * rkl(1)
                    stresre22 = stresre22 + rij(2) * freal * rkl(2)
                    stresre33 = stresre33 + rij(3) * freal * rkl(3)
                    stresre12 = stresre12 + rij(1) * freal * rkl(2)
                    stresre13 = stresre13 + rij(1) * freal * rkl(3)
                    stresre23 = stresre23 + rij(2) * freal * rkl(3)

                    tmp_fch(k,1:3,ithr) = tmp_fch(k,1:3,ithr) + freal * rkl(1:3)
                    tmp_fch(l,1:3,ithr) = tmp_fch(l,1:3,ithr) - freal * rkl(1:3)

                end select
            end if
        end if
    end do
end do
!$OMP end parallel do

! ################################################################################3

pot = pot - vself
pot_INTER = pot

stresre11 = stresre11 * factor3
stresre22 = stresre22 * factor3
stresre33 = stresre33 * factor3
stresre12 = stresre12 * factor3
stresre13 = stresre13 * factor3
stresre23 = stresre23 * factor3

stressr(1,1) = stressr11; stressr(2,2) = stressr22
stressr(3,3) = stressr33; stressr(1,2) = stressr12
stressr(1,3) = stressr13; stressr(2,3) = stressr23
stresre(1,1) = stresre11; stresre(2,2) = stresre22
stresre(3,3) = stresre33; stresre(1,2) = stresre12
stresre(1,3) = stresre13; stresre(2,3) = stresre23
stressr(2,1) = stressr(1,2); stressr(3,1) = stressr(1,3)
stressr(3,2) = stressr(2,3); stresre(2,1) = stresre(1,2)
stresre(3,1) = stresre(1,3); stresre(3,2) = stresre(2,3)

! force units = J/mts = Newtons ...
do i = 1, MM % N_of_atoms
    do k = 1 , numthr
        atom(i) % fsr(1:3) = atom(i) % fsr(1:3) + tmp_fsr(i,1:3,k)
        atom(i) % fch(1:3) = atom(i) % fch(1:3) + tmp_fch(i,1:3,k)
    end do
    atom(i) % f_MM(1:3) = ( atom(i) % fsr(1:3) + atom(i) % fch(1:3) ) * Angs_2_mts
end do

deallocate ( tmp_fsr , tmp_fch , erfkr )

end subroutine FORCEINTER
!
!
!
!===================
 function ERFC ( X )
!===================
 real*8 :: ERFC
 real*8 :: A1, A2, A3, A4, A5, P, T, X, XSQ, TP
 parameter ( A1 = 0.254829592d0, A2 = -0.284496736d0 ) 
 parameter ( A3 = 1.421413741d0, A4 = -1.453122027d0 ) 
 parameter ( A5 = 1.061405429d0, P  =  0.3275911d0   ) 

 T    = 1.0d0 / ( 1.0d0 + P * X )
 XSQ  = X * X
 TP   = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
 ERFC = TP * EXP ( -XSQ )

end function ERFC
!
!
end module F_inter_m

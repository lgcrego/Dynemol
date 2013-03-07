module F_inter_m

    use constants_m
    use omp_lib
    use for_force 
    use MD_read_m   , only : atom , MM , molecule
    use project     , only : MM_system , MM_molecular , MM_atomic

    public :: FORCEINTER
    
    ! module variables ...
    real*8  , public  , save :: stressr(3,3), stresre(3,3)

contains
!
!
!
!====================
subroutine FORCEINTER
!====================
implicit none

! local variables ...
 integer :: i, j, k, l, atk, atl, j1, j2
 integer :: OMP_get_thread_num, ithr, numthr, nresid, nresidl, nresidk
 real*8  :: rjkq, rklq, rjksq, rklsq, tmp, pikap, erfkrq, chrgk, chrgl, vsr, vreal, freal, sr2, sr6, sr12, fs, KRIJ, expar     
 real*8  :: stressr11, stressr22, stressr33, stressr12, stressr13, stressr23
 real*8  :: stresre11, stresre22, stresre33, stresre12, stresre13, stresre23 
 real*8, dimension (3) :: rij, rjk, rkl

 stressr(:,:) = 0.d0; stresre(:,:) = 0.d0
 stressr11 = 0.d0; stressr22 = 0.d0; stressr33 = 0.d0
 stressr12 = 0.d0; stressr13 = 0.d0; stressr23 = 0.d0
 stresre11 = 0.d0; stresre22 = 0.d0; stresre33 = 0.d0
 stresre12 = 0.d0; stresre13 = 0.d0; stresre23 = 0.d0

 numthr = OMP_get_max_threads()

 allocate ( tmp_fsr(MM % N_of_atoms,3,numthr)      , source=0.d0 )
 allocate ( tmp_fch(MM % N_of_atoms,3,numthr)      , source=0.d0 )
 allocate ( erfkr(MM % N_of_atoms,MM % N_of_atoms) , source=0.d0 )

 do i = 1, MM % N_of_atoms 
    atom(i) % fsr(:)    = 0.0d0
    atom(i) % fch(:)    = 0.0d0
 end do

 pot    = 0.d0
 vself  = 0.d0 
 ecoul  = 0.d0
 evdw   = 0.d0
 eintra = 0.d0

! ##################################################################
! vself part of the Coulomb calculation

 do i = 1, MM % N_of_atoms 
    nresid = atom(i) % nresid
    if ( molecule(nresid) % N_of_atoms > 1 ) then
       j1 = ( nresid - 1 ) * molecule(nresid) % N_of_atoms + 1 
       j2 = nresid * molecule(nresid) % N_of_atoms
       do j =  j1 , j2  
          if ( i /= j ) then
             rjk(:) = atom(i) % xyz(:) - atom(j) % xyz(:)
             rjk(:) = rjk(:) - MM % box(:) * ANINT( rjk(:) * MM % ibox(:) )
             rjkq   = sum( rjk(:) * rjk(:) )
             rjksq  = SQRT(rjkq)
             tmp = KAPPA * rjksq
             erfkr(i,j) = ( 1.d0 - ERFC(tmp) ) / rjksq
             erfkr(i,j) = erfkr(i,j) * coulomb * 1.d-20
          end if
       end do
    end if
 end do

 pikap = 0.5d0 * vrecut + rsqpi * KAPPA * coulomb * 1.d-20

 do i = 1, MM % N_of_atoms
    vself = vself + pikap * atom(i) % charge * atom(i) % charge
    nresid = atom(i) % nresid 
    if ( molecule(nresid) % N_of_atoms > 1 ) then
       j1 = ( nresid - 1 ) * molecule(nresid) % N_of_atoms + 1
       j2 = nresid * molecule(nresid) % N_of_atoms
       do j =  j1 , j2
          if ( i /= j ) then
             erfkrq = 0.5d0 * ( erfkr(i,j) + vrecut )
             vself  = vself + atom(i) % charge * atom(j) % charge * erfkrq
          endif
       end do
    endif
 end do
 eintra = eintra + vself

!##############################################################################
!$OMP parallel DO schedule(GUIDED,200)                                                                                              &
!$OMP private (k, l, rklq, rklsq, chrgk, chrgl, sr2, sr6, sr12, KRIJ, rij, rkl, fs, vsr, vreal, expar, freal, nresidk, nresidl)     &
!$OMP reduction (+ : pot, ecoul, stressr11, stressr22, stressr33, stressr12, stressr13, stressr23,                                  &
!$OMP                            stresre11, stresre22, stresre33, stresre12, stresre13, stresre23)

! LJ and Coulomb calculation

 do k = 1 , MM % N_of_atoms - 1
    do l = k , MM % N_of_atoms
       if ( atom(k) % nresid /= atom(l) % nresid ) then

          ithr    = OMP_get_thread_num() + 1

          nresidk = atom(k) % nresid 
          nresidl = atom(l) % nresid             
          rij(:)  = molecule(nresidk) % cm(:) - molecule(nresidl) % cm(:)
          rij(:)  = rij(:) - MM % box * ANINT( rij(:) * MM % ibox(:) )
          chrgk   = atom(k) % charge
          chrgl   = atom(l) % charge
          rkl(:)  = atom(k) % xyz(:) - atom(l) % xyz(:)
          rkl(:)  = rkl(:) - MM % box(:) * ANINT( rkl(:) * MM % ibox(:) )
          rklq    = sum( rkl(:) * rkl(:) )

          if ( rklq < rcutsq ) then

             select case(forcefield)

                case( 1 )
                ! Born-Mayer ; short range ...

                case( 2 )
                ! Lennard Jones ; short range ...
                atk   = (k + molecule(nresidk) % N_of_atoms) - nresidk * molecule(nresidk) % N_of_atoms
                atl   = (l + molecule(nresidl) % N_of_atoms) - nresidl * molecule(nresidl) % N_of_atoms

                rklsq = SQRT(rklq)
                select case ( MM % CombinationRule )
                     case (2) 
                     ! AMBER FF :: GMX COMB-RULE 2
                     sr2   = ( ( atom(k) % sig + atom(l) % sig ) * ( atom(k) % sig + atom(l) % sig ) ) / rklq
                     case (3)
                     ! OPLS  FF :: GMX COMB-RULE 3
                     sr2   = ( ( atom(k) % sig * atom(l) % sig ) * ( atom(k) % sig * atom(l) % sig ) ) / rklq
                end select
                sr6   = sr2 * sr2 * sr2
                sr12  = sr6 * sr6
                fs    = 24.d0 * ( atom(k) % eps * atom(l) % eps * 1.d-20 ) * ( 2.d0 * sr12 - sr6 )
                fs    = fs / rklq - fscut(atk,atl) / rklsq
                vsr   = 4.d0 * ( atom(k) % eps * atom(l) % eps * 1.d-20 ) * ( sr12 - sr6 )
                vsr   = vsr - vscut(atk,atl) + fscut(atk,atl) * ( rklsq - rcut )
                pot   = pot + vsr
                evdw  = evdw + vsr

                stressr11 = stressr11 + rij(1) * fs * rkl(1)
                stressr22 = stressr22 + rij(2) * fs * rkl(2)
                stressr33 = stressr33 + rij(3) * fs * rkl(3)
                stressr12 = stressr12 + rij(1) * fs * rkl(2)
                stressr13 = stressr13 + rij(1) * fs * rkl(3)
                stressr23 = stressr23 + rij(2) * fs * rkl(3)

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
                vreal = coulomb*1.d-20 * chrgk * chrgl * ERFC(KRIJ)/rklsq
                vreal = vreal - vrecut * chrgk * chrgl + frecut * chrgk * chrgl * ( rklsq-rcut ) * 1.d-20

                pot   = pot + vreal
                ecoul = ecoul + vreal
         
                stresre11 = stresre11 + rij(1) * freal * rkl(1) * 1.d-20
                stresre22 = stresre22 + rij(2) * freal * rkl(2) * 1.d-20
                stresre33 = stresre33 + rij(3) * freal * rkl(3) * 1.d-20
                stresre12 = stresre12 + rij(1) * freal * rkl(2) * 1.d-20
                stresre13 = stresre13 + rij(1) * freal * rkl(3) * 1.d-20
                stresre23 = stresre23 + rij(2) * freal * rkl(3) * 1.d-20

                tmp_fch(k,1:3,ithr) = tmp_fch(k,1:3,ithr) + freal * rkl(1:3)
                tmp_fch(l,1:3,ithr) = tmp_fch(l,1:3,ithr) - freal * rkl(1:3)

             end select
          endif
       end if
    end do
 end do

!$OMP end parallel do
! ################################################################################3

 pot = pot - vself

 stressr(1,1) = stressr11; stressr(2,2) = stressr22
 stressr(3,3) = stressr33; stressr(1,2) = stressr12
 stressr(1,3) = stressr13; stressr(2,3) = stressr23
 stresre(1,1) = stresre11; stresre(2,2) = stresre22
 stresre(3,3) = stresre33; stresre(1,2) = stresre12
 stresre(1,3) = stresre13; stresre(2,3) = stresre23
 stressr(2,1) = stressr(1,2); stressr(3,1) = stressr(1,3)
 stressr(3,2) = stressr(2,3); stresre(2,1) = stresre(1,2)
 stresre(3,1) = stresre(1,3); stresre(3,2) = stresre(2,3)

 do i = 1, MM % N_of_atoms
    do k = 1 , numthr
        atom(i) % fsr(1:3) = atom(i) % fsr(1:3) + tmp_fsr(i,1:3,k)
        atom(i) % fch(1:3) = atom(i) % fch(1:3) + tmp_fch(i,1:3,k)
    end do
    atom(i) % ftotal(1:3) = ( atom(i) % fsr(1:3) + atom(i) % fch(1:3) ) * 1.d-10
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

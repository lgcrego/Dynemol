module FF_cutoff

    use constants_m
    use type_m     , only : warning
    use MD_read_m  , only : MM , atom , species , FF , FF_SP_mtx 
    use gmx2mdflex , only : SpecialPairs
    use Build_DWFF , only : HOH => HOH_diss_parms
    use for_force  , only : rcut, vrecut, frecut, rcutsq, vscut, fscut, KAPPA, vself

    public :: FF_cutoff_sphere 

contains
!
!
!
!===========================
 subroutine FF_cutoff_sphere
!===========================
    implicit none
    
    ! local variables
    integer :: i, j, atmax
    real*8  :: expar, ERFC, KRIJ, total_q2
   
    rcutsq = rcut**2
   
    atmax = sum( species(:) % N_of_atoms )                 
   
    If(.not. allocated(fscut)) allocate ( fscut(atmax,atmax) , source = D_zero )
    If(.not. allocated(vscut)) allocate ( vscut(atmax,atmax) , source = D_zero )
   
    do i = 1, atmax
    do j = 1, atmax
   
            select case ( FF_SP_mtx(i,j) )
   
                   case(0) ! <== not a SpecialPair
                           if( FF(i)%LJ .AND. FF(j)%LJ ) then
                               call Lennard_Jones( i , j )
                           elseif( FF(i)%Buck .AND. FF(j)%Buck ) then
                               call Buckingham( i , j )
                           endif
                   case(1) 
                           call Lennard_Jones( i , j )
                   case(2)
                           call Buckingham( i , j )
                   case(3)
                           call DWFF( i , j )
                   case default
                           CALL warning("unknown non-bonding special pair code in FF_SP_mtx")
                           STOP
            end select
   
    end do
    end do 
   
    KRIJ   = KAPPA * rcut
    vrecut = coulomb * ERFC(KRIJ) / rcut
    expar  = exp(-KRIJ**2)
    frecut = coulomb * ( ERFC(KRIJ) + TWO*irsqPI*KAPPA*rcut*expar ) / rcutsq
   
    ! vself part of the Coulomb calculation
    total_q2 = sum( atom(:)%charge**2 )
    vself = (HALF*vrecut + irsqPI*KAPPA*coulomb) * total_q2
    vself = vself*factor3

end subroutine FF_cutoff_sphere
!
!
!
!
!=================================
 subroutine Lennard_Jones( k , l )
!=================================
    implicit none
    integer , intent(in)  :: k 
    integer , intent(in)  :: l 
    
    ! local variables ...
    integer :: n 
    real*8  :: sr2 , sr6 , sr12 , eps
    logical :: flag1 , flag2 
    
    ! Lennard Jones ...
    
    if( FF_SP_mtx(k,l) == 0 ) then 
          select case ( MM % CombinationRule )
          
              case (2) 
                  ! AMBER FF :: GMX COMB-RULE 2
                  sr2 = (FF(k)%sig + FF(l)%sig)**2 / rcutsq
          
              case (3)
                  ! OPLS  FF :: GMX COMB-RULE 3
                  sr2 = (FF(k)%sig * FF(l)%sig )**2 / rcutsq
          
          end select
          eps = FF(k)%eps * FF(l)%eps
    
    else
          
          n_loop: do  n = 1, size(SpecialPairs)
          
             flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( FF(k) % MMSymbol ) ) .AND. &
                     ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( FF(l) % MMSymbol ) )
             flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( FF(k) % MMSymbol ) ) .AND. &
                     ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( FF(l) % MMSymbol ) )
          
             if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
                sr2 = SpecialPairs(n)%Parms(1)**2 / rcutsq
                eps = SpecialPairs(n)%Parms(2) 
                exit n_loop
             end if
          
          end do n_loop
    endif
    
    sr6  = sr2 * sr2 * sr2
    sr12 = sr6 * sr6
    
    !Energy at spherical surface of radius rcut
    vscut(k,l) = 4.d0 * eps * (sr12 - sr6)
    
    ! Force
    fscut(k,l) = 24.d0 * eps * (2.d0*sr12 - sr6)
    fscut(k,l) = fscut(k,l) / rcut       

end subroutine Lennard_Jones
!
!
!
!==============================
 subroutine Buckingham( k , l )
!==============================
    implicit none
    integer , intent(in)  :: k 
    integer , intent(in)  :: l 
    
    ! local variables ...
    integer :: n 
    real*8  :: sr2 , sr6 
    real*8  :: A , B , C 
    logical :: flag1 , flag2 
    
    ! Bukingham Potential and Forces ...
    
    if( FF_SP_mtx(k,l) == 0 ) then 
    
          ! Combination Rules
          A = FF(k)% BuckA * FF(l)% BuckA
          B = FF(k)% BuckB + FF(l)% BuckB
          C = FF(k)% BuckC * FF(l)% BuckC
          
    else
          
          n_loop: do  n = 1, size(SpecialPairs)
          
             flag1 = ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( FF(k) % MMSymbol ) ) .AND. &
                     ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( FF(l) % MMSymbol ) )
             flag2 = ( adjustl( SpecialPairs(n) % MMSymbols(2) ) == adjustl( FF(k) % MMSymbol ) ) .AND. &
                     ( adjustl( SpecialPairs(n) % MMSymbols(1) ) == adjustl( FF(l) % MMSymbol ) )
          
             if ( flag1 .OR. flag2 ) then      ! <== apply SpecialPair parms ... 
                A = SpecialPairs(n)% Parms(1) 
                B = SpecialPairs(n)% Parms(2)
                C = SpecialPairs(n)% Parms(3)
                exit n_loop
             end if
          
          end do n_loop
    end if
    
    sr2 = 1.d0 / rcutsq
    sr6 = sr2 * sr2 * sr2
    
    !Energy at spherical surface of radius rcut
    vscut(k,l) = A*exp(-B*rcut) - C*sr6 
    
    !Force
    fscut(k,l) = A*B*exp(-B*rcut) - SIX*C*sr6/rcut

end subroutine Buckingham
!
!
!
!
!========================
 subroutine DWFF( i , j )
!========================
    implicit none
    integer , intent(in)  :: i
    integer , intent(in)  :: j
    
    ! local parameters ...
    real*8, parameter :: a1 = 1.1283791671d0  ! <== 2/sqrt(PI)
    
    ! local variables ...
    real*8 :: k , irkl , ir2 , ir6 , ir7 , rkl
    real*8 :: zeta , erfc_zeta , arg , exp_arg2
    real*8 :: a2, a3, U0 , Ecoul , Fcoul , f_sr , E_sr
    real*8 :: A, B, C
    character(len=2) :: type1, type2
   
    !--------------------------------------------
    ! identify (i,j) pair
    type1 = atom(i)% MMSymbol
    type2 = atom(j)% MMSymbol
    select case (trim(type1)//'-'//trim(type2))
    case ('HX-HX')
        k = 3
    case ('OX-OX')
        k = 2
    case ('HX-OX' , 'OX-HX')
        k = 1
    end select
    !--------------------------------------------
 
    rkl  = rcut
    irkl = D_one / rcut
    ir2  = D_one / rcutsq
    
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
        ir7 = ir6 * irkl
        
        ! SR Energy
        E_sr = A*erfc_zeta - C*ir6
        ! SR Force
        f_sr = A*( erfc_zeta + a1*exp(-zeta**2) )*irkl - SIX*C*ir7
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
    Fcoul = coulomb * (U0*ir2 - a3*exp_arg2*irkl)
    
    !----------------------------------------------------
    ! total: intent(out)
    !----------------------------------------------------
    vscut(i,j) = E_sr + Ecoul
    fscut(i,j) = f_sr + Fcoul

end subroutine DWFF
!
!
!
end module FF_cutoff

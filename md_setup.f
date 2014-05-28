module setup_m

    use constants_m
    use atomicmass
    use parameters_m    , only : PBC
    use MD_read_m       , only : MM , atom , molecule , species , FF
    use gmx2mdflex      , only : SpecialPairs

    public :: setup , Molecular_CM , move_to_box_CM , offset

contains
!
!
!
!================
 subroutine SETUP
!================
 use for_force , only: rcut, vrecut, frecut, rcutsq, vscut, fscut, KAPPA, forcefield
 implicit none
 
! local variables
 integer :: i, j, k, atmax
 real*8  :: sr2, sr6, sr12, expar, ERFC, KRIJ, eps, sig
 logical :: flag1 , flag2 

 rcutsq  = rcut * rcut

 fscut = 0.d0
 vscut = 0.d0
 
!################################################
 atmax = sum( species(:) % N_of_atoms )                 
 allocate ( fscut(atmax,atmax), vscut(atmax,atmax) )

 if (forcefield == 1) then
  ! Born-Mayer
 else
  ! Lennard-Jones
 do i = 1, atmax
    do j = 1, atmax
       select case ( MM % CombinationRule )
            case (2)
            ! AMBER FF :: GMX COMB-RULE 2 
            sig = FF(i) % sig + FF(j) % sig
            case (3)
            ! OPLS  FF :: GMX COMB-RULE 3  
            sig = FF(i) % sig * FF(j) % sig
       end select
       eps  = FF(i) % eps * FF(j) % eps

       ! Nbond_parms directive on ...
       read_loop: do  k = 1, size(SpecialPairs)
            flag1 = ( adjustl( SpecialPairs(k) % MMSymbols(1) ) == adjustl( FF(i) % MMSymbol ) ) .AND. &
                    ( adjustl( SpecialPairs(k) % MMSymbols(2) ) == adjustl( FF(j) % MMSymbol ) )
            flag2 = ( adjustl( SpecialPairs(k) % MMSymbols(2) ) == adjustl( FF(i) % MMSymbol ) ) .AND. &
                    ( adjustl( SpecialPairs(k) % MMSymbols(1) ) == adjustl( FF(j) % MMSymbol ) )
            if ( flag1 .OR. flag2 ) then
                sig = SpecialPairs(k) % Parms(1) * SpecialPairs(k) % Parms(1)
                eps = SpecialPairs(k) % Parms(2) * SpecialPairs(k) % Parms(2) 
                exit read_loop
            end if
       cycle  read_loop
       end do read_loop

       sr2  = ( sig * sig ) / rcutsq
       sr6  = sr2 * sr2 * sr2
       sr12 = sr6 * sr6
       vscut(i,j) = 4.d0  * ( eps * factor3 ) * (sr12 - sr6)
       fscut(i,j) = 24.d0 * ( eps * factor3 ) * (2.d0 * sr12 - sr6)
       fscut(i,j) = fscut(i,j) / rcut       
     end do
   end do 
 endif

!###########################################################
 sr2    = 1.d0/rcutsq
 KRIJ   = KAPPA * rcut
 vrecut = coulomb * factor3 * ERFC(KRIJ) / rcut
 expar  = exp( - (KRIJ * KRIJ) )
 frecut = coulomb * sr2 * ( ERFC(KRIJ) + 2. * rsqpi * KAPPA * rcut * expar )

end subroutine SETUP
!
!
!
!======================
 subroutine offset( a )
!======================
implicit none
integer , allocatable , intent(out) :: a(:)

! local variables ...
integer :: i

allocate( a(size(species)) )

a = [ (sum( species(1:i-1)%N_of_atoms ) , i=1,size(species)) ]

end subroutine offset
!
!
!
!========================
 subroutine Molecular_CM
!========================
 implicit none

! local variables ...
 integer              :: i, xyz, k, l 
 real*8               :: massa 
 real*8, dimension(3) :: t0, t, t1, dr
 logical              :: SmallMolecule

! calculates the center of mass of molecule i ...

 l = 1
 do i = 1 , MM % N_of_molecules

      t0(1:3) = atom(l) % xyz(1:3) 
      t(1:3)  = atom(l) % xyz(1:3) * Atomic_mass( atom(l) % AtNo ) ! massa do primeiro átomo da molécula i

      SmallMolecule = NINT(float(molecule(i)%N_of_atoms) / float(MM%N_of_atoms)) == 0

      do k = 1 , molecule(i) % N_of_atoms - 1
           t1(1:3) = atom(l+k) % xyz(1:3)
           dr(1:3) = t1(1:3) - t0(1:3)
           do xyz = 1 , 3
                  if ( abs( dr(xyz) ) > MM % box(xyz) * HALF .AND. SmallMolecule ) then
                     ! fix devided molecule i ...
                     t1(xyz) = t1(xyz) - sign( MM % box(xyz) , dr(xyz) ) 
                  endif
           end do
           massa = Atomic_mass( atom(l+k) % AtNo ) ! massa do átomo k da molécula i
           t(1:3) = t(1:3) + massa * t1(1:3)
      end do  

      molecule(i) % cm(1:3) = imol * t(1:3) / molecule(i) % mass
      l = l + molecule(i) % N_of_atoms

 end do

end subroutine Molecular_CM
!
!
!
!==========================
 subroutine move_to_box_CM
!==========================
implicit none

! local varibales ...
integer :: i, j, l
real*8  :: massa, masstot
real*8, dimension(3) :: t, p, rcm, vcm

! determines atomic Center of Mass and its velocity for the box of atoms ... 
p(:) = 0.d0
t(:) = 0.d0
masstot = 0.d0

l = 1
do i = 1 , MM % N_of_molecules 
    do j = l , l + molecule(i) % N_of_atoms - 1
       massa = Atomic_mass( atom(j) % AtNo )
       p(:) = p(:) + massa * atom(j) % xyz(:)
       t(:) = t(:) + massa * atom(j) % vel(:)
       masstot = masstot + massa
    end do
    l = l + molecule(i) % N_of_atoms
end do
rcm(:) = p(:) / masstot  ! <== center of mass of the box of atoms
vcm(:) = t(:) / masstot  ! <== velocity of the center of mass

! transform atomic coordinates and velocities to CM frame ...
l = 1
do i = 1 , MM % N_of_molecules
    do j = l , l + molecule(i) % N_of_atoms -1
       atom(j) % xyz(:) = atom(j) % xyz(:) - rcm(:)  ! <== atomic coordinates with the origin at CM
       atom(j) % vel(:) = atom(j) % vel(:) - vcm(:)  ! <== atomic velocities measured with respect with vcm
    end do
    l = l + molecule(i) % N_of_atoms
end do

end subroutine move_to_box_CM
!
!
!
!===================
 function ERFC ( X )
!===================
 implicit none
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
end module setup_m

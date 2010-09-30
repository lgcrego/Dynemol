module Multipole_Core

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Semi_Empirical_Parms
    use Structure_Builder 

    type(R3_vector) , allocatable , public , protected :: DP_matrix_AO(:,:)

    public :: Dipole_Matrix 
    public :: rotationmultipoles
    public :: multipole_messages
    public :: multipoles1c
    public :: multipoles2c
    public :: Util_Multipoles

    private

contains
!
!
!
!=====================================================================
 subroutine Dipole_Matrix( system , basis , L_vec , R_vec , Total_DP )
!=====================================================================
implicit none
type(structure)             , intent(inout) :: system
type(STO_basis)             , intent(in)    :: basis(:)
complex*16                  , intent(in)    :: L_vec(:,:) , R_vec(:,:)
real*8          , optional  , intent(out)   :: Total_DP(3) 

! local variables ...
real*8  :: Sparsity(3)
integer :: NonZero(3) , M_size , i

If( verbose ) Print 153
!----------------------------------------------------------
!       initialize DIPOLE MATRIX M(i,j)[x,y,z]

 CALL Util_Multipoles

! size of M matrix ...
M_size = sum(atom(system%AtNo)%DOS)

If( allocated(DP_matrix_AO) ) deallocate( DP_matrix_AO )
allocate(DP_matrix_AO(M_size,M_size))

CALL Build_DIPOLE_Matrix(system,basis)

forall(i=1:3) NonZero(i) = count(DP_matrix_AO(:,:)%dp(i) /= 0.d0)

Sparsity(:) = dfloat(NonZero(:))/dfloat((M_size**2))

If( verbose ) Print 73, Sparsity  

if ( DP_Moment ) CALL Dipole_Moment( system , basis , L_vec , R_vec , Total_DP )

!----------------------------------------------------------
If( verbose ) then
    Print*, '>> Dipole Moment done <<'
    Print 155
end If

 include 'formats.h'

end subroutine Dipole_Matrix
!
!
!
!=====================================================================
 subroutine Dipole_Moment( system , basis , L_vec , R_vec , DP_total )
!=====================================================================
implicit none
type(structure)             , intent(inout) :: system
type(STO_basis)             , intent(in)    :: basis(:)
complex*16                  , intent(in)    :: L_vec(:,:) , R_vec(:,:)
real*8          , optional  , intent(out)   :: DP_total(3) 

! local variables ...
integer                       :: i, j, states, xyz, n_basis, Fermi_state
real*8                        :: Nuclear_DP(3), Electronic_DP(3), Total_DP(3)
real*8          , allocatable :: R_vector(:,:)
complex*16      , allocatable :: a(:,:), b(:,:)
logical         , allocatable :: mask(:)
type(R3_vector) , allocatable :: origin_Dependent(:), origin_Independent(:)

! local parameters ...
real*8          , parameter   :: Debye_unit = 4.803204d0

! define system for DP_Moment calculation ...
allocate( mask(size(basis)) , source = .true. )
mask = merge( mask , basis%FMO , count(basis%FMO) == I_zero )

CALL Center_of_Charge( system , R_vector )

! Nuclear dipole ; if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0)
Nuclear_DP = D_zero 

! Electronic dipole 
 n_basis      =  size(basis)
 Fermi_state  =  sum( system%Nvalen ) / two
 
 allocate( a(n_basis,n_basis)              , source = C_zero )
 allocate( b(n_basis,n_basis)              , source = C_zero )
 allocate( origin_Dependent(Fermi_state)   )
 allocate( origin_Independent(Fermi_state) )

 do xyz = 1 , 3

!   origin dependent DP = sum{C_dagger * vec{R} * S_ij * C}

    forall(states=1:Fermi_state)

        forall( i=1:n_basis ) a(states,i) = L_vec(states,i) * R_vector(basis(i)%atom,xyz)

        origin_Dependent(states)%DP(xyz) = 2.d0 * real( sum( a(states,:)*R_vec(:,states) , mask ) )

    end forall    

!   origin independent DP = sum{C_dagger * vec{DP_matrix_AO(i,j)} * C}

    b = DP_matrix_AO%DP(xyz)

    CALL gemm( L_vec , b , a , 'N' , 'N' , C_one , C_zero )    

    forall( states=1:Fermi_state ) origin_Independent(states)%DP(xyz) = 2.d0 * real( sum( a(states,:)*L_vec(states,:) , mask ) )

 end do

 forall(xyz=1:3) Electronic_DP(xyz) = sum( origin_Dependent%DP(xyz) + origin_Independent%DP(xyz) )

 Total_DP = ( Nuclear_DP - Electronic_DP ) * Debye_unit

 If( present(DP_total) ) DP_total = Total_DP

 If( verbose ) Print 154, Total_DP, dsqrt(sum(Total_DP*Total_DP))

 deallocate(R_vector,a,b)
 deallocate(origin_Dependent)
 deallocate(origin_Independent)

 include 'formats.h'

end subroutine Dipole_Moment
!
!
!
!===========================================
 subroutine Center_of_Charge( a , R_vector )
!===========================================
implicit none
type(structure)                 , intent(inout) :: a
real*8          , allocatable   , intent(out)   :: R_vector(:,:)

! local variables ...
integer               :: i , j
real*8                :: total_valence
real*8  , allocatable :: Qi_Ri(:,:) 
logical , allocatable :: mask(:)

! define system for DP_Moment calculation ...
allocate( mask(a%atoms) , source = .true. )
mask = merge( mask , a%FMO , count(a%FMO) == I_zero )

! sum_i = (q_i * vec{r}_i) / sum_i q_i ...

allocate( Qi_Ri(a%atoms,3) , source = D_zero )

forall( j=1:3 , i=1:a%atoms , mask(i) ) Qi_Ri(i,j) = a%Nvalen(i) * a%coord(i,j)

total_valence = sum( a%Nvalen , mask )

forall(j=1:3) a%Center_of_Charge(j) = sum( Qi_Ri(:,j) , mask ) / total_valence

! atomic positions measured from the Center of Charge
allocate( R_vector(a%atoms,3) , source = D_zero )
forall( j=1:3 , i=1:a%atoms , mask(i) ) R_vector(i,j) = a%coord(i,j) - a%Center_of_Charge(j)

deallocate( Qi_Ri , mask )

end subroutine Center_of_Charge
!
!
!
!============================================
subroutine Build_DIPOLE_Matrix(system, basis)
!============================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)

! local variables
real*8  :: expa, expb, xab , yab , zab , Rab 
integer :: AtNo_a , AtNo_b
integer :: a , b , ia , ib , ja , jb 
integer :: na , la , ma 
integer :: nb , lb , mb
integer :: lmult , i , j

real*8  , parameter :: tol = 1.d-10 
integer , parameter :: mxl = 5 , mxmult = 3 , mxlsup = max(mxl,mxmult)
real*8  , parameter :: cutoff_Angs = 10.d0

real*8 , dimension((mxmult+1)**2,-mxl:mxl,-mxl:mxl)        :: qlm
real*8 , dimension(-mxlsup:mxlsup,-mxlsup:mxlsup,0:mxlsup) :: rl , rl2

lmult = 1 ! <== DIPOLE MOMENT

forall(i=1:3) DP_matrix_AO(:,:)%dp(i) = 0.d0

do ib = 1  , system%atoms
do ia = 1 , system%atoms  

! calculate rotation matrix for the highest l

    call RotationMultipoles(system,ia,ib,xab,yab,zab,Rab,lmult,rl,rl2)

    If(Rab*a_Bohr > cutoff_Angs) goto 10

    do jb = 1 , atom(system%AtNo(ib))%DOS  ;  b = system%BasisPointer(ib) + jb
    do ja = 1 , atom(system%AtNo(ia))%DOS  ;  a = system%BasisPointer(ia) + ja

        na = basis(a)%n ;   la = basis(a)%l ;   ma = basis(a)%m
        nb = basis(b)%n ;   lb = basis(b)%l ;   mb = basis(b)%m

        CALL Multipole_Messages(na,nb,la,lb)

!---------------------------------------------------------------------------------------------------- 
!       sum over zeta coefficients
        do i = 1 , basis(a)%Nzeta
        do j = 1 , basis(b)%Nzeta
   
            expa = basis(a)%zeta(i)
            expb = basis(b)%zeta(j)

            if( ia==ib ) then

!               CALLS THE SUBROUTINE FOR THE MULTIPOLES OF ONE-CENTER DISTRIBUTIONS

                qlm = 0.d0   ! check this !!!!

                call multipoles1c(na, la, expa, nb, lb, expb, lmult, qlm)

            else 

!               CALLS THE SUBROUTINE FOR THE MULTIPOLES OF TWO-CENTER DISTRIBUTIONS

                qlm = 0.d0   ! check this !!!!!

                call multipoles2c(na, la, expa, nb, lb, expb, xab, yab, zab, Rab, lmult, rl, rl2, qlm)

            end if

!           p_x(a,b) 
            DP_matrix_AO(a,b)%dp(1) = DP_matrix_AO(a,b)%dp(1) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(4,ma,mb)
!           p_y(a,b)
            DP_matrix_AO(a,b)%dp(2) = DP_matrix_AO(a,b)%dp(2) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(2,ma,mb)
!           p_z(a,b)
            DP_matrix_AO(a,b)%dp(3) = DP_matrix_AO(a,b)%dp(3) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(3,ma,mb)

        end do
        end do
!---------------------------------------------------------------------------------------------------- 

    enddo
    enddo
10 end do
end do

end subroutine Build_DIPOLE_Matrix
!
!
!
!=======================================================================
subroutine RotationMultipoles(system,ia,ib,xab,yab,zab,Rab,lmult,rl,rl2)
!=======================================================================
implicit none
type(structure)            , intent(in)  :: system
integer                    , intent(in)  :: ia , ib
real*8                     , intent(out) :: xab, yab, zab, Rab
integer                    , intent(in)  :: lmult
real*8  , dimension(:,:,:) , intent(out) :: rl , rl2

! local variables ...
real*8  :: xa, ya, za, xb, yb, zb , xy , sinal , cosal , sinbet , cosbet , singa , cosga
integer :: AtNo_a , AtNo_b
integer :: la_max , lb_max , lmax

real*8  , parameter :: tol = 1.d-10
integer , parameter :: mxl = 5 , mxmult = 3 , mxlsup = max(mxl,mxmult)

! INPUT DATA
! na, la, expa, xa, ya, za:  N, L, EXPONENT and cartesian coordinates of function on A
! nb, lb, expb, xb, yb, zb:  N', L', EXPONENT' and cartesian coordinates of function on B
!
! coordinates must be in a.u. 

AtNo_a = system%AtNo(ia)
AtNo_b = system%AtNo(ib)

xa = system%coord(ia,1) / a_Bohr
ya = system%coord(ia,2) / a_Bohr
za = system%coord(ia,3) / a_Bohr
xb = system%coord(ib,1) / a_Bohr
yb = system%coord(ib,2) / a_Bohr
zb = system%coord(ib,3) / a_Bohr

! Loads the matrices for the rotation of (normalized) spherical harmonics
! from the lined-up system to the molecular frame
! Rotation matrices are stored in rl(m,m',l) with  (-l <= m,m' <= l)

la_max = atom(AtNo_a)%AngMax
lb_max = atom(AtNo_b)%AngMax

! -----------------------------------------------------------------
! THIS BLOCK SHOULD BE RUN ONLY ONCE FOR A GIVEN PAIR OF CENTERS
! SINCE THE ROTATION MATRICES ARE THE SAME FOR ALL THE PAIRS OF 
! STO FUNCTIONS BELONGING TO A PAIR OF CENTERS AB. subroutine 
! rotar should then be called with lmax equal to the maximum value
! of the three following: 
!     highest l quantum number in the basis set of center A 
!     highest l quantum number in the basis set of center B
!     highest order of the multipoles to be computed (lmult)
!
! Euler angles for rotating from the molecular frame to a lined-up system  
! (Z axes coincident, X and Y axes parallel)

xab = xb - xa
yab = yb - ya
zab = zb - za
Rab = dsqrt(xab*xab + yab*yab + zab*zab)
xy  = dsqrt(xab*xab + yab*yab) 
if (xy .gt. tol) then
    sinal = yab / xy
    cosal = xab / xy
else
    sinal = 0.d0
    cosal = 1.d0
endif
sinbet = xy / Rab
cosbet = zab / Rab
singa = 0.d0
cosga = 1.d0
lmax = max(lmult, la_max, lb_max)

call rotar(lmax, mxlsup, cosal, sinal, cosbet, sinbet, cosga, singa, rl2, rl)

! END OF BLOCK
! -----------------------------------------------------------------
end subroutine RotationMultipoles
!
!
!
subroutine Util_Multipoles

implicit real*8 (a-h,o-z)

integer , parameter :: mxl = 5 , mxmult = 3 
integer , parameter :: mxltot = mxl + mxmult , mxlsup = max(mxl,mxmult) 
integer , parameter :: mxemes = mxltot 
integer , parameter :: mxreal = 1000 , mxfact = 150 , mxind = 200 , mxroot=50

common/indcomn/  indk12((mxltot+1)**2,(mxltot+1)**2)
common/const/    re(0:mxreal), reali(0:mxreal), fact(0:mxfact), facti(0:mxfact), ang((mxl+1)*(mxl+2)/2), &
                 root(0:mxroot), rooti(mxroot), ind(0:mxind)
common/emescom/  ssv(-mxemes:mxemes,-mxemes:mxemes),sdv(-mxemes:mxemes,-mxemes:mxemes), &
                 msv(-mxemes:mxemes,-mxemes:mxemes),mdv(-mxemes:mxemes,-mxemes:mxemes)

! ---------------------------------------------------------------
! THIS BLOCK SHOULD BE EXECUTED ONLY ONCE IN THE MAIN PROGRAM
! Computes some auliary vectors: factorials, integers to real, ...
! to save computational time

fact(0) = 1.d0
facti(0) = 1.d0
do i = 1, mxfact
    fact(i) = fact(i-1) * i       !  i!
    facti(i) = 1.d0 / fact(i)     !  1.d0 / i!
enddo
re(0) = 0.d0
reali(0) = 0.d0 ! Just a trick to avoid problems in some recursion
do i = 1, mxreal
    re(i) = re(i-1) + 1.d0        ! dfloat(i)
    reali(i) = 1.d0 / re(i)       ! 1.d0 / dfloat(i)
enddo
ind(0) = 0
do i = 1, mxind
    ind(i) = ind(i-1) + i         !  i*(i+1)/2
enddo
root(0) = 0.d0
do i = 1, mxroot
    root(i) = dsqrt(re(i))        !  dsqrt(i)
    rooti(i) = 1.d0 / root(i)     !  1.d0 / dsqrt(i)
enddo
! ang(l*(l+1)/2+m+1) = dsqrt( (2*l+1) * fact(l-m) / (2 * pi * (1 + delta(m,0)) * fact(l+m)) )
ang(1) = 0.282094791773878d0
raiz2 = dsqrt(2.d0)
lm = 1
do l = 1, mxl
    lm = lm + 1
    ang(lm) = ang(1) * dsqrt(re(2*l+1))
    aux = ang(lm) * raiz2
    do m = 1, l
        lm = lm + 1
        aux = aux / dsqrt(re(l-m+1)*re(l+m))
        ang(lm) = aux
    enddo
enddo
! Tabulates the coefficients for decomposing the products of two 
! associated Legendre functions into associated Legendre functions
    call acof(mxlsup)
    call bcof(mxlsup)
! Tabulates some auxiliary indices for locating the previous coefficients
do l2 = 0 , mxltot
    do l1 = 0 , mxltot
        do m2 = -l2, l2
            do m1 = -l1, l1
                l1l1 = ind(l1)
                l2l2 = ind(l2)
                m1a = iabs(m1)
                m2a = iabs(m2)
                if ( l1.eq.l2 ) then
                    k1 = l1l1 + max(m1a,m2a)
                    k12 = ind(k1) + l1l1 + min(m1a,m2a)
                elseif (l1.gt.l2) then
                    k1 = l1l1 + m1a 
                    k12 = ind(k1) + l2l2 + m2a
                else 
                    k1 = l2l2 + m2a 
                    k12 = ind(k1) + l1l1 + m1a
                endif
                indk12(l1*(l1+1)+m1+1,l2*(l2+1)+m2+1)=k12 
            end do
        end do
    end do
end do
! Tabulates the coefficients for the decomposition of products
! of two functions depending on phi (sin (m*phi), cos (m*phi))
! into functions of the same type
do m2 = -mxltot, mxltot
    do m1 = -mxltot, mxltot
        call emes ( m1, m2, ms, md, ss, sd )
        msv(m1,m2) = ms
        mdv(m1,m2) = md
        ssv(m1,m2) = ss
        sdv(m1,m2) = sd 
    enddo
enddo
!     END OF BLOCK
!     -------------------------------------------------------------

end subroutine Util_Multipoles
!
!
!
subroutine Multipole_Messages(na,nb,la,lb)

integer , intent(in) :: na , nb , la , lb

integer , parameter :: mxl = 5 , mxmult = 3 , mxn = 5+2*mxmult

      if (na .gt. mxn) then
         print*, 'Error. Maximum value of na exceeded'
         print*, 'Current value: ', na
         print*, 'Maximum allowed: ', mxn
         print*, 'Change parameter mxn and recompile'
         print*, 'Stop'
         Stop
      endif

      if (nb .gt. mxn) then
         print*, 'Error. Maximum value of nb exceeded'
         print*, 'Current value: ', nb
         print*, 'Maximum allowed: ', mxn
         print*, 'Change parameter mxn and recompile'
         print*, 'Stop'
         Stop
      endif

      if (la .gt. mxl) then
         print*, 'Error. Maximum value of na exceeded'
         print*, 'Current value: ', la
         print*, 'Maximum allowed: ', mxl
         print*, 'Change parameter mxl and recompile'
         print*, 'Stop'
         Stop
      endif

      if (lb .gt. mxn) then
         print*, 'Error. Maximum value of nb exceeded'
         print*, 'Current value: ', lb
         print*, 'Maximum allowed: ', mxl
         print*, 'Change parameter mxl and recompile'
         print*, 'Stop'
         Stop
      endif

end subroutine Multipole_Messages
!
!
!
!
!
!
!     ***************************************************************************
!     Subrutine for computing multipolar moments of one-center STO distributions
!
!     Multipolar moments are defined as:
!
!        Q(l,m) = ( (2-delta(m,0)) * (l-|m|)! / (l+|m|)! ) * 
!                 Integrate[ z(l,m,r_A) * chi_A(N,L,M,r_A) * chi_B(N',L',M',r_B) ]
!
!        zlm(l,m,r_A) being the regular real solid harmonics:
!
!        zlm(l,m,r_A) = r_A**l * zlm(theta_A,phi_A)
!
!     where the integral extends over the whole space (R^3) and 
!     chi_A(N,L,M,r_A), chi_B(N',L',M',r_B) are STO functions centered at A 
!     with the usual meaning of the quntum numbers N, L, M
!
!     With this definition, the long-range potential generated by the charge 
!     distribution  chi_A(N,L,M,r_A) chi_B(N',L',M',r_B)   results:
!
!        Vlong(r_A) = Sum[ Sum[ Qlm(l,m) * zlm(l,m,r_A) / r_A**(2*l+1)
!              , {m, -l, l}], {l, 0, L+L'}]
!
!     where only values of l with the same parity as L+L' appear, and 
!     two  m  components at most are different from zero for each distribution.
!
!     See notebook: comprueba_multipolos.nb
!
!     Coded by Rafael Lopez (rafael.lopez@uam.es)  October 6th 2008.
!
subroutine multipoles1c(na, la, exa, nb, lb, exb, lmult, qlm)

implicit real*8 (a-h,o-z)

integer , intent(in)    :: na, la, nb, lb, lmult
real*8  , intent(in)    :: exa, exb

parameter (mxl = 5, mxmult = 3, mxn = 5+2*mxmult)
parameter (mxltot = mxl + mxmult, mxlsup = max(mxl,mxmult) )
parameter ( mxemes = mxltot )
parameter ( mxlcof = mxlsup*(mxlsup+3)/2 )
parameter ( mxkcof = mxlcof*(mxlcof+3)/2 )
parameter (mxreal = 1000, mxfact = 150, mxind = 200, mxroot=50)
parameter ( tol = 1.d-10 )

common/abpp/     app(0:mxkcof,0:2*mxl+1), bpp(0:mxkcof,0:2*mxl+1)
common/const/    re(0:mxreal), reali(0:mxreal), fact(0:mxfact), facti(0:mxfact), ang((mxl+1)*(mxl+2)/2), &
                 root(0:mxroot), rooti(mxroot), ind(0:mxind)
common/emescom/  ssv(-mxemes:mxemes,-mxemes:mxemes),sdv(-mxemes:mxemes,-mxemes:mxemes), &
                 msv(-mxemes:mxemes,-mxemes:mxemes),mdv(-mxemes:mxemes,-mxemes:mxemes)
common/indcomn/  indk12((mxltot+1)**2,(mxltot+1)**2)

real*8 :: qlm((mxmult+1)**2, -mxl:mxl, -mxl:mxl)
real*8 :: auxvec(0:mxmult)

      do mb = 0,lb
      do ma = 0, la
      do lm = 1, (lmult+1)*(lmult+2)/2
         qlm(lm,ma,mb) = 0.d0
      enddo
      enddo
      enddo
      aux = 1.d0 / (exa+exb)
      auxvec(0) = pi4 * fact(na+nb) * aux**(na+nb+1)
      do l = 1, lmult
         auxvec(l) = auxvec(l-1) * re(na+nb+l) * aux
      enddo
      rnor = dsqrt((exa+exa)**(na+na+1) * (exb+exb)**(nb+nb+1) / ( fact(na+na) * fact(nb+nb) ) ) 
      lma = la * (la+1)/2 + 1
      lmb = lb * (lb+1)/2 + 1
      ldf = iabs(la-lb)
      do ma = -la, la
         do mb = -lb, lb 
            anor = rnor * ang(lma+iabs(ma)) * ang(lmb+iabs(mb)) 
            ms = msv(ma,mb)
            md = mdv(ma,mb)
            ss = ssv(ma,mb)
            sd = sdv(ma,mb)
            k12 = indk12(lb*(lb+1)+mb+1,la*(la+1)+ma+1)
            msabs = iabs(ms)
            if (msabs .le. lmult .and. dabs(ss) .ge. 1.d-5) then
               do l = la+lb, max(ldf,msabs), -2
                  if (l .le. lmult) then
                     qlm(l*(l+1)+ms+1,ma,mb) = ss * app(k12,l) * anor * auxvec(l) * reali(l+l+1) 
                  endif
               enddo 
            endif
            mdabs = iabs(md)
            if (mdabs .le. lmult .and. dabs(sd) .ge. 1.d-5) then
               do l = la+lb, max(ldf,mdabs), -2
                  if (l .le. lmult) then
                     qlm(l*(l+1)+md+1,ma,mb) = sd * bpp(k12,l) * anor * auxvec(l) * reali(l+l+1) 
                  endif
               enddo 
            endif
         enddo
      enddo
      return
end subroutine multipoles1c

!     ***************************************************************************
!     Subrutine for computing multipolar moments of two-center STO distributions
!
!     Multipolar moments are defined as:
!
!        Q(l,m) = ( (2-delta(m,0)) * (l-|m|)! / (l+|m|)! ) * 
!                 Integrate[ z(l,m,r_A) * chi_A(N,L,M,r_A) * chi_B(N',L',M',r_B) ]
!
!        zlm(l,m,r_A) being the regular real solid harmonics:
!
!        zlm(l,m,r_A) = r_A**l * zlm(theta_A,phi_A)
!
!     where the integral extends over the whole space (R^3) and 
!     chi_A(N,L,M,r_A), chi_B(N',L',M',r_B) are STO functions centered at A and B
!     with the usual meaning of the quntum numbers N, L, M
!
!     With this definition, the long-range potential generated by the charge 
!     distribution  chi_A(N,L,M,r_A) chi_B(N',L',M',r_B)   results:
!
!        Vlong(r_A) = Sum[ Sum[ Qlm(l,m) * zlm(l,m,r_A) / r_A**(2*l+1)
!              , {m, -l, l}], {l, 0, Infinity}]
!
!     Multipolar moments are computed in terms of overlap integrals in a lined-up axis system:
!        < N L M | N' L' M >    M >= 0
!
!     See notebook: comprueba_multipolos.nb
!
!     Coded by Rafael Lopez (rafael.lopez@uam.es)  October 6th 2008.
!
subroutine multipoles2c(na, la, exa, nb, lb, exb, xab, yab, zab, Rab, lmult, rl, rl2, qlm)

implicit real*8 (a-h,o-z)

parameter (mxl = 5, mxmult = 3, mxn = 5+2*mxmult)
parameter (mxltot = mxl + mxmult, mxlsup = max(mxl,mxmult) )
parameter ( mxemes = mxltot )
parameter ( mxlcof = mxlsup*(mxlsup+3)/2 )
parameter ( mxkcof = mxlcof*(mxlcof+3)/2 )
parameter (mxreal = 1000, mxfact = 150, mxind = 200, mxroot=50)
parameter ( tol = 1.d-10 )

common/abpp/     app(0:mxkcof,0:2*mxl+1), bpp(0:mxkcof,0:2*mxl+1)
common/const/    re(0:mxreal), reali(0:mxreal), fact(0:mxfact), facti(0:mxfact), ang((mxl+1)*(mxl+2)/2), &
                 root(0:mxroot), rooti(mxroot), ind(0:mxind)
common/emescom/  ssv(-mxemes:mxemes,-mxemes:mxemes),sdv(-mxemes:mxemes,-mxemes:mxemes), &
                 msv(-mxemes:mxemes,-mxemes:mxemes),mdv(-mxemes:mxemes,-mxemes:mxemes)
common/indcomn/  indk12((mxltot+1)**2,(mxltot+1)**2)

real*8 :: sol(0:mxl, 0:mxmult+mxl, 0:(mxl+mxmult)/2)
real*8 :: qlm((mxmult+1)**2, -mxl:mxl, -mxl:mxl)
real*8 :: auxvec(-mxl:mxl)
real*8 :: rl(-mxlsup:mxlsup,-mxlsup:mxlsup,0:mxlsup)
real*8 :: rl2(-mxlsup:mxlsup,-mxlsup:mxlsup,0:mxlsup)

!     computes the overlap integrals on a lined-up system
      call solap(na, la, exa, nb, lb, exb, Rab, lmult, sol)
!     computes the multipolar integrals on a lined-up system
      do mb = 0,lb
      do ma = 0, la
      do lm = 1, (lmult+1)*(lmult+2)/2
         qlm(lm,ma,mb) = 0.d0
      enddo
      enddo
      enddo
      rnor = dsqrt((exa+exa)**(na+na+1) * (exb+exb)**(nb+nb+1) / ( fact(na+na) * fact(nb+nb) ) ) 
      lma = la * (la+1)/2 + 1
      lmb = lb * (lb+1)/2 + 1
      lm = 0
      do l = 0, lmult
         ldf = iabs(l-la)
         do m = -l, l
            lm = lm + 1
            mabs = iabs(m)
            do ma = -la, la 
               anor = rnor * ang(lma+iabs(ma)) * ang(ind(l)+mabs+1) 
               ms = msv(m,ma)
               md = mdv(m,ma)
               ss = ssv(m,ma)
               sd = sdv(m,ma)
               k12 = indk12(l*(l+1)+m+1,la*(la+1)+ma+1)
               msabs = iabs(ms)
               if (msabs .le. lb .and. dabs(ss) .ge. 1.d-5) then
                  k = 0
                  bux = 0.d0
                  do lp = l+la, max(ldf,msabs), -2
                     bux = bux + ss*app(k12,lp) * sol(msabs,l,k)
                     k = k + 1
                  enddo
                  qlm(lm,ma,ms) = bux * anor * ang(lmb+msabs) 
               endif
               mdabs = iabs(md)
               if (mdabs .le. lb .and. dabs(sd) .ge. 1.d-5) then
                  k = 0
                  bux = 0.d0
                  do lp = l+la, max(ldf,mdabs), -2
                     bux = bux + sd*bpp(k12,lp) * sol(mdabs,l,k)
                     k = k + 1
                  enddo
                  qlm(lm,ma,md) = bux * anor * ang(lmb+mdabs) 
               endif
            enddo
         enddo
      enddo
!     rotates the multipolar integrals to the molecular frame
!     rotates function on B
      do ma = -la, la
         do lm = 1, (lmult+1)**2
            do mb = -lb, lb
               auxvec(mb) = qlm(lm,ma,mb)
            enddo
            do m = -lb, lb
               aux = 0.d0
               do mp = -lb, lb
                  aux = aux + rl(m,mp,lb) * auxvec(mp)
               enddo
               qlm(lm,ma,m) = aux
            enddo
         enddo
      enddo
      do mb = -lb, lb
      do ma = -la, la
         lm = 0
         do l = 0, lmult
            do m = -l, l
               lm = lm + 1
            enddo
         enddo
      enddo
      enddo
      do mb = -lb, lb
         do lm = 1, (lmult+1)**2
            do ma = -la, la
               auxvec(ma) = qlm(lm,ma,mb)
            enddo
            do m = -la, la
               aux = 0.d0
               do mp = -la, la
                  aux = aux + rl(m,mp,la) * auxvec(mp)
               enddo
               qlm(lm,m,mb) = aux
            enddo
         enddo
      enddo
      do mb = -lb, lb
      do ma = -la, la
         lm = 0
         do l = 0, lmult
            do m = -l, l
               lm = lm + 1
            enddo
         enddo
      enddo
      enddo
      do mb = -lb, lb
         do ma = -la, la
            lm0 = 0
            do l = 0, lmult
               lm = lm0
               do m = -l, l
                  lm = lm + 1
                  auxvec(m) = qlm(lm,ma,mb)
               enddo
               lm = lm0
               bux = pi4 * reali(l+l+1)
               do m = -l, l
                  aux = 0.d0
                  do mp = -l, l
                     aux = aux + rl(m,mp,l) * auxvec(mp)
                  enddo
                  lm = lm + 1
                  qlm(lm,ma,mb) = bux * aux * ang(ind(l)+iabs(m)+1)
               enddo
               lm0 = lm
            enddo
         enddo
      enddo
      return
end subroutine multipoles2c

!**********************************************************************
!     Subroutine for computing two-center overlap integrals of STO on a lined-up system
!     that are computed by recursion starting from the basic overlap integrals:
!        < N 0 0 | N' 0 0 >
!     Recurrence relations:
!
!     For increasing  M :
!
!     < N M+1 M+1 | N' M+1 M+1 > = ((2M+1)**2/(4R**2)) * (
!             2 < N+1 M M | N'+1 M M > + 2 R**2 < N+1 M M | N'-1 M M >
!             + 2 R**2 < N-1 M M | N'+1 M M > - < N+3 M M | N'-1 M M >
!             - < N-1 M M | N'+3 M M > - R**4 < N-1 M M | N'-1 M M > )
!
!
!     For increasing  L:
!
!     < N L+1 M | N' L' M > = ((2L+1)/(2 (L-M+1) R) ) * (
!             < N+1 L M | N' L' M >  + R**2 < N-1 L M | N' L' M >
!             - < N-1 L M | N'+2 L' M > )
!          - ( (L+M)/(L-M+1) ) < N L-1 M | N' L' M >
!
!
!     For increasing  L':
!
!     < N L M | N' L'+1 M > = -((2L'+1)/(2 (L'-M+1) R) ) * (
!             < N L M | N'+1 L' M >  + R**2 < N L M | N'-1 L' M >
!             - < N+2 L M | N'-1 L' M > )
!          - ( (L'+M)/(L'-M+1) ) < N L M | N' L'-1 M >
!
!     Starting integrals for the recursion are computed by:
!
!     < N 0 0 | N' 0 0 > = ( 4 pi Exp[-z*R] N! N'! / (z+z')**(N+N'+1) )
!         * Sum[ Sum[
!            (N+N'-j-j')! [(z+z')*R]**(j+j') * 1F1[j'+1,j+j'+2,(z-z')R]
!            / ( (N-j)! (N'-j')! (j+j'+1)!
!         , {j',0,N'}] , {j,0,N}]
!
!     where  1F1[a,b,z] is the corresponding confluent hypergeometric function.
!
!     See notebook: Langhoff.nb
!
subroutine solap(na, la, exa, nb, lb, exb, R, lmult, sol)
      implicit real * 8 (a-h,o-z)
      parameter (pi4 = 12.56637061435917d0)
      parameter (mxl = 5, mxmult = 3, mxn = 5+2*mxmult)
      parameter (mxltot = mxl + mxmult )
      parameter (mxh = mxn+2*mxltot+4*mxl)
      parameter (mxreal = 1000, mxfact = 150, mxind = 200, mxroot=50)
      common /const/ re(0:mxreal),reali(0:mxreal),fact(0:mxfact),facti(0:mxfact),ang((mxl+1)*(mxl+2)/2),&
                     root(0:mxroot), rooti(mxroot), ind(0:mxind)
      dimension h1f1(0:mxh+1,0:2*(mxh+1))
      dimension sol(0:mxl, 0:mxmult+mxl, 0:(mxl+mxmult)/2)
      dimension smat(0:mxh,0:mxh), smataux(0:mxh,0:mxh)
      dimension smatm(0:mxh,0:mxh)
      lmaux = min(lmult, max(la,lb))
      if (exa .ge. exb) then
         nmin = na - la
         nmax = na + la + lmult + lmult + lb + lb + lmaux
         npmin = nb - lb
         npmax = nb + lb + la + la + lmult + lmult + lmaux
         z = exa
         zp = exb
      else
         nmin = nb - lb
         nmax = nb + lb + la + la + lmult + lmult + lmaux
         npmin = na - la
         npmax = na + la + lb + lb + lmult + lmult + lmaux
         z = exb
         zp = exa
      endif
      expzR = dexp(-z*R)
!     Computes the Confluent Hypergeometric functions (h1f1)
      ngmax = nmax+npmax+2   ! Top value of gamma parameter
      namax = npmax+1        ! Top value of alpha parameter
      arg = (z-zp) * R       ! Argument of the hypergeometric
      iz = z
      ia = (ngmax-iz) / 2
      if (ia .lt. 1) ia = 1
      if (ia .ge. namax) ia = namax-1

!     Computes 1F1[ia,ngmax,z]  and  1F1[ia+1,ngmax,z] by their series

      auxa = arg * re(ia) * reali(ngmax)
      f11a = 1.d0
      auxb = arg * re(ia+1) * reali(ngmax)
      f11b = 1.d0
      imax = 1000
      do i = 1, imax
         f11a = f11a + auxa
         f11b = f11b + auxb
         auxa = auxa * re(ia+i) * arg * reali(ngmax+i) * reali(i+1)
         auxb = auxb * re(ia+i+1) * arg * reali(ngmax+i) * reali(i+1)
         if ( auxb .lt. 1.d-17*f11b) go to 100
      enddo
      write(6,1000) ia+1, ngmax, arg, imax, f11b, auxb
 1000 format (1x,'Convergence in hypergeometric function 1F1[', &
              i3, ',', i3, ',', e22.15, '] not reached ', &
              /1x, 'after ', i5, ' terms in the series', /1x, &
              'Series value = ', e22.15, 5x, 'Last summand = ', e22.15)
  100 continue
      h1f1(ia,ngmax) = f11a
      h1f1(ia+1,ngmax) = f11b

!     Computes the column with gamma = ngmax by the recursion:
!        a * 1F1[a+1,g,x] = (x+2a-g) * 1F1[a,g,x] + (g-a) * 1F1[a-1,g,x] 
!     which is stable upwards for  a > (g - x) / 2
!               and downwards for a < (g - x) / 2

!     Upwards recursion:
      aux = arg-re(ngmax)
      do i = ia+1, namax-1
         h1f1(i+1,ngmax) = ((aux+re(i+i)) * h1f1(i,ngmax) + re(ngmax-i) * h1f1(i-1,ngmax)) * reali(i)
      enddo
!     Downwards recursion:
      do i = ia, 2, -1
         h1f1(i-1,ngmax) = ( re(i) * h1f1(i+1,ngmax) - (aux + re(i+i)) * h1f1(i,ngmax) ) * reali(ngmax-i)
      enddo
!     Computes the remaining columns (gamma < ngmax) by the recursion:
!        g * 1F1[a,g,x] = a * 1F1[a+1,g+1,x] + (g-a) * 1F1[a,g+1,x]
!     for the rows: a = 1,  and
!        1F1[a,g,x]= (x/g) * 1F1[a,g+1,x] + 1F1[a-1,g,x]
!     for the rows: 2 <= a < ngmax
!     which are stable as they are used
      do j = ngmax-1, 2, -1
         h1f1(1,j) = ( h1f1(2,j+1) + re(j-1) * h1f1(1,j+1)) * reali(j)
         aux = arg * reali(j)
         do i = 1, min(namax-1,j-2)
            h1f1(i+1,j) = aux * h1f1(i+1,j+1) + h1f1(i,j)
         enddo
      enddo
!     Computes the basic integrals  < N 0 0 | N' 0 0 >
      if (exa .ge. exb) then
         zzpR = (z + zp) * R
         zzpi = 1.d0 / (z+zp)
         zzpin = zzpi**(nmin+npmin+1)
         do n = nmin, nmax
            zzpinp = zzpin
            do np = npmin, npmax+nmin-n
               auxk = fact(n+np)*zzpinp
               sumk = 0.d0
               do k = 0, n+np
                  sumj = 0.d0
                  do j = max(0,k-n), min(k,np)
                     sumj = sumj + h1f1(j+1,k+2) * facti(n+j-k) * facti(np-j)
                  enddo
                  sumk = sumk + sumj * auxk
                  auxk = auxk * zzpR * reali(n+np-k) * reali(k+2)
               enddo
               smat(n,np) = pi4 * expzR * fact(n) * fact(np) * sumk
               smataux(n,np) = smat(n,np)
               smatm(n,np) = smat(n,np)
               zzpinp = zzpinp * zzpi
            enddo
            zzpin = zzpin * zzpi
         enddo
      else
         zzpR = (z + zp) * R
         zzpi = 1.d0 / (z+zp)
         zzpin = zzpi**(nmin+npmin+1)
         do n = nmin, nmax
            zzpinp = zzpin
            do np = npmin, npmax+nmin-n
               auxk = fact(n+np)*zzpinp
               sumk = 0.d0
               do k = 0, n+np
                  sumj = 0.d0
                  do j = max(0,k-n), min(k,np)
                     sumj = sumj + h1f1(j+1,k+2) * facti(n+j-k) * facti(np-j)
                  enddo
                  sumk = sumk + sumj * auxk
                  auxk = auxk * zzpR * reali(n+np-k) * reali(k+2)
               enddo
               smat(np,n) = pi4 * expzR * fact(n) * fact(np) * sumk
               smataux(np,n) = smat(np,n)
               smatm(np,n) = smat(np,n)
               zzpinp = zzpinp * zzpi
            enddo
            zzpin = zzpin * zzpi
         enddo                               
 1122 format(5i4,3x,e22.15)
      endif
!     Redefines the limits of the indices taking into account that
!     matrix smat is stored as  smat(na,nb)
      nmin = na - la
      nmax = na + la + lb + lb + lmult + lmult + lmaux
      npmin = nb - lb
      npmax = nb + lb + la + la + lmult + lmult + lmaux
      Rinv = 1.d0 / R
      R2 = R * R
      R4 = R2 * R2
      updltm0i = .5d0
      mtop = min(lmult+la,lb)
      do m = 0, mtop

!     Recurrence on LB

         if (lb .gt. m) then
            do n = nmin, nmax-2
               saux = smat(n,npmin)
               do np = npmin+1, npmax-n+nmin-1
                  sbux = -re(m+m+1) * .5d0 * Rinv * ( smat(n,np+1) + R2 * saux - smat(n+2,np-1) )
                  saux = smat(n,np)
                  smataux(n,np) = saux
                  smat(n,np) = sbux
               enddo
            enddo
            do ll = m+1, lb-1
               do n = nmin, nmax+m-2*(ll+1-m)
                  saux = smat(n,npmin+ll-m)
                  do np = npmin+1+ll-m,npmax-n+nmin-m-ll+m-1
                     sbux = -re(ll+ll+1)*.5d0*reali(ll-m+1)*Rinv*(smat(n,np+1)&
                            +R2*saux-smat(n+2,np-1))-re(ll+m)*reali(ll-m+1)*smataux(n,np)
                     saux = smat(n,np)
                     smataux(n,np) = saux
                     smat(n,np) = sbux
                  enddo
               enddo
            enddo
         endif

!     Recurrence on LA

         do k = max(0,(la+1-m)/2), (la+lmult-m)/2
            sol(m,m+2*k-la,k) = smat(na+m+2*k-la,nb)   ! l = 2k-la+m
         enddo
         if (la+lmult .gt. m) then
            do np = npmin+lb-m, npmax-lb+m-2
               saux = smat(nmin,np)
               do n = nmin+1, nmax-2*(lb-m)-np+npmin+lb-m-1
                  sbux = re(m+m+1) * .5d0 * Rinv * ( smat(n+1,np) + R2 * saux - smat(n-1,np+2) )
                  saux = smat(n,np)
                  smataux(n,np) = saux
                  smat(n,np) = sbux
               enddo
            enddo
            do k = max(0,(la-m)/2), (la+lmult-1)/2
               sol(m,m+2*k+1-la,k) = smat(na+m+2*k+1-la,nb)
            enddo
            do ll = m+1, la+lmult-1
               do np = npmin+lb-m, npmax-lb+m-2*(ll+1-m)
                  saux = smat(nmin+ll-m,np)
                  do n = nmin+1+ll-m, nmax-2*(lb-m)-np+npmin+lb-m-ll+m-1
                     sbux = re(ll+ll+1) * .5d0 * reali(ll-m+1) * Rinv * ( smat(n+1,np) + R2 * saux &
                          - smat(n-1,np+2) ) - re(ll+m) * reali(ll-m+1) * smataux(n,np)
                     saux = smat(n,np)
                     smataux(n,np) = saux
                     smat(n,np) = sbux
                  enddo
               enddo
               do k = max(0,(la-ll)/2), (la+lmult-ll-1)/2
                  sol(m,ll+1+2*k-la,k) = smat(na+ll+1+2*k-la,nb)
                enddo
            enddo
         endif

!     Recurrence on M

         if (m .eq. mtop) go to 200
         nmin = nmin + 1
         nmax = nmax - 3
         npmin = npmin + 1
         npmax = npmax - 3
         aux = updltm0i * .25d0 * re(m+m+1) * re(m+m+1) * Rinv * Rinv 
         do n = nmin, nmax
            do np = npmin, npmax-n+nmin
                smataux(n,np) = aux * (2.d0 * (smatm(n+1,np+1) + R2 * (smatm(n+1,np-1) + smatm(n-1,np+1)))&
                              - smatm(n+3,np-1) - smatm(n-1,np+3) - R4 * smatm(n-1,np-1) )
            enddo
         enddo
         do n = nmin, nmax
            do np = npmin, npmax-n+nmin
               smat(n,np) = smataux(n,np)
               smatm(n,np) = smataux(n,np)
            enddo
         enddo
         updltm0i = 1.d0
      enddo
  200 continue
      return
end subroutine solap
!
!   *******************************************************************
!
subroutine acof( lmax )
      implicit real * 8 ( a-h,o-z )
      parameter (mxl = 5, mxmult = 3)
      parameter (mxltot = mxl + mxmult )
      parameter ( mxlcof = mxl*(mxl+3)/2 )
      parameter ( mxkcof = mxlcof*(mxlcof+3)/2 )
      common /abpp/ app(0:mxkcof,0:2*mxl+1), bpp(0:mxkcof,0:2*mxl+1)

      if ( lmax .gt. mxl ) then
         write(6,*) ' lmax greater than mxl in acof'
         stop 
      endif

      do j = 0, 2*mxl+1
      do i = 0, mxkcof
         app(i,j) = 0.d0
      enddo
      enddo
!
!   starting elements app(lm,00)(n) = delta(l,n)
!
      k1 = 0
      do 10 l = 0 , lmax
      do 10 m = 0 , l
         kk = k1*(k1+1) / 2
         app(kk,l) = 1.d0
         k1 = k1 + 1
   10 continue
!
!   elements app(lm,m'm')(n)
!
      do 20 mp = 1 , lmax
         k2 = mp*(mp+1)/2 + mp
         k20 = (mp-1)*mp/2 + mp-1
         do 20 l = mp , lmax
         if ( l.eq.mp ) then
            m1 = mp
         else
            m1 = 0
         endif
         do 20 m = m1 , l
            k1 = l*(l+1)/2 + m
            kk = k1*(k1+1)/2 + k2
            kk0 = k1*(k1+1)/2 + k20
            do 20 n = l-mp , l+mp , 2
               if ( n.ge.m+mp) then
                  app(kk,n) = (2*mp-1) * ( app(kk0,n-1)/dfloat(n+n-1) - app(kk0,n+1)/dfloat(n+n+3) )
               endif
   20 continue
!
!   elements app(lm,l'm')(n)
!
      do 30 mp = 0 , lmax
      do 30 lp = mp+1 , lmax
      k2 = lp*(lp+1)/2 + mp
      k20 = (lp-1)*lp/2 + mp
      k200 = (lp-2)*(lp-1)/2 + mp
      do 30 l = lp , lmax
      if ( l.eq.lp ) then
         m1 = mp
      else
         m1 = 0
      endif
      do 30 m = m1 , l
      k1 = l*(l+1)/2 + m
      kk = k1*(k1+1)/2 + k2
      kk0 = k1*(k1+1)/2 + k20
      kk00 = k1*(k1+1)/2 + k200
      do 30 n = l-lp , l+lp , 2
         if ( n.ge.m+mp) then
            aux = app(kk0,n+1) * (n+m+mp+1) / dfloat(n+n+3)
            if ( n.gt.m+mp ) aux = aux + app(kk0,n-1) * (n-m-mp) / dfloat(n+n-1)
            aux = aux * ( lp+lp-1 )
            if ( lp.gt.mp+1 ) aux = aux - (lp+mp-1) * app(kk00,n)
            app(kk,n) = aux / dfloat(lp-mp)
         endif
   30 continue
      return
end subroutine acof
!
!   *******************************************************************
!
subroutine bcof( lmax )
      implicit real * 8 ( a-h,o-z )
      parameter (mxl = 5, mxmult = 3)
      parameter (mxltot = mxl + mxmult )
      parameter ( mxlcof = mxl*(mxl+3)/2 )
      parameter ( mxkcof = mxlcof*(mxlcof+3)/2 )
      common /abpp/ app(0:mxkcof,0:2*mxl+1), bpp(0:mxkcof,0:2*mxl+1)

      if ( lmax .gt. mxl ) then
         write(6,*) ' lmax greater than mxl in bcof'
         stop 
      endif

      do j = 0, 2*mxl+1
      do i = 0, mxkcof
         bpp(i,j) = 0.d0
      enddo
      enddo

!
!   starting elements bpp(lm,00)(n) = delta(l,n)
!
      k1 = 0
      do 10 l = 0 , lmax
      do 10 m = 0 , l
         kk = k1*(k1+1) / 2
         bpp(kk,l) = 1.d0
         k1 = k1 + 1
   10 continue
!
!   elements bpp(lm,m'm')(n)
!
      do 20 mp = 1 , lmax
         k2 = mp*(mp+1)/2 + mp
         k20 = (mp-1)*mp/2 + mp-1
         do 20 l = mp , lmax
         if ( l.eq.mp ) then
            m1 = mp
         else
            m1 = 0
         endif
         do 20 m = m1 , l
            k1 = l*(l+1)/2 + m
            kk = k1*(k1+1)/2 + k2
            kk0 = k1*(k1+1)/2 + k20
            do 20 n = l-mp , l+mp , 2
               if ( mp.gt.m ) then
                  t1 = 1.d0
                  t2 = 1.d0
               else
                  t1 = -(n-(m-mp+1))*(n-(m-mp+1)+1)
                  t2 = -(n+(m-mp+1))*(n+(m-mp+1)+1)
               endif

               if ( n.ge.abs(m-mp)) then
                  if (n.eq.0) then
                      bux=0.d0
                  else
                      bux=t1*bpp(kk0,n-1)/dfloat(n+n-1)
                  endif
                  bpp(kk,n) = (2*mp-1) * ( bux - t2*bpp(kk0,n+1)/dfloat(n+n+3) )
               endif

   20 continue
!
!   elements bpp(lm,l'm')(n)
!
      do 30 mp = 0 , lmax
      do 30 lp = mp+1 , lmax
      k2 = lp*(lp+1)/2 + mp
      k20 = (lp-1)*lp/2 + mp
      k200 = (lp-2)*(lp-1)/2 + mp
      do 30 l = lp , lmax
      if ( l.eq.lp ) then
         m1 = mp
      else
         m1 = 0
      endif
      do 30 m = m1 , l
      k1 = l*(l+1)/2 + m
      kk = k1*(k1+1)/2 + k2
      kk0 = k1*(k1+1)/2 + k20
      kk00 = k1*(k1+1)/2 + k200
      do 30 n = l-lp , l+lp , 2
         mmp = abs(m-mp)
         if ( n.ge.mmp) then
            aux = bpp(kk0,n+1) * (n+mmp+1) / dfloat(n+n+3)
            if ( n.gt.mmp ) aux = aux + bpp(kk0,n-1) * (n-mmp) / dfloat(n+n-1)
            aux = aux * ( lp+lp-1 )
            if ( lp.gt.mp+1 ) aux = aux - (lp+mp-1) * bpp(kk00,n)
            bpp(kk,n) = aux / dfloat(lp-mp)
         endif
   30 continue
      return
end subroutine bcof
!
!   ******************************************************************
!
subroutine emes ( m1, m2, ms, md, ss, sd )
      implicit real * 8 (a-h,o-z)
      parameter ( pt5 = 0.5d0 )
      s1 = sign(1,m1)
      s2 = sign(1,m2)
      s12 = s1 * s2
      m1a = iabs(m1)
      m2a = iabs(m2)
      ms = s12 * ( m1a + m2a )
      md = s12 * iabs( m1a - m2a )
      if ( ms.eq.md ) then
         ss = 1.d0
         sd = 0.d0
         return
      endif
      if ( m1.lt.0 .and. m2.lt.0 ) then
         ss = -pt5
      else
         ss = pt5
      endif
      if ( s12.gt.0.d0 ) then
         sd = pt5
      elseif ( md.eq.0 ) then
         sd = 0.d0
         elseif ( sign(1,m1a-m2a) .eq. s1 ) then
            sd = - pt5
         else
            sd = pt5
      endif
      return
end subroutine emes 
!***********************************************************************
!                                                                      *
!   subroutine rotar                                                   *
!                                                                      *
!   this subroutine yields the rotation matrices r(m',m;l) that are    *
!   necessary to perform a coordinate transformation used to align     *
!   two sets of normalized real spherical harmonics centered           *
!   at different points(a and b).                                      *
!   this transformation converts each original real spherical harmonic *
!   in a linear combination of the real spherical harmonics with the   *
!   same l and different m.                                            *
!   the maximum value for the orbital quantum number l is 12, to extend*
!   this program to greater values of l it is necessary to extend the  *
!   common sqroot (needed in the subroutine dlmn) with the values of   *
!   the square roots of the first 2*ltot+1 integers and their          *
!   reciprocals.                                                       *
!                                                                      *
!***********************************************************************
subroutine rotar(lmax,ltot,cosal,sinal,cosbet,sinbet,cosga,singa,dl,rl)
      implicit real*8(a-h,o-z)
      parameter (mxl = 5)
      parameter (mxreal = 1000, mxfact = 150, mxind = 200, mxroot=50)
      dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
      dimension dl(-ltot:ltot,-ltot:ltot,0:ltot)
      data zero/0.0d0/,one/1.0d0/
      common /const/ re(0:mxreal), reali(0:mxreal), fact(0:mxfact)&
         , facti(0:mxfact), ang((mxl+1)*(mxl+2)/2)&
         , root(0:mxroot), rooti(mxroot), ind(0:mxind)
!
!     computation of the initial matrices d0, r0, d1 and r1
!
      dl(0,0,0)  = one
      rl(0,0,0)  = one
      if(lmax.eq.0) go to 201
      dl(1,1,1)  = (one+cosbet)/two
      dl(1,0,1)  =-sinbet*rooti(2)
      dl(1,-1,1) = (one-cosbet)/two
      dl(0,1,1)  =-dl(1,0,1)
      dl(0,0,1)  = dl(1,1,1)-dl(1,-1,1)
      dl(0,-1,1) = dl(1,0,1)
      dl(-1,1,1) = dl(1,-1,1)
      dl(-1,0,1) = dl(0,1,1)
      dl(-1,-1,1)= dl(1,1,1)
      cosag  = cosal * cosga - sinal * singa
      cosamg = cosal * cosga + sinal * singa
      sinag  = sinal * cosga + cosal * singa
      sinamg = sinal * cosga - cosal * singa
      rl(0,0,1)  = dl(0,0,1)
      rl(1,0,1)  = root(2) * dl(0,1,1) * cosal
      rl(-1,0,1) = root(2) * dl(0,1,1) * sinal
      rl(0,1,1)  = root(2) * dl(1,0,1) * cosga
      rl(0,-1,1) =-root(2) * dl(1,0,1) * singa
      rl(1,1,1)  = dl(1,1,1) * cosag - dl(1,-1,1) * cosamg
      rl(1,-1,1) =-dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
      rl(-1,1,1) = dl(1,1,1) * sinag - dl(1,-1,1) * sinamg
      rl(-1,-1,1)= dl(1,1,1) * cosag + dl(1,-1,1) * cosamg
!
!     the remaining matrices are calculated using symmetry and
!     recurrence relations by means of the subroutine dlmn.
!
      if ( dabs(sinbet) .lt. 1.d-14 ) then
          tgbet2 = zero
          sinbet = zero
          cosbet = dsign(one,cosbet)
      elseif ( dabs(sinbet) .lt. 1.d-10 ) then
          tgbet2 = zero
          print*, 'WARNING in ROTAR: sinbet = ', sinbet, ' takes  0'
          sinbet = zero
          cosbet = dsign(one,cosbet)
      else
         tgbet2 = ( one - cosbet ) / sinbet
      endif
      do 10 l = 2, lmax
         l1 = l
         call dlmn(l1,ltot,sinal,cosal,cosbet,tgbet2,singa,cosga,dl,rl)
   10 continue
  201 continue
      return
end subroutine rotar
!***********************************************************************
!                                                                      *
!   subroutine dlmn                                                    *
!                                                                      *
!   this subroutine generates the matrices dl(m',m;l) for a fixed value*
!   of the orbital quantum number l, and it needs the dl(l-2;m',m) and *
!   dl(l-1;m',m) matrices. this subroutine uses symmetry and recurrence*
!   relations. the matrices dl(m',m;l) are the rotation matrices for   *
!   complex spherical harmonics.                                       *
!                                                                      *
!***********************************************************************
subroutine dlmn(l,ltot,sinal,cosal,cosbet,tgbet2,singa,cosga,dl,rl)
      implicit real*8 (a-h,o-z)
      parameter (mxl = 5)
      parameter (mxreal = 1000, mxfact = 150, mxind = 200, mxroot=50)
      dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
      dimension dl(-ltot:ltot,-ltot:ltot,0:ltot)
      data one/1.0d0/
      common /const/ re(0:mxreal), reali(0:mxreal), fact(0:mxfact) &
         , facti(0:mxfact), ang((mxl+1)*(mxl+2)/2) &
         , root(0:mxroot), rooti(mxroot), ind(0:mxind)
      iinf=1-l
      isup=-iinf
!
!     computation of the dl(m',m;l) matrix, mp is m' and m is m.
!
!
!     first row by recurrence: see equations 19 and 20 of reference (6)
!
      dl(l,l,l)=dl(isup,isup,l-1)*(one+cosbet)/two
      dl(l,-l,l)=dl(isup,-isup,l-1)*(one-cosbet)/two
      do 10 m=isup,iinf,-1
         dl(l,m,l)=-tgbet2 * root(l+m+1) * rooti(l-m) * dl(l,m+1,l)
   10 continue
!
!     the rows of the upper quarter triangle of the dl(m',m;l) matrix
!     see equation 21 of reference (6)
!
      al=l
      al1= al-one
      tal1= al+al1
      ali=one / al1
      cosaux = cosbet*al*al1
      do 11 mp=l-1,0,-1
         amp=mp
         laux=l+mp
         lbux=l-mp
         aux= rooti(laux) * rooti(lbux) * ali
         cux= root(laux-1) * root(lbux-1) * al
         do 12 m=isup,iinf,-1
            am=m
            lauz=l+m
            lbuz=l-m
            auz= rooti(lauz) * rooti(lbuz)
            factor= aux * auz
            term=tal1*(cosaux-am*amp)*dl(mp,m,l-1)
            if(lbuz.ne.1.and.lbux.ne.1) then
               cuz= root(lauz-1) * root(lbuz-1)
               term=term-dl(mp,m,l-2)*cux*cuz
            endif
            dl(mp,m,l)=factor*term
   12       continue
         iinf=iinf+1
         isup=isup-1
   11 continue
!
!     the remaining elements of the dl(m',m;l) matrix are calculated
!     using the corresponding symmetry relations:
!     reflection ---> ((-1)**(m-m')) dl(m,m';l) = dl(m',m;l), m'<=m
!     inversion ---> ((-1)**(m-m')) dl(-m',-m;l) = dl(m',m;l)
!
!
!     reflection
!
      sign=one
      iinf=-l
      isup=l-1
      do 13 m=l,1,-1
         do 14 mp=iinf,isup
            dl(mp,m,l)=sign*dl(m,mp,l)
            sign=-sign
   14    continue
         iinf=iinf+1
         isup=isup-1
   13 continue
!
!     inversion
!
      iinf=-l
      isup=iinf
      do 15 m=l-1,-l,-1
         sign=-one
         do 16 mp=isup,iinf,-1
            dl(mp,m,l)=sign*dl(-mp,-m,l)
            sign=-sign
   16    continue
         isup=isup+1
   15 continue
!
!     computation of the rotation matrices rl(m',m;l) for real spherical
!     harmonics using the matrices dl(m',m;l) for complex spherical
!     harmonics: see equations 10 to 18 of reference (6)
!
      rl(0,0,l)=dl(0,0,l)
      cosmal = cosal
      sinmal = sinal
      sign = - one
      do 17 mp = 1, l
         cosmga = cosga
         sinmga = singa
         aux = root(2) * dl(0,mp,l)
         rl(mp,0,l) = aux * cosmal
         rl(-mp,0,l)= aux * sinmal
         do 18 m = 1, l
            aux = root(2) * dl(m,0,l)
            rl(0,m,l) = aux * cosmga
            rl(0,-m,l)=-aux * sinmga
            d1 = dl(-mp,-m,l)
            d2 = sign * dl(mp,-m,l)
            cosag = cosmal * cosmga - sinmal * sinmga
            cosagm= cosmal * cosmga + sinmal * sinmga
            sinag = sinmal * cosmga + cosmal * sinmga
            sinagm= sinmal * cosmga - cosmal * sinmga
            rl(mp,m,l)  = d1 * cosag + d2 * cosagm
            rl(mp,-m,l) =-d1 * sinag + d2 * sinagm
            rl(-mp,m,l) = d1 * sinag + d2 * sinagm
            rl(-mp,-m,l)= d1 * cosag - d2 * cosagm
            aux    = cosmga * cosga - sinmga * singa
            sinmga = sinmga * cosga + cosmga * singa
            cosmga = aux
   18    continue
         sign = - sign
         aux    = cosmal * cosal - sinmal * sinal
         sinmal = sinmal * cosal + cosmal * sinal
         cosmal = aux
   17 continue
      return
end subroutine dlmn

!***********************************************************************
!                                                                      *
!   subroutine m2matprt                                                *
!                                                                      *
!   this subroutine is used to print out the matrices rl(l;m',m) with
!   an appropriate format.                                             *
!                                                                      *
!***********************************************************************
subroutine m2cmatprt(lmax,ltot,rl)
      implicit real*8 (a-h,o-z)
      dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
!
 1000 format(12x,3(i3,15x),i3)
 1002 format(1x,i3,2x,4d18.10)
 1003 format(/,2x,'rotation matrix for l = ',i2)
!
      do 2000 l=0,lmax
         write(6,1003) l
         ind=1
         ilow=-l
    1    iupp=ilow+3
         if(iupp-l)3,2,2
    2    iupp=l
         ind=0
    3    write(6,1000)(j,j=ilow,iupp)
         do 4 i=-l,l
            write(6,1002)i,(rl(i,j,l),j=ilow,iupp)
    4    continue
         if(ind.eq.0) go to 5
         ilow=ilow+4
         go to 1
    5    continue
 2000 continue
      return
end subroutine m2cmatprt
!
!
end module Multipole_Core

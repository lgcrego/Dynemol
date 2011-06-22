module Coulomb_SMILES_m
    use constants_m             
    use type_m                  
    use util_m                  , only : fact
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : Extended_Cell              , &
                                         Generate_Structure         , &
                                         Basis_Builder              , &
                                         ExCell_basis


    public  :: Int_Coulomb

    private

    ! module parameters ...
    real*8  , parameter :: umbral = 1.d-10                                      ! Umbral para almacenamiento de las integrales
    integer , parameter :: mxl = 3 , mxbuff = 100000 , mxequiv = 1000
    integer , parameter :: idmrot = (mxl+1)*(2*mxl+1)*(2*mxl+3)/3
    integer , parameter :: idmrotmx = mxequiv * idmrot

!! attention : using atomic units for length  ###########
contains
!
!
!
!=======================
subroutine Int_Coulomb
!=======================
implicit none


! local variables ...
integer :: ia , ib , ind_A , ind_B

CALL Generate_Structure(1)

CALL Basis_Builder( Extended_Cell , ExCell_basis )

! atom indices ...
ia = 1
ib = 2
! basis index ...
ind_A = 3
ind_B = 8

CALL Coulomb( Extended_Cell , ExCell_basis , ind_A , ind_B )

end subroutine Int_Coulomb
!
!
!========================================================
!          the order orbitals are stored
! 
!       S      -->  1   --> l = 0  ,  m =  0           
!       Py     -->  2   --> l = 1  ,  m = -1    
!       Pz     -->  3   --> l = 1  ,  m =  0         
!       Px     -->  4   --> l = 1  ,  m = +1
!       Dxy    -->  5   --> l = 2  ,  m = -2      
!       Dyz    -->  6   --> l = 2  ,  m = -1
!       Dz2    -->  7   --> l = 2  ,  m =  0     
!       Dxz    -->  8   --> l = 2  ,  m = +1        
!       Dx2y2  -->  9   --> l = 2  ,  m = +2        
!========================================================
!
!
!=============================================================
 subroutine Coulomb( system , basis , ind_A , ind_B )
!=============================================================
implicit none
type(structure) , intent(in)  :: system 
type(STO_basis) , intent(in)  :: basis(:)
integer         , intent(in)  :: ind_A
integer         , intent(in)  :: ind_B

! local arrays ...
real*8  , allocatable                                       :: coul(:,:,:,:) , coul_tmp(:,:,:,:) 
real*8  , dimension (-mxl:mxl,-mxl:mxl,0:mxl)               :: rl , rl2
real*8  , dimension (idmrotmx)                              :: rotmat
real*8  , dimension (-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl)   :: v
real*8  , dimension (-mxl:mxl)                              :: vaux

real*8  , allocatable :: V_coul_El(:,:) , V_coul_Hl(:,:)

! local parameters ...
integer :: spdf_indx(0:3) = [1,2,5,10]

! local variables ...
real*8  :: x1 , x2 , x3 , x4 , rn1 , rn2 , rn3 , rn4 , Rab , soma
integer :: i , j , k , l , ij , icaso, indx1 , indx2, indx3, indx4
integer :: na_1 , na_2 , nb_1 , nb_2 , la_1 , la_2 , lb_1 , lb_2, ma_1, ma_2, mb_1, mb_2
integer :: ia, ja1, a1, ja2, a2
integer :: ib, jb1, b1, jb2, b2

!-------------------------------------------------


allocate( V_coul_El(size(basis),size(basis)) , source=D_zero )
allocate( V_coul_Hl(size(basis),size(basis)) , source=D_zero )




CALL consta

! calculate Coulomb potential ...
allocate( Coul     (-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl) , source=0.d0 )
allocate( Coul_tmp (-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl) , source=0.d0 )

do ia = 1 , system % atoms
do ib = 1 , system % atoms

    ! Coulomb rotation matrix ...
    CALL Rotation_Matrix( system , rotmat , icaso , ia , ib , Rab )

    If( Rab < low_prec ) goto 100

    do ja1 = 1 , atom( system%AtNo(ia) )%DOS , 3   ;   a1 = system%BasisPointer(ia) + ja1
    do ja2 = 1 , atom( system%AtNo(ia) )%DOS , 3   ;   a2 = system%BasisPointer(ia) + ja2

        na_1 = basis(a1)% n     ;       la_1 = basis(a1)% l         
        na_2 = basis(a2)% n     ;       la_2 = basis(a2)% l         

        do jb1 = 1 , atom( system%AtNo(ib) )%DOS  , 3  ;   b1 = system%BasisPointer(ib) + jb1
        do jb2 = 1 , atom( system%AtNo(ib) )%DOS  , 3  ;   b2 = system%BasisPointer(ib) + jb2

            Coul = D_zero

            nb_1 = basis(b1)% n     ;       lb_1 = basis(b1)% l         
            nb_2 = basis(b2)% n     ;       lb_2 = basis(b2)% l         

            do i = 1 , basis(a1)% Nzeta
            do j = 1 , basis(a2)% Nzeta
            do k = 1 , basis(b1)% Nzeta
            do l = 1 , basis(b2)% Nzeta

                x1 = basis(a1)%  zeta(i)    ;   rn1 = sqrt( (x1+x1)**(na_1+na_1+1)/fact(na_1+na_1) )
                x2 = basis(a2)%  zeta(j)    ;   rn2 = sqrt( (x2+x2)**(na_2+na_2+1)/fact(na_2+na_2) )
                x3 = basis(b1)%  zeta(k)    ;   rn3 = sqrt( (x3+x3)**(nb_1+nb_1+1)/fact(nb_1+nb_1) )
                x4 = basis(b2)%  zeta(l)    ;   rn4 = sqrt( (x4+x4)**(nb_2+nb_2+1)/fact(nb_2+nb_2) )

                Coul_tmp = D_zero

                ! ===============================================================
                !   Integrales (AA|BB) : < i2(A) j2(B) | Coul | i1(A) j1(B) >  
                ! ===============================================================

                ! lined-up frame ...
                CALL Coul0sim( na_1 , la_1 , x1 , na_2 , la_2 , x2 , nb_1 , lb_1 , x3 , nb_2 , lb_2 , x4 , rn1 , rn2 , rn3 , rn4 , Rab , Coul_tmp )

                ! rotate Coulomb matrix to crystal frame ...
                CALL Rotate_Coulomb( la_1 , la_2 , lb_1 , lb_2 , rotmat , Coul_tmp )

                Coul = Coul + basis(a1)%coef(i) * basis(a2)%coef(j) * basis(b1)%coef(k) * basis(b2)%coef(l) * Coul_tmp

            end do  !   l
            end do  !   k
            end do  !   j
            end do  !   i

            ! ===============================================================================================
            ! build ELECTRON potential ... 
            do ma_2 = -la_2 , la_2      ;    indx2 = la_2 + ma_2 + system%BasisPointer(ia) + spdf_indx(la_2)
            do ma_1 = -la_1 , la_1      ;    indx1 = la_1 + ma_1 + system%BasisPointer(ia) + spdf_indx(la_1)

                V_coul_El(indx1,indx2) = D_zero

                do mb_2 = -lb_2 , lb_2      ;    indx4 = lb_2 + mb_2 + system%BasisPointer(ib) + spdf_indx(lb_2)
                do mb_1 = -lb_1 , lb_1      ;    indx3 = lb_1 + mb_1 + system%BasisPointer(ib) + spdf_indx(lb_1)

                    V_coul_El(indx1,indx2) = V_coul_El(indx1,indx2) + Coul( ma_1 , ma_2 , mb_1 , mb_2 )

                end do
                end do

            end do
            end do

            ! ===============================================================================================
            ! build HOLE potential ... 
            do mb_2 = -lb_2 , lb_2      ;    indx4 = lb_2 + mb_2 + system%BasisPointer(ib) + spdf_indx(lb_2)
            do mb_1 = -lb_1 , lb_1      ;    indx3 = lb_1 + mb_1 + system%BasisPointer(ib) + spdf_indx(lb_1)

                V_coul_Hl(indx3,indx4) = D_zero

                do ma_2 = -la_2 , la_2      ;    indx2 = la_2 + ma_2 + system%BasisPointer(ia) + spdf_indx(la_2)
                do ma_1 = -la_1 , la_1      ;    indx1 = la_1 + ma_1 + system%BasisPointer(ia) + spdf_indx(la_1)

                    V_coul_Hl(indx3,indx4) = V_coul_Hl(indx3,indx4) + Coul( ma_1 , ma_2 , mb_1 , mb_2 )

                end do
                end do

            end do
            end do
            ! ===============================================================================================
           
        end do
        end do

    end do
    end do

100 continue

end do    
end do  

do i = 1 , size(basis)
do j = i , size(basis)
    write(8,*) i,j,V_coul_El(i,j),V_coul_El(j,i)
    write(9,*) i,j,V_coul_Hl(i,j),V_coul_Hl(j,i)
end do
end do













deallocate( Coul , Coul_tmp )

stop


end subroutine Coulomb
!
!
!
!==============================================================
 subroutine Rotate_Coulomb( l1 , l2 , l3 , l4 , rotmat , Coul )
!==============================================================
implicit none
integer , intent(in)    :: l1
integer , intent(in)    :: l2
integer , intent(in)    :: l3
integer , intent(in)    :: l4
real*8  , intent(in)    :: rotmat(:)
real*8  , intent(inout) :: Coul(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl) 

! local arrays ...
real*8  , dimension (-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl)   :: v
real*8  , dimension (-mxl:mxl)                              :: vaux

! local variables ...
real*8  :: soma
integer :: i , j , k , l , ij 
integer :: icaso, knt, m1, m2, m3, m4, lmax


knt = 0

! Carga las integrales semilla ...
v(-l1:l1,-l2:l2,-l3:l3,-l4:l4) = Coul(-l1:l1,-l2:l2,-l3:l3,-l4:l4)

!-------------------------------------------------------------------
! Rota las integrales si es preciso ...

lmax = max(l1,l2,l3,l4)
if( (icaso >= 0) .and. (lmax > 0) ) then

    if (l1 > 0) then
        knt = (2*l1-1)*(2*l1+1)*l1/3
        do 60 m2 = -l2, l2
        do 60 m3 = -l3, l3
        do 60 m4 = -l4, l4
            vaux(-l1:l1) = v(-l1:l1,m2,m3,m4)
            ij = knt
            do 70 i = -l1, l1
                soma = 0.d0
                do 80 j = -l1, l1
                    ij = ij + 1
                    soma = soma + rotmat(ij) * vaux(j)
80              continue
                v(i,m2,m3,m4) = soma
70          continue
60      continue
    endif

    if(l2 > 0) then
        knt = (2*l2-1)*(2*l2+1)*l2/3
        do 160 m1 = -l1, l1
        do 160 m3 = -l3, l3
        do 160 m4 = -l4, l4
            vaux(-l2:l2) = v(m1,-l2:l2,m3,m4)
            ij = knt
            do 170 i = -l2, l2
                soma = 0.d0
                do 180 j = -l2, l2
                    ij = ij + 1
                    soma = soma + rotmat(ij) * vaux(j)
180             continue
                v(m1,i,m3,m4) = soma
170         continue
160     continue
    endif

    if(l3 > 0) then
        knt = (2*l3-1)*(2*l3+1)*l3/3
        do 260 m1 = -l1, l1
        do 260 m2 = -l2, l2
        do 260 m4 = -l4, l4
            vaux(-l3:l3) = v(m1,m2,-l3:l3,m4)
            ij = knt
            do 270 i = -l3, l3
                soma = 0.d0
                do 280 j = -l3, l3
                    ij = ij + 1
                    soma = soma + rotmat(ij) * vaux(j)
280             continue
                v(m1,m2,i,m4) = soma
270         continue
260     continue
    endif

    if(l4 > 0) then
        knt = (2*l4-1)*(2*l4+1)*l4/3
        do 360 m1 = -l1, l1
        do 360 m2 = -l2, l2
        do 360 m3 = -l3, l3
            vaux(-l4:l4) = v(m1,m2,m3,-l4:l4)
            ij = knt
            do 370 i = -l4, l4
                soma = 0.d0
                do 380 j = -l4, l4
                    ij = ij + 1
                    soma = soma + rotmat(ij) * vaux(j)
380             continue
                v(m1,m2,m3,i) = soma
370         continue
360     continue
    endif

end if

! Fin de la rotacion de las integrales  ...

Coul(-l1:l1,-l2:l2,-l3:l3,-l4:l4) = v(-l1:l1,-l2:l2,-l3:l3,-l4:l4)

end subroutine Rotate_Coulomb
!
!
!***********************************************************************
!   subroutine rotar                                                   *
!                                                                      *
!   this subroutine yields the rotation matrices r(m',m;l) that are    *
!   necessary to perform a coordinate transformation used to align     *
!   two sets of real spherical harmonics centered at different points  *
!   (a and b).                                                         *
!   this transformation converts each original real spherical harmonic *
!   in a linear combination of the real spherical harmonics with the   *
!   same l and different m.                                            *
!   the maximum value for the orbital quantum number l is 12, to extend*
!   this program to greater values of l it is necessary to extend the  *
!   common sqroot (needed in the subroutine dlmn) with the values of   *
!   the square roots of the first 2*ltot+1 integers and their          *
!   reciprocals.                                                       *
!                                                                      * 
!   Modificada para calcular las matrices de rotacion al menos hasta   *
!   l = mxmlt (el orden del multipolo mas alto) ...                    *
!=======================================================================
 subroutine Rotation_Matrix( system , rotmat , icaso , ia , ib , Rab )
!=======================================================================
implicit none
type(structure) , intent(in)    :: system
real*8          , intent(out)   :: rotmat(:)
integer         , intent(out)   :: icaso
integer         , intent(in)    :: ia , ib
real*8          , intent(out)   :: Rab

! local arrays ...
real*8  , dimension (-mxl:mxl,-mxl:mxl,0:mxl)   :: rl , rl2

! local variables ...
real*8  :: xaq, yaq, zaq, xbq, ybq, zbq, xa, ya, za, xb, yb, zb, xab, yab, zab, xy 
real*8  :: sinal, cosal, sinbet, cosbet, singa, cosga, cosalf, cosgam, sinalf, singam, den
integer :: i, j, k, l, ma, mb
integer :: ltot, knt
integer :: AtNo_a , AtNo_b
integer :: la_max , lb_max , lmax

! local parameters ...
real*8  , parameter :: tol = 1.d-10

!------------------------------------------------------------------
! Pre-process information for a pair of centers ...

AtNo_a = system%AtNo (ia)
AtNo_b = system%AtNo (ib)

! coordinates must be in a.u. 
xa  = system%coord (ia,1) !/ a_Bohr
ya  = system%coord (ia,2) !/ a_Bohr
za  = system%coord (ia,3) !/ a_Bohr
xb  = system%coord (ib,1) !/ a_Bohr
yb  = system%coord (ib,2) !/ a_Bohr
zb  = system%coord (ib,3) !/ a_Bohr
xab = xb - xa
yab = yb - ya
zab = zb - za
xy  = dsqrt(xab*xab + yab*yab)
Rab = dsqrt(xab*xab + yab*yab + zab*zab) 

cosbet = zab / Rab
cosalf = 1.d0
cosgam = 1.d0

if( abs(cosbet) > (1.d0+tol) ) then
    print*, '>>> problem in Rotation_Matrix <<<'
    stop
elseif( (abs(cosbet) > 1.d0) .or. (1.d0-abs(cosbet)) < tol ) then
        cosbet = sign(1.d0,cosbet)
        sinbet = 0.d0
    else
        sinbet = sqrt(1.d0 - cosbet*cosbet)
endif

if( abs(1.d0-abs(cosbet)) < tol ) then
    if( cosbet > 0.d0 ) then
        ! Caso -1: los puntos estan alineados y sobre el eje Z y rab esta orientado hacia Z positivo ...
        icaso = -1
    else
        ! Caso  0: los puntos estan alineados y sobre el eje Z pero rab esta orientado hacia Z negativo ...
        icaso = 0
        sinbet = 0.d0
        sinalf = 0.d0
        singam = 0.d0
    endif
else
    !  Caso 1: los puntos estan alineados sobre un eje distinto del Z ...
    icaso = 1
    xab = xb - xa
    yab = yb - ya
    den = 1.d0 / sqrt(xab*xab+yab*yab)
    cosalf = xab * den
    sinalf = yab * den
    singam = 0.d0
endif

!  Mira si los centros coinciden con los de la integral original ...
if ( (abs(1.d0 - cosalf) < tol) .and. (abs(1.d0 - cosbet) < tol) .and. (abs(1.d0 - cosgam) < tol) ) then
    ! Caso -2: los tres puntos coinciden con los de la integral original ...
    icaso = -2
endif

!-------------------------------------------------------------------
! Llama a ROTcARt si es necesario rotar las integrales

lmax = max( atom(AtNo_a)%AngMax ,  atom(AtNo_b)%AngMax )

if( (icaso >= 0) .and. (lmax /= 0) ) then

    ltot = mxl
    call rotar_local( LMAX , LTOT , COSALf , SiNALf , COSBET , SiNBET , COSGAm , SiNGAm , rl2 , rl )

    knt = 0
    do l = 0, lmax
        do ma = -l, l
            do mb = -l, l
                knt = knt + 1
                rotmat(knt) = rl(ma,mb,l)
            enddo
        enddo
    enddo

endif
! rotation matrix for pair (ia,ib) finished ...

end subroutine Rotation_Matrix
!
!
!=================================================================================================
 subroutine rotar_local( lmax , ltot , cosal , sinal , cosbet , sinbet , cosga , singa , dl , rl )
!=================================================================================================
implicit none
integer , intent(in)    :: lmax
integer , intent(in)    :: ltot
real*8  , intent(inout) :: cosal
real*8  , intent(inout) :: sinal
real*8  , intent(inout) :: cosbet
real*8  , intent(inout) :: sinbet
real*8  , intent(inout) :: cosga
real*8  , intent(inout) :: singa
real*8  , intent(out)   :: rl(-ltot:ltot,-ltot:ltot,0:ltot)
real*8  , intent(out)   :: dl(-ltot:ltot,-ltot:ltot,0:ltot)

! local variables ...
real*8  :: cosag , sinag , sinamg , tgbet , cosamg , tgbet2
integer :: i , j , k , l , l1 , l2 , l3 , l4

! local parameters ...
real*8  , parameter :: root2 = dsqrt(2.0d0)
real*8  , parameter :: zero  = 0.0d0
real*8  , parameter :: one   = 1.0d0
real*8  , parameter :: two   = 2.0d0
!
!     computation of the initial matrices d0, r0, d1 and r1
!
      dl(0,0,0)  = one
      rl(0,0,0)  = one
      if(lmax.eq.0) go to 201
      dl(1,1,1)  = (one+cosbet)/two
      dl(1,0,1)  =-sinbet/root2
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
      rl(1,0,1)  = root2 * dl(0,1,1) * cosal
      rl(-1,0,1) = root2 * dl(0,1,1) * sinal
      rl(0,1,1)  = root2 * dl(1,0,1) * cosga
      rl(0,-1,1) =-root2 * dl(1,0,1) * singa
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

end subroutine rotar_local
!
!
!
end module Coulomb_SMILES_m

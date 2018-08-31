c  Modificada en julio de 2002 para calcular las matrices de rotacion
c  al menos hasta  l = mxmlt (el orden del multipolo mas alto)
c
c  difiere del m2caux2.f en la subrutina stoel1. he cambiado la forma
c  de almacenamiento de las matrices stoa y stob. ahora funciona para
c  cualquier valor de mxn y mxl
c
c  modificada en septiembre de 2006 para los nuevos algoritmos de las
c  integrales tricentricas
c
c---------------------------------------------------------------------
c
      subroutine subrot (rl,rl2,rlvec)
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      include 'm2cdat.inc'
      dimension rlvec(mxrotdm)
     &        , rl(-mxl:mxl,-mxl:mxl,0:mxl)
     &        , rl2(-mxl:mxl,-mxl:mxl,0:mxl)

      do i = 1, mxrotdm
         rlvec(i) = 0.d0
      enddo

c
c     calcula las matrices de rotacion y las almacena en el vector
c     rlvec (por columnas)
c     calcula tambien las distancias entre pares de centros
c
      ienlac = 0
      irl = 0
      do 50 k = 2, ncen
         lb = lmaxc(k)
         do 60 i = 1, k-1
            lmax = max(1,lmaxc(i),lb,mxmlt)
c
c  el "1" en el max(...) es debido a que es necesario en algunas
c  ocasiones (3c1e) rotar coordenadas.
c
            ienlac = ienlac + 1
            xab = r(k,1) - r(i,1)
            yab = r(k,2) - r(i,2)
            zab = r(k,3) - r(i,3)
            xy = dsqrt(xab*xab + yab*yab)
            a = distab(ienlac) 
            if (xy .gt. tol) then
               sinal = yab / xy
               cosal = xab / xy
            else
               sinal = 0.d0
               cosal = 1.d0
            endif
            sinbet = xy / a
            cosbet = zab / a
            singa = 0.d0
            cosga = 1.d0
            call rotar (lmax, mxl, cosal, sinal, cosbet, sinbet, cosga,
     ,               singa, rl2, rl)
            irlind(ienlac) = irl + 1
            do 62 l = 0, lmax
               do 64 mp = -l,l
               do 64 m = -l,l
                  irl = irl+1
                  rlvec(irl) = rl(m,mp,l)
   64          continue
   62       continue
   60    continue
   50 continue
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine leerl (j, i, rlvec, rl)
      implicit real * 8 (a-h,o-z)
c
c     importante:  j: indice del centro "a"
c                  i: indice del centro "b"
c
      include 'comun.inc'
      include 'm2cdat.inc'
      dimension rlvec(mxrotdm), rl(-mxl:mxl,-mxl:mxl,0:mxl)
c
c       recupera en  rl  las matrices de rotacion del centro  i  al  j
c

      if (i .eq. j) stop 'error en leerl:  i = j '
      lmax = max( 1 , lmaxc(i), lmaxc(j), mxmlt )
c
c  ver comentario en subrot acerca del "1" en el max(...)
c
      if (i .gt. j) then
         k = irlind(ind(i-2)+j)
         do 10 l = 0, lmax
         do 10 mp = -l,l
         do 10 m = -l,l
            rl(m,mp,l) = rlvec(k)
            k = k + 1
   10    continue
      else
         k = irlind(ind(j-2)+i)
         sl = 1.d0
         do 20 l = 0, lmax
           sm = -1.d0
           do 30 mp = -l,-1
             do 40 m = -l,l
               rl(m,mp,l) = sm * rlvec(k)
               k = k + 1
   40        continue
             sm = -sm
   30      continue
           sm = sl
           do 50 mp = 0,l
             do 60 m = -l,l
               rl(m,mp,l) = sm * rlvec(k)
               k = k + 1
   60        continue
             sm = -sm
   50      continue
           sl = -sl
   20    continue
      endif
  902 format(2x,'the rotation matrices for lmax up to ',i2,' are:')
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine leerl2 (i, j, rlvec, rl, rl2)
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      include 'm2cdat.inc'
      dimension rlvec(mxrotdm), rl(-mxl:mxl,-mxl:mxl,0:mxl)
     &        , rl2(-mxl:mxl,-mxl:mxl,0:mxl)
c
c       recupera en  rl  las matrices de rotacion del centro  i  al  j
c             y en  rl2  las matrices de rotacion del centro  j  al  i
c
      if (i .eq. j) stop 'error en leerl:  i = j '
      lmax = max( 1, lmaxc(i), lmaxc(j), mxmlt )
c
c  ver comentario en subrot acerca del "1" en el max(...)
c
      k = irlind( ind( max(i,j)-2 )+min(i,j) )
      sl = 1.d0
      do 20 l = 0, lmax
        sm = -1.d0
        do 30 mp = -l,-1
          do 40 m = -l,l
            rl(m,mp,l) = rlvec(k)
            rl2(m,mp,l) = sm * rlvec(k)
            k = k + 1
   40     continue
          sm = -sm
   30   continue
        sm = sl
        do 50 mp = 0,l
          do 60 m = -l,l
               rl(m,mp,l) = rlvec(k)
               rl2(m,mp,l) = sm * rlvec(k)
               k = k + 1
   60     continue
          sm = -sm
   50   continue
        sl = -sl
   20 continue

  902 format(2x,'the rotation matrices for lmax up to ',i2,' are:')

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine consta
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      PARAMETER (PI = 3.14159265358979312D0)
      PARAMETER (dospi = 2.d0 * PI)
      parameter ( raiz2 = 1.414213562373095d0 )
      Common / extrapolacion / extrapol(6,2)
      Common / cuadratura5 / ugausspat0(5), wgausspat0(5)
      Common / cuadratura11 / ugausspat1(11), wgausspat1(11)
      Common / cuadratura23 / ugausspat2(23), wgausspat2(23)
      Common / cuadratura47 / ugausspat3(47), wgausspat3(47)
      Common / cuadratura95 / ugausspat4(95), wgausspat4(95)
      Common / cuadratura191 / ugausspat5(191), wgausspat5(191)
      Common / cuadratura383 / ugausspat6(383), wgausspat6(383)
      common /emecom/ msv(-lmxmnach:lmxmnach,-lmxmnach:lmxmnach),
     &                mdv(-lmxmnach:lmxmnach,-lmxmnach:lmxmnach),
     &                ssv(-lmxmnach:lmxmnach,-lmxmnach:lmxmnach),
     &                sdv(-lmxmnach:lmxmnach,-lmxmnach:lmxmnach)

      ind(0) = 0
    1 continue
      do 10 i = 1, mxind
         ind(i) = ind(i-1) + i
   10 continue
      root(0 ) = 0.d0
      do 15 i = 1, mxroot
         root(i) = dsqrt(dfloat(i))
         rooti(i) = 1.d0 / root(i)
   15 continue
      fact(0) = 1.d0
      facti(0) = 1.d0
      facts(-1) = 0.177245385090552d1
      factsi(-1) = 1.D0/facts(-1)                                       
      facts(0) = facts(-1) * .5d0
      factsi(0) = 1.D0/facts(0)
      do 20 i = 1, mxfact
         fact(i) = fact(i-1) * dfloat(i)
         facts(i) = facts(i-1) * (2*i+1) * .5d0
         facti(i) = 1.d0 / fact(i)
         factsi(i) = 1.d0 / facts(i)                                    
   20 continue
      real(0) = 0.d0
      dosl1(0) = 1.d0
      dosl1p(0) = 1.d0
      rll1(0) = 0.d0
      ri(0) = 1.d250
      do 30 l = 1, mxreal
         real(l) = dfloat(l)
         real(-l) = -real(l)                                            
         dosl1(l) = dfloat(2*l+1)
         dosl1(-l) = dfloat(-2*l+1)                                     
         dosl1p(l) = 1 / dosl1(l)
         dosl1p(-l) = 1 / dosl1(-l)                                     
         ri(l) = 1.d0 / real(l)
         ri(-l) = -ri(l)                                                
         rll1(l) = dfloat(l*(l+1))
   30 continue
      c1r(0) = real(2) * ri(3)
      c2r(0) = c1r(0)
      c3r(0) = 0.d0
      do 40 l = 1, mxc123
         c1r(l) = rll1(l+1) * dosl1p(l) * dosl1p(l+1)
         c2r(l) = 1.d0 - real(l+1)**2 * dosl1p(l) * dosl1p(l+1)
     &          - real(l) ** 2 * dosl1p(l-1) * dosl1p(l)
         c3r(l) = rll1(l-1) * dosl1p(l-1) * dosl1p(l)
   40 continue
      ang(0,0) = 0.282094791773878d0
      do 50 l = 1, mxangl
         ang(l,0) = ang(0,0) * dsqrt(dfloat(2*l+1))
         aux = ang(l,0) * root(2)
         do 50 m = 1, min(l,mxangm)
            aux = aux / dsqrt(dfloat((l-m+1)*(l+m)))
            ang(l,m) = aux
   50 continue

      do 60 i = 0, linter
        faux(i,0) = dosl1(i)
        do 60 m = 1, min(i,2*mxl)
           faux(i,m) = dosl1(i) * (fact(i-m) / fact(i+m))**2
   60 continue
      rl2l1p(0) = 0.d0
      rl12l1(0) = 1.d0
      do 61 i = 1, linter
        rl2l1p(i) = real(i) * dosl1p(i)
        rl12l1(i) = real(i+1) * dosl1p(i)
   61 continue

      lmax = mxlcf

      call acof(lmax)
      call bcof(lmax)
      call stoel1
      call zzzhv
      call cpcomg
      call cgen
c     para las integrales tricentricas con armonicos irregulares
c     cuadratura xknum 0-1
      call abcpes5(ugausspat0,wgausspat0,0.D0)
      call abcpes11(ugausspat1,wgausspat1,0.D0)
      call abcpes23(ugausspat2,wgausspat2,0.D0)
      call abcpes47(ugausspat3,wgausspat3,0.D0)
      call abcpes95(ugausspat4,wgausspat4,0.D0)
      call abcpes191(ugausspat5,wgausspat5,0.D0)
      call abcpes383(ugausspat6,wgausspat6,0.D0)
c     extrapolacion
      extrapol(1,1)=6.30957D-2
      extrapol(1,2)=1.D-2
      extrapol(2,1)=1.58489D-2
      extrapol(2,2)=1.D-3
      extrapol(3,1)=3.98107D-3
      extrapol(3,2)=1.D-4
      extrapol(4,1)=1.D-3
      extrapol(4,2)=1.D-5
      extrapol(5,1)=2.51189D-4
      extrapol(5,2)=1.D-6
      extrapol(6,1)=6.30957D-5
      extrapol(6,2)=1.D-7
      do m2 = -lmxmnach, lmxmnach
      do m1 = -lmxmnach, lmxmnach
         call emes ( m1, m2, ms, md, ss, sd )                           
         msv(m1,m2) = ms                                                
         mdv(m1,m2) = md                                                
         ssv(m1,m2) = ss                                                
         sdv(m1,m2) = sd 
      enddo
      enddo
      do l2 = 0,lmxmnach                                                
      do l1 = 0,lmxmnach                                                
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
            il=-1
            do l = l1+l2,Iabs(l1-l2),-2                                 
            il=il+1
            av((lmxmnach+1)*k12+il+1)=app(k12,l)
            bv((lmxmnach+1)*k12+il+1)=bpp(k12,l)
            end do                                                      
      end do                                                            
      end do                                                            
      end do                                                            
      end do                                                            
      dm(0) = 1.d0                          
      unom(0) = 1.d0                          
      dpf(0) = 1.d0                          
      p5(0) = 1.d0                          
      do i = 1,40        
         dm(i) = dm(i-1)*(-2.D0)    
         p5(i) = p5(i-1)*0.5D0     
         dpf(i) = 2.D0*dpf(i-1)/dfloat(i)
         unom(i) = unom(i-1)*(-1.D0)      
      enddo                                   
      Do l=0,10                     
      f=Dfloat(2*l+1)*1.591549430918953D-1          
      xarm(l*(l+1)+1)=Dsqrt(0.5D0*f)                   
      Do m=1,l                                         
      xarm(l*(l+1)+m+1)=Dsqrt(f*fact(l-m)*facti(l+m))  
      xarm(l*(l+1)-m+1)=xarm(l*(l+1)+m+1)              
      enddo                                            
      enddo                                            
cip                                                                     

      kalg1 = 0
      kalg2 = 0

      do m2 = -2*mxl, 2*mxl
      do m1 = -2*mxl, 2*mxl
         call emes ( m1, m2, ms, md, ss, sd )
         msvec(m1,m2) = ms
         mdvec(m1,m2) = md
         ssvec(m1,m2) = ss
         sdvec(m1,m2) = sd
      enddo
      enddo

      lma = 0
      do la = 0, mxl
      do ma = -la, la
         mma = iabs(ma)
         lma = lma + 1
         if (ma .eq. 0) then
            aux = raiz2 * dospi
         else
            aux = dospi
         endif
         lmb = 0
         do lb = 0, mxl
         do mb = -lb, lb
            mmb = iabs(mb)
            lmb = lmb + 1
            if (mb .eq. 0) then
               bux = raiz2 * aux
            else
               bux = aux
            endif
            auxlnor(lma,lmb) = bux * dsqrt(fact(la+mma) * fact(lb+mmb)
     &                    / (fact(la-mma) * fact(lb-mmb)) )
         enddo
         enddo
      enddo
      enddo

*     Los coeficientes  cfbk0  y  cfbk1  se usan en la relacion de
*     recurrencia de las funciones K de Bessel redefinidas como:
*             bk(n) = Exp[z] * BesselK[n+1/2,z] / (n+1/2)!
*     cfbk0(i) = 1 / ((2i+1) * (2i+3))
*     cfbk1(i) = (2i+1) / (2i+3)

      do i = 0, mxk
         cfbk0(i) = ri(i+i+3) * ri(i+i+1)
         cfbk1(i) = real(i+i+1) * ri(i+i+3)
      enddo

*     tabula los coeficientes   ckplm   para los factores radiales
*     en la traslacion de distribuciones bicentricas

      call subckplma(mxk,mxk,ckplm)
      call subckplmb(mxk,mxk,ckplm)
      call subckplmc(mxk,mxk,ckplm)

*     tabula los coeficientes   cffk   para los factores radiales
*     de las funciones 2s  (para acelerar el calculo)

      do k = 1, mxk
         cffk21(k) = -2.d0*real(k)
         cffk22(k) = real(2*k+1)
         cffk23(k) = -2.d0*real(k+1)
      enddo

      cffk21(0) = 0.d0
      cffk22(0) = 1.d0
      cffk23(0) = -2.d0

      return
      end

c***********************************************************************
c                                                                      *
c   subroutine rotar                                                   *
c                                                                      *
c   this subroutine yields the rotation matrices r(m',m;l) that are    *
c   necessary to perform a coordinate transformation used to align     *
c   two sets of real spherical harmonics centered at different points  *
c   (a and b).                                                         *
c   this transformation converts each original real spherical harmonic *
c   in a linear combination of the real spherical harmonics with the   *
c   same l and different m.                                            *
c   the maximum value for the orbital quantum number l is 12, to extend*
c   this program to greater values of l it is necessary to extend the  *
c   common sqroot (needed in the subroutine dlmn) with the values of   *
c   the square roots of the first 2*ltot+1 integers and their          *
c   reciprocals.                                                       *
c                                                                      *
c***********************************************************************
      subroutine rotar(lmax,ltot,cosal,sinal,cosbet,sinbet,cosga,singa,
     &                 dl,rl)
      implicit real*8(a-h,o-z)
      dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
      dimension dl(-ltot:ltot,-ltot:ltot,0:ltot)
      data root2/1.4142135623730950488016887242096970d0/,
     -     zero/0.0d0/,one/1.0d0/,two/2.0d0/
c
c     computation of the initial matrices d0, r0, d1 and r1
c
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
c
c     the remaining matrices are calculated using symmetry and
c     recurrence relations by means of the subroutine dlmn.
c
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
      end
c***********************************************************************
c                                                                      *
c   subroutine dlmn                                                    *
c                                                                      *
c   this subroutine generates the matrices dl(m',m;l) for a fixed value*
c   of the orbital quantum number l, and it needs the dl(l-2;m',m) and *
c   dl(l-1;m',m) matrices. this subroutine uses symmetry and recurrence*
c   relations. the matrices dl(m',m;l) are the rotation matrices for   *
c   complex spherical harmonics.                                       *
c                                                                      *
c***********************************************************************
      subroutine dlmn(l,ltot,sinal,cosal,cosbet,tgbet2,singa,cosga,
     &                dl,rl)
      implicit real*8 (a-h,o-z)
      dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
      dimension dl(-ltot:ltot,-ltot:ltot,0:ltot)
      include 'comun.inc'
      data root2/1.4142135623730950488016887242096970d0/,
     -     one/1.0d0/,two/2.0d0/
      iinf=1-l
      isup=-iinf
c
c     computation of the dl(m',m;l) matrix, mp is m' and m is m.
c
c
c     first row by recurrence: see equations 19 and 20 of reference (6)
c
      dl(l,l,l)=dl(isup,isup,l-1)*(one+cosbet)/two
      dl(l,-l,l)=dl(isup,-isup,l-1)*(one-cosbet)/two
      do 10 m=isup,iinf,-1
         dl(l,m,l)=-tgbet2 * root(l+m+1) * rooti(l-m) * dl(l,m+1,l)
   10 continue
c
c     the rows of the upper quarter triangle of the dl(m',m;l) matrix
c     see equation 21 of reference (6)
c
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
c
c     the remaining elements of the dl(m',m;l) matrix are calculated
c     using the corresponding symmetry relations:
c     reflexion ---> ((-1)**(m-m')) dl(m,m';l) = dl(m',m;l), m'<=m
c     inversion ---> ((-1)**(m-m')) dl(-m',-m;l) = dl(m',m;l)
c
c
c     reflexion
c
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
c
c     inversion
c
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
c
c     computation of the rotation matrices rl(m',m;l) for real spherical
c     harmonics using the matrices dl(m',m;l) for complex spherical
c     harmonics: see equations 10 to 18 of reference (6)
c
      rl(0,0,l)=dl(0,0,l)
      cosmal = cosal
      sinmal = sinal
      sign = - one
      do 17 mp = 1, l
         cosmga = cosga
         sinmga = singa
         aux = root2 * dl(0,mp,l)
         rl(mp,0,l) = aux * cosmal
         rl(-mp,0,l)= aux * sinmal
         do 18 m = 1, l
            aux = root2 * dl(m,0,l)
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
      end

************************************************************************
*                                                                      *
*   subroutine matprt                                                  *
*                                                                      *
*   this subroutine is used to print out the matrices rl(l;m',m) with
*   an appropriate format.                                             *
*                                                                      *
************************************************************************
*     subroutine matprt(lmax,ltot,rl)
      subroutine m2cmatprt(lmax,ltot,rl)
      implicit real*8 (a-h,o-z)
      dimension rl(-ltot:ltot,-ltot:ltot,0:ltot)
c
 1000 format(12x,3(i3,15x),i3)
 1002 format(1x,i3,2x,4d18.10)
 1003 format(/,2x,'rotation matrix for l = ',i2)
c
c     do 2000 l=1,lmax
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
      end

c
c   *******************************************************************
c
      subroutine acof( lmax )
      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'

      do j = 0, 2*mxlcf+1
      do i = 0, mxkcof
         app(i,j) = 0.d0
      enddo
      enddo

      if ( lmax .gt. mxlcf ) stop ' lmax mayor que mxlcf en acof'
c
c   elementos de partida app(lm,00)(n) = delta(l,n)
c
      k1 = 0
      do 10 l = 0 , lmax
      do 10 m = 0 , l
         kk = k1*(k1+1) / 2
         app(kk,l) = 1.d0
         k1 = k1 + 1
   10 continue
c
c   elementos del tipo app(lm,m'm')(n)
c
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
                  app(kk,n) = (2*mp-1) * ( app(kk0,n-1)/dfloat(n+n-1)
     &                                 - app(kk0,n+1)/dfloat(n+n+3) )
               endif
   20 continue
c
c   elementos del tipo app(lm,l'm')(n)
c
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
            if ( n.gt.m+mp )
     &        aux = aux + app(kk0,n-1) * (n-m-mp) / dfloat(n+n-1)
            aux = aux * ( lp+lp-1 )
            if ( lp.gt.mp+1 ) aux = aux - (lp+mp-1) * app(kk00,n)
            app(kk,n) = aux / dfloat(lp-mp)
         endif
   30 continue
      return
      end
c
c   *******************************************************************
c
      subroutine bcof( lmax )
      implicit real * 8 ( a-h,o-z )
      include 'comun.inc'

      if ( lmax .gt. mxlcf ) stop ' lmax mayor que mxlcf en bcof'

      do j = 0, 2*mxlcf+1
      do i = 0, mxkcof
         bpp(i,j) = 0.d0
      enddo
      enddo

c
c   elementos de partida bpp(lm,00)(n) = delta(l,n)
c
      k1 = 0
      do 10 l = 0 , lmax
      do 10 m = 0 , l
         kk = k1*(k1+1) / 2
         bpp(kk,l) = 1.d0
         k1 = k1 + 1
   10 continue
c
c   elementos del tipo bpp(lm,m'm')(n)
c
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
                  bpp(kk,n) = (2*mp-1) * ( bux
     &                                 - t2*bpp(kk0,n+1)/dfloat(n+n+3) )
               endif

   20 continue
c
c   elementos del tipo bpp(lm,l'm')(n)
c
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
            if ( n.gt.mmp )
     &        aux = aux + bpp(kk0,n-1) * (n-mmp) / dfloat(n+n-1)
            aux = aux * ( lp+lp-1 )
            if ( lp.gt.mp+1 ) aux = aux - (lp+mp-1) * bpp(kk00,n)
            bpp(kk,n) = aux / dfloat(lp-mp)
         endif
   30 continue
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine indab ( l1 , m1 , l2 , m2 , l1l1 , l2l2 , k12 )
c
c  devuelve la fila, k, en las matrices a y b, correspondiente a la
c  descomposicion de dos polinomios de legendre
c
      implicit real*8 (a-h,o-z)
      include 'comun.inc'

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
      return
      end
c
c   ******************************************************************
c
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
      end
c***********************************************************************
c                                                                      *
c   subroutine matpr                                                   *
c                                                                      *
c   this subroutine is used to print out                               *
c   matrices with an appropriate format.                               *
c                                                                      *
c                                                                      *
c***********************************************************************
      subroutine matpr(idim, jdim, imax, jmax, sm)
      implicit real*8 (a-h,o-z)
      dimension sm(idim,jdim)
c
 1000 format(12x,3(i3,15x),i3)
 1002 format(1x,i3,2x,4d18.10)
c
         indi=1
         jlow=1
    1    jupp=jlow+3
         if(jupp-jmax)3,2,2
    2    jupp=jmax
         indi=0
    3    write(6,1000)(j,j=jlow,jupp)
         do 4 i= 1, imax
            write(6,1002)i,(sm(i,j),j=jlow,jupp)
    4    continue
         if(indi.eq.0) go to 5
         jlow=jlow+4
         go to 1
    5    continue
      return
      end
c
c  *******************************************************************
c
      subroutine ro1 ( la, lb, lc, ld, v, r )
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      dimension r(-mxl:mxl,-mxl:mxl,0:mxl)
      dimension v(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl),aux(-mxl:mxl)
      if ( la.eq.0 ) return
      do 3 md =-ld,ld
      do 3 mc =-lc,lc
      do 3 mb =-lb,lb
         do 2 ma =-la,la
            bux = 0.d0
            do 1 m = -la,la
               bux = bux + r(ma,m,la) * v(m,mb,mc,md)
    1       continue
            aux(ma) = bux
    2   continue
        do 4 ma = -la,la
           v(ma,mb,mc,md) = aux(ma)
    4   continue
    3 continue
      return
      end
c
c  *******************************************************************
c
      subroutine ro2 ( la, lb, lc, ld, v, r )
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      dimension r(-mxl:mxl,-mxl:mxl,0:mxl)
      dimension v(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl),aux(-mxl:mxl)
      if ( lb.eq.0 ) return
      do 3 md =-ld,ld
      do 3 mc =-lc,lc
      do 3 ma =-la,la
         do 2 mb =-lb,lb
            bux = 0.d0
            do 1 m = -lb,lb
               bux = bux + r(mb,m,lb) * v(ma,m,mc,md)
    1       continue
            aux(mb) = bux
    2   continue
        do 4 mb = -lb,lb
           v(ma,mb,mc,md) = aux(mb)
    4   continue
    3 continue
      return
      end
c
c  *******************************************************************
c
      subroutine ro3 ( la, lb, lc, ld, v, r )
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      dimension r(-mxl:mxl,-mxl:mxl,0:mxl)
      dimension v(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl),aux(-mxl:mxl)
      if ( lc.eq.0 ) return
      do 3 md =-ld,ld
      do 3 mb =-lb,lb
      do 3 ma =-la,la
         do 2 mc =-lc,lc
            bux = 0.d0
            do 1 m = -lc,lc
               bux = bux + r(mc,m,lc) * v(ma,mb,m,md)
    1       continue
            aux(mc) = bux
    2   continue
        do 4 mc = -lc,lc
           v(ma,mb,mc,md) = aux(mc)
    4   continue
    3 continue
      return
      end
c
c  *******************************************************************
c
      subroutine ro4 ( la, lb, lc, ld, v, r )
      implicit real * 8 (a-h,o-z)
      include 'comun.inc'
      dimension r(-mxl:mxl,-mxl:mxl,0:mxl)
      dimension v(-mxl:mxl,-mxl:mxl,-mxl:mxl,-mxl:mxl),aux(-mxl:mxl)
      if ( ld.eq.0 ) return
      do 3 mc =-lc,lc
      do 3 mb =-lb,lb
      do 3 ma =-la,la
         do 2 md =-ld,ld
            bux = 0.d0
            do 1 m = -ld,ld
               bux = bux + r(md,m,ld) * v(ma,mb,mc,m)
    1       continue
            aux(md) = bux
    2   continue
        do 4 md = -ld,ld
           v(ma,mb,mc,md) = aux(md)
    4   continue
    3 continue
      return
      end
c
c     ******************************************************************
c
      subroutine dfexi(x,ex,ier)
c     
c     exponential integral  
c     
      implicit real * 8 (a-h,o-z)
      dimension cp(8)
      xa   =dabs(x)
      if (xa .lt. 1.d-40) then
         ex   =1.d75
         write (6,*) 'argumento de la exponencial integral: ', xa
         write (6,*) 'devuelve como valor de la funcion: ', ex
         return
      endif
      if (xa .ge. 170.d0) then
         ex   =0.d0
         return
      endif
      rx   =1/xa
      if (x .gt. 0.d0) then

*        argumento positivo

         if (xa .ge. 4.d0) then
           ex   =dexp(-xa)*rx*(0.5d0+(0.5d0-rx*((((((((
     1     +1.4899084997294817d1)*rx  +2.1060773714263329d3)*rx
     2     +6.5760969874802118d3)*rx  +6.1688521005547635d3)*rx
     3     +2.3960194324749054d3)*rx  +4.2753267120198854d2)*rx
     4     +3.4406199500668489d1)*rx  +9.999999999999734d-1)/((((((((
     5     +2.3781389910216022d3)*rx  +1.1149775287109662d4)*rx
     6     +1.6324145355778350d4)*rx  +1.0337075308584098d4)*rx
     7     +3.1902723748954330d3)*rx  +4.9434507020990365d2)*rx
     8     +3.6406199500645980d1)*rx  +1.d0)))
           return
         endif
         if (xa .ge. 1.d0) then
           ex   =dexp(-xa)*((((((((
     1     +4.0137766494066472d-1)*rx  +8.836718088038439d+0 )*rx
     2     +4.2520203476884078d+1)*rx  +6.9927945129100302d+1)*rx
     3     +4.5593064425338982d+1)*rx  +1.1848310555494584d+1)*rx
     4     +9.999955193013903d-1 )*rx  +8.677459548384437d-8)/((((((((
     5     +9.088045691888692d-2 )*rx  +3.7955900376212224d+0)*rx
     6     +3.1497184917044075d+1)*rx  +8.973110971252898d+1 )*rx
     7     +1.0664518376991388d+2)*rx  +5.6443356956180320d+1)*rx
     8     +1.2848193537915665d+1)*rx  +1.d0)
           return
         endif
         ex   =-dlog(xa)+((((((
     1     +5.0341618409756839d0)*xa  +1.0109380616190553d2)*xa
     2     +1.8421108866799964d3)*xa  +9.938313889620368d3 )*xa
     3     +5.7721724713944435d4)*xa  -4.4178547172821655d4)/(((((xa
     4     +3.7229835283332743d1)*xa  +6.3541941837838170d2)*xa
     5     +6.1061079424575950d3)*xa  +3.2597188129027473d4)*xa
     6     +7.653733233376136d4)
         return
      else

*     argumento negativo

         if (xa .ge. 174.d0) then
            ex   =-1.d75
            write (6,*) 'argumento de la exponencial integral: ', xa
            write (6,*) 'devuelve como valor de la funcion: ', ex
            return
         endif
         if (xa .gt. 24.d0) then
           ex   =-dexp(xa)*rx*(1.d0+rx*(1.0000000000000049d0
     1     +1.9999999999904810d0/(-3.0000000032098127d0 +xa
     2     -2.9999989404032496d0/(-5.0000664041313100d0 +xa
     3     -7.992435957763397d0 /(-7.0681097789502936d0 +xa
     4     -1.2018776354715474d1/(-1.5285662363692964d1 +xa
     5     +7.0483184718042468d1/(-7.631477016202536d+0 +xa
     6     +1.1717922050208646d2/(-2.7979852862430539d1 +xa
     7     +1.3779039023574800d2/(-1.8194966492986891d1 +xa
     8     +3.9727710910041452d0/(-2.2312767077763241d2 +xa
     9     +3.9784597716741472d4/(+1.7533880126546597d2 +xa)))))))))))
           return
         endif
         if (xa .ge. 12.d0) then
           ex   =-dexp(xa)*rx*(9.999933106160569d-1
     1     +1.0015338520453427d0/(-1.8450862323912787d0 +xa
     2     -1.0935561953910912d1/(+2.6525758184527998d1 +xa
     3     +1.9910044708177425d2/(+2.4954877304020594d1 +xa
     4     +1.1928324239686010d3/(-3.3236125793439623d1 +xa
     5     +4.4294131783379284d1/(-9.134835699998743d-1 +xa
     6     +2.5388193156307080d2/(-2.1057407995480405d1 +xa
     7     +5.9949323256674074d1/(-1.0006419139892848d1 +xa
     8     +6.4038004053524156d1/(-1.8600921217264376d1 +xa
     9     +9.792403599217290d+1/(-1.6477211724634631d0 +xa))))))))))
           return
         endif
         if (xa .ge. 6.d0) then
           ex   =-dexp(xa)*rx*(9.989576665165517d-1
     1     +1.1462525324901619d0/(+5.7311670574450802d0 +xa
     2     -1.9914960023123516d2/(+4.1810242256285662d0 +xa
     3     +3.4136521252437554d2/(+5.8865824075328111d0 +xa
     4     +5.2316556873455861d1/(-1.9413296751443070d1 +xa
     5     +3.1727948925436933d2/(+7.894722092944572d0  +xa
     6     -8.387670841896407d0 /(+2.3273023383903914d1 +xa
     7     +9.654052174292803d2 /(-3.6778311347831146d1 +xa
     8     +2.6398300731802459d0/(-2.4694098344836127d0 +xa)))))))))
           return
         endif
         rx   =xa/3.d0-1.d0
         cp(1)=+4.1262666724891194d2
         cp(2)=+2.8544688181364702d4
         cp(3)=+7.371477901846574d4
         cp(4)=+8.986832916437583d6
         cp(5)=-3.3689956420159190d6
         cp(6)=+7.0592160959005675d8
         cp(7)=-2.5038999488635136d8
         cp(8)=+5.9956994689237001d9
         t1   =+1.1063954724163958d1
         t2   =0.d0
         do 20 k=1,8
            t0   =rx*t1
            t    =t0-t2+t0
            t2   =t1
            t1   =t+cp(k)
   20    continue
         ex   =t1-t0
         cp(1)=-3.0252368223822741d3
         cp(2)=+8.537710001807491d4
         cp(3)=-1.4957545720255922d6
         cp(4)=+1.7715830801079988d7
         cp(5)=-1.4498071439302388d8
         cp(6)=+8.028277829469565d8
         cp(7)=-2.7967335112298459d9
         cp(8)=+2.5592649760761635d9
         t1   =+5.1257812500000000d1
         t2   =0.d0
         do 30 k=1,8
            t0   =rx*t1
            t    =t0-t2+t0
            t2   =t1
            t1   =t+cp(k)
   30    continue
         ex   =ex/(t1-t0)
         x01 =  0.372507410781366621d0
         x02 =  0.131401834109754090d-16
         rx   =(xa-x01)-x02
         if (dabs(rx) .le. .37d-1) then
           t0   =rx/x01
           t    =t0*((((((
     1     +3.6546222413236843d-4)*t0  +4.8908925378927915d-2)*t0
     2     +5.1049927962321940d-1)*t0  +1.7027705960680930d+0)*t0
     3     +2.2306993766689975d+0)*t0  +1.d0)/((((((
     4     +1.3111415119497771d-2)*t0  +2.2881793399052641d-1)*t0
     5     +1.2176596296015153d+0)*t0  +2.7347869510692584d+0)*t0
     6     +2.7306993766689975d+0)*t0  +1.d0)
         endif
         t    =dlog(xa/x01)
         ex   =-t-rx*ex
      endif
      return
      end
c
c  *****************************************************************
c
      subroutine stoel1
c
c  subrutina que calcula la descomposicion en coordenadas elip. de 1 sto
c  r**n * (-1)**m * p(l,m), colocado en a ( matriz stoa ) y en b
c  ( matriz stob ). se factoriza el producto
c  (a/2)**n * ( (xi**2-1)(1-eta**2) )**(m/2)
c
c  nueva version. 19-10-94
c
      implicit real*8 (a-h,o-z)
      include 'comun.inc'
      dimension dfact(0:10)
c     dfact(i) = (2i-1)!!
      dfact(0) = 1.d0
      do 2 i = 1 , 10
    2   dfact(i) = (i+i-1) * dfact(i-1)
      nmax = mxn
      lmax = mxl
      nmax1 = nmax+1
      lmax1 = lmax+1
      lmax12 = lmax*lmax1/2
      ln = nmax1*lmax1 - lmax12
      do 3 m = 0 , lmax
      im = m*( 6*ln + (m-1)*( m-2-3*nmax1) ) / 6
      do 3 l = m , lmax
         if(l.eq.m) then
            k = im
            stoa(k,0,0) = dfact(m)
            stob(k,0,0)= dfact(m)
         else
            k = im + ( (l-m) * ( 2*nmax1 + 1 -l - m ) ) / 2
            k1= im + ( (l-1-m) * ( 2*nmax1 + 1 -(l-1) - m ) ) / 2
            if (l.gt.m+1)
     &      k2= im + ( (l-2-m) * ( 2*nmax1 + 1 -(l-2) - m ) ) / 2
            do 5 i = 0 , l-m
            do 5 j = 0 , l-m
               aux = stoa(k1,i,j)
               aux2= -stob(k1,i,j)
               if(i.gt.0 .and. j.gt.0) then
                  aux = aux + stoa(k1,i-1,j-1)
                  aux2= aux2+ stob(k1,i-1,j-1)
               endif
               bux = 0.d0
               bux2= 0.d0
               if (l.gt.m+1) then
                  if (i.gt.1) then
                     bux = bux + stoa(k2,i-2,j)
                     bux2= bux2+ stob(k2,i-2,j)
                  endif
                  if (j.gt.1) then
                     bux = bux + stoa(k2,i,j-2)
                     bux2= bux2+ stob(k2,i,j-2)
                  endif
                  if (i.gt.0 .and. j.gt.0) then
                     bux = bux + 2*stoa(k2,i-1,j-1)
                     bux2= bux2- 2*stob(k2,i-1,j-1)
                  endif
               endif
               stoa(k,i,j)=( (2*l-1)*aux - (l+m-1)*bux ) / dfloat(l-m)
               stob(k,i,j)=( (2*l-1)*aux2 - (l+m-1)*bux2 ) / dfloat(l-m)
    5       continue
         endif
            do 6 n = l+1, nmax
               k1 = k
               k = k+1
               do 6  i = 0 , n-m
               do 6  j = 0 , n-m
                  aux = 0.d0
                  aux2= 0.d0
                  if (i.gt.0) then
                     aux = aux + stoa(k1,i-1,j)
                     aux2= aux2+ stob(k1,i-1,j)
                  endif
                  if (j.gt.0) then
                     aux = aux + stoa(k1,i,j-1)
                     aux2= aux2- stob(k1,i,j-1)
                  endif
                  stoa(k,i,j) = aux
                  stob(k,i,j)= aux2
    6       continue
    3 continue

      return
      end
c
c   ********************************************************************
c
      function indzzz(m,la,lb,k)
      common / zzzdat / lmax , mxlb
      l1 = max(la,lb)
      l2 = min(la,lb)
      mcom = (m*(16 + 14*lmax - m - 4*m**2 - 2*lmax*m**2 + m**3
     &     + 20*mxlb+ 36*lmax*mxlb- 24*mxlb**2 + 24*lmax*mxlb**2
     &     - 32*mxlb** 3))/12
      lpcom = ((l2 - m)*(5 + 3*lmax + 3*l2 + 3*lmax*l2 - 2*l2**2
     &      + 3*m + 3*lmax*m - 2*l2*m - 2*m**2))/6
      lcom = (l1-l2)*(1 + l2)
      kcom = 1 + ( k - (l1-l2) ) / 2
      indzzz = mcom + lpcom + lcom + kcom
      return
      end
c
c     ******************************************************************
c
      subroutine cpcomg
      IMPLICIT REAL * 8 (A-H,O-Z)
      include 'comun.inc'
      PARAMETER (MDIM = mxl)
      PARAMETER (PI = 3.14159265358979312D0)
      PARAMETER (RAIZPI = 1.77245385090551602D0)
      PARAMETER (ltot = MDIM)
      COMMON/CPCOM/ CPAM((MDIM+1)*(MDIM+2)/2,0:MDIM,0:MDIM),
     &              CNAM((MDIM+1)*(MDIM+2)/2,0:MDIM,0:MDIM)
      DIMENSION CPA(0:MDIM,0:MDIM)
      DIMENSION CNA(0:MDIM,0:MDIM)

      lm = 0
      DO 30 L = 0, MDIM
      DO 30 M = 0, l
         lm = lm + 1
         DO 40 I = 0, MDIM
         DO 40 J = 0, MDIM
            CPA(I,J) = 0.D0
            CNA(I,J) = 0.D0
  40     CONTINUE
         CALL RECUR (L, M, CPA, CNA)
         DO 50 I = 0, MDIM
         DO 50 J = 0, MDIM
           CPAM(lM,I,J) = CPA(I,J)
           CNAM(lM,I,J) = CNA(I,J)
  50    CONTINUE
  30  CONTINUE

      return
      end
c
c    *****************************************************************
c
          SUBROUTINE RECUR (LA, MA, CPA, CNA)
          IMPLICIT REAL * 8 (A-H,O-Z)
          include 'comun.inc'
          PARAMETER (MDIM = mxl)
          DIMENSION CP1(0:MDIM,0:MDIM), CN1(0:MDIM,0:MDIM),
     ,            CP2(0:MDIM,0:MDIM), CN2(0:MDIM,0:MDIM),
     ,            CPA(0:MDIM,0:MDIM), CNA(0:MDIM,0:MDIM)


*       intializes the matrices to zeroes

          DO 10 I = 0, MDIM
          DO 10 J = 0, MDIM
             CP1(I,J) = 0.D0
             CN1(I,J) = 0.D0
             CP2(I,J) = 0.D0
             CN2(I,J) = 0.D0
   10   CONTINUE

          CP1(0,0) = 1.D0
          MAA = IABS(MA)
          MAA1= MAA+1
          INDg= 0

*       RECURsion on MA

          IF (MAA .GT. 0) THEN
             DO 20 M = 0, MAA-1
                IF (INDg.EQ. 0) THEN
                     CALL REQM (M, CP1, CN1, CP2, CN2)
                     INDg= 1
                ELSE
                     CALL REQM (M, CP2, CN2, CP1, CN1)
                     INDg= 0
                ENDIF
   20      CONTINUE
          ENDIF

*       RECURsion on LA

          IF (LA .GT. MAA) THEN
             AUX = DFLOAT(MAA+MAA1)
             IF (INDg.EQ. 0) THEN
                DO 30 I = 0, MAA
                DO 30 J = 0, I
                     CP2(I,J) = AUX * CP1(I,J)
                     CP2(J,I) = AUX * CP1(J,I)
                     CN2(I,J) = AUX * CN1(I,J)
                     CN2(J,I) = AUX * CN1(J,I)
   30           CONTINUE
                DO 35 I = 0, MAA1
                     CP2(I,MAA1) = 0.D0
                     CP2(MAA1,I) = 0.D0
                     CN2(I,MAA1) = 0.D0
                     CN2(MAA1,I) = 0.D0
   35           CONTINUE
                INDg= 1
             ELSE
                DO 40 I = 0, MAA
                DO 40 J = 0, I
                     CP1(I,J) = AUX * CP2(I,J)
                     CP1(J,I) = AUX * CP2(J,I)
                     CN1(I,J) = AUX * CN2(I,J)
                     CN1(J,I) = AUX * CN2(J,I)
   40         CONTINUE
                DO 45 I = 0, MAA1
                     CP1(I,MAA1) = 0.D0
                     CP1(MAA1,I) = 0.D0
                     CN1(I,MAA1) = 0.D0
                     CN1(MAA1,I) = 0.D0
   45           CONTINUE
                INDg= 0
             ENDIF
             IF (LA .GT. MAA+1) THEN
                DO 50 L = MAA+1, LA-1
                     IF (INDg.EQ. 0) THEN
                        CALL REQL (L, MAA, CP1, CP2)
                        CALL REQL (L, MAA, CN1, CN2)
                        INDg= 1
                     ELSE
                        CALL REQL (L, MAA, CP2, CP1)
                        CALL REQL (L, MAA, CN2, CN1)
                        INDg= 0
                     ENDIF
   50         CONTINUE
             ENDIF
          ENDIF


*       Charges the final results oN CPA and CNA

          IF (INDg.EQ. 0) THEN
             DO 60 I = 0, LA
             DO 60 J = 0, I
                CPA(I,J) = CP1(I,J)
                CPA(J,I) = CP1(J,I)
                CNA(I,J) = CN1(I,J)
                CNA(J,I) = CN1(J,I)
   60      CONTINUE
             INDg= 1
          ELSE
             DO 70 I = 0, LA
             DO 70 J = 0, I
                CPA(I,J) = CP2(I,J)
                CPA(J,I) = CP2(J,I)
                CNA(I,J) = CN2(I,J)
                CNA(J,I) = CN2(J,I)
   70      CONTINUE
             INDg= 0
          ENDIF

          RETURN
          END

c
c  ****************************************************************
c
          SUBROUTINE REQM (M, CP1, CN1, CP2, CN2)
          IMPLICIT REAL * 8 (A-H,O-Z)
          include 'comun.inc'
          PARAMETER (MDIM = mxl)
          DIMENSION CP1(0:MDIM,0:MDIM), CN1(0:MDIM,0:MDIM),
     ,            CP2(0:MDIM,0:MDIM), CN2(0:MDIM,0:MDIM)

          AUX = DFLOAT(2*M+1)
          M1 = M+1
          CP2(0,0) = 0.D0

          DO 5 I = 0, M1
             CP1(I,M1) = 0.D0
             CP1(M1,I) = 0.D0
             CN1(I,M1) = 0.D0
             CN1(M1,I) = 0.D0
    5     CONTINUE

          DO 10 I = 1, M1
             CP2(I,0) = AUX *  CP1(I-1,0)
             CN2(I,0) = AUX *  CN1(I-1,0)
   10   CONTINUE

          DO 20 J = 1, M1
             CP2(0,J) = - AUX *  CN1(0,J-1)
             CN2(0,J) = AUX *  CP1(0,J-1)
   20   CONTINUE

          DO 30 I = 1, M1
          DO 30 J = 1, M1
             CP2(I,J) = AUX * ( CP1(I-1,J) - CN1(I,J-1) )
             CN2(I,J) = AUX * ( CP1(I,J-1) + CN1(I-1,J) )
   30   CONTINUE

          RETURN
          END
c
c  ****************************************************************
c
          SUBROUTINE REQL (L, M, CP1, CP2)
          IMPLICIT REAL * 8 (A-H,O-Z)
          include 'comun.inc'
          PARAMETER (MDIM = mxl)
          DIMENSION CP1(0:MDIM,0:MDIM), CP2(0:MDIM,0:MDIM)

          BUX = 1.D0 / DFLOAT(L-M+1)
          AUX = DFLOAT(2*L+1) * BUX
          BUX = DFLOAT(L+M) * BUX
          L1 = L+1

          DO 10 I = 0, L1
             CP2(I,L1) = 0.D0
             CP2(I,L) = 0.D0
             CP2(L1,I) = 0.D0
             CP2(L,I) = 0.D0
             CP1(I,L1) = 0.D0
             CP1(L1,I) = 0.D0
   10     CONTINUE

          DO 30 I = L1, 2, -1
          DO 30 J = L1, 2, -1
             CP2(I,J) = AUX *  CP1(I,J) - BUX * ( CP2(I,J)
     +              +  CP2(I,J-2) + CP2(I-2,J) )
   30   CONTINUE

          DO 40 I = L1, 2, -1
             CP2(I,1) = AUX *  CP1(I,1) - BUX * ( CP2(I,1)
     +              + CP2(I-2,1) )
             CP2(I,0) = AUX *  CP1(I,0) - BUX * ( CP2(I,0)
     +              + CP2(I-2,0) )
   40   CONTINUE

          DO 50 J = L1, 2, -1
             CP2(1,J) = AUX *  CP1(1,J) - BUX * ( CP2(1,J)
     +              + CP2(1,J-2) )
             CP2(0,J) = AUX *  CP1(0,J) - BUX * ( CP2(0,J)
     +              + CP2(0,J-2) )
   50   CONTINUE

          CP2(1,1) = AUX * CP1(1,1) - BUX * CP2(1,1)
          CP2(1,0) = AUX * CP1(1,0) - BUX * CP2(1,0)
          CP2(0,1) = AUX * CP1(0,1) - BUX * CP2(0,1)
          CP2(0,0) = AUX * CP1(0,0) - BUX * CP2(0,0)

          RETURN
          END
C***********************************************************************
C                                                                      *
C     SUBROUTINE CGEN                                                  *
C                                                                      *
C                                                                      *
C***********************************************************************
      SUBROUTINE CGEN
      IMPLICIT REAL* 8 (A-H,O-Z)
      include 'comun.inc'
      PARAMETER ( NDIM = ndim3c )
      dimension C3c(0:NDIM3c,0:NDIM3c,0:NDIM3c+NDIM3c)
C
C     SUBRUTINA PARA GENERAR LOS COEFICIENTES C(N,NP,K)
C
      KONT = 1
      DO 10 NP = 0, NDIM
         C3c(0,NP,0) = (-1.D0)**NP
         VC3c(KONT) = C3c(0,NP,0)
         KONT = KONT + 1
         DO 20 K = 0, NP-1
            C3c(0,NP,K+1) = C3c(0,NP,K)*DFLOAT(K-NP)/DFLOAT(K+1)
            VC3c(KONT) = C3c(0,NP,K+1)
            KONT = KONT + 1
   20    CONTINUE
         DO 30 N = 1, NDIM
            C3c(N,NP,0) = (-1.D0)**NP
            VC3c(KONT) = C3c(N,NP,0)
            KONT = KONT + 1
            DO 40 K = 1, N+NP-1
               C3c(N,NP,K) = C3c(N-1,NP,K-1)+C3c(N-1,NP,K)
               VC3c(KONT) = C3c(N,NP,K)
            KONT = KONT + 1
   40       CONTINUE
            C3c(N,NP,N+NP) = 1.D0
            VC3c(KONT) = C3c(N,NP,N+NP)
            KONT = KONT + 1
   30    CONTINUE
   10 CONTINUE
      RETURN
      END

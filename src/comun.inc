*
*   Version de enero de 2008
*

      logical lsdisf

      real*4 tarray, tresult

*                                     PARAMETERS
*                                     ==========

*     Los siguientes  "parameter"  se usan fundamentalmente para
*     dimensionar el programa:
*
*         mxcap: numero maximo de capas (grupos de funciones que solo
*                  difieren en el numero cuantico "m")
*         mxfun: numero maximo de funciones de base
*         mxcen: numero maximo de centros (atomos+nucleos desnudos+
*                                          centros de base sin nucleo)
c                IMPORTANTE!!!  Para evitar posibles problemas en la 
c                parte de aprovechamiento de la simetria local,
c                conviene tomar siempre   mxcen > 6
*         mxl:   maximo valor del numero cuantico "l" en las funciones
*                de la base
*         mxn:   maximo valor del numero cuantico "n" en las funciones
*                de la base
*         ngmax: maximo numero de gaussianas para la expansion STO-nG
*         nenlac: maximo numero de pares de centros
*         mxequiv: maximo numero de pares, trios o tetradas de centros
*                equivalentes
*         mxrotdm: tamagno maximo del vector para almacenar las matrices
*                 de rotacion entre todos los pares de centros
*         isimloc: dimension de los vectores asociados con la simetria
*                local 
*         mxroot: entero maximo cuya raiz cuadrada se almacena en
*                 un vector auxiliar (para acelerar partes criticas)
*         mxreal: entero maximo cuyo valor en doble prec. se almacena en
*                 un vector auxiliar (para acelerar partes criticas)
*         mxfact: entero maximo cuyo factorial (dbl prec.)se almacena en
*                 un vector auxiliar (para acelerar partes criticas)

      parameter (mxcap = 1, mxfun = 511, mxcen = 2, mxl = 3, mxn = 6)
      parameter (ngmax = 30)
      parameter (nenlac = mxcen * (mxcen-1)/2 )
      parameter (mxequiv = 100)
      parameter (idmrot = (mxl+1)*(2*mxl+1)*(2*mxl+3)/3)
      parameter (idmrotmx = mxequiv * idmrot)
      parameter (mxrotdm = nenlac*(mxl+1)*(4*(mxl+1)**2-1)/3)
      parameter (isimloc = mxcen*(mxcen-1)*(mxcen-2)*(mxcen-3)/8)
      parameter (mxroot = 25, mxreal = 2000, mxfact = 170
     &  , mxc123 = 150, mxangl = 15, mxangm = 2*mxl
     &  , mxind = 3000 )

      parameter ( mxbuff = 100000 ) ! tamanho del buffer de integrales

      parameter (mxmlt = 3)  ! maximo  l  de los multipolos

* para las matrices auxiliares de las integrales de intercambio
      parameter ( linter = 64 )

* el siguiente parameter reserva espacio para todas las matrices
* de la subrutina inte0sim (elip7.f)
      parameter ( mxwlk  = ((linter+1)**2) * ((2*mxn+1)**2) 
     &                   + (linter+1) * ((2*mxn+1)**2) * (2*mxl+1)
     &                   + 2 * (2*mxl+1) * (linter+1) * 2 * (2*mxn+1)
     &                   + 4 * ((2*mxl+1)**2) * (linter+1)* 2 *(2*mxn+1)
     &                   + 4 * ((mxl+1)**2) *(2*mxn+1)**2
     &)

      parameter ( umbrzn = 1.d-10 )  ! Umbral de comparacion para la
                                     ! carga nuclear

c con el parameter isto ajusto stoa(b) (common elipab) de forma general
      parameter(isto = mxl*(mxl-1)*(2*mxl+5)/6 + mxn*(mxl+1)*(mxl+2)/2)

c dimension para el data de las hibridas
      parameter ( izzzh =
     &  ((1+3*mxl+2*mxl**2)*(3-mxl-6*mxl**2+6*mxn+8*mxl*mxn))/3 )

c para las tricentricas de atraccion nuclear
      PARAMETER ( NDIM3c = mxn + 3*mxl )

c para las integrales (AB|CC)
      parameter (mxlexp = 10)                                           
      parameter (mxldst = mxl+mxl)                                      
      parameter (mxltot = mxldst + mxlexp)                              
      parameter ( mxl2 = mxl+mxl )                                      
      parameter ( mxlcf = mxl+mxlexp )                                  
      parameter ( nlegABCC = 64, nlagABCC = 20)
      parameter ( nqdABCC = 4*nlegABCC+nlagABCC)

* para las matrices de descomposicion del producto de armonicos
      parameter ( mxlcof = mxlcf*(mxlcf+3)/2 )                          
      parameter ( mxkcof = mxlcof*(mxlcof+3)/2 )                        

c para los programas de integrales tricentricas con armonicos irregulares
      parameter (nmxnach=6, lmxnach=5, lmxmnach=10)  ! Limitaciones del programa
                                                     ! de Nacho
      parameter ( mxemes = max(lmxmnach, 2*mxl) )    ! Dimensiones de los
                                                     ! componentes del common
                                                     ! emescom
      parameter ( mxauxil = 40 )
      parameter ( mxxarm = 121 )
      parameter ( mxlcofn = lmxmnach*(lmxmnach+3)/2 )                   
      parameter ( mxkcofn = mxlcofn*(mxlcofn+3)/2 )                     
      parameter ( mxmul = 14641 )                                       

c maxima longitud de las expansiones en la traslacion de distribuciones
      parameter ( mxk = 650)

*                        COMMONS
*                        =======

      common / abcof / av((mxkcofn+1)*(lmxmnach+1)),
     &                 bv((mxkcofn+1)*(lmxmnach+1))
                                                                        
      common / abpp / app(0:mxkcof,0:2*mxlcf+1)
     &  , bpp(0:mxkcof,0:2*mxlcf+1)             

      common / angcom / ang(0:mxangl,0:mxangm)

      common / an3c / VC3c((2*NDIM3c+1)*(NDIM3c+1)**2)

      Common / auxiliares / dpf(0:mxauxil), p5(0:mxauxil)
     &       , unom(0:mxauxil) , dm(0:mxauxil)

      common / auxlnorcom / auxlnor((mxl+1)**2,(mxl+1)**2)

      common / cfkcom / cfbk1(0:mxk) , cfbk0(0:mxk) , ckplm(3311,0:mxk)
     &       , cffk21(0:mxk), cffk22(0:mxk), cffk23(0:mxk)
                                                                        
      common / const / ri(-mxreal:mxreal), real(-mxreal:mxreal)         
     ,     , fact(0:mxfact), facts(-1:mxfact) , rll1(0:mxreal)          
     ,     , dosl1(-mxreal:mxreal), dosl1p(-mxreal:mxreal)              
     ,     , facti(0:mxfact), factsi(-1:mxfact)                         
     ,     , c1r(0:mxc123),c2r(0:mxc123),c3r(0:mxc123)

      common / elipab / stoa(0:isto,0:mxn,0:mxn)
     &                 ,stob(0:isto,0:mxn,0:mxn)

      common /emescom/ ssvec(-mxemes:mxemes,-mxemes:mxemes)
     & ,sdvec(-mxemes:mxemes,-mxemes:mxemes)
     & ,msvec(-mxemes:mxemes,-mxemes:mxemes)
     & ,mdvec(-mxemes:mxemes,-mxemes:mxemes)

      common / indcom / ind(0:mxind)

      common /indcomn/ indk12((lmxmnach+1)**2,(lmxmnach+1)**2)

      common / intercambio / faux( 0:linter, 0:2*mxl)
     &               , rl2l1p(0:linter),rl12l1(0:linter)

      common / kontcom / kalg1, kalg2

      common / logsdisf / lsdisf(mxcap,mxcap)

      common / simetria / nseed, npares, nternas, ncuatern

      common / sqroot / root(0:mxroot), rooti(mxroot)

      Common / xnorarm / xarm(mxxarm)


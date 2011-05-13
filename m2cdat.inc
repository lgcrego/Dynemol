*     Version de marzo de 1999
*       incluye un indice  indbs  que identifica la base asociada a cada
*       centro
      logical lsolo1e, langstrom, lnosimetr, l4cent, lgeomfix, lno34c2e
      character*20 projectname

      common /datos/ tol, umbral, rep
     &   , xx(mxcap), r(mxcen,3) , rnor(mxcap), zn(mxcen)
     &   , distab(nenlac)
     &   , ncen, nbas, ncaps
     &   , ibscap(mxcap), nf(mxcap), ngini(mxcen+1), ngfin(mxcen)
     &   , nn(mxcap) , ll(mxcap), lmaxc(mxcen), indbs(mxcen)
     &   , lsolo1e, langstrom, lnosimetr, l4cent, lgeomfix, lno34c2e
     &   , projectname

      common / rotac / irlind(nenlac)

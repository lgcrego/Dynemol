module CG_class_m

    use type_m
    use constants_m             , only : a_Bohr
    use Semi_Empirical_Parms    , only : atom
    use parameters_m            , only : driver, profiling, &
                                         DP_Moment ,        &
                                         Alpha_Tensor 
    use Structure_Builder       , only : Extended_Cell 
    use GA_QCModel_m            , only : GA_eigen ,         &
                                         GA_DP_Analysis ,   &
                                         AlphaPolar 
    use cost_tuning_m           , only : evaluate_cost                              

    private

    public :: CG_OPT

    type , extends(OPT) :: CG_OPT
        real*8          , allocatable   :: p(:)
        integer                         :: ITMAX = 50               ! <== 50-200 is a good compromise of accuracy and safety
        real*8                          :: BracketSize = 1.d-4      ! <== this value may vary between 1.0d-3 and 1.0d-5
        type(STO_basis) , allocatable   :: basis(:)
        character (len=11)              :: driver
        logical                         :: profiling = .false.
    contains
        procedure :: cost
        procedure :: cost_variation
        procedure :: output => Dump_OPT_parameters
        procedure , private :: modify_EHT_parameters
        procedure , private :: normalization
    end type CG_OPT

    interface CG_OPT
        module procedure   preprocess 
    end interface

    ! module variables ...
    type(STO_basis) , allocatable :: CG_basis(:)  
    type(R_eigen)                 :: CG_UNI  

contains
!
!
!
!=================================================
 function preprocess( GA_basis , GA ) result( me )
!=================================================
implicit none
type(STO_basis) , intent(in)  :: GA_basis(:)
type(OPT)       , intent(in)  :: GA

type(CG_OPT) :: me

me % driver     = driver
me % profiling  = profiling

me % DP = GA % DP

allocate( me % basis    (size(GA_basis)                     ) , source=GA_basis    )
allocate( me % erg      (size(GA%erg)                       ) , source=GA%erg      )
allocate( me % EHSymbol (size(GA%EHSymbol)                  ) , source=GA%EHSymbol )
allocate( me % key      (size(GA%key(:,1)),size(GA%key(1,:))) , source=GA%key      )

allocate( me % p(GA%GeneSize) , source=D_zero ) 

end function preprocess
!
!
!
!
!===================
 function cost( me )
!===================
implicit none
class(CG_OPT) , intent(inout)  :: me
real*8                         :: cost

!local variables ...
integer     :: info
real*8      :: CG_DP(3) , CG_Alpha_ii(3)

If( .NOT. allocated(CG_basis) ) allocate( CG_basis(size(me%basis)) , source = me%basis)

call me % modify_EHT_parameters( )

info = 0
CALL GA_eigen( Extended_Cell , me%basis , CG_UNI , info )

If( DP_Moment ) CALL GA_DP_Analysis( Extended_Cell , me%basis , CG_UNI%L ,CG_UNI%R , CG_DP )

If( Alpha_Tensor ) CALL AlphaPolar( Extended_Cell , me%basis , CG_Alpha_ii )

cost = evaluate_cost( CG_UNI , me%basis , CG_DP , CG_Alpha_ii )

end function cost
!
!
!
!====================================
 subroutine cost_variation( me , df )
!====================================
implicit none
class(CG_OPT)   , intent(in)    :: me
real*8          , intent(out)   :: df(:)

! local parameters ...
real*8  , parameter :: small = 1.d-8

! local variables ...
integer         :: i , GeneSize
type(CG_OPT)    :: before , after


GeneSize = size(me%p(:))

do i = 1 , GeneSize

    after  = me
    before = me

    after  % p(i) = me % p(i) + small
    before % p(i) = me % p(i) - small

    df(i) = ( after%cost() - before%cost() ) / (two*small)

end do

end subroutine cost_variation
!
!
!
!======================================
 subroutine modify_EHT_parameters( me )
!======================================
implicit none
class(CG_OPT)   , intent(inout)    :: me

! local variables ...
integer :: L , gene , EHS , N_of_EHSymbol 
integer :: indx(size(CG_basis)) , k , i
real*8  :: zeta(2) , coef(2)

! -----------------------------------------------
!       changing basis: editting functions ...
! -----------------------------------------------
!                CG%key storage
!
!                     EHSymbol    --->
!               | 1 - S
!               | 2 - P
!               V 3 - D
!                 4 - IP
!                 5 - zeta1
!                 6 - zeta2
!                 7 - k_WH
! -----------------------------------------------

indx = [ ( i , i=1,size(CG_basis) ) ]

N_of_EHSymbol = size(me%EHSymbol)

gene = 0

do EHS = 1 , N_of_EHSymbol

!   S , P , D  orbitals ...
do  L = 0 , 2

    If( me%key(L+1,EHS) == 1 ) then

        ! changes VSIP ...
        gene = gene + me%key(4,EHS)
        If( me%key(4,EHS) == 1 ) where( (CG_basis%EHSymbol == me%EHSymbol(EHS)) .AND. (CG_basis%l == L) ) me%basis%IP = me%basis%IP + me%p(gene)

        ! single STO orbitals ...
        gene = gene + me%key(5,EHS) - me%key(6,EHS)
        If( (me%key(5,EHS) == 1) .AND. (me%Key(6,EHS) == 0) ) &
        where( (CG_basis%EHSymbol == me%EHSymbol(EHS)) .AND. (CG_basis%l == L) ) me%basis%zeta(1) = me%basis%zeta(1) + me%p(gene)  

        ! double STO orbitals ...
        If( (me%key(5,EHS) == 1) .AND. (me%key(6,EHS) ==1) ) then

            ! finds the first EHT atom ...
            k = minloc( indx , dim=1 , MASK = (CG_basis%EHSymbol == me%EHSymbol(EHS)) .AND. (CG_basis%l == L) ) 

            ! calculate  coef(1)[ zeta(1) , zeta(2) , coef(2) ] ...
            gene    = gene + 1
            zeta(1) = abs( me%p(gene) + me%basis(k)%zeta(1) )
            gene    = gene + 1
            zeta(2) = abs( me%p(gene) + me%basis(k)%zeta(2) )
            coef(2) = abs( me%p(gene) )
            CALL me % normalization( zeta , coef , CG_basis(k)%n , k , 1 , 2 )
            where( (CG_basis%EHSymbol == me%EHSymbol(EHS)) .AND. (CG_basis%l == L) ) 
                me % basis % zeta(1) = zeta(1) 
                me % basis % zeta(2) = zeta(2) 
                me % basis % coef(1) = coef(1)
                me % basis % coef(2) = coef(2)
            end where

        End If

        ! changes k_WH ...
        gene = gene + me%key(7,EHS)
        If( me%key(7,EHS) == 1 ) where( (CG_basis%EHSymbol == me%EHSymbol(EHS)) .AND. (CG_basis%l == L) ) me%basis%k_WH = me%basis%k_WH + me%p(gene) 

    end If

end do
end do

end subroutine modify_EHT_parameters
!
!
!
!============================================================
 subroutine normalization( me , zeta , coef , n , k , i , j )
!============================================================
implicit none
class(CG_OPT)   , intent(in)    :: me
real*8          , intent(inout) :: zeta(:)
real*8          , intent(inout) :: coef(:)
integer         , intent(in)    :: n
integer         , intent(in)    :: k
integer         , intent(in)    :: i
integer         , intent(in)    :: j

! local variables ...
real*8  :: zeta_tmp(size(zeta)) , coef_tmp(size(coef))
real*8  :: alpha , prod , soma

zeta_tmp = zeta
coef_tmp = coef

prod = zeta_tmp(1) * zeta_tmp(2)
soma = zeta_tmp(1) + zeta_tmp(2)

alpha = ( four*prod )**n * two*sqrt( prod )
alpha = alpha / soma**(two*n + 1)

coef_tmp(i) = - coef_tmp(j)*alpha + sqrt( 1.d0 + coef_tmp(j)*coef_tmp(j)*(alpha*alpha - 1.d0) )

! if coef > 1 go back to original non-optimized STO parameters ...
If( coef_tmp(i) >= 1.d0 ) then
    zeta(1) = me%basis(k)%zeta(1)
    zeta(2) = me%basis(k)%zeta(2)
    coef(1) = me%basis(k)%coef(1)
    coef(2) = me%basis(k)%coef(2)
else
    coef(i) = coef_tmp(i)
end If

end subroutine normalization
!
!
!
!
!===========================================
 subroutine Dump_OPT_parameters( me , dumb )
!===========================================
implicit none
class(CG_OPT)               , intent(in) :: me
integer         , optional  , intent(in) :: dumb

! local variables ...
integer :: i , j , L , AngMax ,n_EHS , N_of_EHSymbol
integer , allocatable   :: indx_EHS(:)

! local parameters ...
character(1)    , parameter :: Lquant(0:3) = ["s","p","d","f"]
integer         , parameter :: DOS   (0:3) = [ 1 , 4 , 9 , 16]

N_of_EHSymbol = size( me % EHSymbol )

allocate( indx_EHS(N_of_EHSymbol) )

! locate position of the first appearance of EHS-atoms in OPT_basis
indx_EHS = [ ( minloc(me % basis % EHSymbol , 1 , me % basis % EHSymbol == me % EHSymbol(i)) , i=1,N_of_EHSymbol ) ] 

!==================================================================================
! creating file CG_OPT_eht_parameters.dat with temporary optimized parameters ...
!==================================================================================

! print heading ...
write(42,48)

do n_EHS = 1 , N_of_EHSymbol

    i = indx_EHS(n_EHS)

    AngMax = atom(me % basis(i) % AtNo) % AngMax

    do L = 0 , AngMax

        j = (i-1) + DOS(L)
    
        write(42,17)    me % basis(j) % Symbol          ,   &
                        me % basis(j) % EHSymbol        ,   &
                        me % basis(j) % AtNo            ,   &
                   atom(me % basis(j) % AtNo) % Nvalen  ,   &
                        me % basis(j) % Nzeta           ,   &
                        me % basis(j) % n               ,   &
                 Lquant(me % basis(j) % l)              ,   &
                        me % basis(j) % IP              ,   &
                        me % basis(j) % zeta(1)*a_Bohr  ,   &      ! <== zetas of OPT_eht_parameters.output.dat are written in units of a0^{-1} ...
                        me % basis(j) % zeta(2)*a_Bohr  ,   &      ! <== zetas of OPT_eht_parameters.output.dat are written in units of a0^{-1} ...
                        me % basis(j) % coef(1)         ,   &
                        me % basis(j) % coef(2)         ,   &
                        me % basis(j) % k_WH
    end do

enddo

rewind(42)

17 format(t1,A2,t13,A3,t25,I3,t33,I3,t45,I3,t53,I3,t60,A3,t68,F9.5,t78,F9.6,t88,F9.6,t98,F9.6,t108,F9.6,t118,F9.6)

include 'formats.h'

end subroutine Dump_OPT_parameters
!
!
!
end module CG_class_m

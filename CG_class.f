module CG_class_m

    use type_m
    use constants_m
    use parameters_m            , only : DP_Moment 
    use Structure_Builder       , only : Extended_Cell 
    use GA_QCModel_m            , only : GA_eigen ,         &
                                         GA_DP_Analysis 
    use cost_tuning_m           , only : evaluate_cost                              

    private

    public :: CG_OPT

    type , extends(OPT) :: CG_OPT
        real*8                        :: InitialCost
        real*8          , allocatable :: p(:)
        type(STO_basis) , allocatable :: basis(:)
    contains
        procedure :: cost
        procedure :: cost_variation
        procedure :: modify_EHT_parameters
        procedure :: normalization
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
real*8      :: CG_DP(3)

If( .NOT. allocated(CG_basis) ) allocate( CG_basis(size(me%basis)) , source = me%basis)

call me % modify_EHT_parameters( )

info = 0
CALL GA_eigen( Extended_Cell , me%basis , CG_UNI , info )

If( DP_Moment ) CALL GA_DP_Analysis( Extended_Cell , me%basis , CG_UNI%L ,CG_UNI%R , CG_DP )

cost = evaluate_cost( CG_UNI , me%basis , CG_DP )

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
end module CG_class_m

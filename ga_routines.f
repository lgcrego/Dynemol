module GA_m

    use type_m
    use constants_m
    use parameters_m            , only : DP_Moment , F_ , T_ , CG_
    use Semi_Empirical_Parms    , only : element => atom 
    use Structure_Builder       , only : Extended_Cell 
    use GA_QCModel_m            , only : GA_eigen ,         &
                                         GA_DP_Analysis ,   &
                                         Mulliken
    use cost_tuning_m           , only : evaluate_cost                                         
    use CG_driver_m             , only : CG_driver


    public :: Genetic_Algorithm , Mulliken

    private 

    integer , parameter :: Pop_Size       =   80         
    integer , parameter :: N_generations  =   100
    integer , parameter :: Top_Selection  =   20           ! <== top selection < Pop_Size
    real*8  , parameter :: Pop_range      =   0.5d0        ! <== range of variation of parameters
    real*8  , parameter :: Mutation_rate  =   0.2           
    logical , parameter :: Mutate_Cross   =   T_           ! <== false -> pure Genetic Algorithm ; prefer false for fine tunning !

    type(OPT) :: GA

contains
!
!
!
!==========================================================
 subroutine modify_EHT_parameters( basis , GA_basis , Pop )
!==========================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
type(STO_basis) , intent(inout) :: GA_basis(:)
real*8          , intent(in)    :: Pop(:) 

! local variables ...
integer :: L , gene , EHS , N_of_EHSymbol 
integer :: indx(size(basis)) , k , i
real*8  :: zeta(2) , coef(2)

! -----------------------------------------------
!       changing basis: editting functions ...
! -----------------------------------------------
!                GA%key storage
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

indx = [ ( i , i=1,size(basis) ) ]

N_of_EHSymbol = size(GA%EHSymbol)

gene = 0

do EHS = 1 , N_of_EHSymbol

!   S , P , D  orbitals ...
do  L = 0 , 2

    If( GA%key(L+1,EHS) == 1 ) then

        ! changes VSIP ...
        gene = gene + GA%key(4,EHS)
        If( GA%key(4,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%IP = Pop(gene)/2.d0 + basis%IP

        ! single STO orbitals ...
        gene = gene + GA%key(5,EHS) - GA%key(6,EHS)
        If( (GA%key(5,EHS) == 1) .AND. (GA%Key(6,EHS) == 0) ) &
        where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%zeta(1) = abs( Pop(gene) + basis%zeta(1) )

        ! double STO orbitals ...
        If( (GA%key(5,EHS) == 1) .AND. (Ga%key(6,EHS) ==1) ) then

            ! finds the first EHT atom ...
            k = minloc( indx , dim=1 , MASK = (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) 

            ! calculate  coef(1)[ zeta(1) , zeta(2) , coef(2) ] ...
            gene    = gene + 1
            zeta(1) = abs( Pop(gene) + basis(k)%zeta(1) )
            gene    = gene + 1
            zeta(2) = abs( Pop(gene) + basis(k)%zeta(2) )
            coef(2) = abs( Pop(gene) )
            CALL normalization( basis , zeta , coef , GA_basis(k)%n , k , 1 , 2 )
            where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) 
                GA_basis % zeta(1) = zeta(1) 
                GA_basis % zeta(2) = zeta(2) 
                GA_basis % coef(1) = coef(1)
                GA_basis % coef(2) = coef(2)
            end where

        End If

        ! changes k_WH ...
        gene = gene + GA%key(7,EHS)
        If( GA%key(7,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%k_WH = Pop(gene) + basis%k_WH

    end If

end do
end do

end subroutine modify_EHT_parameters
!
!
!
!===================================================== 
subroutine Genetic_Algorithm(system, basis, OPT_basis)
!=====================================================	
implicit none
type(structure)                 , intent(in)  :: system
type(STO_basis)                 , intent(in)  :: basis(:)
type(STO_basis) , allocatable   , intent(out) :: OPT_basis(:)

! local variables ...
real*8          , allocatable   :: Pop(:,:) , Old_Pop(:,:) , a(:,:) , semente(:,:) , pot(:,:) , Custo(:) 
real*8                          :: GA_DP(3)
integer         , allocatable   :: indx(:)
integer                         :: i , j , l , generation , info , Pop_start
type(R_eigen)                   :: GA_UNI
type(STO_basis) , allocatable   :: CG_basis(:) , GA_basis(:) , GA_Selection(:,:)

! reading input-GA key ...
CALL Read_GA_key

!-----------------------------------------------
!        SETTING-UP initial populations
!-----------------------------------------------

! Initial Populations ...
allocate( Pop     (Pop_Size , GA%GeneSize) )
allocate( Old_Pop (Pop_Size , GA%GeneSize) )
allocate( indx    (Pop_Size)               )

CALL random_seed

Pop_start = 1
CALL generate_RND_Pop( Pop_start , Pop )       

indx = [ ( i , i=1,Pop_Size ) ]

!-----------------------------------------------

! create new basis ...
allocate( GA_basis (size(basis)) )
GA_basis = basis

allocate( custo(Pop_size) )

do generation = 1 , N_generations

    do i = Pop_start , Pop_Size

        ! intent(in):basis ; intent(inout):GA_basis ...
        CALL modify_EHT_parameters( basis , GA_basis , Pop(i,:) ) 

        info = 0
        CALL  GA_eigen( Extended_Cell , GA_basis , GA_UNI , info )

        If (info /= 0) then 
            custo(i) = 1.d14
            continue
        end if

        If( DP_Moment ) CALL GA_DP_Analysis( Extended_Cell , GA_basis , GA_UNI%L , GA_UNI%R , GA_DP )

!       gather data and evaluate population cost ...
        custo(i) =  evaluate_cost( GA_UNI , GA_basis , GA_DP )

    end do

!   evolve populations ...    
    CALL sort2(custo,indx)

    Old_Pop = Pop
    Pop( 1:Pop_Size , : ) = Old_pop( indx(1:Pop_Size) , : )

    Pop_start = Top_Selection + 1

!   Mutation_&_Crossing preserves the top-selections ...
    If( Mutate_Cross) then
        CALL Mutation_and_Crossing( Pop )
    else
        CALL generate_RND_Pop( Pop_start , Pop )       
    end If

    indx = [ ( i , i=1,Pop_Size ) ]

    Print 160 , generation , N_generations

end do

!----------------------------------------------------------------
! Prepare grid of parameters for CG fine tuning optimization ...
!----------------------------------------------------------------
If( CG_ ) then

    allocate( GA_Selection( size(basis) , Top_Selection ) )

    do i = 1 , Top_Selection 
        ! optimized parameters by GA method : intent(in):basis ; intent(inout):GA_basis ...    
        CALL modify_EHT_parameters( basis , GA_basis , Pop(i,:) )

        GA_Selection(:,i) = GA_basis
    end do

    CALL CG_driver( Extended_Cell , GA , GA_Selection , CG_basis )

    ! create OPT basis ...
    allocate( OPT_basis (size(basis)) )
    OPT_basis = CG_basis

    deallocate( GA_basis , CG_basis , GA_Selection )

else

    ! optimized parameters by GA method : intent(in):basis ; intent(inout):GA_basis ...    
    CALL modify_EHT_parameters( basis , GA_basis , Pop(1,:) )

    ! create OPT basis ...
    allocate( OPT_basis (size(basis)) )
    OPT_basis = GA_basis

    deallocate( GA_basis )

end if

! saving the optimized parameters ...
CALL Dump_OPT_parameters( OPT_basis )

deallocate( GA_UNI%L , GA_UNI%R , GA_UNI%erg )
deallocate( Pop , indx , Old_Pop ) 

include 'formats.h'

end subroutine Genetic_Algorithm
!
!
!
!==============================================
 subroutine generate_RND_Pop( Pop_start , Pop )
!==============================================
implicit none
integer               , intent(in)    :: Pop_start
real*8  , allocatable , intent(inout) :: Pop(:,:) 

! local variables ...
integer               :: i , j
real*8  , allocatable :: a(:,:) , seed(:,:) , pot(:,:) 

!-----------------------------------------------
!           SETTING-UP populations
!-----------------------------------------------

allocate( a    (Pop_Size , GA%GeneSize) )
allocate( seed (Pop_Size , GA%GeneSize) )
allocate( pot  (Pop_Size , GA%GeneSize) )

CALL random_seed
        
do i = Pop_start , Pop_size
    do j = 1 , GA%GeneSize

        CALL random_number( a   (i,j) )
        CALL random_number( seed(i,j) )

        pot(i,j) = int( 2*seed(i,j) )
        Pop(i,j) = ((-1)**pot(i,j)) * a(i,j) * Pop_range

    end do
end do

deallocate( a , seed , pot )

!-----------------------------------------------

end subroutine generate_RND_Pop
!
!
!
!=======================================
 subroutine Mutation_and_Crossing( Pop )
!=======================================
implicit none
real*8  , intent(inout) :: Pop(:,:)

! local variables ...
real*8  , allocatable   :: mat2(:,:) , ma(:) , za(:) , sem(:) , po(:)
real*8                  :: rn , rp
integer , allocatable   :: p(:) , n(:)
integer                 :: i , j , pt , pm

! local parameters ...
integer , parameter     :: N_crossings      = Pop_Size - Top_Selection 
integer , parameter     :: Ponto_cruzamento = N_crossings / 2 

! Random number for crossing ...
allocate( n(N_crossings) )
do i = 1 , N_crossings
    call random_number( rn )
    n(i) = int(Pop_Size*rn) + 1
end do

allocate( p(Ponto_Cruzamento) )
do i = 1 , Ponto_cruzamento
    call random_number( rp )
    p(i) = min( int(GA%GeneSize*rp) + 1 , GA%GeneSize-1 )
end do

! Crossing ...
allocate( mat2(Pop_Size,GA%GeneSize) )
mat2 = Pop
do i = 1 , 2
    do j = 1 , Ponto_Cruzamento
        Pop( Top_Selection+i+(2*j-2) , 1      : p(j)        )  =  mat2( n(i+(2*j-2)) , 1      : p(j)        )
        Pop( Top_Selection+i+(2*j-2) , p(j)+1 : GA%GeneSize )  =  mat2( n((2*j+1)-i) , p(j)+1 : GA%GeneSize )
    end do
end do

! Mutation ...
pt = 2*int( N_crossings * GA%GeneSize * Mutation_rate )
allocate( ma(pt) )

do i = 1 , pt
    call random_number( ma (i) )
end do

pm = pt / 2
allocate( za  (pm) )
allocate( sem (pm) )
allocate( po  (pm) )
do i = 1 , pm
    call random_number( za (i) )
    call random_number( sem(i) )

    po(i) = int( 2 * sem(i) )
    Pop(int(N_crossings*ma(i))+Top_Selection+1,int(GA%GeneSize*ma(i+pm))+1) = ((-1)**po(i)) * za(i)
end do

! deallocating ...
deallocate( n , p , mat2 , ma , za , sem , po )

end subroutine Mutation_and_Crossing
!
!
!
!======================
 subroutine Read_GA_key
!======================
implicit none

! local variables ...
integer :: i , j , ioerr , nr , N_of_EHSymbol
character(1) :: dumb

OPEN(unit=3,file='input-GA.dat',status='old',iostat=ioerr,err=10)
nr = 0
do 
    read(3,*,IOSTAT=ioerr) dumb
    if(ioerr < 0) EXIT
    nr = nr + 1
end do    

N_of_EHSymbol = nr - 1

! allocatting EH_keys ...
allocate( GA%EHSymbol    ( N_of_EHSymbol) )
allocate( GA%key      (7 , N_of_EHSymbol) )

! read the input-GA ...
rewind 3
read(3,*) dumb

Print 40 
Print 41
do j = 1 , N_of_EHSymbol

    read(3,42)  GA%EHSymbol(j) , ( GA%key(i,j) , i=1,7 )

    write(*,42) GA%EHSymbol(j) , ( GA%key(i,j) , i=1,7 )

end do

CLOSE(3)

Print 43

GA%GeneSize = sum( [ ( count(GA%key(1:3,j)==1) * count(GA%key(4:7,j)==1) , j=1,N_of_EHSymbol ) ] )

10 if( ioerr > 0 ) stop "input-GA.dat file not found; terminating execution"

include 'formats.h'

end subroutine Read_GA_key
!
!
!
!===========================================
 subroutine Dump_OPT_parameters( OPT_basis )
!===========================================
implicit none
type(STO_basis) , intent(inout) :: OPT_basis(:)

! local variables ...
integer :: i , j , L , AngMax ,n_EHS , N_of_EHSymbol
integer , allocatable   :: indx_EHS(:)

! local parameters ...
character(1)    , parameter :: Lquant(0:3) = ["s","p","d","f"]
integer         , parameter :: DOS   (0:3) = [ 1 , 4 , 9 , 16]

N_of_EHSymbol = size( GA%EHSymbol )

allocate( indx_EHS(N_of_EHSymbol) )

! locate position of the first appearance of EHS-atoms in OPT_basis
indx_EHS = [ ( minloc(OPT_basis%EHSymbol , 1 , OPT_basis%EHSymbol == GA%EHSymbol(i)) , i=1,N_of_EHSymbol ) ] 

! creating file OPT_eht_parameters.output.dat with the optimized parameters ...
open( unit=13, file='OPT_eht_parameters.output.dat', status='unknown' )

! print heading ...
write(13,48)

do n_EHS = 1 , N_of_EHSymbol

    i = indx_EHS(n_EHS)

    AngMax = element(OPT_basis(i)%AtNo)%AngMax

    do L = 0 , AngMax

        j = (i-1) + DOS(L)
    
        write(13,17)    OPT_basis(j)%Symbol          ,   &
                        OPT_basis(j)%EHSymbol        ,   &
                        OPT_basis(j)%AtNo            ,   &
                element(OPT_basis(j)%AtNo)%Nvalen    ,   &
                        OPT_basis(j)%Nzeta           ,   &
                        OPT_basis(j)%n               ,   &
                 Lquant(OPT_basis(j)%l)              ,   &
                        OPT_basis(j)%IP              ,   &
                        OPT_basis(j)%zeta(1)         ,   &
                        OPT_basis(j)%zeta(2)         ,   &
                        OPT_basis(j)%coef(1)         ,   &
                        OPT_basis(j)%coef(2)         ,   &
                        OPT_basis(j)%k_WH
    end do

enddo
close(13)

17 format(t1,A2,t13,A3,t25,I3,t33,I3,t45,I3,t53,I3,t60,A3,t68,F9.5,t78,F9.6,t88,F9.6,t98,F9.6,t108,F9.6,t118,F9.6)

include 'formats.h'

end subroutine Dump_OPT_parameters
!
!
!
!===============================================================
 subroutine normalization( basis , zeta , coef , n , k , i , j )
!===============================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
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
    zeta(1) = basis(k)%zeta(1)
    zeta(2) = basis(k)%zeta(2)
    coef(1) = basis(k)%coef(1)
    coef(2) = basis(k)%coef(2)
else
    coef(i) = coef_tmp(i)
end If

end subroutine normalization
!
!
!
!=========================
 subroutine  sort2(ra,ira)
!=========================
implicit none
real*8  , intent(inout) :: ra(:)
integer , intent(inout) :: ira(:)

! local variables ...
integer  :: irra, l, n, ir, i, j
real*8   :: rra

!----------------------------------------------------------
!  SORT A(I) , SO THAT THE ELEMENTS IA(I) FOLLOW TOGETHER
!----------------------------------------------------------
      n = size(ra)
      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         rra  = ra(l)
         irra = ira(l)
      else
         rra  = ra(ir)
         irra = ira(ir)
         ra(ir)  = ra(1)
         ira(ir) = ira(1)
         ir = ir - 1
         if(ir .eq. 1) then
             ra(1)  = rra
             ira(1) = irra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(ra(j) .lt. ra(j+1)) j = j + 1
        endif
      if(rra .lt. ra(j)) then
        ra(i)  = ra(j)
        ira(i) = ira(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      ra(i)  = rra
      ira(i) = irra
      goto 10

end subroutine sort2
!
!
!
end module GA_m

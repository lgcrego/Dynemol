module GA_m

    use type_m
    use constants_m
    use MM_types                , only : LogicalKey
    use parameters_m            , only : DP_Moment , F_ , T_ , CG_ ,    &
                                         Pop_size , N_generations ,     &
                                         Top_Selection , Pop_range ,    &
                                         Mutation_rate , Mutate_Cross , &
                                         Alpha_Tensor , OPT_parms ,     &
                                         Adaptive_
    use Semi_Empirical_Parms    , only : atom 
    use Structure_Builder       , only : Extended_Cell 
    use OPT_Parent_class_m      , only : GA_OPT
    use FF_OPT_class_m          , only : FF_OPT
    use EH_CG_driver_m          , only : CG_driver
    use GA_QCModel_m            , only : GA_eigen ,                     &
                                         GA_DP_Analysis ,               &
                                         AlphaPolar ,                   &
                                         Adaptive_GA  
    use cost_EH                 , only : evaluate_cost                                         
    use cost_MM                 , only : SetKeys ,                      &
                                         KeyHolder


    public :: Genetic_Algorithm , Dump_OPT_parameters

    interface Genetic_Algorithm
        module procedure Genetic_Algorithm_EH
        module procedure Genetic_Algorithm_MM
    end interface

    private 

    type(GA_OPT) :: GA

contains
!
!
!
!================================================= 
subroutine Genetic_Algorithm_EH( basis, OPT_basis)
!=================================================	
implicit none
type(STO_basis)                 , intent(inout) :: basis(:)
type(STO_basis) , allocatable   , intent(out)   :: OPT_basis(:)

! local variables ...
real*8          , allocatable   :: Pop(:,:) , Old_Pop(:,:) , cost(:) 
real*8                          :: GA_DP(3) , Alpha_ii(3)
integer         , allocatable   :: indx(:)
integer                         :: i , generation , info , Pop_start , GeneSize
type(R_eigen)                   :: GA_UNI
type(STO_basis) , allocatable   :: CG_basis(:) , GA_basis(:) , GA_Selection(:,:)

! reading input-GA key ...
CALL Read_GA_key( basis )

!-----------------------------------------------
!        SETTING-UP initial populations
!-----------------------------------------------

GeneSize = GA%GeneSize

! Initial Populations ...
allocate( Pop     (Pop_Size , GeneSize) )
allocate( Old_Pop (Pop_Size , GeneSize) )
allocate( indx    (Pop_Size)               )

CALL random_seed

Pop_start = 1
CALL generate_RND_Pop( Pop_start , Pop )       

! this keeps the input EHT parameters in the population ...
Pop(1,:) = D_zero

indx = [ ( i , i=1,Pop_Size ) ]

!-----------------------------------------------

! create new basis ...
allocate( GA_basis (size(basis)) )
GA_basis = basis

allocate( cost(Pop_size) )

! enable on the fly evaluation cost ...
Adaptive_GA%mode = Adaptive_

do generation = 1 , N_generations

    ! sharing these variables with ga_QCModel ...
    Adaptive_GA%gen = generation ; Adaptive_GA%Ngen = N_generations

    !do i = Pop_start , Pop_Size
     do i = 1 , Pop_Size

        ! intent(in):basis ; intent(inout):GA_basis ...
        CALL modify_EHT_parameters( basis , GA_basis , Pop(i,:) ) 

        CALL  GA_eigen( Extended_Cell , GA_basis , GA_UNI )

        If( DP_Moment )    CALL GA_DP_Analysis( Extended_Cell , GA_basis , GA_UNI%L , GA_UNI%R , GA_DP )

        If( Alpha_Tensor ) CALL AlphaPolar( Extended_Cell , GA_basis , Alpha_ii )

!       gather data and evaluate population cost ...
        cost(i) =  evaluate_cost( Extended_Cell , GA_UNI , GA_basis , GA_DP , Alpha_ii )

    end do

!   select the fittest ...    
    CALL sort2(cost,indx)

    Old_Pop = Pop
    Pop( 1:Pop_Size , : ) = Old_pop( indx(1:Pop_Size) , : )

    Pop_start = Pop_size/2 + 1

!   Mutation_&_Crossing preserves the top-selections ...
    If( Mutate_Cross .AND. (mod(generation,4) /= 0) ) then
        CALL Mutation_and_Crossing( Pop )
    else
        CALL generate_RND_Pop( Pop_start , Pop )       
    end If

    indx = [ ( i , i=1,Pop_Size ) ]

    If( Adaptive_ ) then
        Print 159 , generation , N_generations
    else
        Print 160 , generation , N_generations
    EndIf
    Print*, cost(1)
    write(23,*) generation , cost(1)

    ! saving the temporary optimized parameters ...
    ! intent(in):basis ; intent(inout):GA_basis ...
    CALL modify_EHT_parameters( basis , GA_basis , Pop(1,:) ) 
    CALL Dump_OPT_parameters( GA_basis , output = "tmp" )

end do

close(23)

! switch-off on the fly evaluation cost ...
Adaptive_GA%mode = .false.

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

    CALL CG_driver( GA , GA_Selection , CG_basis )

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

end subroutine Genetic_Algorithm_EH
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
integer               :: i , j , GeneSize
real*8  , allocatable :: a(:,:) , seed(:,:) , pot(:,:) 

GeneSize = size(Pop(1,:))

!-----------------------------------------------
!           SETTING-UP populations
!-----------------------------------------------

allocate( a    (Pop_Size , GeneSize) )
allocate( seed (Pop_Size , GeneSize) )
allocate( pot  (Pop_Size , GeneSize) )

CALL random_seed ! <== distribution within the range 0 <= x < 1.
        
do i = Pop_start , Pop_size
    do j = 1 , GeneSize

        CALL random_number( a   (i,j) )
        CALL random_number( seed(i,j) )

        pot(i,j) = int( 2*seed(i,j) )  ! <== bimodal function (-1)^pot = -1 , +1 
        Pop(i,j) = ((-1)**pot(i,j)) * a(i,j) * Pop_range

    end do
end do

! truncate variations to 1.d-5 ...
Pop = Pop * 1.d5 ; Pop = int(Pop) ; Pop = Pop * 1.d-5

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
real*8  , allocatable   :: aux(:,:), a(:), b(:), seed(:), pot(:)
real*8                  :: rn, rp
integer , allocatable   :: p(:), n(:)
integer                 :: i, j, MT, HALF_MT, NX, XP, GeneSize, odd, even, ind1, ind2

GeneSize = size(Pop(1,:))

!N_crossings and XingPoint ...
NX = Pop_Size / 2
XP = NX / 2 

! Start Xing ...
!---------------------------------------------------------------------------
! random population pointer ...
allocate( n(NX) )
do i = 1 , NX
    call random_number( rn ) ! <== distribution within the range 0 <= x < 1.
    n(i) = int(Pop_Size*rn) + 1
end do

! random gene pointer ...
allocate( p(XP) )
do i = 1 , XP
    call random_number( rp )
    p(i) = min( int(GeneSize*rp) + 1 , GeneSize-1 )
end do

allocate( aux(Pop_Size,GeneSize) , source=Pop )

j=0
do odd = 1 , 2*XP-1 , 2

   even = odd + 1
   j=j+1

   Pop( NX + odd ,      1:p(j)     ) = aux( n(odd)  ,       1:p(j)      )
   Pop( NX + odd , p(j)+1:GeneSize ) = aux( n(even) , p(j)+1 : GeneSize )

   Pop( NX + even,      1:p(j)     ) = aux( n(even) ,       1:p(j)      )
   Pop( NX + even, p(j)+1:GeneSize ) = aux( n(odd)  , p(j)+1 : GeneSize )

end do

deallocate( n , p , aux )
!---------------------------------------------------------------------------
! End Xing...

! Start Mutation ...
!---------------------------------------------------------------------------
! integer # of genes to mutate 
MT = int( NX * GeneSize * Mutation_rate )
HALF_MT = MT / 2

allocate( b(MT), a(HALF_MT), seed(HALF_MT), pot(HALF_MT) )

do i = 1 , MT
    call random_number( b(i) ) ! <== distribution within the range 0 <= x < 1.
end do

do i = 1 , HALF_MT

    call random_number( a (i)   )
    call random_number( seed(i) )

    pot(i) = int( two * seed(i) )  ! <== bimodal function (-1)^pot = -1 , +1 
    ind1   = int( NX  * b(i)    )
    ind2   = int( GeneSize * b(HALF_MT+i) )

    Pop( NX+1+ind1 , ind2+1) = (-1)**pot(i) * a(i) * Pop_range

end do

! truncate variations to 1.d-5 ...
Pop = Pop * 1.d5 ; Pop = int(Pop) ; Pop = Pop * 1.d-5

deallocate( b , a , seed , pot )
!---------------------------------------------------------------------------
! End Mutation ...

end subroutine Mutation_and_Crossing
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
        If( GA%key(4,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%IP = Pop(gene) + basis%IP

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
!===============================
 subroutine Read_GA_key( basis )
!===============================
implicit none
type(STO_basis) , intent(inout) :: basis(:)

! local variables ...
integer :: i , j , ioerr , err , n , N_of_EHSymbol
character(1) :: dumb

OPEN(unit=3,file='input-GA.dat',status='old',iostat=ioerr,err=10)
n = 0
do 
    read(3,*,IOSTAT=ioerr) dumb
    if(ioerr < 0) EXIT
    n = n + 1
end do    

N_of_EHSymbol = n - 1

! allocatting EH_keys: [s,p,d,IP,zeta,coef,k_WH] ...
allocate( GA%EHSymbol    ( N_of_EHSymbol) )
allocate( GA%key      (7 , N_of_EHSymbol) )

! read the input-GA ...
rewind 3
read(3,*) dumb

Print 40 
Print 41
do j = 1 , N_of_EHSymbol

    read(3,42)   GA%EHSymbol(j) , ( GA%key(i,j) , i=1,7 )

    write(*,421) GA%EHSymbol(j) , ( GA%key(i,j) , i=1,7 )

end do

CLOSE(3)

Print 43

GA%GeneSize = sum( [ ( count(GA%key(1:3,j)==1) * count(GA%key(4:7,j)==1) , j=1,N_of_EHSymbol ) ] )

do j = 1 , N_of_EHSymbol

    If( GA%key(1,j) /= 0 ) &   ! <== optimizing s orbital ...
    where( adjustl(basis% EHSymbol) == adjustl(GA% EHSymbol(j)) .AND. basis%L == 0 ) basis%Nzeta = max( GA% key(5,j)+GA% key(6,j) , basis%Nzeta )

    If( GA%key(2,j) /= 0 ) &   ! <== optimizing p orbital ...
    where( adjustl(basis% EHSymbol) == adjustl(GA% EHSymbol(j)) .AND. basis%L == 1 ) basis%Nzeta = max( GA% key(5,j)+GA% key(6,j) , basis%Nzeta )

    If( GA%key(3,j) /= 0 ) &   ! <== optimizing d orbital ...
    where( adjustl(basis% EHSymbol) == adjustl(GA% EHSymbol(j)) .AND. basis%L == 2 ) basis%Nzeta = max( GA% key(5,j)+GA% key(6,j) , basis%Nzeta )

end do

If( OPT_parms ) then
     Print*, ">> OPT_parms being used as input <<"
else
     Print*, ">> OPT_parms were not used <<"
end if

10 if( ioerr > 0 ) stop "input-GA.dat file not found; terminating execution"

include 'formats.h'

end subroutine Read_GA_key
!
!
!
!====================================================
 subroutine Dump_OPT_parameters( OPT_basis , output )
!====================================================
implicit none
type(STO_basis)            , intent(inout) :: OPT_basis(:)
character(len=*), optional , intent(in)    :: output

! local variables ...
integer               :: i , j , L , AngMax ,n_EHS , N_of_EHSymbol
integer , allocatable :: indx_EHS(:)
integer               :: unit_tag

! local parameters ...
character(1)    , parameter :: Lquant(0:3) = ["s","p","d","f"]
integer         , parameter :: DOS   (0:3) = [ 1 , 4 , 9 , 16]

N_of_EHSymbol = size( GA%EHSymbol )

allocate( indx_EHS(N_of_EHSymbol) )

! locate position of the first appearance of EHS-atoms in OPT_basis
indx_EHS = [ ( minloc(OPT_basis%EHSymbol , 1 , OPT_basis%EHSymbol == GA%EHSymbol(i)) , i=1,N_of_EHSymbol ) ] 

If( present(output) .AND. output=="STDOUT" ) then

    Print*,""
    Print*,""
    unit_tag = 6

elseIf( present(output) .AND. output=="tmp" ) then

    ! creating file opt_eht_parms.output with the optimized parameters ...
    open( unit=13, file='opt.trunk/opt_eht_parms.output', status='unknown' )
    unit_tag = 13

else

    ! creating file opt_eht_parms.output with the optimized parameters ...
    open( unit=13, file='opt_eht_parms.output', status='unknown' )
    unit_tag = 13

end If

! print heading ...
write(unit_tag,48)

do n_EHS = 1 , N_of_EHSymbol

    i = indx_EHS(n_EHS)

    AngMax = atom(OPT_basis(i)%AtNo)%AngMax

    do L = 0 , AngMax

        j = (i-1) + DOS(L)
    
  write(unit_tag,17)    OPT_basis(j)%Symbol          ,   &
                        OPT_basis(j)%EHSymbol        ,   &
                        OPT_basis(j)%residue         ,   &
                        OPT_basis(j)%AtNo            ,   &
                   atom(OPT_basis(j)%AtNo)%Nvalen    ,   &
                        OPT_basis(j)%Nzeta           ,   &
                        OPT_basis(j)%n               ,   &
                 Lquant(OPT_basis(j)%l)              ,   &
                        OPT_basis(j)%IP              ,   &
                        OPT_basis(j)%zeta(1)*a_Bohr  ,   &      ! <== zetas of opt_eht_parms.output are written in units of a0^{-1} ...
                        OPT_basis(j)%zeta(2)*a_Bohr  ,   &      ! <== zetas of opt_eht_parms.output are written in units of a0^{-1} ...
                        OPT_basis(j)%coef(1)         ,   &
                        OPT_basis(j)%coef(2)         ,   &
                        OPT_basis(j)%k_WH
    end do

enddo
If( unit_tag == 13 ) close(13)

17 format(t1,A2,t13,A3,t26,A3,t36,I3,t45,I3,t57,I3,t65,I3,t72,A3,t80,F9.5,t90,F9.6,t100,F9.6,t110,F9.6,t120,F9.6,t130,F9.6)

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
!
!=========================================================
!
!
!
!
!=======================================================================
 subroutine Genetic_Algorithm_MM( MM_parms , GA_Selection , directives )
!=======================================================================	
implicit none
type(FF_OPT)                  , intent(inout) :: MM_parms
real*8          , allocatable , intent(out)   :: GA_Selection(:,:)
character(*)                  , intent(in)    :: directives

! local variables ...
real*8          , allocatable   :: Pop(:,:) , Old_Pop(:,:) , cost(:) , p0(:)
integer         , allocatable   :: indx(:)
integer                         :: i , generation , Pop_start ,  GeneSize
type(LogicalKey)                :: key

!-----------------------------------------------
!        SETTING-UP MM_parms object
!-----------------------------------------------

CALL SetKeys

key = KeyHolder( size(KeyHolder) )

MM_parms = FF_OPT( key , kernel = "NormalModes" , directives = directives )

GeneSize = MM_parms % N_of_freedom 

! creating reference copy of FF vector ...
allocate( p0(GeneSize) , source = MM_parms % p )

!-----------------------------------------------
!        SETTING-UP initial populations
!-----------------------------------------------

! Initial Populations ...
allocate( Pop     (Pop_Size , GeneSize) )
allocate( Old_Pop (Pop_Size , GeneSize) )
allocate( indx    (Pop_Size)            )

CALL random_seed

Pop_start = 1
CALL generate_RND_Pop( Pop_start , Pop )       

! the input parameters constitute one of the genes of Pop ...
Pop(1,:) = D_zero

indx = [ ( i , i=1,Pop_Size ) ]

!-----------------------------------------------
! setting normal mode data in FF_OPT_class ...

allocate( cost(Pop_size) )

do generation = 1 , N_generations

    do i = Pop_start , Pop_Size

        ! virtual displacements in the FF parameter space ...
        MM_parms % p(:) = p0(:) * (D_one + Pop(i,:))

        ! evaluate  Pop(i,:)'s  cost ...
        cost(i) = MM_parms % cost()

    end do

!   evolve populations ...    
    CALL sort2(cost,indx)

    Old_Pop = Pop
    Pop( 1:Pop_Size , : ) = Old_pop( indx(1:Pop_Size) , : )

    Pop_start = Top_Selection + 1

!   Mutation_&_Crossing preserves the top-selections ...
    If( Mutate_Cross) then

        CALL Mutation_and_Crossing( Pop )

    else

        CALL generate_RND_Pop( Pop_start , Pop )       

        If( generation < N_generations) forall(i = Pop_Start:Pop_Size) Pop(i,:) = Pop(1,:) - Pop(i,:)

    end If

    indx = [ ( i , i=1,Pop_Size ) ]

    Print 161 , generation , N_generations , cost(1) , directives 

end do

allocate( GA_Selection( size(p0) , Top_Selection ) )

forall( i=1:Top_Selection ) GA_Selection(:,i) = p0(:) * (D_one + Pop(i,:))

MM_parms%p = p0(:) * (D_one + Pop(1,:)) 

deallocate( Pop , indx , Old_Pop ) 

include 'formats.h'

end subroutine Genetic_Algorithm_MM
!
!
!
end module GA_m

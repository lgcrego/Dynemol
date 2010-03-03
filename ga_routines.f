module GA_m

    use type_m
    use constants_m
    use Semi_Empirical_Parms    , only : atom  
    use Structure_Builder       , only : Extended_Cell , Basis_Builder
    use QCModel_Huckel          , only : EigenSystem
    use Multipole_Core          , only : Dipole_Matrix

    public :: Genetic_Algorithm

    private 

    integer , parameter :: Pop_Size       =   40         
    integer , parameter :: N_generations  =   80         
    integer , parameter :: Top_Selection  =   15        
    real*8  , parameter :: Pop_range      =   1.0d-1        ! <== range of variation of parameters
    real*8  , parameter :: Mutation_rate  =   0.3           
    logical , parameter :: Mutate_Cross   =  .false.        ! <== false -> pure Genetic Algorithm ; prefer false for fine tunning !

    type(OPT) :: GA

contains

!===============================================
 function evaluate_cost( system , REF )
!===============================================
implicit none
type(structure) , intent(in)  :: system
type(OPT)       , intent(in)  :: REF
real*8                        :: evaluate_cost

! local variables ...
real*8   :: chi(10) , weight(10)
integer  :: k , HOMO , LUMO

! general definitions ...
HOMO = system%N_of_electrons / 2
LUMO = HOMO + 1

chi(:) = 0.d0   ;   weight(:) = 1.d0

! HOMO-LUMO gaps ...
chi(1) = ( GA%erg(9) - GA%erg(8) ) - 8.0    ; weight(1) = 1.0

chi(2) = ( GA%erg(9) - GA%erg(7) ) - 8.0    ; weight(1) = 1.0 

chi(3) = ( GA%erg(8) - GA%erg(6) ) - 0.1    ; weight(1) = 5.0 

chi(4) = ( GA%erg(7) - GA%erg(6) ) - 0.1    ; weight(1) = 5.0 

! Total DIPOLE moment ...
chi(5) = dot_product( GA%DP , GA%DP ) - dot_product( REF%DP , REF%DP )   

evaluate_cost = sqrt( dot_product(chi,chi) )

end function evaluate_cost
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

N_of_EHSymbol = size(GA%EHSymbol)

gene = 0

do EHS = 1 , N_of_EHSymbol

!   S , P , D  orbitals ...
do  L = 0 , 2

    If( GA%key(L+1,EHS) == 1 ) then

        gene = gene + GA%key(4,EHS)
        If( GA%key(4,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%IP = Pop(gene) + basis%IP

        gene = gene + GA%key(5,EHS)
        If( GA%key(5,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%zeta(1) = Pop(gene) + basis%zeta(1)

        gene = gene + GA%key(6,EHS)
        If( GA%key(6,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%zeta(2) = Pop(gene) + basis%zeta(2)

        gene = gene + GA%key(7,EHS)
        If( GA%key(7,EHS) == 1 ) where( (GA_basis%EHSymbol == GA%EHSymbol(EHS)) .AND. (GA_basis%l == L) ) GA_basis%k_WH = Pop(gene) + basis%k_WH

    end If
end do
end do

end subroutine modify_EHT_parameters
!
!
!
!========================================================= 
subroutine Genetic_Algorithm(system, basis, REF, GA_basis)
!========================================================= 	
implicit none
type(structure)                 , intent(in)  :: system
type(STO_basis)                 , intent(in)  :: basis(:)
type(OPT)                       , intent(in)  :: REF
type(STO_basis) , allocatable   , intent(out) :: GA_basis(:)

! local variables ...
real*8          , allocatable   :: Pop(:,:) , Old_Pop(:,:) , a(:,:) , semente(:,:) , pot(:,:) , Custo(:) 
real*8                          :: GA_DP(3)
integer         , allocatable   :: indx(:)
integer                         :: i , j , l , generation , info 
type(eigen)                     :: GA_UNI

! reading input-GA key ...
CALL Read_GA_key

!-----------------------------------------------
!           SETTING-UP populations
!-----------------------------------------------

! Initial Populations ...
allocate( Pop     (Pop_Size , GA%GeneSize) )
allocate( Old_Pop (Pop_Size , GA%GeneSize) )
allocate( a       (Pop_Size , GA%GeneSize) )
allocate( semente (Pop_Size , GA%GeneSize) )
allocate( pot     (Pop_Size , GA%GeneSize) )
allocate( indx    (Pop_Size)               )

CALL random_seed
        
do i = 1 , Pop_size
    do j = 1 , GA%GeneSize
        CALL random_number( a(i,j) )
        CALL random_number( semente(i,j) )
        pot(i,j) = int(2 * semente(i,j))
        Pop(i,j) = ((-1)**pot(i,j)) * a(i,j) * Pop_range
    end do
end do
indx = [ ( i , i=1,Pop_Size ) ]

deallocate( a , semente , pot )

!-----------------------------------------------

! create new basis ...
allocate( GA_basis  (size(basis)) )
allocate( GA%erg    (size(basis)) )
GA_basis = basis

allocate( custo(Pop_size) )

do generation = 1 , N_generations

    do i = 1 , Pop_Size

        CALL modify_EHT_parameters( basis , GA_basis , Pop(i,:) )

        info = 0
        CALL EigenSystem( Extended_Cell, GA_basis, GA_UNI , info )
                
        If (info /= 0) then 
            custo(i) = 1.d14
            continue
        end if

        CALL Dipole_Matrix( Extended_Cell, GA_basis, GA_UNI%L, GA_UNI%R, GA_DP )

        GA%erg  =  GA_UNI%erg
        GA%DP   =  GA_DP
                
!       evaluate population cost ...
        custo(i) = evaluate_cost( system , REF )

    end do

!   evolve populations ...    
    CALL sort2(custo,indx)
        
    Old_Pop = Pop
    Pop( 1:Top_Selection , : ) = Old_pop( indx(1:Top_Selection) , : )
        
    If( Mutate_Cross) CALL Mutation_and_Crossing( Pop )

    indx = [ ( i , i=1,Pop_Size ) ]

    Print 160 , generation , N_generations

end do

!-----------------------------------------------

! optimized Huckel parameters are ...
CALL modify_EHT_parameters( basis , GA_basis , Pop(1,:) )

! saving the optimized parameters ...
CALL Dump_OPT_parameters( GA_basis )

deallocate( GA%erg , GA_UNI%L , GA_UNI%R , GA_UNI%erg)
deallocate( Pop , indx , Old_Pop ) 

include 'formats.h'

end subroutine Genetic_Algorithm
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
    p(i) = int(GA%GeneSize*rp) + 1
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
pt = N_crossings * GA%GeneSize * Mutation_rate 
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

do j = 1 , N_of_EHSymbol

    read(3,11) GA%EHSymbol(j) , ( GA%key(i,j) , i=1,7 )

end do

CLOSE(3)

GA%GeneSize = sum( [ ( count(GA%key(1:3,j)==1) * count(GA%key(4:7,j)==1) , j=1,N_of_EHSymbol ) ] )

10 if( ioerr > 0 ) stop "input-GA.dat file not found; terminating execution"

11 FORMAT(A3,t17,I1,t25,I1,t33,I1,t41,I1,t49,I1,t61,I1,t73,I1)

end subroutine Read_GA_key
!
!
!
!==========================================
 subroutine Dump_OPT_parameters( GA_basis )
!==========================================
implicit none
type(STO_basis) , intent(inout) :: GA_basis(:)

! local variables ...
integer :: i , j , L , AngMax ,n_EHS , N_of_EHSymbol
integer , allocatable   :: indx_EHS(:)

! local parameters ...
character(1)    , parameter :: Lquant(0:3) = ["s","p","d","f"]

N_of_EHSymbol = size( GA%EHSymbol )

allocate( indx_EHS(N_of_EHSymbol) )

indx_EHS = [ ( minloc(GA_basis%EHSymbol , 1 , GA_basis%EHSymbol == GA%EHSymbol(i)),i=1,N_of_EHSymbol ) ] 

! creating file OPT_eht_parameters.output.dat with the optimized parameters ...
open( unit=13, file='OPT_eht_parameters.output.dat', status='unknown' )

do n_EHS = 1 , N_of_EHSymbol

    i = indx_EHS(n_EHS)

    AngMax = atom(GA_basis(i)%AtNo)%AngMax 

    do L = 0 , AngMax
        j = i+L
        write(13,17)    GA_basis(j)%Symbol          ,   &
                        GA_basis(j)%EHSymbol        ,   &
                        GA_basis(j)%AtNo            ,   &
                   atom(GA_basis(j)%AtNo)%Nvalen    ,   &
                        GA_basis(j)%Nzeta           ,   &
                        GA_basis(j)%n               ,   &
                 Lquant(GA_basis(j)%l)              ,   &
                        GA_basis(j)%IP              ,   &
                        GA_basis(j)%zeta(1)         ,   &
                        GA_basis(j)%zeta(2)         ,   &
                        GA_basis(j)%coef(1)         ,   &
                        GA_basis(j)%coef(2)         ,   &
                        GA_basis(j)%k_WH
    end do

enddo
close(13)

17 format(A3,t6,A3,t10,I3,t17,I3,t24,I3,t30,I3,t36,A3,t43,F8.3,t52,F8.4,t61,F8.4,t70,F8.4,t79,F8.4,t88,F8.4)

end subroutine Dump_OPT_parameters
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
end module GA_m

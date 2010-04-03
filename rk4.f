module RK_m

    use type_m
    use constants_m
    use QOptics_m

    public :: RK4_dynamics

    private

 contains
!
!
!
!----------------------------------------------------
 subroutine RK4_dynamics( system , basis , QM , FMO )
!----------------------------------------------------
 type(structure) , intent(in) :: system
 type(STO_basis) , intent(in) :: basis(:)
 type(C_eigen)   , intent(in) :: QM
 type(C_eigen)   , intent(in) :: FMO

! . local variables ...
 integer                    :: nstp , n_OK , n_BAD , dim_f
 complex*16 , allocatable   :: f(:,:) , dfdt(:,:) , fscal(:,:)

 real*8                     :: t , dt , tdid , step 
 real*8     , parameter     :: eps = 1.d-3


 Print 165, 'Starting Runge-Kutta' 

 OPEN(unit=27,file='RK.log',status='new') 
 CLOSE(27)

!-------------------------------------------------------- 
! . transition dipole matrices ... 

 CALL RK_setup( system , basis , QM , FMO , f , dim_f )

!=========================================================
!                   DYNAMIC DRIVER
!=========================================================

 t = t_i

 dt = (t_f - t_i) / float(n_t)

 n_OK  = 0
 n_BAD = 0

 allocate( dfdt (dim_f,dim_f) )
 allocate( fscal(dim_f,dim_f) )

 do nstp = 1 , 10000

    CALL evolution( t , f , dfdt )

    fscal = cdabs(f) + cdabs(dt*dfdt)

    CALL rkqs(f , dfdt , t , dt , eps , fscal , tdid , step )

    if( tdid == dt ) then
        n_OK = n_OK + 1
    else
        n_BAD = n_BAD + 1
    end if

    CALL RK_log( f , t , nstp , dt , tdid , step )

print*, t , real( sum( (/( f(i,i) , i=1,10 )/) ) ) , real( sum( (/( f(i,i) , i=1,100 )/) ) )  
write(11,*) t , real( sum( (/( f(i,i) , i=1,10 )/) ) ) , real( sum( (/( f(i,i) , i=1,100 )/) ) )  

    if( t >= t_f ) return

    dt = step
    t  = t + tdid

 end do    

 pause 'too many steps in RK'

 include 'formats.h'

 return

 end subroutine RK4_dynamics
!
!
!
!-------------------------------------------------------------------
 SUBROUTINE rkqs( y , dydx , x , htry , eps , yscal , hdid , hnext )
!-------------------------------------------------------------------
    IMPLICIT NONE
    complex*16 , DIMENSION(:,:) , INTENT(INOUT) :: y
    complex*16 , DIMENSION(:,:) , INTENT(IN)    :: dydx,yscal
	REAL*8                      , INTENT(INOUT) :: x
	REAL*8                      , INTENT(IN)    :: htry,eps
	REAL*8                      , INTENT(OUT)   :: hdid,hnext

! . local variables	...
	REAL*8                                                  :: errmax , h , htemp , xnew
    complex*16 , DIMENSION( size(y(:,1)) , size(y(1,:)) )   :: yerr , ytemp

! . local parameters ...
	real*8      , parameter     :: pgrow    =  -2.d-1 
	real*8      , parameter     :: pshrnk   =  -2.5d-1 
	real*8      , parameter     :: safety   =  9.d-1 
	real*8      , parameter     :: errcon   =  6.d-4  !1.89d-4
    real*8      , parameter     :: tenth    =  0.1d0

    h=htry
    do
        call rkck( y , dydx , x , h , ytemp , yerr )

        errmax = maxval( cdabs ( yerr(:,:)/yscal(:,:) ) ) / eps

        if (errmax <= 1.0) exit

        htemp=SAFETY*h*(errmax**PSHRNK)

        h = sign( max( dabs(htemp) , tenth*dabs(h) ) , h )

        xnew = x + h

        if (xnew == x) pause 'stepsize underflow in rkqs'
    end do

    if (errmax > ERRCON) then
        hnext=SAFETY*h*(errmax**PGROW)
    else
        hnext=five * h
    end if

    hdid=h
    x=x+h

    y(:,:)=ytemp(:,:)

END SUBROUTINE rkqs
!
!
!
!------------------------------------------------
SUBROUTINE rkck( y , dydx , x , h , yout , yerr )
!------------------------------------------------
    IMPLICIT NONE
    complex*16  , INTENT(IN)    :: y(:,:) , dydx(:,:)
    real*8      , INTENT(IN)    :: x , h
    complex*16  , INTENT(OUT)   :: yout(:,:) , yerr(:,:)

! . local variables ...
    complex*16  , allocatable   :: ytemp(:,:) , ak2(:,:) , ak3(:,:) , ak4(:,:) , ak5(:,:) , ak6(:,:)
    integer                     :: dim_y

! . local parameters ...
    real*8      , PARAMETER     :: A2=0.2d0 , A3=0.3d0 , A4=0.6d0 , A5=1.0d0 , A6=0.875d0 , B21=0.2d0 , B31=3.0d0/40.0d0 , B32=9.0d0/40.0d0     , &
                                   B41=0.3d0 , B42=-0.9d0 , B43=1.2d0 , B51=-11.0d0/54.0d0 , B52=2.5d0 , B53=-70.0d0/27.0d0 , B54=35.0d0/27.0d0 , &
                                   B61=1631.0d0/55296.0d0 , B62=175.0d0/512.0d0 , B63=575.0d0/13824.0d0 , B64=44275.0d0/110592.0d0              , &
                                   B65=253.0d0/4096.0d0 , C1=37.0d0/378.0d0 , C3=250.0d0/621.0d0 , C4=125.0d0/594.0d0 , C6=512.0d0/1771.0d0     , &
                                   DC1=C1-2825.0d0/27648.0d0 , DC3=C3-18575.0d0/48384.0d0 , DC4=C4-13525.0d0/55296.0d0 , DC5=-277.0d0/14336.0d0 , &
                                   DC6=C6-0.25d0

 dim_y = size(y(:,1))
 CALL allocation( ytemp , ak2 , ak3 , ak4 , ak5 , ak6 , dim_y )

    ytemp = y + B21*h*dydx

    CALL evolution( x+A2*h , ytemp , ak2 )

    ytemp = y + h*(B31*dydx+B32*ak2)

    CALL evolution( x+A3*h , ytemp , ak3 )

    ytemp = y + h*(B41*dydx+B42*ak2+B43*ak3)

    CALL evolution( x+A4*h , ytemp , ak4 )

    ytemp = y + h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)

    CALL evolution( x+A5*h , ytemp , ak5 )

    ytemp = y + h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)

    CALL evolution( x+A6*h , ytemp , ak6 )

    yout = y + h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)

    yerr = h * (DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

 deallocate( ytemp , ak2 , ak3 , ak4 , ak5 , ak6 )

END SUBROUTINE rkck
!
!
!
!
!---------------------------------------------------------------------
 subroutine Allocation( ytemp , ak2 , ak3 , ak4 , ak5 , ak6 , dim_y )
!---------------------------------------------------------------------
complex*16 , allocatable , intent(out) :: ytemp(:,:)
complex*16 , allocatable , intent(out) :: ak2(:,:)
complex*16 , allocatable , intent(out) :: ak3(:,:)
complex*16 , allocatable , intent(out) :: ak4(:,:)
complex*16 , allocatable , intent(out) :: ak5(:,:)
complex*16 , allocatable , intent(out) :: ak6(:,:)
integer                                :: dim_y

allocate( ytemp (dim_y,dim_y) ) 
allocate( ak2   (dim_y,dim_y) ) 
allocate( ak3   (dim_y,dim_y) ) 
allocate( ak4   (dim_y,dim_y) ) 
allocate( ak5   (dim_y,dim_y) ) 
allocate( ak6   (dim_y,dim_y) ) 

end subroutine Allocation
!
!
!
!---------------------------------------------------
subroutine RK_log( f , t , nstp , dt , tdid , step )
!---------------------------------------------------
complex*16  , intent(in) :: f(:,:)
real*8      , intent(in) :: t
integer     , intent(in) :: nstp
real*8      , intent(in) :: dt
real*8      , intent(in) :: tdid
real*8      , intent(in) :: step

! . local variable ...
integer :: dim_f

OPEN(unit=27,file='RK.log',status='old',access='append') 

dim_f = size(f(:,1))

if( tdid == dt ) then
    write(27,180) nstp , t , 'OK' 
else
    write(27,180) nstp , t , 'BAD'
end if

write(27,182) dt
write(27,183) tdid
write(27,184) step
write(27,185) 1.d0 - real( sum( (/( f(i,i) , i=1,dim_f )/) ) )

close(27)

include 'formats.h'

end subroutine RK_log


end module RK_m

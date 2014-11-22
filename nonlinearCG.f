module NonlinearCG_m

        use constants_m             , only: real_large , prec => mid_prec 
        use CG_class_m              , only: CG_OPT

        implicit none

        public :: Fletcher_Reeves_Polak_Ribiere_minimization

        private

        ! module parameters ...
            integer   , parameter :: ITMAX = 3!50     ! <== 100-300 is a good compromise of accuracy and safety

        ! module types ...
        type f1com
            integer               :: ncom
            real*8  , allocatable :: xicom(:)
            real*8  , allocatable :: pcom(:)
        end type f1com 

        ! module variables ...
        integer     :: NMAX  
        real*8      :: BracketSize = 1.d-4      ! <== this value may vary between 1.0d-3 and 1.0d-5
        type(f1com) :: f1

contains      

!=================================================================================
 subroutine Fletcher_Reeves_Polak_Ribiere_minimization( this , n , local_minimum ) 
!=================================================================================
type(CG_OPT) , intent(inout) :: this
integer      , intent(in)    :: n
real*8       , intent(out)   :: local_minimum

!local variables ...
INTEGER :: its,j,iter
real*8  :: dgg,fp,gam,gg
real*8  :: g(N),h(N),xi(N)

! Initializations ...

        NMAX = n
        allocate( f1%xicom(NMAX) , f1%pcom(NMAX) )

        fp = this % cost()
        call this % cost_variation(xi)

        do j=1,n
           g(j)  = -xi(j)
           h(j)  = g(j)
           xi(j) = h(j)
        end do 

        do its=1,ITMAX                                                                          ! Loop over iterations.
           iter=its

           call Linear_Minimization( this , xi , n , local_minimum )                            ! Next statement is the normal return:

           If( this % profiling ) then
               Print*, its , local_minimum
               write(32,*) its , local_minimum 
               if (this % driver == "MM_Optimize" ) call this%output( iter )
           end If

           if( local_minimum > fp ) then
               local_minimum = real_large
               goto 100
           end if

           if( (2.d0*abs(local_minimum-fp)) <= prec*(abs(local_minimum)+abs(fp)+prec) )  goto 100

           fp=local_minimum
           call this % cost_variation(xi)

           gg=0.0d0
           dgg=0.0d0
           do j=1,n
              gg=gg+g(j)**2
!              dgg=dgg+xi(j)**2                ! This statement for Fletcher-Reeves.
              dgg=dgg+(xi(j)+g(j))*xi(j)       ! This statement for Polak-Ribiere.
           end do 
           if(gg.eq.0.)return                                                                   ! Unlikely. If gradient is exactly zero then we are al-
           gam=dgg/gg                                                                           ! ready done.
           do j=1,n
              g(j)=-xi(j)
              h(j)=g(j)+gam*h(j)
              xi(j)=h(j)
           end do 
        enddo 

100     deallocate( f1%xicom , f1%pcom )
       
end subroutine Fletcher_Reeves_Polak_Ribiere_minimization
!        
!
!
!======================================================
subroutine Linear_Minimization(this,xi,n,local_minimum)
!======================================================
type(CG_OPT) , intent(inout) :: this
real*8       , intent(inout) :: xi(:)
integer      , intent(in)    :: n
real*8       , intent(out)   :: local_minimum

! local variables ...
INTEGER     :: j
real*8      :: ax,bx,fa,fb,fx,xmin,xx,p(n)

 p = this%p

 f1%ncom=n  
 do j=1,f1%ncom
    f1%pcom(j)  = p(j)
    f1%xicom(j) = xi(j)
 end do

! Initial guess for brackets.
 ax=0.0d0                                            
 xx=BracketSize

 ! call mnbrak(this,ax,xx,bx,fa,fx,fb)                  ! <== this is not stable for the present optimization case

 local_minimum = brent(this,ax,xx,bx,xmin)

 If( this % driver == "Genetic_Alg" .AND. this % profiling ) CALL this % output

 do j=1,n                                               ! Construct the vector results to return.
    xi(j)=xmin*xi(j)
    p(j)=p(j)+xi(j)
 end do 

 this%p = p

end subroutine Linear_Minimization
!
!
!
!========================================
SUBROUTINE mnbrak(this,ax,bx,cx,fa,fb,fc)
!========================================
type(CG_OPT) , intent(inout) :: this
real*8       , intent(inout) :: ax , bx , cx
real*8       , intent(inout) :: fa , fb , fc

!local parameters ...
real*8  , PARAMETER :: GOLD=1.618034d0 , GLIMIT=2000.d0 , WEE=1.d-20

! local variables ...
real*8 :: dum,fu,q,r,u,ulim

        fa=f1dim(this,ax)
        fb=f1dim(this,bx)
        if(fb.gt.fa)then                                                    ! Switch roles of a and b so that we can go downhill in the
           dum=ax                                                           ! direction from a to b.
           ax=bx
           bx=dum
           dum=fb
           fb=fa
           fa=dum
        end if
        cx=bx+GOLD*(bx-ax)                                                  ! First guess for c.
        fc=f1dim(this,cx)
 1      if(fb.ge.fc)then                                                    ! “do while”: keep returning here until we bracket.
        r=(bx-ax)*(fb-fc)                                                   ! Compute u by parabolic extrapolation from a, b, c. WEE
        q=(bx-cx)*(fb-fa)                                                   ! is used to prevent any possible division by zero.
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),WEE),q-r))
        ulim=bx+GLIMIT*(cx-bx)                                              ! We won’t go farther than this. Test various possibilities:
        if((bx-u)*(u-cx).gt.0.d0)then                                       ! Parabolic u is between b and c: try it.
           fu=f1dim(this,u)
           if(fu.lt.fc)then                                                 ! Got a minimum between b and c.
              ax=bx
              fa=fb
              bx=u
              fb=fu
              return
           else if(fu.gt.fb)then                                            ! Got a minimum between between a and u.
              cx=u
              fc=fu
              return
           end if
           u=cx+GOLD*(cx-bx)                                                ! Parabolic fit was no use. Use default magnification.
           fu=f1dim(this,u)
        else if((cx-u)*(u-ulim).gt.0.d0)then                                ! Parabolic fit is between c and its allowed
           fu=f1dim(this,u)                                                 ! limit.
           if(fu.lt.fc)then
              bx=cx
              cx=u
              u=cx+GOLD*(cx-bx)
              fb=fc
              fc=fu
              fu=f1dim(this,u)
           end if
        else if((u-ulim)*(ulim-cx).ge.0.d0)then                             ! Limit parabolic u to maximum allowed
           u=ulim                                                           ! value.
           fu=f1dim(this,u)
        else                                                                ! Reject parabolic u, use default magnification.
           u=cx+GOLD*(cx-bx)
           fu=f1dim(this,u)
        endif
        ax=bx                                                               ! Eliminate oldest point and continue.
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
        endif
        
end subroutine mnbrak
! 
!
!
!======================
 function f1dim(this,x)
!======================
type(CG_OPT) , intent(inout) :: this
real*8       , intent(in)    :: x
real*8                       :: f1dim

! local variables ...
INTEGER     :: j
real*8      :: xt(NMAX)

        do j = 1 , f1%ncom
           xt(j) = f1%pcom(j)+x*f1%xicom(j)
        end do 

        this%p(:) = xt(:)

        f1dim = this % cost()
      
end function f1dim
!
!
!
!=================================
function brent(this,ax,bx,cx,xmin)
!=================================
type(CG_OPT) , intent(inout) :: this
real*8       , intent(in)    :: ax , bx , cx
real*8       , intent(out)   :: xmin
real*8                       :: brent

! local parameters ...  
real*8 , PARAMETER :: CGOLD=0.3819660d0

! local variables ...
INTEGER     :: iter
real*8      :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

        a=min(ax,cx)                    ! a and b must be in ascending order, though the input
        b=max(ax,cx)                    ! abscissas need not be.
        v=bx                            ! Initializations...
        w=v
        x=v
        e=0.0d0                         ! This will be the distance moved on the step before last.
        fx=f1dim(this,x)
        fv=fx
        fw=fx
        do iter=1,ITMAX                 !Main program loop.
           xm=0.5*(a+b)
           tol1=prec*abs(x)+prec
           tol2=2.*tol1
           if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3    ! Test for done here.
           if(abs(e).gt.tol1) then                    ! Construct a trial parabolic fit.
           r=(x-w)*(fx-fv)
           q=(x-v)*(fx-fw)
           p=(x-v)*q-(x-w)*r
           q=2.0*(q-r)
           if(q.gt.0.) p=-p
           q=abs(q)
           etemp=e
           e=d
           if( (abs(p)>=abs(.5*q*etemp)) .or. (p<=q*(a-x)) .or. (p>=q*(b-x)) ) goto 1
! The above conditions determine the acceptability of the parabolic fit. Here it is o.k.:
           d=p/q                                                            ! Take the parabolic step.
           u=x+d
           if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
           goto 2                                                           ! Skip over the golden section step.
           endif
 1         if(x.ge.xm) then                                                 ! We arrive here for a golden section step, which we take
              e=a-x                                                         ! into the larger of the two segments.
           else
              e=b-x
           end if
           d=CGOLD*e                                                        ! Take the golden section step.
 2         if(abs(d).ge.tol1) then                                          ! Arrive here with d computed either from parabolic fit, or
              u=x+d                                                         ! else from golden section.
           else
              u=x+sign(tol1,d)
           end if
           fu=f1dim(this,u)                                                 ! This is the one function evaluation per iteration,
           if(fu.le.fx) then                                                ! and now we have to decide what to do with our function
              if(u.ge.x) then                                               ! evaluation. Housekeeping follows:
                 a=x
              else
                 b=x
              end if
              v=w
              fv=fw
              w=x
              fw=fx
              x=u
              fx=fu
           else
              if(u.lt.x) then
                 a=u
              else
                 b=u
              end if
              if(fu.le.fw .or. w.eq.x) then
                 v=w
                 fv=fw
                 w=u
                 fw=fu
              else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                 v=u
                 fv=fu
              end if
           endif                                      ! Done with housekeeping. Back for another iteration.
        end do
 3      xmin=x                                        ! Arrive here ready to exit with best values.
        brent=fx
 
end function brent
!
!
!
end module NonlinearCG_m

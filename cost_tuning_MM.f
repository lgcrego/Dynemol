module cost_MM

    use type_m
    use constants_m
    use parameters_m  , only: T_ , F_ , N_of_AdSteps 
    use MM_types      , only: MMOPT_Control, LogicalKey

    public :: evaluate_cost , nmd_REF_erg , nmd_NOPT_ERG , SetKeys , KeyHolder 

    private 

    ! module variables ...
    integer                        :: AStep = 0
    integer          , allocatable :: nmd_sorted_indx(:) , nmd_deg_indx(:,:)
    real*8           , allocatable :: nmd_REF_erg(:) , nmd_NOPT_ERG(:) , AdiabaticStep(:)
    type(LogicalKey) , allocatable :: KeyHolder(:)

contains
!
!
!
!==================
 subroutine SetKeys
!==================
implicit none

If( .not. allocated(KeyHolder) ) allocate( KeyHolder(1) )


    KeyHolder(1)%comment = "==> optimize all"
    KeyHolder(1)%bonds   = [T_,F_,F_]
    KeyHolder(1)%angs    = [T_,F_]
    KeyHolder(1)%diheds  = [F_,F_,T_,F_,F_,F_,F_]



end subroutine SetKeys
!
!
!
!========================================================
 function evaluate_cost( Hesse_erg , nmd_indx , control )
!========================================================
implicit none
real*8              , optional               , intent(in)    :: Hesse_erg(:)
integer             , optional , allocatable , intent(inout) :: nmd_indx(:)
type(MMOPT_Control)                          , intent(in)    :: control
real*8              :: evaluate_cost

! local variables ...
integer      :: i 
real*8       :: order_cost, split_cost
real*8       :: chi(100)    = D_zero 
real*8       :: weight(100) = D_zero
character(6) :: prepare_environment

If( control% new_adiabat ) then
    AStep = AStep + 1 
    evaluate_cost = real_large
    return
end If

If( control% preprocess ) then

    evaluate_cost = real_large

    ! reading command line arguments for initiating parameterization ...
    If( COMMAND_ARGUMENT_COUNT() == 0 ) pause "Quit and choose option: newOPT , repeat , resume." 

    CALL GET_COMMAND_ARGUMENT( 1 , prepare_environment )

    select case( prepare_environment )
         case( "repeat" ) ;  CALL SYSTEM( "mv OPT_nmd_indx.old OPT_nmd_indx.inpt" )
         case( "resume" ) ;  CALL SYSTEM( "mv OPT_nmd_indx.out OPT_nmd_indx.inpt" )
    end select

    select case( prepare_environment )

         case( "repeat" , "resume" ) 
              OPEN( unit=14 , file='OPT_nmd_indx.inpt' )
                   read(14,*) i
                   allocate( nmd_indx(i)  )    
                   read(14,*) ( nmd_indx(i) , i=1,size(nmd_indx) )
              close(14) 

         case( "newOPT" )

         allocate( nmd_indx , source = [7,8,10,9,11,13,12,15,14,16,17,19,21,20,18,22,23,28,25,24,26,27,30,29, &
                                        31,37,34,36,32,38,35,33,39,40,41,43,42,48,46,44,54,53,59,57,51,45,47, &
                                        50,49,52,58,55] )
         case default 
              pause "Quit and choose option: newOPT , repeat , resume." 

    end select

end If

!------------------------
! NMD frequencies (cm-1})
! ascending energy order
!------------------------

chi(1)  = Hesse_erg(nmd_indx(1))  - 26.9d0                             ; weight(1)  = 1.0d0

chi(2)  = Hesse_erg(nmd_indx(2))  - 94.2d0                             ; weight(2)  = 1.0d0

chi(3)  = Hesse_erg(nmd_indx(3))  - 125.5d0                            ; weight(3)  = 1.0d0

chi(4)  = Hesse_erg(nmd_indx(4))  - 176.8d0                            ; weight(4)  = 1.0d0

chi(5)  = Hesse_erg(nmd_indx(5))  - 203.1d0                            ; weight(5)  = 1.0d0

chi(6)  = Hesse_erg(nmd_indx(6))  - 231.2d0                            ; weight(6)  = 5.0d0

chi(7)  = Hesse_erg(nmd_indx(7))  - 254.2d0                            ; weight(7)  = 1.0d0

chi(8)  = Hesse_erg(nmd_indx(8))  - 293.3d0                            ; weight(8)  = 1.0d0

chi(9)  = Hesse_erg(nmd_indx(9))  - 352.4d0                            ; weight(9)  = 1.0d0

chi(10) = Hesse_erg(nmd_indx(10)) - 360.7d0                            ; weight(10) = 1.0d0

chi(11) = Hesse_erg(nmd_indx(11)) - 416.8d0                            ; weight(11) = 1.0d0

chi(12) = Hesse_erg(nmd_indx(12)) - 429.3d0                            ; weight(12) = 1.0d0

chi(13) = Hesse_erg(nmd_indx(13)) - 448.0d0                            ; weight(13) = 1.0d0

chi(14) = Hesse_erg(nmd_indx(14)) - 461.9d0                            ; weight(14) = 5.0d0

chi(15) = Hesse_erg(nmd_indx(15)) - 466.6d0                            ; weight(15) = 1.0d0

chi(16) = Hesse_erg(nmd_indx(16)) - 520.4d0                            ; weight(16) = 1.0d0

chi(17) = Hesse_erg(nmd_indx(17)) - 534.1d0                            ; weight(17) = 1.0d0

chi(18) = Hesse_erg(nmd_indx(18)) - 537.2d0                            ; weight(18) = 1.0d0

chi(19) = Hesse_erg(nmd_indx(19)) - 545.4d0                            ; weight(19) = 1.0d0

chi(20) = Hesse_erg(nmd_indx(20)) - 547.3d0                            ; weight(20) = 5.0d0

chi(21) = Hesse_erg(nmd_indx(21)) - 580.3d0                            ; weight(21) = 1.0d0

chi(22) = Hesse_erg(nmd_indx(22)) - 618.3d0                            ; weight(22) = 1.0d0

chi(23) = Hesse_erg(nmd_indx(23))  - 625.3d0                           ; weight(23)= 1.0d0

chi(24) = Hesse_erg(nmd_indx(24))  - 639.5d0                           ; weight(24)= 1.0d0

chi(25) = Hesse_erg(nmd_indx(25))  - 643.2d0                           ; weight(25)= 1.0d0

chi(26) = Hesse_erg(nmd_indx(26))  - 743.0d0                           ; weight(26)= 1.0d0

chi(27) = Hesse_erg(nmd_indx(27))  - 745.0d0                           ; weight(27)= 1.0d0

chi(28) = Hesse_erg(nmd_indx(28))  - 757.8d0                           ; weight(28)= 1.0d0

chi(29) = Hesse_erg(nmd_indx(29))  - 773.3d0                           ; weight(29)= 1.0d0

chi(30) = Hesse_erg(nmd_indx(30))  - 777.6d0                           ; weight(30)= 1.0d0

chi(31) = Hesse_erg(nmd_indx(31))  - 791.4d0                           ; weight(31)= 1.0d0

chi(32) = Hesse_erg(nmd_indx(32))  - 796.9d0                           ; weight(32)= 1.0d0

chi(33) = Hesse_erg(nmd_indx(33))  - 819.1d0                           ; weight(33)= 1.0d0

chi(34) = Hesse_erg(nmd_indx(34))  - 819.8d0                           ; weight(34)= 1.0d0

chi(35) = Hesse_erg(nmd_indx(35))  - 848.9d0                           ; weight(35)= 1.0d0

chi(36) = Hesse_erg(nmd_indx(36))  - 869.1d0                           ; weight(36)= 1.0d0

chi(37) = Hesse_erg(nmd_indx(37))  - 885.9d0                           ; weight(37)= 1.0d0

chi(38) = Hesse_erg(nmd_indx(38))  - 890.8d0                           ; weight(38)= 1.0d0

chi(39) = Hesse_erg(nmd_indx(39))  - 904.4d0                           ; weight(39)= 1.0d0

chi(40) = Hesse_erg(nmd_indx(40))  - 932.6d0                           ; weight(40)= 1.0d0

chi(41) = Hesse_erg(nmd_indx(41))  - 950.6d0                           ; weight(41)= 1.0d0

chi(42) = Hesse_erg(nmd_indx(42))  - 955.3d0                           ; weight(42)= 1.0d0

chi(43) = Hesse_erg(nmd_indx(43))  - 960.7d0                           ; weight(43)= 1.0d0

chi(44) = Hesse_erg(nmd_indx(44))  - 964.8d0                           ; weight(44)= 1.0d0

chi(45) = Hesse_erg(nmd_indx(45))  -  985.5d0                          ; weight(45)= 5.0d0

chi(46) = Hesse_erg(nmd_indx(46))  - 1037.9d0                          ; weight(46)= 1.0d0

chi(47) = Hesse_erg(nmd_indx(47))  - 1053.0d0                          ; weight(47)= 1.0d0

chi(48) = Hesse_erg(nmd_indx(48))  - 1089.0d0                          ; weight(48)= 1.0d0

chi(49) = Hesse_erg(nmd_indx(49))  - 1103.5d0                          ; weight(49)= 1.0d0

chi(50) = Hesse_erg(nmd_indx(50))  - 1127.4d0                          ; weight(50)= 1.0d0

chi(51) = Hesse_erg(nmd_indx(51))  - 1142.8d0                          ; weight(51)= 1.0d0

chi(52) = Hesse_erg(nmd_indx(52))  - 1144.2d0                          ; weight(52)= 1.0d0

!--------------------------------------------------------------------

If( control% preprocess ) then
    allocate( nmd_REF_erg     , source = Hesse_erg(nmd_indx)-chi(1:size(nmd_indx))         )
    allocate( nmd_NOPT_erg    , source = Hesse_erg(nmd_indx)                               )
    allocate( nmd_deg_indx    , source = IdentifyDegenerates( Hesse_erg , chi , nmd_indx ) )
    allocate( nmd_sorted_indx , source = sort( nmd_indx )                                  )
    allocate( AdiabaticStep   , source = chi/float(N_of_AdSteps)                           )
end If

! include out of order cost ...
order_cost = D_zero
If( control% LineUpCost ) order_cost = sum([( abs(nmd_indx(i)-nmd_sorted_indx(i)) , i=1,size(nmd_indx) )]) 

! include degenerate splitting cost ...
split_cost = sum(abs(chi(nmd_deg_indx(:,1))-chi(nmd_deg_indx(:,2))))

! finally apply weight on chi and evaluate cost ...
If( control% use_no_weights) weight = D_one
If( control% adiabatic_OPT ) chi    = chi - AdiabaticStep*(N_of_AdSteps - AStep)

chi = chi * weight
evaluate_cost = sqrt(dot_product(chi,chi)) + order_cost + split_cost  

end function evaluate_cost
!
!
!
!=====================================================================
 function IdentifyDegenerates( Hesse_erg , chi , nmd_indx ) result(me)
!=====================================================================
implicit none
real*8        , intent(in) :: Hesse_erg(:)
real*8        , intent(in) :: chi(:)
integer       , intent(in) :: nmd_indx(:)

! local variables ...
integer               :: i , j , k
integer , allocatable :: me(:,:) , tmp(:,:)

allocate( tmp(size(nmd_indx),2) )
k = 1
do i = 1 , size(nmd_indx)
    do j = i+1 , size(nmd_indx)
         If( abs(chi(i)-Hesse_erg(nmd_indx(i))) == abs(chi(j)-Hesse_erg(nmd_indx(j))) ) then

              tmp(k,1) = i
              tmp(k,2) = j

              k = k + 1

         end If
    end do
end do

allocate  ( me , source=tmp(1:k-1,:) )
deallocate( tmp )

end function IdentifyDegenerates
!
!
!
!===========================
 function  sort(b) result(a)
!===========================
implicit none
integer , intent(in) :: b(:)

! local variables ...
integer :: ra, l, n, ir, i, j
integer :: a(size(b))

!---------------------
!      SORT A(I) 
!---------------------
      a = b
      n = size(a)
      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         ra  = a(l)
      else
         ra = a(ir)
         a(ir) = a(1)
         ir = ir - 1
         if(ir .eq. 1) then
             a(1) = ra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(a(j) .lt. a(j+1)) j = j + 1
        endif
      if(ra .lt. a(j)) then
        a(i) = a(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      a(i) = ra
      goto 10

end function sort
!
!
!
end module cost_MM

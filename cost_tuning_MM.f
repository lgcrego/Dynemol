module cost_MM

    use type_m
    use constants_m
    use parameters_m  , only: T_ , F_ , N_of_CGSteps 
    use MM_types      , only: MMOPT_Control, LogicalKey

    public :: evaluate_cost , nmd_REF_erg , nmd_NOPT_ERG , SetKeys , KeyHolder , overweight , chi

    private 

    ! module variables ...
    integer                        :: AStep = 0
    integer          , allocatable :: nmd_sorted_indx(:) , nmd_deg_indx(:,:)
    real*8           , allocatable :: nmd_REF_erg(:) , nmd_NOPT_ERG(:) , AdiabaticStep(:) , overweight(:)
    real*8                         :: chi(100) = D_zero 
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
real*8       :: weight(100) = D_zero
real*8       :: order_cost, split_cost
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

         allocate( nmd_indx , source = [7,8,9,10,12,11,13,14,19,20,24,16,15,17,18,23,21,22,25,30,26,27,28,29, &
                                        36,34,35,31,32,33] )

         case default 
              pause "Quit and choose option: newOPT , repeat , resume." 

    end select

end If

!------------------------
! NMD frequencies (cm-1})
! ascending energy order
!------------------------

chi(1) = Hesse_erg(nmd_indx(1))  - 395.d0                             ; weight(1) = 1.0d0

chi(2) = Hesse_erg(nmd_indx(2))  - 395.d0                             ; weight(2) = 1.0d0

chi(3) = Hesse_erg(nmd_indx(3))  - 599.d0                             ; weight(3) = 1.0d0

chi(4) = Hesse_erg(nmd_indx(4))  - 599.d0                             ; weight(4) = 1.0d0

chi(5) = Hesse_erg(nmd_indx(5))  - 677.d0                             ; weight(5) = 1.0d0

chi(6) = Hesse_erg(nmd_indx(6))  - 700.d0                             ; weight(6) = 1.0d0

chi(7) = Hesse_erg(nmd_indx(7))  - 838.d0                             ; weight(7) = 1.0d0

chi(8) = Hesse_erg(nmd_indx(8))  - 838.d0                             ; weight(8) = 1.0d0

chi(9) = Hesse_erg(nmd_indx(9))  - 954.d0                             ; weight(9) = 1.0d0

chi(10)= Hesse_erg(nmd_indx(10))  - 954.d0                            ; weight(10)= 1.0d0

chi(11)= Hesse_erg(nmd_indx(11))  - 981.d0                            ; weight(11)= 1.0d0

chi(12)= Hesse_erg(nmd_indx(12))  - 994.d0                            ; weight(12)= 1.0d0

chi(13)= Hesse_erg(nmd_indx(13))  - 1014.d0                           ; weight(13)= 1.0d0

chi(14)= Hesse_erg(nmd_indx(14))  - 1048.d0                           ; weight(14)= 1.0d0

chi(15)= Hesse_erg(nmd_indx(15))  - 1048.d0                           ; weight(15)= 1.0d0

chi(16)= Hesse_erg(nmd_indx(16))  - 1137.d0                           ; weight(16)= 1.0d0

chi(17)= Hesse_erg(nmd_indx(17))  - 1162.d0                           ; weight(17)= 1.0d0

chi(18)= Hesse_erg(nmd_indx(18))  - 1162.d0                           ; weight(18)= 1.0d0

chi(19)= Hesse_erg(nmd_indx(19))  - 1332.d0                           ; weight(19)= 1.0d0

chi(20)= Hesse_erg(nmd_indx(20))  - 1376.d0                           ; weight(20)= 1.0d0

chi(21)= Hesse_erg(nmd_indx(21))  - 1477.d0                           ; weight(21)= 1.0d0

chi(22)= Hesse_erg(nmd_indx(22))  - 1477.d0                           ; weight(22)= 1.0d0

chi(23)= Hesse_erg(nmd_indx(23))  - 1612.d0                           ; weight(23)= 1.0d0

chi(24)= Hesse_erg(nmd_indx(24))  - 1612.d0                           ; weight(24)= 1.0d0

chi(25)= Hesse_erg(nmd_indx(25))  - 3080.d0                           ; weight(25)= 1.0d0

chi(26)= Hesse_erg(nmd_indx(26))  - 3092.d0                           ; weight(26)= 1.0d0

chi(27)= Hesse_erg(nmd_indx(27))  - 3092.d0                           ; weight(27)= 1.0d0

chi(28)= Hesse_erg(nmd_indx(28))  - 3110.d0                           ; weight(28)= 1.0d0

chi(29)= Hesse_erg(nmd_indx(29))  - 3110.d0                           ; weight(29)= 1.0d0

chi(30)= Hesse_erg(nmd_indx(30))  - 3121.d0                           ; weight(30)= 1.0d0

!--------------------------------------------------------------------

If( control% preprocess ) then
    allocate( nmd_REF_erg     , source = Hesse_erg(nmd_indx)-chi(1:size(nmd_indx))         )
    allocate( nmd_NOPT_erg    , source = Hesse_erg(nmd_indx)                               )
    allocate( nmd_deg_indx    , source = IdentifyDegenerates( Hesse_erg , chi , nmd_indx ) )
    allocate( nmd_sorted_indx , source = sort( nmd_indx )                                  )
    allocate( AdiabaticStep   , source = chi/float(N_of_CGSteps)                           )

    allocate( overweight(size(nmd_REF_erg)) )
    forall(i=1:size(nmd_REF_erg)) overweight(i) = overweight(i) + abs(chi(i)/nmd_REF_erg(i))
end If

! include out of order cost ...
order_cost = D_zero
If( control% LineUpCost ) order_cost = sum([( FOUR*abs(nmd_indx(i)-nmd_sorted_indx(i)) , i=1,size(nmd_indx) )]) 

! include degenerate splitting cost ...
split_cost = TWO*sum(abs(chi(nmd_deg_indx(:,1))-chi(nmd_deg_indx(:,2))))

! finally apply weight on chi and evaluate cost ...
If( control% use_no_weights) weight = D_one
If( control% use_overweight) forall(i=1:size(overweight)) weight(i) = weight(i) + overweight(i)
If( control% adiabatic_OPT ) chi = chi - AdiabaticStep*(N_of_CGSteps - AStep)

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

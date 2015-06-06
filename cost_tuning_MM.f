module cost_MM

    use type_m
    use constants_m
    use parameters_m  , only: T_ , F_

    public :: evaluate_cost , nmd_REF_erg , nmd_NOPT_ERG , SetKeys , KeyHolder , LogicalKey 

    private 

    ! module types ...
    type LogicalKey
        logical       :: bonds(3)  = F_
        logical       :: angs(2)   = F_
        logical       :: diheds(7) = F_
        character(20) :: comment
    end type

    ! module variables ...
    integer          , allocatable :: nmd_deg_indx(:,:)
    real*8           , allocatable :: nmd_REF_erg(:) , nmd_NOPT_ERG(:) 
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
!=========================================================
 function evaluate_cost( Hesse_erg , nmd_indx , instance )
!=========================================================
implicit none
real*8           , optional               , intent(in)     :: Hesse_erg(:)
integer          , optional , allocatable , intent(inout)  :: nmd_indx(:)
character(*)     , optional               , intent(in)     :: instance
real*8           :: evaluate_cost

! local variables ...
integer      :: i 
real*8       :: chi(100)    = D_zero
real*8       :: weight(100) = D_zero
logical      :: preprocess
character(6) :: prepare_environment

preprocess = present(instance) .AND. instance == "preprocess"

If( preprocess ) then

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

chi(1)  = Hesse_erg(nmd_indx(1))  - 26.9d0                             ; weight(1)  = 4.0d0

chi(2)  = Hesse_erg(nmd_indx(2))  - 94.2d0                             ; weight(2)  = 1.0d0

chi(3)  = Hesse_erg(nmd_indx(3))  - 125.5d0                            ; weight(3)  = 3.0d0

chi(4)  = Hesse_erg(nmd_indx(4))  - 176.8d0                            ; weight(4)  = 1.0d0

chi(5)  = Hesse_erg(nmd_indx(5))  - 203.1d0                            ; weight(5)  = 1.0d0

chi(6)  = Hesse_erg(nmd_indx(6))  - 231.2d0                            ; weight(6)  = 4.0d0

chi(7)  = Hesse_erg(nmd_indx(7))  - 254.2d0                            ; weight(7)  = 2.0d0

chi(8)  = Hesse_erg(nmd_indx(8))  - 293.3d0                            ; weight(8)  = 1.0d0

chi(9)  = Hesse_erg(nmd_indx(9))  - 352.4d0                            ; weight(9)  = 1.0d0

chi(10) = Hesse_erg(nmd_indx(10)) - 360.7d0                            ; weight(10) = 2.0d0

chi(11) = Hesse_erg(nmd_indx(11)) - 416.8d0                            ; weight(11) = 2.0d0

chi(12) = Hesse_erg(nmd_indx(12)) - 429.3d0                            ; weight(12) = 1.0d0

chi(13) = Hesse_erg(nmd_indx(13)) - 448.0d0                            ; weight(13) = 1.0d0

chi(14) = Hesse_erg(nmd_indx(14)) - 461.9d0                            ; weight(14) = 4.0d0

chi(15) = Hesse_erg(nmd_indx(15)) - 466.6d0                            ; weight(15) = 1.0d0

chi(16) = Hesse_erg(nmd_indx(16)) - 520.4d0                            ; weight(16) = 1.0d0

chi(17) = Hesse_erg(nmd_indx(17)) - 534.1d0                            ; weight(17) = 1.0d0

chi(18) = Hesse_erg(nmd_indx(18)) - 537.2d0                            ; weight(18) = 1.0d0

chi(19) = Hesse_erg(nmd_indx(19)) - 545.4d0                            ; weight(19) = 1.0d0

chi(20) = Hesse_erg(nmd_indx(20)) - 547.3d0                            ; weight(20) = 3.0d0

chi(21) = Hesse_erg(nmd_indx(21)) - 580.3d0                            ; weight(21) = 1.0d0

chi(22) = Hesse_erg(nmd_indx(22)) - 618.3d0                            ; weight(22) = 1.0d0

chi(23) = Hesse_erg(nmd_indx(23))  - 625.3d0                           ; weight(23)= 3.0d0

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

chi(39) = Hesse_erg(nmd_indx(39))  - 904.4d0                           ; weight(39)= 3.0d0

chi(40) = Hesse_erg(nmd_indx(40))  - 932.6d0                           ; weight(40)= 1.0d0

chi(41) = Hesse_erg(nmd_indx(41))  - 950.6d0                           ; weight(41)= 1.0d0

chi(42) = Hesse_erg(nmd_indx(42))  - 955.3d0                           ; weight(42)= 1.5d0

chi(43) = Hesse_erg(nmd_indx(43))  - 960.7d0                           ; weight(43)= 1.5d0

chi(44) = Hesse_erg(nmd_indx(44))  - 964.8d0                           ; weight(44)= 1.5d0

chi(45) = Hesse_erg(nmd_indx(45))  -  985.5d0                          ; weight(45)= 4.0d0

chi(46) = Hesse_erg(nmd_indx(46))  - 1037.9d0                          ; weight(46)= 1.0d0

chi(47) = Hesse_erg(nmd_indx(47))  - 1053.0d0                          ; weight(47)= 1.0d0

chi(48) = Hesse_erg(nmd_indx(48))  - 1089.0d0                          ; weight(48)= 1.0d0

chi(49) = Hesse_erg(nmd_indx(49))  - 1103.5d0                          ; weight(49)= 3.0d0

chi(50) = Hesse_erg(nmd_indx(50))  - 1127.4d0                          ; weight(50)= 1.0d0

chi(51) = Hesse_erg(nmd_indx(51))  - 1142.8d0                          ; weight(51)= 1.0d0

chi(52) = Hesse_erg(nmd_indx(52))  - 1144.2d0                          ; weight(52)= 1.0d0

!--------------------------------------------------------------------

If( preprocess ) then

    allocate( nmd_REF_erg  , source = Hesse_erg(nmd_indx)-chi(1:size(nmd_indx))         )
    allocate( nmd_NOPT_erg , source = Hesse_erg(nmd_indx)                               )
    allocate( nmd_deg_indx , source = IdentifyDegenerates( Hesse_erg , chi , nmd_indx ) )

end If

! finally apply weight on chi and evaluate cost ...
If( present(instance) .AND. instance == "use_no_weights") weight = D_one

chi = chi * weight
evaluate_cost = sqrt( dot_product(chi,chi) ) + sum(abs(chi(nmd_deg_indx(:,1))-chi(nmd_deg_indx(:,2))))

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
end module cost_MM

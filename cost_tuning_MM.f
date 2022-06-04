module cost_MM

    use type_m
    use constants_m
    use MM_input      , only: nmd_window
    use MM_types      , only: MMOPT_Control, LogicalKey

    public :: evaluate_cost , nmd_REF_erg , nmd_NOPT_ERG , SetKeys , KeyHolder , overweight , chi

    private 

    ! module variables ...
    integer          , allocatable :: nmd_sorted_indx(:) , nmd_deg_indx(:,:)
    real*8           , allocatable :: nmd_REF_erg(:) , nmd_NOPT_ERG(:) , overweight(:)
    real*8                         :: chi(100) = D_zero 
    type(LogicalKey) , allocatable :: KeyHolder(:)

    ! module parameters ...
    logical, parameter :: T_ = .true. , F_ = .false.

contains
!
!==================================================================================================================================
! BONDS
! case ('harm' , 1 )                                        case ('Mors' , 3 )
! V = 0.5*k_1*( rijsq - k_2) )^2                            V = k_1*[ 1.0 - ( rijsq - k_2)*exp(-k_3) ]^2
!
! ANGLES
! case ('harm' , 1 )                                        case ('urba' , 5 )
! V = 0.5*k_1*( phi - k_2 )^2                               V = V_harm + 0.5*k_3*( riksq - k_4)^2
!
! DIHEDRALS
! case ('cos' , 1 )    
! V = k_phi * [ 1 + cos( n * phi - phi_s ) ]                                ! Eq. 4.60 (GMX 5.0.5 manual)
! V = k_2 * (1.0 + cos[ int(k_3)*phi - k_1 ]  )
!
! case ('harm' , 2 )   
! V = 1/2.k( xi - xi_0 )²                                                   ! Eq. 4.59 (GMX 5.0.5 manual)
! V = 0.5 * k_2 * ( phi - k_1 )^2
!
! case ('cos3' , 3 )    
! V = C0 + C1*cos(phi - 180) + C2*cos^2(phi - 180) + C3*cos^3(phi - 180) + C4*cos^4(phi - 180) + C5*cos(phi - 180)
! psi = phi-180                                                             ! Eq. 4.61 (GMX 5.0.5 manual)
! V = k_1 + k_2*cos(psi) + k_3*cos(psi)^2 + k_4*cos(psi)^3 + k_5*cos(psi)^4 + k_6*cos(psi)^5
!
! case ('imp' ,  4 )    
! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (improper)                     ! Eq. 4.60 (GMX 5.0.5 manual)
! V = k_2 * ( 1.0 + cos [ int(k_3)*phi - k_1 ] )
!
! case ('chrm' , 9 )   
! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple)                     ! Eq. 4.60 (GMX 5.0.5 manual)
!==================================================================================================================================
!
!==================
 subroutine SetKeys
!==================
implicit none

If( .not. allocated(KeyHolder) ) allocate( KeyHolder(1) )

    KeyHolder(1)%comment = "==> optimize all"
    KeyHolder(1)%bonds        = [T_,F_,F_]
    KeyHolder(1)%angs         = [T_,F_,F_,F_]
    KeyHolder(1)%diheds(1:7)  = [F_,F_,T_,F_,F_,F_,F_]
    KeyHolder(1)%dihedtype    = 3 

! NOTICE: KeyHolder is defined as a VECTOR to allow for sequential optimization routines with different keys ...
! however, this feature is not implemented and KeyHolder behaves as a SCALAR ...
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

         allocate( nmd_indx , source = [7,8,9,10,12,11,13,14,19,20,23,15,16,17,18,24,21,22,25,30,26,27,28,29,31,34,35,32,33,36] )

         case default 
              pause "Quit and choose option: newOPT , repeat , resume." 

    end select

end If


!------------------------
! NMD frequencies (cm-1})
! ASCENDING energy order
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

!--------------------------------------------------------------------

If( control% preprocess ) then

    allocate( nmd_REF_erg     , source = Hesse_erg(nmd_indx)-chi(1:size(nmd_indx))         )
    allocate( nmd_NOPT_erg    , source = Hesse_erg(nmd_indx)                               )
    allocate( nmd_deg_indx    , source = IdentifyDegenerates( Hesse_erg , chi , nmd_indx ) )
    allocate( nmd_sorted_indx , source = sort( nmd_indx )                                  )

    allocate( overweight(size(nmd_REF_erg)) , source = D_zero )
    forall(i=1:size(nmd_REF_erg)) overweight(i) = abs(chi(i)/nmd_REF_erg(i))   

end If

! include out of order cost ...
order_cost = D_zero
If( control% LineUpCost ) order_cost = sum([( FOUR*abs(nmd_indx(i)-nmd_sorted_indx(i)) , i=1,size(nmd_indx) )]) 

! include splitting cost, for modes that should be degenerate but are not ...
split_cost = TWO*sum(abs(chi(nmd_deg_indx(:,1))-chi(nmd_deg_indx(:,2))))

! finally apply weight on chi and evaluate cost ...
If( control% use_no_weights) weight = D_one
If( control% use_overweight) forall(i=1:size(overweight)) weight(i) = weight(i) + overweight(i)

If( nmd_window% inicio/=0 .AND. nmd_window% fim/=0) then
    weight(:nmd_window%inicio-1 ) = D_zero
    weight( nmd_window%fim   +1:) = D_zero
EndIf

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
         ! identifies (i,j) <== (nmd_REF_erg(i) = nmd_REF_erg(j)) ...
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

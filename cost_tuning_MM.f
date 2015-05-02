module cost_MM

    use type_m
    use constants_m
    use parameters_m     , only : read_nmd_indx_

    public :: evaluate_cost , nmd_REF_erg , nmd_NOPT_ERG , SetKeys , KeyHolder , LogicalKey 

    private 

    ! module types ...
    type LogicalKey
        logical       :: bonds(3)
        logical       :: angs(2)
        logical       :: diheds(7)
        character(30) :: comment
    end type

    ! module variables ...
    real*8           , allocatable :: nmd_REF_erg(:) , nmd_NOPT_ERG(:)
    type(LogicalKey) , allocatable :: KeyHolder(:)

contains
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
integer  :: i 
real*8   :: chi(30)    = D_zero
real*8   :: weight(30) = D_zero


select case (instance)

    case("preprocess")

        evaluate_cost = real_large

        If( read_nmd_indx_ ) then 

            OPEN( unit=14 , file='OPT_nmd_indx.inpt' )
                read(14,*) i
                allocate( nmd_indx(i)  )    
                read(14,*) ( nmd_indx(i) , i=1,size(nmd_indx) )
            close(14) 
        else
            
            allocate( nmd_indx , source = [7,8,9,10,12,11,13,14,19,20,24,15,16,17,18,23,21,22,25,30,26,27,29,28] )
           !allocate( nmd_indx , source = [7,8,9,10,12,11,13,14,19,20,24,15,16,17,18,23,21,22,25,30,26,27,29,28,31] )

        end If

    case default

!------------------------
! NMD frequencies (cm-1})
!------------------------

chi(1) = Hesse_erg(nmd_indx(1))  - 395.d0                             ; weight(1) = 2.0d0

chi(2) = Hesse_erg(nmd_indx(2))  - 395.d0                             ; weight(2) = 2.0d0

chi(3) = Hesse_erg(nmd_indx(3))  - 599.d0                             ; weight(3) = 4.0d0

chi(4) = Hesse_erg(nmd_indx(4))  - 599.d0                             ; weight(4) = 4.0d0

chi(5) = Hesse_erg(nmd_indx(5))  - 677.d0                             ; weight(5) = 1.0d0

chi(6) = Hesse_erg(nmd_indx(6))  - 700.d0                             ; weight(6) = 4.0d0

chi(7) = Hesse_erg(nmd_indx(7))  - 838.d0                             ; weight(7) = 1.0d0

chi(8) = Hesse_erg(nmd_indx(8))  - 838.d0                             ; weight(8) = 1.0d0

chi(9) = Hesse_erg(nmd_indx(9))  - 954.d0                             ; weight(9) = 1.0d0

chi(10)= Hesse_erg(nmd_indx(10))  - 954.d0                            ; weight(10)= 1.0d0

chi(11)= Hesse_erg(nmd_indx(11))  - 981.d0                            ; weight(11)= 1.0d0

chi(12)= Hesse_erg(nmd_indx(12))  - 994.d0                            ; weight(12)= 3.0d0

chi(13)= Hesse_erg(nmd_indx(13))  -  1014.d0                          ; weight(13)= 4.0d0

chi(14)= Hesse_erg(nmd_indx(14))  -  1048.d0                          ; weight(14)= 3.0d0

chi(15)= Hesse_erg(nmd_indx(15))  -  1048.d0                          ; weight(15)= 3.0d0

chi(16)= Hesse_erg(nmd_indx(16))  - 1137.d0                           ; weight(16)= 1.0d0

chi(17)= Hesse_erg(nmd_indx(17))  - 1162.d0                           ; weight(17)= 2.0d0

chi(18)= Hesse_erg(nmd_indx(18))  - 1162.d0                           ; weight(18)= 2.0d0

chi(19)= Hesse_erg(nmd_indx(19))  - 1332.d0                           ; weight(19)= 1.0d0

chi(20)= Hesse_erg(nmd_indx(20))  - 1376.d0                           ; weight(20)= 1.0d0

chi(21)= Hesse_erg(nmd_indx(21))  - 1477.d0                           ; weight(21)= 1.0d0

chi(22)= Hesse_erg(nmd_indx(22))  - 1477.d0                           ; weight(22)= 1.0d0

chi(23)= Hesse_erg(nmd_indx(23))  - 1612.d0                           ; weight(23)= 1.0d0

chi(24)= Hesse_erg(nmd_indx(24))  - 1612.d0                           ; weight(24)= 1.0d0

!chi(25)= Hesse_erg(nmd_indx(25))  - 3080.d0                           ; weight(25)= 1.0d0

!--------------------------------------------------------------------

If( .NOT. allocated(nmd_REF_erg ) ) allocate( nmd_REF_erg  , source = Hesse_erg(nmd_indx)-chi(1:size(nmd_indx)) )
If( .NOT. allocated(nmd_NOPT_erg) ) allocate( nmd_NOPT_erg , source = Hesse_erg(nmd_indx) )

! apply weight on chi and evaluate cost ...

chi = chi * weight
evaluate_cost = sqrt( dot_product(chi,chi) )

end select

end function evaluate_cost
!
!
!
!==================
 subroutine SetKeys
!==================
implicit none

! local variables ...
logical :: F_ = .false. , T_ = .true. 

If( .not. allocated(KeyHolder) ) allocate( KeyHolder(1) )

    KeyHolder(1)%comment = "==> optimize all"
    KeyHolder(1)%bonds   = [T_,F_,F_]
    KeyHolder(1)%angs    = [T_,F_]
    KeyHolder(1)%diheds  = [F_,F_,T_,F_,F_,F_,F_]

end subroutine SetKeys
!
!
!
end module cost_MM

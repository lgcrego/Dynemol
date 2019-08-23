module diagnosis_m

use types_m     , only : universe

contains
!
!
!==========================
subroutine diagnosis( sys )
!==========================
implicit none
type(universe)  , intent(in)    :: sys

! local variables
integer         :: i , j
character(2)    :: yn , option

Write(*,'(/a)') sys%Surface_Characteristics

write(*,'(/a)',advance='no') ">>> details of the input file ? (y/n) "
read (*,'(a)') yn

if ( yn == "y" ) then

    do
        write(*,'(/a)', advance='no') '>>>  Symbol (1)  ;  MMSymbol  (2)  ;  Fragment (3)  ;  resid (4)  ;  nresid (5)  ;   &
                                            copy (6)  ;  delete (7)  ;  translate  (8)  ;  rotate (9)  ;  T/F (10)  '
        read (*,'(a)') option

        select case ( option )

            case ('1')
                Write(*,10) sys % atom % Symbol 

            case ('2')
                Write(*,20) sys % atom % MMSymbol 

            case ('3')
                Write(*,30) sys % atom % Fragment

            case ('4')
                Write(*,40) sys % atom % resid

            case ('5')
                Write(*,50) sys % atom % nresid

            case ('6')
                Write(*,70) sys % atom % copy

            case ('7')
                Write(*,70) sys % atom % delete

            case ('8')
                Write(*,70) sys % atom % translate

            case ('9')
                Write(*,70) sys % atom % rotate

            case ('10')
                Write(*,80) ( ( sys % atom(i) % TorF(j) , j=1,3 ) , " " , i=1,sys%N_of_atoms ) 

            case default
                exit

        end select

    end do

end if

10 Format(50a3)
20 Format(37a4)
30 Format(50a3)
40 Format(37a4)
50 Format(50i3)
60 Format(50a3)
70 Format(50L3)
80 Format(152a1)

end subroutine diagnosis

end module diagnosis_m
!
!
!
!

module diagnosis_m

use ansi_colors
use types_m     , only : universe

public :: diagnosis

private

contains
!
!
!==========================
subroutine diagnosis( sys )
!==========================
implicit none
type(universe)  , intent(in)  :: sys

! local variables
integer      :: i, j
character(2) :: yn, option

write(*,'(/a)') bold // cyan // sys%System_Characteristics // reset

write(*,'(/a)', advance='no') bold // yellow // '>>> Details of the input file? (y/n): ' // reset
read (*,'(a)') yn

if (yn == 'y') then

    do
        write(*,'(/a)', advance='no') bold // cyan // &
            '>>>  Symbol (1); MMSymbol (2); Fragment (3); resid (4); nresid (5); ' // &
            'copy (6); delete (7); translate (8); rotate (9); T/F (10): ' // reset
        read (*,'(a)') option

        select case (option)

            case ('1')
                write(*,'(a)') green // '>> Symbol:' // reset
                write(*,10) sys%atom%Symbol

            case ('2')
                write(*,'(a)') green // '>> MMSymbol:' // reset
                write(*,20) sys%atom%MMSymbol

            case ('3')
                write(*,'(a)') green // '>> Fragment:' // reset
                write(*,30) sys%atom%Fragment

            case ('4')
                write(*,'(a)') green // '>> Resid:' // reset
                write(*,40) sys%atom%resid

            case ('5')
                write(*,'(a)') green // '>> nresid:' // reset
                write(*,50) sys%atom%nresid

            case ('6')
                write(*,'(a)') green // '>> Copy flag:' // reset
                write(*,70) sys%atom%copy

            case ('7')
                write(*,'(a)') red // '>> Delete flag:' // reset
                write(*,70) sys%atom%delete

            case ('8')
                write(*,'(a)') green // '>> Translate:' // reset
                write(*,70) sys%atom%translate

            case ('9')
                write(*,'(a)') green // '>> Rotate:' // reset
                write(*,70) sys%atom%rotate

            case ('10')
                write(*,'(a)') green // '>> True/False flags:' // reset
                write(*,80) ( ( sys%atom(i)%TorF(j), j=1,3 ), ' ', i=1, sys%N_of_atoms )

            case default
                exit

        end select
    end do

end if

10 Format(50a3)
20 Format(37a4)
30 Format(50a3)
40 Format(37a4)
50 Format(30i5)
60 Format(50a3)
70 Format(50L3)
80 Format(152a1)

end subroutine diagnosis

end module diagnosis_m
!
!
!
!

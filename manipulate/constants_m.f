MODULE constants_m
    INTEGER    , PARAMETER :: ABOVE=+100000000
    INTEGER    , PARAMETER :: BELOW=-100000000
	REAL*8     , PARAMETER :: PI=3.141592653589793238462643383279502884197d0 
	REAL*8     , PARAMETER :: PIO2=1.57079632679489661923132169163975144209858d0 
	REAL*8     , PARAMETER :: TWOPI=6.283185307179586476925286766559005768394d0 
	REAL*8     , PARAMETER :: FOURPI=12.56637061435917295385d0
    REAL*8     , PARAMETER :: PI4 = FOURPI
    REAL*8     , PARAMETER :: RADIAN = 180.D0 / PI
    real*8     , PARAMETER :: rad_2_deg = 180.0d0 / PI
    real*8     , PARAMETER :: deg_2_rad = PI / 180.0d0
	REAL*8     , PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967d0 
    REAL*8     , PARAMETER :: SQRT2PI=2.50662827463100050242d0
	REAL*8     , PARAMETER :: EULER=0.5772156649015328606065120900824024310422d0 
	REAL*8     , PARAMETER :: D_zero= 0.d0
	REAL*8     , PARAMETER :: D_one = 1.d0
    REAL*8     , PARAMETER :: HALF  = 5.d-1
    REAL*8     , PARAMETER :: THIRD = 1.d0/3.d0
    REAL*8     , PARAMETER :: TWO   = 2.d0
    REAL*8     , PARAMETER :: THREE = 3.d0
    REAL*8     , PARAMETER :: FOUR  = 4.d0
    REAL*8     , PARAMETER :: FIVE  = 5.d0
    REAL*8     , PARAMETER :: SIX   = 6.d0
	REAL*8     , PARAMETER :: a_Bohr=0.52917720859d0
    REAL*8     , PARAMETER :: h_bar=6.58264d-4 ! <== Planck's  constant  (eV * ps)
    REAL*8     , PARAMETER :: u_mass=1.660538782d-24 ! Atomic Mass Unit <== (gram/particle)
    REAL*8     , PARAMETER :: J_2_meV=6.24181d21 ! Joule to meV converter
    REAL*8     , PARAMETER :: J_2_kcalmol=1.44d20 ! Joule to kcal/mol converter
    REAL*8     , PARAMETER :: low_prec=1.d-8
    REAL*8     , PARAMETER :: high_prec=1.d-14
    REAL*8     , PARAMETER :: infty=1.d10
    complex*16 , parameter :: zi=(0.d0,1.d0)
END MODULE constants_m

module ansi_colors
! ANSI (American National Standards Institute)
! ANSI Escape Codes are sequences of characters used to control the 
! display of a terminal or command-line interface
    implicit none
    character(len=*), parameter :: esc = char(27)

    character(len=*), parameter :: reset = esc // '[0m'
    character(len=*), parameter :: bold  = esc // '[1m'

    character(len=*), parameter :: red     = esc // '[31m'
    character(len=*), parameter :: green   = esc // '[32m'
    character(len=*), parameter :: yellow  = esc // '[33m'
    character(len=*), parameter :: blue    = esc // '[34m'
    character(len=*), parameter :: magenta = esc // '[35m'
    character(len=*), parameter :: cyan    = esc // '[36m'
    character(len=*), parameter :: white   = esc // '[37m'
    character(len=*), parameter :: orange  = esc // '[38;5;208m'

    character(len=*), parameter :: bg_red     = esc // '[31;47m'
    character(len=*), parameter :: bg_green   = esc // '[32;47m'
    character(len=*), parameter :: bg_blue    = esc // '[34;47m'
    character(len=*), parameter :: bg_magenta = esc // '[35;47m'

end module ansi_colors
!
!
!
MODULE color_funcs

use ansi_colors 
implicit none

    public:: red_, green_, yellow_, blue_, magenta_, orange_
    public:: red_bg, green_bg, blue_bg, magenta_bg

    private

    contains

    !================================================
    function green_(string) 
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: green_
        green_ = green//string//reset
    end function green_
!   
    !===============================================
    function green_bg(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: green_bg
        green_bg = bg_green//string//reset
    end function green_bg
!
    !================================================
    function red_(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: red_
        red_ = red//string//reset
    end function red_
!   
    !===============================================
    function red_bg(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: red_bg
        red_bg = bg_red//string//reset
    end function red_bg
!   
    !================================================
    function yellow_(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: yellow_
        yellow_ = yellow//string//reset
    end function yellow_
!   
    !================================================
    function orange_(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: orange_
        orange_ = orange//string//reset
    end function orange_
!   
    !================================================
    function blue_(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: blue_
        blue_ = blue//string//reset
    end function blue_
!   
    !===============================================
    function blue_bg(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: blue_bg
        blue_bg = bg_blue//string//reset
    end function blue_bg
!   
    !================================================
    function magenta_(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: magenta_
        magenta_ = magenta//string//reset
    end function magenta_
!   
    !===============================================
    function magenta_bg(string)
        character(len=*), intent(in) :: string

        !local variable
        character(len=:), allocatable :: magenta_bg
        magenta_bg = bg_magenta//string//reset
    end function magenta_bg

END MODULE color_funcs

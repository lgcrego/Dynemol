module util_m

    use constants_m

    public:: fact , binomial
    public:: TO_UPPER_CASE , count_lines , split_line , seq_in_range , read_CSV_file

private

contains
!
!
!
!-------------------------
elemental function fact(n)
!-------------------------
real*8               :: fact
integer , intent(in) :: n

integer :: i

fact = 1.d0
if(n > 1) then
    do i = 2 , n
        fact = fact * dfloat(i)
    end do
end if

end function fact
!
!
!
!-------------------------------
elemental function binomial(n,k)
!-------------------------------
integer , intent(in) :: n,k
real*8               :: binomial

binomial = fact(n)/(fact(k)*fact(n-k))

end function binomial
!
!
!
!======================================
 pure FUNCTION TO_UPPER_CASE ( STRING )
!======================================
 implicit none
 CHARACTER ( LEN = * )              , INTENT(IN)    :: STRING
 CHARACTER ( LEN = LEN ( STRING ) )                 :: TO_UPPER_CASE

! Local parameters ...
INTEGER, PARAMETER :: BIG_A = ICHAR ( "A" ), LITTLE_A = ICHAR ( "a" ), LITTLE_Z = ICHAR ( "z" )

! Local scalars ...
INTEGER :: I, ICHR

! Loop over the characters in the string ...
DO I = 1,LEN ( STRING )

!   Get the ASCII order for the character to be converted ...
    ICHR = ICHAR ( STRING(I:I) )

!   Use the order to change the case of the character ...
    IF ( ( ICHR >= LITTLE_A ) .AND. ( ICHR <= LITTLE_Z ) ) THEN
        TO_UPPER_CASE(I:I) = CHAR ( ICHR + BIG_A - LITTLE_A )
    ELSE
        TO_UPPER_CASE(I:I) = STRING(I:I)
    END IF
END DO

END FUNCTION TO_UPPER_CASE
!
!
!
!================================================
 FUNCTION count_lines (f_name) result(n_of_lines)
!================================================
implicit none
character(len=*) , intent(in) :: f_name

! Local variables ...
INTEGER      :: i , ioerr , n_of_lines
character(1) :: dumb

! counts the number of lines of a file
                                                                                                                                               
OPEN( unit=3 , file=f_name , status='old' , action="read" )
i = 0
do 
    read(3,*,IOSTAT=ioerr) dumb
    if(ioerr < 0) EXIT
    i = i + 1
end do    
close(3)

n_of_lines = i

end FUNCTION count_lines
!
!
!
!====================================================
 FUNCTION count_columns (f_name) result(n_of_columns)
!====================================================
implicit none
character(len=*) , intent(in) :: f_name

! Local variables ...
INTEGER      :: i , ioerr , n_of_columns
character(80) :: line

! counts the number of columns of a file

OPEN( unit=3 , file=f_name , status='old' , action="read" )
read(3,'(a)') line
n_of_columns = size( parse_this(line) )
close(3)

end FUNCTION count_columns
!
!
!
!=========================================
 FUNCTION split_line (line) result(tokens)
!=========================================
implicit none
character(len=*) :: line

! Local variables ...
integer :: i , size_line , n_tokens , ioerr
character(len=3) , allocatable :: tokens(:) , aux(:)

! split line in tokens 

size_line = len(line)
allocate ( aux(size_line) )                                                                                                                   
read(line,*,iostat=ioerr) (aux(i) , i=1,size_line)
n_tokens = i-1

allocate ( tokens(n_tokens) ) 
forall(i=1:n_tokens) tokens(i) = trim(aux(i))
deallocate(aux)

end FUNCTION split_line
!
!
!
!========================================
function seq_in_range(string) result(seq)
!========================================
implicit none
character(len=*) :: string

! Local variables ...
integer               :: i , i1 , i2 , colon
integer , allocatable :: seq(:)

! returns a sequence of numbers for string = "i1:i2"

colon = index(string,":")
read(string(1:colon-1),*) i1
read(string(colon+1: ),*) i2

seq = [ ( i , i=i1,i2 ) ]

end function seq_in_range
!
!
!       
!=====================================
function parse_this(line) result(indx)
!=====================================
implicit none 
character(len=*)        , intent(in)  :: line

! local variables ...
integer                        :: i , num
character(len=9) , allocatable :: tokens(:)
integer          , allocatable :: indx(:)

! parse a string containing a list of numbers separated by spaces, ex.: 3 5 7 12 
! OR
! returns a sequence of numbers for string = "i1:i2" 

allocate( tokens , source = split_line(line) ) 

allocate(indx(0))

do i = 1 , size(tokens)
    if( scan(tokens(i),":") /= 0 ) &
    then
        indx = [ indx , seq_in_range(tokens(i)) ]
    else
        read(tokens(i),*) num 
        indx = [ indx , num ] 
    end if
end do

deallocate( tokens )

end function parse_this
!
!
!
!
!=========================================
subroutine read_CSV_file( f_name , dados )
!=========================================
implicit none
character(*)                   , intent(in)  :: f_name
character(len=9) , allocatable , intent(out) :: dados(:,:)

! local variables ...
integer                        :: i , nl , f_unit , row_length , iostat
character(len=10000)           :: line
character(len=9) , allocatable :: tokens(:)

nl = count_lines(f_name)

OPEN( file=f_name , status='old' , action="read" , newunit=f_unit)

! find number of records in the row ...
       read( f_unit , '(a)' , IOSTAT=iostat )  line
       backspace(f_unit)
       tokens = split_line(line)
       row_length = size(tokens)
       
       ! matrix of strings
       allocate( dados(nl,row_length) )

       do i = 1 , nl
            read( f_unit , '(a)' , IOSTAT=iostat )  line
            if (iostat < 0) stop "Unexpected EoF"
            dados(i,:) = split_line(line)
            end do
       
close(f_unit)

end subroutine read_CSV_file
!
!
!
end module util_m
!!!!!!!!!!!!!!!!!!!!!!!!!  ZIGGURAT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Copyright (C) 2000  George Marsaglia and Wai Wan Tsang
!    Copyright (C) 2013-2019  Jason M. Wood, Montana State University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> The ziggurat algorithm is a rejection sampling algorithm for pseudo-random
!> number sampling. The original C version was created by George Marsaglia and
!> Wai Wan Tsang [1].  The Fortran 90 version compatible with OpenMP was
!> created by Jason M. Wood.
!>
!> #### Reference
!>  [1] George Marsaglia, Wai Wan Tsang.  2000.  The Ziggurat Method for
!>        Generating Random Variables.  'Journal of Statistical Software'.
!>        Vol. 5, Issue 8.  http://www.jstatsoft.org/v05/i08/
!>
!> @author George Marsaglia
!> @author Wai Wan Tsang
!> @author Jason M. Wood
!> @copyright GNU General Public License
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ziggurat
  ! Load intrinsic modules.
  use, intrinsic :: iso_fortran_env

  implicit none

  private

  ! Declare public methods.
  public :: ziggurat_t
  public :: ziggurat_rexp
  public :: ziggurat_rnor
  public :: ziggurat_seed
  public :: ziggurat_shr3
  public :: ziggurat_uni

  ! The state variables for the ziggurat algorithm.
  type :: ziggurat_t
    integer(kind = int32) :: hz
    integer(kind = int32) :: iz
    integer(kind = int32) :: jz
    integer(kind = int32) :: jsr
    integer(kind = int32) :: ke(0:255)
    integer(kind = int32) :: kn(0:127)
    real(kind = real32)   :: fe(0:255)
    real(kind = real32)   :: fn(0:127)
    real(kind = real32)   :: we(0:255)
    real(kind = real32)   :: wn(0:127)
  end type ziggurat_t

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This procedure sets the seed and creates the RNOR and REXP tables for
  !> the ziggurat algorithm.
  !>
  !> state: The state variables for the ziggurat algorithm.
  !> iii: The seed for the ziggurat algorithm.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ziggurat_seed (state, iii)
    type(ziggurat_t), intent(inout)   :: state
    integer(kind = int32), intent(in) :: iii
    ! Local parameters.
    real(kind = real64), parameter :: m1 = 2147483648.0d0
    real(kind = real64), parameter :: m2 = 4294967296.0d0
    real(kind = real64), parameter :: ve = 3.949659822581572d-3
    real(kind = real64), parameter :: vn = 9.91256303526217d-3
    ! Local variables.
    integer(kind = int32) :: i
    real(kind = real64)   :: de
    real(kind = real64)   :: dn
    real(kind = real64)   :: te
    real(kind = real64)   :: tn
    real(kind = real64)   :: q
    dn = 3.442619855899
    de = 7.697117470131487
    tn = dn
    te = de
    state%jsr = iii
    ! Tables for RNOR:
    q = vn / exp (-0.5 * dn * dn)
    state%kn(0) = int (dn / q * m1, kind = int32)
    state%kn(1) = 0
    state%wn(0) = real (q / m1, kind = real32)
    state%wn(127) = real (dn / m1, kind = real32)
    state%fn(0) = 1.0
    state%fn(127) = real (exp (-0.5 * dn * dn), kind = real32)
    do i = 126, 1, -1
      dn = sqrt (-2.0 * log (vn / dn + exp (-0.5 * dn * dn)))
      state%kn(i + 1) = int (dn / tn * m1, kind = int32)
      tn = dn
      state%fn(i) = real (exp (-0.5 * dn * dn), kind = real32)
      state%wn(i) = real (dn / m1, kind = real32)
    end do
    ! Tables for REXP:
    q = ve / exp (-de)
    state%ke(0) = int (de / q * m2, kind = int32)
    state%ke(1) = 0
    state%we(0) = real (q / m2, kind = real32)
    state%we(255) = real (de / m2, kind = real32)
    state%fe(0) = 1.0
    state%fe(255) = real (exp (-de), kind = real32)
    do i = 254, 1, -1
      de = -log (ve / de + exp (-de))
      state%ke(i + 1) = int (de / te * m2, kind = int32)
      te = de
      state%fe(i) = real (exp (-de), kind = real32)
      state%we(i) = real (de / m2, kind = real32)
    end do
  end subroutine ziggurat_seed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates pseudo-random 32-bit integers.
  !>
  !> state: The state variables for the ziggurat algorithm.
  !> return: A pseudo-random 32-bit integer.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_shr3 (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    integer(kind = int32)           :: return_value
    state%jz = state%jsr
    state%jsr = ieor (state%jsr, ishft (state%jsr, 13))
    state%jsr = ieor (state%jsr, ishft (state%jsr, -17))
    state%jsr = ieor (state%jsr, ishft (state%jsr, 5))
    return_value = int (state%jz + state%jsr, kind = int32)
    return
  end function ziggurat_shr3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates a UNIFORMLY DISTRIBUTED pseudo-random value in the range [0,1).
  !>
  !> state: The state variables for the ziggurat algorithm.
  !> return: A uniformly distributed pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_uni (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    return_value = 0.5 + ziggurat_shr3(state) * 0.2328306e-9
    return
  end function ziggurat_uni

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates a NORMALLY DISTRIBUTED pseudo-random value with a mean of 0 and
  !> a variance of 1.
  !>
  !> state: The state variables for the ziggurat algorithm.
  !> return: A normally distributed pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_rnor (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    ! Local parameters.
    real(kind = real32), parameter :: r = 3.442620      ! Start of the right tail.
    real(kind = real32), parameter :: rinv = 0.2904764  ! 0.2904764 is 1/r.
    ! Local variables.
    real(kind = real32) :: x
    real(kind = real32) :: y
    real(kind = real32) :: z
    state%hz = ziggurat_shr3(state)
    state%iz = iand (state%hz, 127)
    if (abs (state%hz) .lt. state%kn(state%iz)) then
      return_value = state%hz * state%wn(state%iz)
    else
      ! RNOR rejection occurred, generate variates from the residue.
      do
        x = state%hz * state%wn(state%iz)
        ! Handle the base strip.
        if (state%iz .eq. 0) then
          do
            x = -log (ziggurat_uni(state)) * rinv
            y = -log (ziggurat_uni(state))
            if (y + y .ge. x * x) exit
          end do
          if (state%hz .gt. 0) then
            return_value = r + x
          else
            return_value = -r - x
          end if
          return
        end if
        ! Handle the wedges of other strips.
        z = state%fn(state%iz) + ziggurat_uni(state) * &
          (state%fn(state%iz - 1) - state%fn(state%iz))
        if (z .lt. exp (-0.5 * x * x)) then
          return_value = x
          return
        end if
        ! Try to exit do loop.
        state%hz = ziggurat_shr3(state)
        state%iz = iand (state%hz, 127)
        if (abs (state%hz) .lt. state%kn(state%iz)) then
          return_value = state%hz * state%wn(state%iz)
          return
        end if
      end do
    end if
    return
  end function ziggurat_rnor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Generates EXPONENTIALLY DISTRIBUTED pseudo-random values.
  !>
  !> state: The state variables for the ziggurat algorithm.
  !> return: A exponential pseudo-random value.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ziggurat_rexp (state) result (return_value)
    type(ziggurat_t), intent(inout) :: state
    real(kind = real32)             :: return_value
    ! Local variables.
    real(kind = real32) :: x
    real(kind = real32) :: y
    state%jz = ziggurat_shr3(state)
    state%iz = iand (state%jz, 255)
    if (abs (state%jz) .lt. state%ke(state%iz)) then
      x = state%jz * real (state%we(state%iz), kind = real32)
      return_value = real (x, kind = real32)
    else
      ! REXP rejection occurred, generate variates from the residue.
      do
        ! Handles the base strip.
        if (state%iz .eq. 0) then
          return_value = 7.69711 - log (ziggurat_uni(state))
          exit
        end if
        ! Handle the wedges of other strips.
        x = state%jz * real (state%we(state%iz), kind = real32)
        y = state%fe(state%iz) + ziggurat_uni(state) * &
          (state%fe(state%iz - 1) - state%fe(state%iz))
        if (y .lt. exp (-x)) then
          return_value = real (x, kind = real32)
          exit
        end if
        ! Try to exit do loop.
        state%jz = ziggurat_shr3(state)
        state%iz = iand (state%jz, 255)
        if (state%jz .lt. state%ke(state%iz)) then
          x = state%jz * real (state%we(state%iz), kind = real32)
          return_value = real (x, kind = real32)
          exit
        end if
      end do
    end if
    return
  end function ziggurat_rexp

end module ziggurat




module util_m

    use constants_m
    use f95_precision
    use blas95
    use lapack95

    public:: fact , binomial , det , Frobenius_norm
    public:: TO_UPPER_CASE, count_lines, split_line, seq_in_range, read_general_file, read_CSV_file, read_file_name, renumber_sequence

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
!
!=================================
function det( A ) result(det_of_A)
!=================================
implicit none 
real*8  , intent(in) :: A(:,:)

!local variables ...
integer :: i
real*8  :: det_of_A
real*8  , allocatable :: aux(:,:)

! calculates the determinant of the matrix A
allocate( aux , source=A )
CALL getrf(aux)

det_of_A = 1.d0
do i = 1 , size(aux(:,1))
    det_of_A = det_of_A * aux(i,i) 
    end do

deallocate( aux )

end function det
!
!
!
!==========================================
function Frobenius_norm( A ) result(F_norm)
!==========================================
implicit none 
real*8  , intent(in) :: A(:,:)

!local variables ...
integer :: i , m
real*8  :: F_norm
real*8  , allocatable :: aux(:,:)

! calculates the Frobenius norm of the matrix A
m = size(A(:,1))
allocate(aux(m,m))

call gemm( A , A , aux , 'N','T')

F_norm = 0.d0
do i = 1 , m
   F_norm = F_norm + aux(i,i)
   end do

deallocate(aux)

end function Frobenius_norm
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
INTEGER :: n_of_columns
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
character(len=9) , allocatable :: tokens(:) , aux(:)

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
!
!====================================
subroutine read_general_file( dados )
!====================================
implicit none
real*8 , allocatable , intent(out) :: dados(:,:)

! local variables ...
integer :: i , N_of_lines , N_of_columns
real*8  :: allofthem(72)
logical :: exist
logical :: done 
integer , allocatable :: columns(:)
character(len=72) :: f_name , command , line

CALL system("clear")

done = .false.
do while (.NOT. done) 
        write(*,'(a)',advance='no') 'name of the file:  '
        read (*,*) f_name
        
        inquire(file=f_name, EXIST=exist)
        if( .NOT. exist ) then
            print*, ' >>> this file  was  not  found !' 
        else
            n_of_lines   = count_lines(f_name)
            n_of_columns = count_columns(f_name)
            done = .true. 
        end if
end do

command = "head "//f_name
CALL system(command)

! choose columns to be read ...
write(*,'(/3a)') 'choose columns to read (2 column numbers separated by space, press Enter at the end):  '
read (*,'(a)') line
allocate( columns(2) )
columns =  parse_this(line)

allocate( dados(N_of_lines,2) )

OPEN( unit=3 , file=f_name , status='old' , action="read" )

do i = 1 , n_of_lines
   read(3,*) allofthem(:n_of_columns)
   dados(i,:) = allofthem(columns(:))
end do

close(3)

end subroutine read_general_file 
!
!
!
!=============================================
subroutine read_file_name( f_name , file_type)
!=============================================
 implicit none
 character(len=30)            , intent(out) :: f_name
 character(len=*) ,  optional , intent(in)  :: file_type

! local variables ...
logical :: exist
logical :: done 

CALL system("clear")

if( present(file_type) ) &
then
    select case(file_type)
           case("pdb")
               write(*,'(/,a)') "ls *.pdb"
               call system("ls *.pdb") 
           case("xyz")
               write(*,'(/,a)') "ls *.xyz"
               call system("ls *.xyz") 
    end select
end if

done = .false.
do while (.NOT. done) 
        write(*,'(/,a)',advance='no') 'name of the file:  '
        read (*,*) f_name
        
        inquire(file=f_name, EXIST=exist)
        if( .NOT. exist ) then
            print*, ' >>> file  not  found.  Try again !' 
        else
            done = .true. 
        end if
end do

end subroutine read_file_name
!
!
!
!=================================================
function renumber_sequence( seq ) result( new_seq)
!=================================================
    implicit none
    integer, intent(in) :: seq(:)

    ! local variables ...
    integer, allocatable :: new_seq(:)
    integer :: i, current_val, group_id, n

    n = size(seq)
    allocate(new_seq(n))

    group_id = 1
    new_seq(1) = group_id
    current_val = seq(1)

    do i = 2, n
        if (seq(i) /= current_val) then
            group_id = group_id + 1
            current_val = seq(i)
        end if
        new_seq(i) = group_id
    end do

end function renumber_sequence
!
!
!
end module util_m


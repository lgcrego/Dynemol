module util_m

    use constants_m
    use f95_precision
    use blas95
    use lapack95

    public:: fact , binomial , det , Frobenius_norm
    public:: TO_UPPER_CASE , count_lines , split_line , seq_in_range

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
 FUNCTION count_lines (f_unit) result(n_of_lines)
!================================================
implicit none
integer          , intent(in) :: f_unit

! Local variables ...
INTEGER      :: i , ioerr , n_of_lines
character(1) :: dumb

i = 0
do 
    read(f_unit,*,IOSTAT=ioerr) dumb
    if(ioerr < 0) EXIT
    i = i + 1
end do    

n_of_lines = i

end FUNCTION count_lines
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

        colon = index(string,":")
        read(string(1:colon-1),*) i1
        read(string(colon+1: ),*) i2

        seq = [ ( i , i=i1,i2 ) ]

end function seq_in_range
!
!
!
end module util_m


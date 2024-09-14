module util_m

    use constants_m

    public:: fact , binomial , TO_UPPER_CASE , count_lines , split_line

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
character(len=3) , allocatable :: tokens(:) , aux(:)

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
end module util_m


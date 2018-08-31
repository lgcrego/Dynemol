module util_m

    use constants_m

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
end module util_m


 module Slater_Type_Orbitals

    use type_m
    use Semi_Empirical_Parms
    use Structure_Builder

    public :: STO_parameters
    public :: s_orb , p_orb , d_x2y2 , d_z2 , d_xyz  

    private

!   sto = 2^(n+1/2) / sqrt((2*n)!) 
    real*8  , parameter , dimension(5) :: sto = [2.0d0 , 1.1547005d0 , 0.42163702d0 , 0.112687234d0 , 0.23756555d0] 

 contains   
!
!
!
!
!========================
 pure function s_orb(r,k)
!========================
 implicit none
 real*8  , intent(in) :: r
 integer , intent(in) :: k

 real*8  :: a , b , coef , s_orb
 integer :: n , i

 real*8 , parameter :: harm = 0.282094792d0    !! <== normalization constant for ang. orbitals

 n = atom(k)%Nquant(0)

 b = 1.d0 
 do i=1,n-1 
    b=b*r 
 end do                         !! <==   b = r^(n-1)

 a = 1.d0
 do i=1,n 
    a=a*atom(k)%zeta(0,1)
 end do 
 a=a*dsqrt(atom(k)%zeta(0,1))   !! <==   a = eta^(n+1/2)

 coef  = atom(k)%coef(0,1) * sto(n) * a * harm  
 s_orb = coef * b * dexp(-atom(k)%zeta(0,1)*r)

 end function s_orb
!
!
!
!==========================
 pure function p_orb(r,x,k)
!==========================
 implicit none
 real*8  , intent(in) :: r , x
 integer , intent(in) :: k

 real*8  :: a , b , coef , p_orb
 integer :: i , n 

 real*8 , parameter :: harm = 0.488602512d0    !! <== normalization constant for ang. orbitals     

 n = atom(k)%Nquant(1) 

 b = 1.d0 
 do i=1,n-2 
    b=b*r 
 end do                           !! <==   b = r^(n-2)

 a = 1.d0
 do i=1,n 
    a=a*atom(k)%zeta(1,1) 
 end do 
 a=a*dsqrt(atom(k)%zeta(1,1))     !! <==   a = eta^(n+1/2)

 coef  = atom(k)%coef(1,1) * sto(n) * a * harm
 p_orb = coef * b * x * dexp(-atom(k)%zeta(1,1)*r)

 end function p_orb
!
!
!
!===============================
 pure function d_x2y2(r,x1,x2,k)
!===============================
 implicit none
 real*8  , intent(in) :: r , x1 , x2
 integer , intent(in) :: k

 integer            :: n , i  
 real*8             :: d_x2y2
 real*8             :: b , a1 , a2 , orb1 , orb2 , coef
 real*8 , parameter :: harm = -0.546274215d0   !! <== normalization constant for ang. orbitals 

 n = atom(k)%Nquant(2)

 b = 1.d0 
 do i=1,n-3
    b=b*r 
 end do                               !! <--   b = r^(n-3)

 a1 = 1.d0
 do i=1,n 
    a1=a1*atom(k)%zeta(2,1) 
 end do 
 a1=a1*dsqrt(atom(k)%zeta(2,1))       !! <--  eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
    a2=a2*atom(k)%zeta(2,2) 
 end do 
 a2=a2*dsqrt(atom(k)%zeta(2,2))       !! <--  eta_2^(n+1/2)
                                    
 coef = sto(n) * harm 
 orb1 = a1 * atom(k)%coef(2,1) * exp(-atom(k)%zeta(2,1)*r)
 orb2 = a2 * atom(k)%coef(2,2) * exp(-atom(k)%zeta(2,2)*r) 

 d_x2y2 = coef * (-x1*x1 + x2*x2) * b * ( orb1 + orb2 )

 end function d_x2y2
!
!
!
!================================
 pure function d_z2(r,x1,x2,x3,k)
!================================
 implicit none
 real*8  , intent(in) :: r , x1 , x2 , x3
 integer , intent(in) :: k

 integer            :: n , i
 real*8             :: d_z2 , b , a1 , a2 , coef , orb1 , orb2
 real*8 , parameter :: harm = -0.315391565d0   !! <== normalization constant for ang. orbitals 
 real*8 , parameter :: two = 2.d0  

 n = atom(k)%Nquant(2)

 b = 1.d0 
 do i=1,n-3 
    b=b*r 
 end do                               !! <--   b = r^(n-3)

 a1 = 1.d0
 do i=1,n 
    a1=a1*atom(k)%zeta(2,1) 
 end do 
 a1=a1*dsqrt(atom(k)%zeta(2,1))       !! <--  eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
     a2=a2*atom(k)%zeta(2,2) 
 end do 
 a2=a2*dsqrt(atom(k)%zeta(2,2))       !! <--  eta_2^(n+1/2)

 coef = sto(n) * harm
 orb1 = a1 * atom(k)%coef(2,1) * dexp(-atom(k)%zeta(2,1)*r)
 orb2 = a2 * atom(k)%coef(2,2) * dexp(-atom(k)%zeta(2,2)*r)

 d_z2 = coef * (-two*x3*x3 + x1*x1 + x2*x2) * b * ( orb1 + orb2 )

 end function d_z2
!
!
!
!==============================
 pure function d_xyz(r,x1,x2,k)
!==============================
 implicit none
 real*8  , intent(in) :: r , x1 , x2
 integer , intent(in) :: k

 integer            :: i , n 
 real*8             :: d_xyz , b , a1 , a2 , orb1 , orb2 , coef
 real*8 , parameter ::  harm = 1.092548431d0    !! <== normalization constant for ang. orbitals 

 n = atom(k)%Nquant(2)
 
 b = 1.d0 
 do i=1,n-3 
    b=b*r 
 end do                               !! <--   b = r^(n-3)
 
 a1 = 1.d0
 do i=1,n 
    a1=a1*atom(k)%zeta(2,1) 
 end do 
 a1=a1*dsqrt(atom(k)%zeta(2,1))       !! <--  eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
    a2=a2*atom(k)%zeta(2,2) 
 end do 
 a2=a2*dsqrt(atom(k)%zeta(2,2))       !! <--  eta_2^(n+1/2)

 coef = sto(n) * harm
 orb1 = a1 * atom(k)%coef(2,1) * exp(-atom(k)%zeta(2,1)*r)
 orb2 = a2 * atom(k)%coef(2,2) * exp(-atom(k)%zeta(2,2)*r)

 d_xyz = coef * x1*x2 * b * ( orb1 + orb2 ) 

 end function d_xyz
!
!
end module Slater_Type_Orbitals

 module Slater_Type_Orbitals

    use type_m
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : Cube_Coef , Cube_Zeta

    public :: s_orb , p_orb , d_x2y2 , d_z2 , d_xyz  

    private

!   sto = 2^(n+1/2) / sqrt((2*n)!) 
    real*8  , parameter , dimension(6) :: sto = [2.0d0 , 1.1547005d0 , 0.42163702d0 , 0.112687234d0 , 0.023756555d0 , 0.004135485d0] 

 contains   
!
!
!
!
!=============================
 pure function s_orb(r,AtNo,k)
!=============================
 implicit none
 real*8  , intent(in) :: r
 integer , intent(in) :: AtNo 
 integer , intent(in) :: k 

 real*8  :: a1 , a2 , b , coef , orb_1 , orb_2 , s_orb
 integer :: n , i

 real*8 , parameter :: harm = 0.282094792d0    !! <== normalization constant for ang. orbitals

 n = atom(AtNo)%Nquant(0)

 b = 1.d0 
 do i=1,n-1 
    b=b*r 
 end do                             !! <==   b = r^(n-1)

 a1 = 1.d0
 do i=1,n 
    a1=a1*Cube_Zeta(k,1)
 end do 
 a1=a1*dsqrt(Cube_Zeta(k,1))        !! <==   a1 = eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
    a2=a2*Cube_Zeta(k,2)
 end do 
 a2=a2*dsqrt(Cube_Zeta(k,2))        !! <==   a2 = eta_2^(n+1/2)

 coef  = sto(n) * harm  
 orb_1 = Cube_Coef(k,1) * a1 * dexp(-Cube_Zeta(k,1)*r)
 orb_2 = Cube_Coef(k,2) * a2 * dexp(-Cube_Zeta(k,2)*r)

 s_orb = coef * b * ( orb_1 + orb_2 )

 end function s_orb
!
!
!
!===============================
 pure function p_orb(r,x,AtNo,k)
!===============================
 implicit none
 real*8  , intent(in) :: r , x
 integer , intent(in) :: AtNo 
 integer , intent(in) :: k

 real*8  :: a1 , a2 , b , coef , orb_1 , orb_2 , p_orb
 integer :: i , n 

 real*8 , parameter :: harm = 0.488602512d0    !! <== normalization constant for ang. orbitals     

 n = atom(AtNo)%Nquant(1) 

 b = 1.d0 
 do i=1,n-2 
    b=b*r 
 end do                             !! <==   b = r^(n-2)

 a1 = 1.d0
 do i=1,n 
    a1=a1*Cube_Zeta(k,1) 
 end do 
 a1=a1*dsqrt(Cube_Zeta(k,1))        !! <==   a1 = eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
    a2=a2*Cube_Zeta(k,2) 
 end do 
 a2=a2*dsqrt(Cube_Zeta(k,2))        !! <==   a2 = eta_2^(n+1/2)

 coef  = sto(n) * harm
 orb_1 = Cube_Coef(k,1) * a1 * dexp(-Cube_Zeta(k,1)*r)
 orb_2 = Cube_Coef(k,2) * a2 * dexp(-Cube_Zeta(k,2)*r)

 p_orb = coef * b * x * ( orb_1 + orb_2 )

 end function p_orb
!
!
!
!====================================
 pure function d_x2y2(r,x1,x2,AtNo,k)
!====================================
 implicit none
 real*8  , intent(in) :: r , x1 , x2
 integer , intent(in) :: AtNo 
 integer , intent(in) :: k

 integer            :: n , i  
 real*8             :: d_x2y2
 real*8             :: b , a1 , a2 , orb1 , orb2 , coef
 real*8 , parameter :: harm = -0.546274215d0   !! <== normalization constant for ang. orbitals 

 n = atom(AtNo)%Nquant(2)

 b = 1.d0 
 do i=1,n-3
    b=b*r 
 end do                                 !! <--   b = r^(n-3)

 a1 = 1.d0
 do i=1,n 
    a1=a1*Cube_Zeta(k,1) 
 end do 
 a1=a1*dsqrt(Cube_Zeta(k,1))            !! <--  eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
    a2=a2*Cube_Zeta(k,2) 
 end do 
 a2=a2*dsqrt(Cube_Zeta(k,2))            !! <--  eta_2^(n+1/2)
                                    
 coef = sto(n) * harm 
 orb1 = a1 * Cube_Coef(k,1) * exp(-Cube_Zeta(k,1)*r)
 orb2 = a2 * Cube_Coef(k,2) * exp(-Cube_Zeta(k,2)*r) 

 d_x2y2 = coef * (-x1*x1 + x2*x2) * b * ( orb1 + orb2 )

 end function d_x2y2
!
!
!
!=====================================
 pure function d_z2(r,x1,x2,x3,AtNo,k)
!=====================================
 implicit none
 real*8  , intent(in) :: r , x1 , x2 , x3
 integer , intent(in) :: AtNo 
 integer , intent(in) :: k

 integer            :: n , i
 real*8             :: d_z2 , b , a1 , a2 , coef , orb1 , orb2
 real*8 , parameter :: harm = -0.315391565d0   !! <== normalization constant for ang. orbitals 
 real*8 , parameter :: two = 2.d0  

 n = atom(AtNo)%Nquant(2)

 b = 1.d0 
 do i=1,n-3 
    b=b*r 
 end do                                 !! <--   b = r^(n-3)

 a1 = 1.d0
 do i=1,n 
    a1=a1*Cube_Zeta(k,1) 
 end do 
 a1=a1*dsqrt(Cube_Zeta(k,1))            !! <--  eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
     a2=a2*Cube_Zeta(k,2) 
 end do 
 a2=a2*dsqrt(Cube_Zeta(k,2))            !! <--  eta_2^(n+1/2)

 coef = sto(n) * harm
 orb1 = a1 * Cube_Coef(k,1) * dexp(-Cube_Zeta(k,1)*r)
 orb2 = a2 * Cube_Coef(k,2) * dexp(-Cube_Zeta(k,2)*r)

 d_z2 = coef * (-two*x3*x3 + x1*x1 + x2*x2) * b * ( orb1 + orb2 )

 end function d_z2
!
!
!
!===================================
 pure function d_xyz(r,x1,x2,AtNo,k)
!===================================
 implicit none
 real*8  , intent(in) :: r , x1 , x2
 integer , intent(in) :: AtNo 
 integer , intent(in) :: k

 integer            :: i , n 
 real*8             :: d_xyz , b , a1 , a2 , orb1 , orb2 , coef
 real*8 , parameter ::  harm = 1.092548431d0    !! <== normalization constant for ang. orbitals 

 n = atom(AtNo)%Nquant(2)
 
 b = 1.d0 
 do i=1,n-3 
    b=b*r 
 end do                                 !! <--   b = r^(n-3)
 
 a1 = 1.d0
 do i=1,n 
    a1=a1*Cube_Zeta(k,1) 
 end do 
 a1=a1*dsqrt(Cube_Zeta(k,1))            !! <--  eta_1^(n+1/2)

 a2 = 1.d0
 do i=1,n 
    a2=a2*Cube_Zeta(k,2) 
 end do 
 a2=a2*dsqrt(Cube_Zeta(k,2))            !! <--  eta_2^(n+1/2)

 coef = sto(n) * harm
 orb1 = a1 * Cube_Coef(k,1) * exp(-Cube_Zeta(k,1)*r)
 orb2 = a2 * Cube_Coef(k,2) * exp(-Cube_Zeta(k,2)*r)

 d_xyz = coef * x1*x2 * b * ( orb1 + orb2 ) 

 end function d_xyz
!
!
end module Slater_Type_Orbitals

module FF_angles

    use type_m   
    use omp_lib
    use constants_m
    use parameters_m , only: PBC 
    use for_force    , only: angpot 
    use MD_read_m    , only: atom , molecule , MM 

    private

    public :: f_angle

contains
!
!
!==================
subroutine f_angle
!==================
implicit none

! local_variables ...
real*8 , dimension (3):: rij, rjk, rik 
real*8  :: rijq, rjkq, rijsq, rjksq, fxyz, riju, riku, rikq, riksq
real*8  :: coephi, coephi0,  rterm, rterm0, phi
integer :: i, j, k, l, n, ati, atj, atk, atl, loop  


do j = 1 , MM % N_of_atoms
   atom(j) % fang(:) = D_zero  ! Bending/Angular
end do
angpot = D_zero

!================================
! Angle - bending potential ...
!             J     K
!              \   /
!               \ / 
!                I 
!================================

do i = 1 , MM % N_of_molecules
    do j = 1 , molecule(i) % Nangs

        ! MIND: the atomic sequence is JIK, with ATOM I IN THE VERTEX
        atj = molecule(i) % angs(j,1)
        ati = molecule(i) % angs(j,2)
        atk = molecule(i) % angs(j,3)

        if( atom(atj) % flex .OR. atom(ati) % flex .OR. atom(atk) % flex ) then

            ! MIND: the atomic sequence is JIK, with ATOM I IN THE VERTEX
            ! rij = r_j - r_i
            rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rijq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq  = SQRT(rijq)
            ! rik = r_k - r_i 
            rik(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
            rik(:) = rik(:) - MM % box(:)*DNINT( rik(:) * MM % ibox(:) ) * PBC(:)
            rikq   = rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3)
            riksq  = SQRT(rikq)

            phi = ACOS( (rij(1)*rik(1) + rij(2)*rik(2) + rij(3)*rik(3) ) / ( rijsq * riksq ) )

            select case ( molecule(i) % angle_type(j) )
           
                case( "harm" , "urba" )
                ! Harmonic and Urey-Bradley potentials ...
 
                coephi = 0.d0
                if (phi < 1.d-12 .OR. abs(pi - phi) < 1.d-12) then
                    coephi = 0.d0
                    rterm  = 0.0d0 
                else 
                    coephi = ( phi - molecule(i) % kang0(j,2) ) / SIN(phi)
                    rterm  = HALF * molecule(i) % kang0(j,1) * ( phi - molecule(i) % kang0(j,2) )**2
                end if
                angpot = angpot + rterm
         
                do l = 1, 3
                    if (l == 1) atl = atj
                    if (l == 2) atl = ati
                    if (l == 3) atl = atk
                    do loop = 1, 3                    ! X,Y,Z axis (n = 1, 2 or 3)
                       riju = rij(loop)
                       riku = rik(loop)
                       fxyz = ( molecule(i) % kang0(j,1) * coephi ) *                     &
                              ( (DEL(atl,atj)-DEL(atl,ati))*riku/(rijsq*riksq) +          &
                              (DEL(atl,atk)-DEL(atl,ati))*riju/(rijsq*riksq) -            &
                              COS(phi)*( (DEL(atl,atj)-DEL(atl,ati))*riju/(rijsq*rijsq) + &
                              (DEL(atl,atk)-DEL(atl,ati))*riku/(riksq*riksq)) )

                       atom(atl) % fang(loop) = atom(atl) % fang(loop) + fxyz

                    end do
                end do

                ! Urey-Bradley bonding ...
                if( molecule(i) % angle_type(j) == "urba" ) then
                    rjk(:) = atom(atk) % xyz(:) - atom(atj) % xyz(:)
                    rjk(:) = rjk(:) - MM % box(:) * DNINT( rjk(:) * MM % ibox(:) ) * PBC(:)
                    rjkq   = rjk(1)*rjk(1) + rjk(2)*rjk(2) + rjk(3)*rjk(3)
                    rjksq  = SQRT(rjkq)

                    coephi0 = molecule(i) % kang0(j,3) * ( rjksq - molecule(i) % kang0(j,4) )/rjksq
                    rterm0  = HALF*molecule(i)%kang0(j,3) * ( rjksq - molecule(i)%kang0(j,4) )**2

                    atom(atk) % fang(:) = atom(atk) % fang(:) - coephi0*rjk(:)
                    atom(atj) % fang(:) = atom(atj) % fang(:) + coephi0*rjk(:)  

                    angpot = angpot + rterm0

                end if
              
            end select 
        end if
    end do
end do

end subroutine f_angle
!                                                                                                                                                                                   
!                                                                                                                                                                                   
!                                                                                                                                                                                   
!====================                                                                                                                                                               
 function DEL ( X,Y )                                                                                                                                                               
!====================                                                                                                                                                               
 real*8  :: DEL                                                                                                                                                                     
 integer :: X, Y                                                                                                                                                                    
                                                                                                                                                                                    
 if (X == Y) then                                                                                                                                                                   
   DEL = 1.                                                                                                                                                                         
 else if (X /= Y) then                                                                                                                                                              
   DEL = 0.                                                                                                                                                                         
 end if                                                                                                                                                                             
                                                                                                                                                                                    
end function DEL   
!
!
end module FF_angles

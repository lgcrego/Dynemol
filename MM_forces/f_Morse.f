module FF_Morse

    use type_m   
    use omp_lib
    use constants_m
    use parameters_m , only: PBC 
    use for_force    , only: Morspot
    use MD_read_m    , only: atom , molecule , MM 
    use gmx2mdflex   , only: MorsePairs

    private

    public :: f_Morse

contains
!
!
!==================
subroutine f_Morse
!==================
implicit none

! local_variables ...
real*8 , dimension (3):: rij 
real*8  :: rij2 
real*8  :: MorsA, MorsB, MorsC, dij
real*8  :: coephi, qterm, qterm0 
integer :: i, j, n, ati, atj 
logical :: flag1 , flag2


do j = 1 , MM % N_of_atoms
   atom(j) % fMorse(:) = D_zero 
end do
Morspot = D_zero

!====================================================================
! Morse Intra/Inter potential for H transfer ...

If( allocated(MorsePairs) ) then

   do i = 1 , MM % N_of_atoms - 1
       do j = i+1 , MM % N_of_atoms
       read_loop2: do  n = 1, size(MorsePairs) 
           flag1 = ( adjustl( MorsePairs(n) % MMSymbols(1) ) == adjustl( atom(i) % MMSymbol ) ) .AND. &
                   ( adjustl( MorsePairs(n) % MMSymbols(2) ) == adjustl( atom(j) % MMSymbol ) )
           flag2 = ( adjustl( MorsePairs(n) % MMSymbols(2) ) == adjustl( atom(i) % MMSymbol ) ) .AND. &
                   ( adjustl( MorsePairs(n) % MMSymbols(1) ) == adjustl( atom(j) % MMSymbol ) ) 
           if ( flag1 .OR. flag2 ) then
               ati = atom(i) % my_id
               atj = atom(j) % my_id
               rij(:) = atom(ati) % xyz(:) - atom(atj) % xyz(:)
               rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
               rij2   = sum( rij(:)**2 )
               dij    = sqrt(rij2)

               MorsA = MorsePairs(n) % Parms(1)
               MorsB = MorsePairs(n) % Parms(2)
               MorsC = MorsePairs(n) % Parms(3)
 
               ! Morse potential ...
               qterm0 = exp( - MorsC*(dij-MorsB) )
               qterm  = MorsA*( d_one - qterm0 )**2
               coephi = TWO * MorsA*MorsC * qterm0 * (d_one - qterm0)

               atom(ati) % fMorse(:) = atom(ati) % fMorse(:) - coephi*rij(:)/dij
               atom(atj) % fMorse(:) = atom(atj) % fMorse(:) + coephi*rij(:)/dij

               Morspot = Morspot + qterm*factor3 

           end if
       end do read_loop2
       end do
   end do

end If

end subroutine f_Morse
!
!
end module FF_Morse

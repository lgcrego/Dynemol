module FF_bonds

    use type_m   
    use omp_lib
    use constants_m
    use parameters_m , only: PBC    
    use for_force    , only: bdpot, harm_bond, morse_bond
    use MD_read_m    , only: atom , molecule , MM 

    private

    public :: f_bond

contains
!
!
!================
subroutine f_bond
!================
implicit none

! local_variables ...
real*8 , allocatable  :: tmp_force(:,:,:)
real*8 , dimension (3):: rij  
real*8  :: rijq, rijsq 
real*8  :: MorsA, MorsB, MorsC
real*8  :: coephi, qterm, qterm0 
integer :: i, j, k, ithr, numthr, ati, atj 

do j = 1 , MM % N_of_atoms
    atom(j) % fbond(:)  = D_zero    ! Stretching/Bonding 
    atom(j) % fMorse(:) = D_zero    ! Non-bonded Morse
end do

numthr = OMP_get_max_threads()
allocate( tmp_force (MM%N_of_atoms,3,numthr) , source = D_zero )

bdpot      = D_zero
harm_bond  = D_zero
morse_bond = D_zero

! Bonding - stretching potential ...

!$OMP parallel DO &
!$OMP default (shared) &
!$OMP private (i , j , ati , atj , rij , rijq , rijsq , qterm , coephi , qterm0 , ithr)  &
!$OMP reduction (+: harm_bond , morse_bond , bdpot)
do i = 1 , MM % N_of_molecules
    do j = 1 ,  molecule(i) % Nbonds

        ithr = OMP_get_thread_num() + 1

        ati = molecule(i) % bonds(j,1)
        atj = molecule(i) % bonds(j,2)

        if( atom(atj)% flex .OR. atom(ati)% flex ) then 

            rij(:)  = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:)  = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rijq    = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq   = SQRT(rijq)

            select case ( molecule(i) % bond_type(j) )

                case ( "harm" )
                ! harmonic potential ...
                qterm   = HALF*molecule(i)%kbond0(j,1) * (rijsq - molecule(i)%kbond0(j,2))**2
                coephi  = molecule(i)%kbond0(j,1)*( rijsq - molecule(i)%kbond0(j,2) ) / rijsq
                harm_bond = harm_bond + qterm

                case ( "Mors" )
                ! Morse potential ...
                MorsA = molecule(i)% kbond0(j,1)
                MorsB = molecule(i)% kbond0(j,2)
                MorsC = molecule(i)% kbond0(j,3)
 
                qterm0 = exp( - MorsC*(rijsq - MorsB) ) 
                qterm  = MorsA*(d_one - qterm0)**2
                coephi = TWO * MorsA*MorsC * qterm0 * (d_one - qterm0) / rijsq 
                morse_bond = morse_bond + qterm
 
            end select

            tmp_force(atj,1:3,ithr) = tmp_force(atj,1:3,ithr) - coephi*rij(:)
            tmp_force(ati,1:3,ithr) = tmp_force(ati,1:3,ithr) + coephi*rij(:)

            bdpot = bdpot + qterm

        end if 
    end do
end do 
!$OMP end parallel do

! force units = J/mts = Newtons ...    
do i = 1, MM % N_of_atoms
   do k=1,numthr              
      atom(i)% fbond(1:3) = atom(i)% fbond(1:3) + tmp_force(i,1:3,k)
   enddo
end do    
deallocate( tmp_force )

end subroutine f_bond
!
!
end module FF_bonds

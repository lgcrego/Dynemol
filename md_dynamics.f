module MD_dynamics_m

    use constants_m
    use syst        ! using all syst
    use MD_read_m   , only: MM , atom , molecule, species
    use for_force   , only: forcefield
    use f_inter_m   , only: stressr , stresre

    public :: VV1 , VV2 , SUMMAT , PRESS_BOUNDARY

    ! local variables ...
    real*8  :: Astres(3,3) , stresvv(3,3)

contains
!
!
!===================
subroutine VV1( dt )
!===================
implicit none
real*8  , intent(in)    :: dt

! local variables ...
real*8  , dimension (3):: ai, newpos, oldpos, distance
real*8                  :: massa , dt_half
integer                 :: i

 dt_half = dt / two

! VV1 ... 
 do i = 1 , MM % N_of_atoms
    if ( atom(i) % free ) then
        oldpos(1:3) = atom(i) % xyz(1:3)
        massa = mol / atom(i) % mass 
        ai(1:3) = atom(i) % ftotal(1:3) * massa
        atom(i) % xyz(1:3) = atom(i) % xyz(1:3) + atom(i) % vel(1:3) * 1.0d10 * dt + dt * dt * ai(1:3) * 0.5d10
        atom(i) % vel(1:3) = atom(i) % vel(1:3) + dt_half * ai(1:3)
        newpos(1:3) = atom(i) % xyz(1:3)
        distance(1:3) = newpos(1:3) - oldpos(1:3)
    end if
 end do

end subroutine VV1
!
!
!
!====================================
 subroutine VV2( Ttrans , step , dt )
!====================================
use atomicmass
implicit none
real*8     , intent(inout) :: Ttrans
integer    , intent(in)    :: step
real*8     , intent(in)    :: dt

! local variables ...
integer :: i, j, j1, j2, nresid
real*8  :: massa, sumtemp, temp, lambda, dt_half
real*8, dimension(3) :: vi

sumtemp = D_zero
stresvv(1:3,1:3) = 0.0d0

! VV2 and thermostat ...
do i = 1 , MM % N_of_molecules
    vi(1:3) = 0.0d0
    nresid = molecule(i) % nr
    j1 = ( nresid - 1 ) * molecule(nresid) % N_of_atoms + 1
    j2 = nresid * molecule(i) % N_of_atoms
    do j = j1 , j2
       if ( atom(j) % free ) then
       massa = atmas( atom(j) % AtNo )
       vi(1:3) = vi(1:3) + massa * atom(j) % vel(1:3)
       end if
    end do   
    massa = molecule(i) % mass
    vi(1:3) = vi(1:3) * imol / massa
    stresvv(1,1) = stresvv(1,1) + massa * vi(1) * vi(1)
    stresvv(2,2) = stresvv(2,2) + massa * vi(2) * vi(2)
    stresvv(3,3) = stresvv(3,3) + massa * vi(3) * vi(3)
    stresvv(1,2) = stresvv(1,2) + massa * vi(1) * vi(2)
    stresvv(1,3) = stresvv(1,3) + massa * vi(1) * vi(3)
    stresvv(2,3) = stresvv(2,3) + massa * vi(2) * vi(3)
    ! 2*kinetic energy of the system ...
    sumtemp = sumtemp + massa * ( vi(1) * vi(1) + vi(2) * vi(2) + vi(3) * vi(3) )
end do
stresvv(2,1) = stresvv(1,2)
stresvv(3,1) = stresvv(1,3)
stresvv(3,2) = stresvv(2,3)

!######################################################
if (sumtemp == 0.d0) then
    temp = 1.0d0
else
    ! instantaneous temperature : E_kin/(3/2*NkB) ... 
    temp = sumtemp * iboltz / real( MM % N_of_molecules )
endif
! Berendsen Thermostat ;  temper = T_ext ...
lambda = (dt / talt) * ( temper / temp - 1.0d0 )
lambda = SQRT(1.0d0 + lambda)

dt_half = dt / two

sumtemp = 0.d0
do i = 1 , MM % N_of_molecules
    vi(1:3) = 0.0d0
    nresid = molecule(i) % nr
    j1 = ( nresid - 1 ) * molecule(nresid) % N_of_atoms + 1
    j2 = nresid * molecule(i) % N_of_atoms
    do j = j1 , j2
        if ( atom(j) % free ) then
            massa = mol / atmas( atom(j) % AtNo )
!            atom(j) % vel(1:3) = atom(j) % vel(1:3) * lambda + ( dt_half * atom(j) % ftotal(1:3) ) * massa
            atom(j) % vel(1:3) = atom(j) % vel(1:3)  + ( dt_half * atom(j) % ftotal(1:3) ) * massa
            vi(1:3) = vi(1:3) + atom(j) % vel(1:3) / massa
        end if
    end do
    massa = molecule(i) % mass
    vi(1:3) = vi(1:3) / massa
    sumtemp = sumtemp + massa * ( sum( vi(:) * vi(:) ) )
end do

! instantaneous temperature of the system after contact with thermostat ...
Ttrans = sumtemp * iboltz / real(MM % N_of_molecules)
Ekin = Ekin + sumtemp
TempToT = TempTot + Ttrans

end subroutine VV2
!
!
!
!===========================
subroutine SUMMAT( density )
!===========================
 implicit none
 real*8 , intent(out) :: density

! local variables ...
 real*8  :: volume, massa
 integer :: i

 volume = MM % box(1) * MM % box(2) * MM % box(3) * 1d-30

 if (forcefield == 1) then
 else
   stresvv(1:3,1:3) = stresvv(1:3,1:3) / (1.0d8 * volume)
   stressr(1:3,1:3) = stressr(1:3,1:3) / (1.0d8 * volume)
   stresre(1:3,1:3) = stresre(1:3,1:3) / (1.0d8 * volume)
   Astres(1:3,1:3) = stresvv(1:3,1:3) + stressr(1:3,1:3) + stresre(1:3,1:3) 
 endif

 massa = 0.0d0
 do i = 1 , MM % N_of_species
    massa = massa + species(i) % N_of_molecules * species(i) % mass
 end do
 
 density = massa / volume
 density = density * 1.0d3 / 1.0d6
 DensTot = DensTot + density

end subroutine SUMMAT
!
!
!
!=========================================
subroutine PRESS_BOUNDARY( pressure , dt )
!=========================================
implicit none
real*8  , intent(inout) :: pressure
real*8  , intent(in)    :: dt

! local variables ...
 integer :: i, j, j1, j2, nresid
 real*8  :: mip

! pressure = instantaneous pressure ;  press = external pressure ... 
 if (forcefield == 1) then
 else
   pressure = ( Astres(1,1)+Astres(2,2)+Astres(3,3) ) / three
 endif

 PressTot = pressure + PressTot

 mip = dt * talp * ( pressure - press )
 mip = (D_one + mip)**(D_one/three)
 MM % box(1:3) = MM % box(1:3) * mip
 
 do i = 1 , MM % N_of_molecules
    nresid = molecule(i) % nr
    j1 = ( nresid - 1 ) * molecule(nresid) % N_of_atoms + 1
    j2 = nresid * molecule(i) % N_of_atoms
    do j = j1 , j2
        if ( atom(j) % free ) then
            atom(j) % xyz(1:3) = atom(j) % xyz(1:3) * mip
            atom(j) % xyz(1:3) = atom(j) % xyz(1:3) - MM % box(1:3) * ANINT( atom(j) % xyz(1:3) * MM % ibox(1:3) )
        end if
    end do
 end do
!
! 
end subroutine PRESS_BOUNDARY
!
!
end module MD_dynamics_m

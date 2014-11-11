module verlet_m

    use constants_m
    use syst        ! using all syst
    use parameters_m , only: PBC 
    use MD_read_m    , only: MM , atom , molecule, species
    use for_force    , only: forcefield
    use f_inter_m    , only: stressr , stresre

    public :: VV1 , VV2 , SUMMAT , PRESS_BOUNDARY

    ! module variables ...
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
real*8  :: ai(3)
real*8  :: massa , dt_half
integer :: i

dt_half = dt / two

! VV1 ... 
do i = 1 , MM % N_of_atoms
    if( atom(i) % flex ) then
        massa = mol / atom(i) % mass 
        ai = atom(i) % ftotal * massa
        atom(i) % xyz = atom(i) % xyz + ( atom(i) % vel*dt + HALF*dt*dt*ai ) * 1.0d10
        atom(i) % vel = atom(i) % vel + dt_half*ai
    end if
end do

end subroutine VV1
!
!
!
!=======================================
 subroutine VV2( Ttrans , kinetic , dt )
!=======================================
implicit none
real*8  , intent(inout) :: Ttrans
real*8  , intent(out)   :: kinetic
real*8  , intent(in)    :: dt

! local variables ...
real*8  :: tmp(3) , V_CM(3) , V_atomic(3)
real*8  :: massa , factor , sumtemp , temp , lambda , dt_half 
integer :: i , j , j1 , j2 , k , nresid , thermostat_type

! select atomic or molecular kinetic energy to calculate the temperature ...
thermostat_type = NINT( float(maxval(molecule%N_of_atoms)) / float(MM%N_of_molecules) )

! VV2 and thermostats ...
sumtemp = D_zero
stresvv = D_zero

select case ( thermostat_type )

    case(0:1)  ! <== molecular ...

        do i = 1 , MM % N_of_molecules
            tmp = D_zero
            nresid = molecule(i) % nr
            j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
            j2 = sum(molecule(1:nresid) % N_of_atoms)
            do j = j1 , j2
                if( atom(j) % flex ) then
                    tmp = tmp + atom(j)%mass * atom(j)%vel
                end if
            end do   

            V_CM = tmp * imol / molecule(i) % mass

            forall(k=1:3) stresvv(k,k) = stresvv(k,k) +  molecule(i) % mass * V_CM(k) * V_CM(k)

            stresvv(1,2) = stresvv(1,2) + molecule(i) % mass * V_CM(1) * V_CM(2)
            stresvv(1,3) = stresvv(1,3) + molecule(i) % mass * V_CM(1) * V_CM(3)
            stresvv(2,3) = stresvv(2,3) + molecule(i) % mass * V_CM(2) * V_CM(3)

            ! 2*kinetic energy (molecular) of the system ...
            sumtemp = sumtemp + molecule(i) % mass * sum( V_CM * V_CM )
        end do

    case (2:)  ! <== atomic ...

        do i = 1 , MM % N_of_atoms 
            If( atom(i) % flex ) then

                V_atomic = atom(i) % vel 

                factor   = imol * atom(i) % mass

                forall( k=1:3) stresvv(k,k) = stresvv(k,k) + factor * V_atomic(k) * V_atomic(k)

                stresvv(1,2) = stresvv(1,2) + factor * V_atomic(1) * V_atomic(2)
                stresvv(1,3) = stresvv(1,3) + factor * V_atomic(1) * V_atomic(3)
                stresvv(2,3) = stresvv(2,3) + factor * V_atomic(2) * V_atomic(3)

                ! 2*kinetic energy (atomic) of the system ...
                sumtemp = sumtemp + factor * sum( V_atomic * V_atomic )

            end if
        end do

end select

stresvv(2,1) = stresvv(1,2)
stresvv(3,1) = stresvv(1,3)
stresvv(3,2) = stresvv(2,3)

!######################################################
if( sumtemp == 0.d0 ) then
    temp = D_ONE
else
    ! instantaneous temperature : E_kin/(3/2*NkB) ... 
    select case (thermostat_type)

        case (0:1) ! <== molecular ...
        temp = sumtemp * iboltz / real( count(molecule%flex) )

        case (2:)  ! <== atomic ...
        temp = sumtemp * iboltz / real( count(atom%flex) )

    end select
endif

! Berendsen Thermostat ; turned off for talt == infty ...
If( talt == infty ) then
    lambda = D_one
else
    lambda = ( dt / (talt*pico_2_sec) ) * ( bath_T / temp - D_ONE )
    lambda = SQRT(D_ONE + lambda)
end If

dt_half = dt / two

sumtemp = 0.d0

select case (thermostat_type)

    case (0:1) ! <== molecular ...

        do i = 1 , MM % N_of_molecules
            tmp = D_zero
            nresid = molecule(i) % nr
            j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
            j2 = sum(molecule(1:nresid) % N_of_atoms)
            do j = j1 , j2
                if ( atom(j) % flex ) then
                    massa = mol / atom(j) % mass
                    atom(j) % vel = atom(j) % vel * lambda + ( dt_half * atom(j) % ftotal ) * massa

                    tmp = tmp + atom(j) % vel / massa

                end if
            end do
            V_CM    = tmp / molecule(i) % mass
            sumtemp = sumtemp + molecule(i) % mass *  sum( V_CM * V_CM ) 
        end do

    case (2:) ! <== atomic ...

        V_atomic = D_zero
        do i = 1 , MM % N_of_atoms
            if( atom(i) % flex ) then
                massa = mol / atom(i) % mass
                atom(i) % vel = atom(i) % vel * lambda + ( dt_half * atom(i) % ftotal ) * massa

                V_atomic = atom(i) % vel 
                factor   = imol * atom(i) % mass
                sumtemp  = sumtemp + factor * sum( V_atomic * V_atomic )

            end if
        end do

end select

! instantaneous temperature of the system after contact with thermostat ...
select case (thermostat_type)
    case (0:1) ! <== molecular ...
    Ttrans =  sumtemp * iboltz / real( count(molecule%flex) ) 

    case (2:)  ! <atomic ...
    Ttrans =  sumtemp * iboltz / real( count(atom%flex) ) 
end select
Ekin    =  Ekin + sumtemp
TempToT =  TempTot + Ttrans


! calculation of the kinetic energy ...
kinetic = 0.d0
do i = 1 , MM % N_of_atoms
    kinetic = kinetic + ( atom(i) % mass ) * sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half
end do
kinetic = kinetic * 1.0d-6 / MM % N_of_Molecules

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
   stresvv = stresvv / (1.0d8 * volume)
   stressr = stressr / (1.0d8 * volume)
   stresre = stresre / (1.0d8 * volume)
   Astres  = stresvv + stressr + stresre
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
    pressure = ( Astres(1,1) + Astres(2,2) + Astres(3,3) ) * third
endif

PressTot = pressure + PressTot

! Pressurestat ; turned off for talp == infty ...
If( talp == infty ) then
    mip = D_one
else
    mip  = dt * ( 107.0d-6 / (talp * pico_2_sec) ) * ( pressure - press )
    mip  = (D_one + mip)**third
end If

MM % box = MM % box * mip
 
do i = 1 , MM % N_of_molecules
    nresid = molecule(i) % nr
    j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
    j2 = sum(molecule(1:nresid) % N_of_atoms)
    do j = j1 , j2
        if ( atom(j) % flex ) then
            atom(j) % xyz(1:3) = atom(j) % xyz(1:3) * mip
            atom(j) % xyz(1:3) = atom(j) % xyz(1:3) - MM % box(1:3) * DNINT( atom(j) % xyz(1:3) * MM % ibox(1:3) ) * PBC(1:3)
        end if
    end do
end do
!
! 
end subroutine PRESS_BOUNDARY
!
!
!
end module verlet_m

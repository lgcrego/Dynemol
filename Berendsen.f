module Berendsen_m

    use constants_m
    use syst                  ! using all syst
    use parameters_m          , only: PBC 
    use MD_read_m             , only: MM , atom , molecule, species
    use f_inter_m             , only: stressr , stresre
    use VV_Parent             , only: VV

    public :: Berendsen , VirialKinetic , barostat

    type, extends(VV) :: Berendsen
    contains
        procedure :: VV1
        procedure :: VV2
    end type  Berendsen

    interface Berendsen
        module procedure  constructor
    end interface

    ! module variables ...
    real*8 :: stresvv(3,3)

contains
!
!
!
!===================================
 function constructor() result( me )
!===================================
implicit none
type(Berendsen) :: me 

!local variable ...

! select atomic or molecular kinetic energy to calculate the temperature ...
me % thermostat_type = NINT( float(maxval(molecule%N_of_atoms)) / float(MM%N_of_molecules) )

end function constructor
!
!
!
!========================
subroutine VV1( me , dt )
!========================
implicit none
class(Berendsen) , intent(inout) :: me
real*8           , intent(in)    :: dt

! local variables ...
real*8  :: ai(3)
real*8  :: massa , dt_HALF , dt2_HALF
integer :: i

dt_HALF  = dt / two
dt2_HALF = dt_HALF * dt

! VV1 ... 
do i = 1 , MM % N_of_atoms
    if( atom(i) % flex ) then
        massa = mol / atom(i) % mass 
        ai = atom(i) % ftotal * massa
        atom(i) % xyz = atom(i) % xyz + ( atom(i) % vel*dt + dt2_HALF*ai ) * mts_2_Angs
        atom(i) % vel = atom(i) % vel + dt_HALF*ai
    end if
end do

end subroutine VV1
!
!
!
!=========================
 subroutine VV2( me , dt )
!=========================
implicit none
class(Berendsen) , intent(inout) :: me
real*8           , intent(in)    :: dt

! local variables ...
integer :: i , j , j1 , j2 , nresid 
real*8  :: massa , factor , sumtemp , temp , lambda , dt_half 
real*8  :: tmp(3) , V_CM(3) , V_atomic(3)

CALL VirialKinetic( me % thermostat_type , sumtemp )

!######################################################
!               Berendsen Thermostat 

if( sumtemp == 0.d0 ) then
    temp = D_ONE
else
    ! instantaneous temperature : E_kin/(3/2*NkB) ... 
    select case (me % thermostat_type)

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

!######################################################

dt_half = dt / two

! VV2 ... 
sumtemp = D_zero
select case (me % thermostat_type)

    case (0:1) ! <== molecular ...

        do i = 1 , MM % N_of_molecules
            tmp = D_zero
            nresid = molecule(i) % nr
            j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
            j2 = sum(molecule(1:nresid) % N_of_atoms)
            do j = j1 , j2
                if ( atom(j) % flex ) then
                    massa = mol / atom(j) % mass
!                    atom(j) % vel = atom(j) % vel * lambda + ( dt_half * atom(j) % ftotal ) * massa
                    atom(j) % vel = atom(j) % vel + ( dt_half * atom(j) % ftotal ) * massa
                    atom(j) % vel = atom(j) % vel * lambda

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
!                atom(i) % vel = atom(i) % vel * lambda + ( dt_half * atom(i) % ftotal ) * massa
                atom(i) % vel = atom(i) % vel + ( dt_half * atom(i) % ftotal ) * massa
                atom(i) % vel = atom(i) % vel * lambda

                V_atomic = atom(i) % vel 
                factor   = imol * atom(i) % mass
                sumtemp  = sumtemp + factor * sum( V_atomic * V_atomic )

            end if
        end do

end select

! instantaneous temperature of the system after contact with thermostat ...
select case (me % thermostat_type)
    case (0:1) ! <== molecular ...
    me % Temperature =  sumtemp * iboltz / real( count(molecule%flex) ) 

    case (2:)  ! <atomic ...
    me % Temperature =  sumtemp * iboltz / real( count(atom%flex) ) 
end select

! calculation of the kinetic energy ...
me % Kinetic = D_zero
do i = 1 , MM % N_of_atoms
    me % Kinetic = me % Kinetic + ( atom(i) % mass ) * sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half
end do
me % Kinetic = me % Kinetic * micro / MM % N_of_Molecules

! calculation of pressure and density ...
CALL barostat( me , dt )

end subroutine VV2
!
!
!
!=====================================================
 subroutine VirialKinetic( thermostat_type , sumtemp )
!=====================================================
implicit none
integer            , intent(in)  :: thermostat_type
real*8  , optional , intent(out) :: sumtemp

! local variables ...
real*8  :: tmp(3) , V_CM(3) , V_atomic(3) , factor , kinetic
integer :: i , j , j1 , j2 , k , nresid 


! VV2 and thermostats ...
kinetic = D_zero
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
            kinetic = kinetic + molecule(i) % mass * sum( V_CM * V_CM )
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
                kinetic = kinetic + factor * sum( V_atomic * V_atomic )

            end if
        end do

end select

stresvv(2,1) = stresvv(1,2)
stresvv(3,1) = stresvv(1,3)
stresvv(3,2) = stresvv(2,3)

If( present(sumtemp) ) sumtemp = kinetic

end subroutine VirialKinetic
!
!
!
!==============================
 subroutine barostat( me , dt )
!==============================
implicit none
class(Berendsen) , intent(inout) :: me
real*8           , intent(in)    :: dt

! local variables ...
integer :: i, j, j1, j2, nresid
real*8  :: mip , volume, massa
real*8  :: Astres(3,3) 

massa   = sum( species(:) % N_of_molecules * species(:) % mass )
volume  = product( MM % box(:) * Angs_2_mts )
me % density = (massa / volume) * milli


stresvv = stresvv / (cm_2_Angs * volume)
stressr = stressr / (cm_2_Angs * volume)
stresre = stresre / (cm_2_Angs * volume)
Astres  = stresvv + stressr + stresre

me % pressure = ( Astres(1,1) + Astres(2,2) + Astres(3,3) ) * third

! Pressurestat ; turned off for talp == infty ...
If( talp == infty ) then
    mip = D_one
else
    mip = dt * ( 107.0d-6 / (talp * pico_2_sec) ) * ( me % pressure - press )
    mip = (D_one + mip)**third
end If

MM % box = MM % box * mip
 
do i = 1 , MM % N_of_molecules
    nresid = molecule(i) % nr
    j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
    j2 = sum(molecule(1:nresid) % N_of_atoms)
    do j = j1 , j2
        if ( atom(j) % flex ) then
            atom(j) % xyz = atom(j) % xyz * mip
            atom(j) % xyz = atom(j) % xyz - MM % box * DNINT( atom(j) % xyz * MM % ibox ) * PBC
        end if
    end do
end do
 
end subroutine barostat
!
!
!
end module Berendsen_m

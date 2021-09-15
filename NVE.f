module NVE_m

    use constants_m
    use syst         ! using all syst
    use MD_read_m    , only: MM , atom , molecule
    use VV_Parent    , only: VV


    public :: NVE

    private

    type, extends(VV) :: NVE
    contains
        procedure :: VV1
        procedure :: VV2
    end type NVE 

    interface NVE
        module procedure  constructor
    end interface

contains
!
!
!
!===================================
 function constructor() result( me )
!===================================
implicit none
type(NVE) :: me 

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
class(NVE) , intent(inout) :: me
real*8     , intent(in)    :: dt

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
class(NVE) , intent(inout) :: me
real*8     , intent(in)    :: dt

! local variables ...
real*8  :: tmp(3) , V_CM(3) , V_atomic(3)
real*8  :: massa , factor , sumtemp , dt_half 
integer :: i , j , j1 , j2 , nresid 

dt_half = dt / two

sumtemp = D_zero

! VV2 ... 
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
                    atom(j) % vel = atom(j) % vel + ( dt_half * atom(j) % ftotal ) * massa

                    tmp = tmp + atom(j) % vel / massa

                end if
            end do
            V_CM    = tmp / molecule(i) % mass
            sumtemp = sumtemp + molecule(i) % mass *  sum( V_CM * V_CM )    ! <== Joule ...
        end do

    case (2:) ! <== atomic ...

        V_atomic = D_zero
        do i = 1 , MM % N_of_atoms
            if( atom(i) % flex ) then
                massa = mol / atom(i) % mass
                atom(i) % vel = atom(i) % vel + ( dt_half * atom(i) % ftotal ) * massa

                V_atomic = atom(i) % vel 
                factor   = imol * atom(i) % mass
                sumtemp  = sumtemp + factor * sum( V_atomic * V_atomic )    ! <== Joule ...

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
me % kinetic = D_zero
do i = 1 , MM % N_of_atoms
    me % kinetic = me % kinetic + ( atom(i) % mass ) * sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half   ! <== J/kmol
end do
me % kinetic = me % kinetic * micro / MM % N_of_Molecules   ! <== kJ/mol

end subroutine VV2
!
!
!
end module NVE_m

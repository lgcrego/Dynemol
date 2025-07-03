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

    ! module parameters ...
    real*8 , parameter :: inv_kilo_mol = imol

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
real*8  :: accelerate(3)
real*8  :: dt_HALF , dt2_HALF
integer :: i

dt_HALF  = dt / two
dt2_HALF = dt_HALF * dt

! VV1 ... 
do i = 1 , MM % N_of_atoms
    if( atom(i) % flex ) then
        accelerate = atom(i)% ftotal / (atom(i)%mass * inv_kilo_mol)
        atom(i) % xyz = atom(i) % xyz + ( atom(i) % vel*dt + dt2_HALF*accelerate ) * mts_2_Angs
        atom(i) % vel = atom(i) % vel + dt_HALF*accelerate
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
real*8  :: E_kinetic , dt_half 
real*8  :: total_Momentum(3) , V_CM(3) , V_atomic(3) , accelerate(3)
integer :: i , j

dt_half = dt / two

E_kinetic = D_zero

! VV2 ... 
select case (me % thermostat_type)

    case (0:1) ! <== molecular ...

        do i = 1 , MM % N_of_molecules
            total_Momentum = D_zero
            do j = molecule(i)%span % inicio , molecule(i)%span % fim
                if ( atom(j) % flex ) then
                    accelerate = atom(j)% ftotal / (atom(j)% mass * inv_kilo_mol)
                    atom(j) % vel = atom(j) % vel + accelerate*dt_half

                    total_Momentum = total_Momentum + (atom(j)% mass * atom(j)%vel) * inv_kilo_mol
                end if
            end do
            V_CM = total_Momentum / molecule(i) % mass
            E_kinetic = E_kinetic + molecule(i) % mass * sum(V_CM*V_CM)    ! <== Joule ...
        end do

    case (2:) ! <== atomic ...
        V_atomic = D_zero
        do i = 1 , MM % N_of_atoms
            if( atom(i) % flex ) then
                accelerate = atom(i)%ftotal / (atom(i)%mass * inv_kilo_mol)
                atom(i) % vel = atom(i) % vel + accelerate*dt_half

                V_atomic = atom(i) % vel 
                E_kinetic = E_kinetic + atom(i)%mass * sum(V_atomic*V_atomic) * inv_kilo_mol    ! <== Joule ...
            end if
        end do

end select

! instantaneous temperature of the system after contact with thermostat ...
select case (me % thermostat_type)
    case (0:1) ! <== molecular ...
    me % Temperature =  E_kinetic * iboltz / real( count(molecule%flex) ) 

    case (2:)  ! <atomic ...
    me % Temperature =  E_kinetic * iboltz / real( count(atom%flex) ) 
end select

! calculation of the kinetic energy ...
me % kinetic = D_zero
do i = 1 , MM % N_of_atoms
    me% kinetic = me% kinetic + ( atom(i)% mass ) * sum(atom(i)% vel(:) * atom(i)% vel(:)) * half   ! <== J/kmol
end do
me % kinetic = me % kinetic * micro / MM % N_of_Molecules   ! <== kJ/mol

end subroutine VV2
!
!
!
end module NVE_m

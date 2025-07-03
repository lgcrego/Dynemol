module Berendsen_m

    use constants_m
    use syst                  ! using all syst
    use MD_read_m             , only: MM , atom , molecule, species
    use Berendsen_barostat    , only: Virial , barostat
    use VV_Parent             , only: VV

    public :: Berendsen 

    private

    type, extends(VV) :: Berendsen
    contains
        procedure :: VV1
        procedure :: VV2
    end type  Berendsen

    interface Berendsen
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
type(Berendsen) :: me 

!local variable ...

! select atomic or molecular kinetic energy to calculate the temperature ...
me % thermostat_type = NINT( float(maxval(molecule%N_of_atoms)) / float(MM%N_of_molecules) )

end function constructor
!
!
!========================
subroutine VV1( me , dt )
!========================
implicit none
class(Berendsen) , intent(inout) :: me
real*8           , intent(in)    :: dt

! local variables ...
real*8  :: acceleration(3)
real*8  :: dt_HALF , dt2_HALF
integer :: i

dt_HALF  = dt / two
dt2_HALF = dt_HALF * dt

! VV1 ... 
do i = 1 , MM % N_of_atoms
    if( atom(i) % flex ) then
        acceleration = atom(i)%ftotal / (atom(i)%mass * inv_kilo_mol)
        atom(i) % xyz = atom(i) % xyz + ( atom(i) % vel*dt + dt2_HALF*acceleration ) * mts_2_Angs
        atom(i) % vel = atom(i) % vel + dt_HALF*acceleration   
    end if
end do

end subroutine VV1
!
!
!=========================
 subroutine VV2( me , dt )
!=========================
implicit none
class(Berendsen) , intent(inout) :: me
real*8           , intent(in)    :: dt

! local variables ...
integer :: i , j 
real*8  :: E_kinetic , temperature , lambda , dt_half 
real*8  :: total_Momentum(3) , V_CM(3) , V_atomic(3) , accelerate(3)

if( using_barostat ) CALL Virial( me % thermostat_type )

E_kinetic= kinetic_erg( me % thermostat_type ) 
!######################################################
!               Berendsen Thermostat 

if( E_kinetic == 0.d0 ) then
    temperature = D_ONE
else
    ! instantaneous temperature : E_kin/(3/2*NkB) ... 
    select case (me % thermostat_type)

        case (0:1) ! <== molecular ...
        temperature = E_kinetic * iboltz / real( count(molecule%flex) )

        case (2:)  ! <== atomic ...
        temperature = E_kinetic * iboltz / real( count(atom%flex) )

    end select
endif

! Berendsen Thermostat ; turned off for talt == infty ...
If( talt == infty ) then
    lambda = D_one
else
    lambda = ( dt / (talt*pico_2_sec) ) * ( bath_T / temperature - D_ONE )
    lambda = SQRT(D_ONE + lambda)
end If

!######################################################

dt_half = dt / two

! VV2 ... 
E_kinetic = D_zero
select case (me % thermostat_type)

    case (0:1) ! <== molecular ...
          do i = 1 , MM % N_of_molecules
              total_Momentum = D_zero
              do j = molecule(i)%span % inicio , molecule(i)%span % fim
                  if ( atom(j) % flex ) then
                      accelerate = atom(j)% ftotal / (atom(j)% mass * inv_kilo_mol)
                      atom(j) % vel = atom(j) % vel + accelerate*dt_half
                      atom(j) % vel = atom(j) % vel * lambda

                      total_Momentum = total_Momentum + (atom(j)%mass * atom(j)%vel) * inv_kilo_mol

                  end if
              end do
              V_CM = total_Momentum / molecule(i) % mass
              E_kinetic = E_kinetic + molecule(i) % mass *  sum( V_CM * V_CM ) 
          end do

    case (2:) ! <== atomic ...
          V_atomic = D_zero
          do i = 1 , MM % N_of_atoms
              if( atom(i) % flex ) then
                  accelerate = atom(i)% ftotal / (atom(i)% mass * inv_kilo_mol)
                  atom(i) % vel = atom(i) % vel + accelerate*dt_half 
                  atom(i) % vel = atom(i) % vel * lambda

                  V_atomic = atom(i) % vel 
                  E_kinetic = E_kinetic + atom(i)%mass * sum( V_atomic * V_atomic ) * inv_kilo_mol
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
me % Kinetic = D_zero
do i = 1 , MM % N_of_atoms
    me % Kinetic = me % Kinetic + ( atom(i) % mass ) * sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half
end do
me % Kinetic = me % Kinetic * micro / MM % N_of_Molecules

me % density = get_density()

! calculation of pressure ...
if(using_barostat) CALL barostat( me % pressure , dt )

end subroutine VV2
!
!
!=======================================================
 function kinetic_erg( thermostat_type ) result(kinetic)
!=======================================================
implicit none
integer , intent(in)  :: thermostat_type

! local variables ...
real*8  :: total_Momentum(3) , V_CM(3) , V_atomic(3) , kinetic
integer :: i , j 

kinetic = D_zero

select case ( thermostat_type )

       case(0:1)  ! <== molecular ...
                  do i = 1 , MM % N_of_molecules
                      total_Momentum = D_zero
                      do j = molecule(i)%span % inicio , molecule(i)%span % fim
                          if( atom(j) % flex ) then
                              total_Momentum = total_Momentum + atom(j)%mass * atom(j)%vel
                          end if
                      end do   
          
                      V_CM = total_Momentum / molecule(i) % mass
          
                      ! 2*kinetic energy (molecular) of the system ...
                      kinetic = kinetic + molecule(i) % mass * sum( V_CM * V_CM )
                  end do
                  kinetic = kinetic * inv_kilo_mol**2
          
       case (2:)  ! <== atomic ...
                  do i = 1 , MM % N_of_atoms 
                      If( atom(i) % flex ) then
                          V_atomic = atom(i) % vel 
                          ! 2*kinetic energy (atomic) of the system ...
                          kinetic = kinetic + atom(i) % mass * sum( V_atomic * V_atomic ) * inv_kilo_mol
                      end if
                  end do

end select

end function kinetic_erg
!
!
!
!====================================
 function get_density result(density)
!====================================
implicit none

! local variables ...
real*8 :: volume, massa, density

massa   = sum( species(:) % N_of_molecules * species(:) % mass )
volume  = product( MM % box(:) * Angs_2_mts )
density = (massa / volume) * milli

end function get_density
!
!
!
end module Berendsen_m

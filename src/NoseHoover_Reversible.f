module NH_Reversible_m

    use constants_m
    use syst         ! using all syst
    use MD_read_m    , only: MM , atom , molecule , species
    use VV_Parent    , only: VV
    use Berendsen_m  , only: Barostat

    public :: NH_Reversible

    private

    type, extends(VV) :: NH_Reversible
    contains
        procedure :: VV1
        procedure :: VV2
    end type NH_Reversible

    interface NH_Reversible
        module procedure  constructor
    end interface

    ! module variables ...
	 real*8 :: Csi = D_zero
         real*8 :: Q

contains
!
!
!
!===================================
 function constructor() result( me )
!===================================
implicit none
type(NH_Reversible) :: me

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
class(NH_Reversible) , intent(inout) :: me
real*8               , intent(in) :: dt

! local variables ...
real*8  :: ai(3) , tmp(3) , V_CM(3)
real*8  :: massa , dt_HALF , dt2_HALF , sigma
integer :: i , j , j1 , j2 , nresid

dt_HALF  = dt / two
dt2_HALF = dt_HALF * dt

! calculation of the kinetic energy and thermostat-related things ...
me % kinetic = D_zero
select case (me % thermostat_type)

    case (0:1) ! <== molecular ...
      sigma = bath_T*boltz*THREE*real( count(molecule%flex) )*HALF
      ! molecular kinetic energy at time t ...
      do i = 1 , MM % N_of_molecules
        tmp = D_zero
        nresid = molecule(i) % nr
        j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
        j2 = sum(molecule(1:nresid) % N_of_atoms)
        do j = j1 , j2
            if ( atom(j) % flex ) then
              massa = mol / atom(j) % mass
              tmp = tmp + atom(j) % vel / massa
            end if
        end do
        V_CM    = tmp / molecule(i) % mass
        me % kinetic = me % kinetic + molecule(i) % mass *  sum( V_CM * V_CM ) * half
      end do

    case (2:) ! <== atomic ...
      sigma = bath_T*boltz*THREE*real(count(atom%flex))*HALF 
      ! atomic kinetic energy at time t ...
      do i = 1 , MM % N_of_atoms 
        me % kinetic = me%kinetic + imol*atom(i)%mass*sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half
      end do

end select

Q   =  D_one/(THREE*real(count(atom%flex))*boltz*bath_T*(talt*pico_2_sec)**2)

Csi =  Csi + Q*(me%kinetic - sigma)*dt

! VV1 ...
do i = 1 , MM % N_of_atoms
    if( atom(i) % flex ) then
        massa = mol / atom(i) % mass
        ai = atom(i) % ftotal * massa
        atom(i) % vel = (atom(i) % vel)*exp(-Csi*dt*HALF) + dt_HALF*ai
        atom(i) % xyz = atom(i) % xyz + (atom(i) % vel*dt)* mts_2_Angs

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
class(NH_Reversible) , intent(inout) :: me
real*8               , intent(in)    :: dt

! local variables ...
real*8  :: tmp(3) , V_CM(3) , V_atomic(3)
real*8  :: massa , factor , sumtemp , dt_half , volume , sigma
integer :: i , j , j1 , j2 , nresid

dt_half = dt / two
sumtemp = D_zero

! calculation of the kinetic energy and thermostat-related things ...
me % kinetic = D_zero
select case (me % thermostat_type)

    case (0:1) ! <== molecular ...
      sigma = bath_T*boltz*THREE*real( count(molecule%flex) )*HALF
      ! molecular kinetic energy at time t ...
      do i = 1 , MM % N_of_molecules
        tmp = D_zero
        nresid = molecule(i) % nr
        j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
        j2 = sum(molecule(1:nresid) % N_of_atoms)
        do j = j1 , j2
            if ( atom(j) % flex ) then
              massa = mol / atom(j) % mass
              tmp = tmp + atom(j) % vel / massa
            end if
        end do
        V_CM    = tmp / molecule(i) % mass
        me % kinetic = me % kinetic + molecule(i) % mass *  sum( V_CM * V_CM ) * half
      end do

    case (2:) ! <== atomic ...
      sigma = bath_T*boltz*THREE*real(count(atom%flex))*HALF 
      ! atomic kinetic energy at time t ...
      do i = 1 , MM % N_of_atoms 
        me % kinetic = me%kinetic + imol*atom(i)%mass*sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half
      end do

end select

Csi = Csi + Q*(me%kinetic - sigma)*dt

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
                    atom(j) % vel = (atom(j) % vel +  dt_half * atom(j) % ftotal * massa)
                    atom(j) % vel = (atom(j) % vel)*exp(-Csi*dt*HALF)
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
                atom(i) % vel = (atom(i) % vel +  dt_half * atom(i) % ftotal * massa)
                atom(i) % vel = (atom(i) % vel)*exp(-Csi*dt*HALF)
                    
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
me % kinetic = D_zero
do i = 1 , MM % N_of_atoms
    me % kinetic = me % kinetic + atom(i) % mass * sum( atom(i) % vel(:) * atom(i) % vel(:) ) * half
end do
me % kinetic = me % kinetic * micro / MM % N_of_Molecules

! calculation of pressure and density ...
massa = sum( species(:) % N_of_molecules * species(:) % mass )
volume = product( MM % box(:) * Angs_2_mts )
me % density = (massa / volume) * milli

end subroutine VV2
!
!
!
end module NH_Reversible_m

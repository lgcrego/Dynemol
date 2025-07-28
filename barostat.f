module Berendsen_Barostat

    use constants_m
    use syst                  ! using all syst
    use parameters_m          , only: PBC 
    use MD_read_m             , only: MM , atom , molecule
    use F_inter_m             , only: virial_tensor

    public :: Ek_tensor , barostat 

    private

    ! module parameters ...
    real*8 , parameter :: inv_kilo_mol = imol

    ! module variables ...
    real*8 , save :: kinetic_tensor(3,3)
 
    contains

!=======================================
 subroutine Ek_Tensor( thermostat_type )
!=======================================
implicit none
integer , intent(in) :: thermostat_type

! local variables ...
real*8  :: total_Momentum(3) , V_CM(3) , V_atomic(3) , massa 
integer :: i , j , k , l

kinetic_tensor = D_zero

select case ( thermostat_type )
    case(0:1)  ! <== molecular ...         
        do i = 1 , MM % N_of_molecules
            total_Momentum = D_zero
            do j = molecule(i)%span % inicio , molecule(i)%span % fim
                if( atom(j) % flex ) then
                    total_Momentum = total_Momentum + atom(j)%mass * atom(j)%vel
                end if
            end do   
            V_CM = total_Momentum/molecule(i)%mass * inv_kilo_mol
            do concurrent (k = 1:3, l = 1:3, l>=k)
                 kinetic_tensor(k,l) = kinetic_tensor(k,l) + molecule(i)%mass*V_CM(k)*V_CM(l)
            end do     
        end do

    case (2:)  ! <== atomic ...
        do i = 1 , MM % N_of_atoms 
            If( atom(i) % flex ) then
                V_atomic = atom(i) % vel 
                massa = atom(i)%mass * inv_kilo_mol
                do concurrent (k = 1:3, l = 1:3, l>=k)
                    kinetic_tensor(k,l) = kinetic_tensor(k,l) + massa*V_atomic(k)*V_atomic(l)
                end do     
            end if
        end do
end select

do concurrent (l = 1:2, k = 1:3, k>l)
  kinetic_tensor(k, l) = kinetic_tensor(l, k)
end do

end subroutine Ek_Tensor
!
!
!
!====================================
 subroutine barostat( pressure , dt )
!====================================
implicit none
real*8 , intent(out) :: pressure
real*8 , intent(in)  :: dt

! local variables ...
integer :: i, j
real*8  :: mip , volume
real*8  :: stress_tensor(3,3) 

volume  = product( MM % box(:) * Angs_2_mts )

kinetic_tensor = kinetic_tensor / (cm_2_Angs * volume)
virial_tensor  = virial_tensor  / (cm_2_Angs * volume)

stress_tensor = kinetic_tensor + virial_tensor 

! scalar pressure, for isotropic systems
pressure = ( stress_tensor(1,1) + stress_tensor(2,2) + stress_tensor(3,3) ) * third

mip = dt * ( 107.0d-6 / (talp * pico_2_sec) ) * ( pressure - press )
mip = (D_one + mip)**third

MM % box = MM % box * mip
 
do i = 1 , MM % N_of_molecules
    do j = molecule(i)%span % inicio , molecule(i)%span % fim
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
end module Berendsen_Barostat

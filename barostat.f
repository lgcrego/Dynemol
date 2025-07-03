module Berendsen_Barostat

    use constants_m
    use syst                  ! using all syst
    use parameters_m          , only: PBC 
    use MD_read_m             , only: MM , atom , molecule

    public :: Virial , barostat , InitializeStressMatrix , BuildStressMatrix , ConcludeStressMatrix

    private

    ! module parameters ...
    real*8 , parameter :: kilo_mol = mol
    real*8 , parameter :: inv_kilo_mol = imol

    ! module variables ...
    real*8 , save :: stresvv(3,3)
    real*8 , save :: stresSR(3,3), stresRE(3,3)

    real*8 , save :: stresSR11 , stresSR22 , stresSR33 , stresSR12 , stresSR13 , stresSR23
    real*8 , save :: stresRE11 , stresRE22 , stresRE33 , stresRE12 , stresRE13 , stresRE23 
 
    contains

!=================================
 subroutine InitializeStressMatrix
!=================================
  implicit none
  !--------------------------------------------------------------
  ! initializing variables for this integration step ...
  !--------------------------------------------------------------
  stresSR(:,:) = D_zero
  stresRE(:,:) = D_zero
  !upper triangle
  stresSR11 = D_zero; stresSR12 = D_zero; stresSR13 = D_zero
  stresSR22 = D_zero; stresSR23 = D_zero; stresSR33 = D_zero
  !upper triangle
  stresRE11 = D_zero; stresRE12 = D_zero; stresRE13 = D_zero
  stresRE22 = D_zero; stresRE23 = D_zero; stresRE33 = D_zero
  !--------------------------------------------------------------
end subroutine InitializeStressMatrix
!
!
!========================================================
 subroutine BuildStressMatrix( fs , Fcoul , cm_kl , rkl )
!========================================================
  implicit none
  real*8 , intent(in) :: fs 
  real*8 , intent(in) :: Fcoul
  real*8 , intent(in) :: cm_kl(:)
  real*8 , intent(in) :: rkl(:) 
  
  stresSR11 = stresSR11 + cm_kl(1) * fs * rkl(1)
  stresSR22 = stresSR22 + cm_kl(2) * fs * rkl(2)
  stresSR33 = stresSR33 + cm_kl(3) * fs * rkl(3)
  stresSR12 = stresSR12 + cm_kl(1) * fs * rkl(2)
  stresSR13 = stresSR13 + cm_kl(1) * fs * rkl(3)
  stresSR23 = stresSR23 + cm_kl(2) * fs * rkl(3)
  
  stresRE11 = stresRE11 + cm_kl(1) * Fcoul * rkl(1)
  stresRE22 = stresRE22 + cm_kl(2) * Fcoul * rkl(2)
  stresRE33 = stresRE33 + cm_kl(3) * Fcoul * rkl(3)
  stresRE12 = stresRE12 + cm_kl(1) * Fcoul * rkl(2)
  stresRE13 = stresRE13 + cm_kl(1) * Fcoul * rkl(3)
  stresRE23 = stresRE23 + cm_kl(2) * Fcoul * rkl(3)
end subroutine BuildStressMatrix
!
!
!===============================
 subroutine ConcludeStressMatrix
!===============================
  implicit none
  !--------------------------------------------------------
  stresSR11 = stresSR11 * factor3
  stresSR22 = stresSR22 * factor3
  stresSR33 = stresSR33 * factor3
  stresSR12 = stresSR12 * factor3
  stresSR13 = stresSR13 * factor3
  stresSR23 = stresSR23 * factor3
  
  stresSR(1,1) = stresSR11; stresSR(2,2) = stresSR22
  stresSR(3,3) = stresSR33; stresSR(1,2) = stresSR12
  stresSR(1,3) = stresSR13; stresSR(2,3) = stresSR23
  
  stresSR(2,1) = stresSR(1,2); stresSR(3,1) = stresSR(1,3)
  stresSR(3,2) = stresSR(2,3) 
  !--------------------------------------------------------
  
  !--------------------------------------------------------
  stresRE11 = stresRE11 * factor3
  stresRE22 = stresRE22 * factor3
  stresRE33 = stresRE33 * factor3
  stresRE12 = stresRE12 * factor3
  stresRE13 = stresRE13 * factor3
  stresRE23 = stresRE23 * factor3
  
  stresRE(1,1) = stresRE11; stresRE(2,2) = stresRE22
  stresRE(3,3) = stresRE33; stresRE(1,2) = stresRE12
  stresRE(1,3) = stresRE13; stresRE(2,3) = stresRE23
  
  stresRE(2,1) = stresRE(1,2); stresRE(3,1) = stresRE(1,3)
  stresRE(3,2) = stresRE(2,3)
  !--------------------------------------------------------
end subroutine ConcludeStressMatrix
!
!
!====================================
 subroutine Virial( thermostat_type )
!====================================
implicit none
integer , intent(in) :: thermostat_type

! local variables ...
real*8  :: total_Momentum(3) , V_CM(3) , V_atomic(3) , massa 
integer :: i , j , k

stresvv = D_zero

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
            forall(k=1:3) stresvv(k,k) = stresvv(k,k) +  molecule(i)%mass * V_CM(k)*V_CM(k)
            stresvv(1,2) = stresvv(1,2) + molecule(i) % mass * V_CM(1) * V_CM(2)
            stresvv(1,3) = stresvv(1,3) + molecule(i) % mass * V_CM(1) * V_CM(3)
            stresvv(2,3) = stresvv(2,3) + molecule(i) % mass * V_CM(2) * V_CM(3)
        end do

    case (2:)  ! <== atomic ...
        do i = 1 , MM % N_of_atoms 
            If( atom(i) % flex ) then
                V_atomic = atom(i) % vel 
                massa = atom(i)%mass * inv_kilo_mol
                forall( k=1:3) stresvv(k,k) = stresvv(k,k) + massa * V_atomic(k) * V_atomic(k)
                stresvv(1,2) = stresvv(1,2) + massa * V_atomic(1) * V_atomic(2)
                stresvv(1,3) = stresvv(1,3) + massa * V_atomic(1) * V_atomic(3)
                stresvv(2,3) = stresvv(2,3) + massa * V_atomic(2) * V_atomic(3)
            end if
        end do

end select

stresvv(2,1) = stresvv(1,2)
stresvv(3,1) = stresvv(1,3)
stresvv(3,2) = stresvv(2,3)

end subroutine Virial
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
real*8  :: Astres(3,3) 

stresvv = stresvv / (cm_2_Angs * volume)
stressr = stressr / (cm_2_Angs * volume)
stresre = stresre / (cm_2_Angs * volume)
Astres  = stresvv + stressr + stresre

pressure = ( Astres(1,1) + Astres(2,2) + Astres(3,3) ) * third

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

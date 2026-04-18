module Build_DWFF

use iso_fortran_env
use constants_m
use MM_types     , only: Dissociative, MM_molecular

    public :: include_DWFF_parameters , HOH_diss_parms

    private

    ! module variables ...
    type(dissociative) :: HOH_diss_parms

    real*8 :: PointCharge_O, PointCharge_H, DiffCharge_O, DiffCharge_H  
    real*8 :: OH_A_repulsive, OH_short_range, OH_q_dsprsion, OH_C6_dsprsion
    real*8 :: OO_A_repulsive, OO_short_range, OO_q_dsprsion, OO_C6_dsprsion
    real*8 :: HH_A_repulsive, HH_short_range, HH_q_dsprsion, HH_C6_dsprsion
    real*8 :: lambda, r0, gama, theta_0

    ! module parameters ...
    real*8 :: temperature = 300.0d0
    real*8 :: pressure    = 1.0d0

    ! module variables ...
    real*8 :: mtxA(0:3,0:5) 

contains
!
!
!===========================================
 subroutine include_DWFF_parameters(species)
!===========================================
implicit none
type(MM_molecular), intent(in):: species(:)

call DWFF_Parameters

call get_DWFF_Parameters

call set_H_pointers(species)

end subroutine include_DWFF_parameters
!
!
!
!
!===========================
 subroutine DWFF_Parameters 
!===========================
implicit none
!-------------------------------------------------------------
!                PARAMETERS FROM THE ARTICLE
!
! Dissociative Water Potential for Molecular Dynamics Simulations
!           T. S. Mahadevan and S. H. Garofalini
!           J. Phys. Chem. B 2007, 111, 8919-8927
!-------------------------------------------------------------

call define_mtxA 

!  Charges on Species
PointCharge_O = -0.904d0
PointCharge_H = +0.452d0

DiffCharge_O  = +0.226d0
DiffCharge_H  = -0.113d0

!  Parameters of the Two-Body Potential
OH_A_repulsive = 2.283d-16  ! <== Joule
OH_short_range = PolynomialFunction( )       ! <== Angs
OH_q_dsprsion  = 24.d0      ! <== Angs
OH_C6_dsprsion = D_zero     ! <== J*Angs^6

OO_A_repulsive = 4.250d-17  ! <== Joule
OO_short_range = 0.610d0    ! <== Angs
OO_q_dsprsion  = 24.d0      ! <== Angs
OO_C6_dsprsion = 4.226d-18  ! <== J*Angs^6

HH_A_repulsive = D_zero     ! <== Joule
HH_short_range = D_zero     ! <== Angs
HH_q_dsprsion  = 24.d0      ! <== Angs
HH_C6_dsprsion = D_zero     ! <== J*Angs^6

! Three-Body Parameters
lambda   = 3.0d-17          ! <== Joule
theta_0  = 100.d0           ! <== degrees
!------------------------------------------
! original parameters
r0       = 1.6d0            ! <== Angs
gama     = 1.3d0            ! <== Angs
!------------------------------------------
! tuned parameters: 
! Phys. Chem. Chem. Phys.,2015, 17, 10934
!------------------------------------------
!r0       = 1.5d0            ! <== Angs
!gama     = 1.089d0          ! <== Angs
!------------------------------------------

!-------------------------------------------------------------
!          CONVERT PARAMETERS TO DYNEMOL UNITS     
! the following values are being converted to Dynemol internal units ...                                                                       
! Dynemol uses Joule and Newton units to evaluate the eqs. of motion ...
!-------------------------------------------------------------

!  Parameters of the Two-Body Potential
OH_A_repulsive = OH_A_repulsive / factor3
OH_C6_dsprsion = OH_C6_dsprsion / factor3

OO_A_repulsive = OO_A_repulsive / factor3
OO_C6_dsprsion = OO_C6_dsprsion / factor3

HH_A_repulsive = HH_A_repulsive / factor3
HH_C6_dsprsion = HH_C6_dsprsion / factor3

! Three-Body Parameters
lambda  = lambda / factor3
theta_0 = theta_0 * deg_2_rad   

end subroutine DWFF_Parameters 
!
!
!
!========================================
function PolynomialFunction ( ) result(f)
!========================================
implicit none

! local variables ...
integer :: m , n
real*8  :: f

f = 0.d0
do m = 0 , 3
   do n = 0 , 5
      f = f + mtxA(m,n) * pressure**m * temperature**n
   end do
end do

end function PolynomialFunction 
!
!
!
!=====================
subroutine define_mtxA
!=====================
implicit none

mtxA(0,0) = 0.655726502d0
mtxA(1,0) = 3.403472d-4
mtxA(2,0) = -4.057853d-8
mtxA(3,0) = 1.657262d-12

mtxA(0,1) = -1.04442689d-2
mtxA(1,1) = -3.986929d-6
mtxA(2,1) = 4.677537d-10
mtxA(3,1) = -1.838785d-14

mtxA(0,2) = 8.31892416d-5
mtxA(1,2) = 1.742261d-8
mtxA(2,2) = -2.007873d-12
mtxA(3,2) = 7.549619d-17

mtxA(0,3) = -3.07929142d-7
mtxA(1,3) = -3.364186d-11
mtxA(2,3) = 3.800411d-15
mtxA(3,3) = -1.355453d-19

mtxA(0,4) = 5.44770929d-10
mtxA(1,4) = 2.419996d-14
mtxA(2,4) = -2.672717d-18
mtxA(3,4) = 8.939302d-23

mtxA(0,5) = -3.73609493d-13
mtxA(1,5) = D_zero
mtxA(2,5) = D_zero
mtxA(3,5) = D_zero

end subroutine define_mtxA 
!
!
!
!==================================
 subroutine set_H_pointers(species)
!==================================
implicit none
type(MM_molecular), intent(in):: species(:)

! local variables
integer :: HOH

! locate HOH species in array species
do HOH = 1 , size(species)
   if ( species(HOH)%DWFF ) exit
end do

! for HOH species, do the following
associate( MMSymbol => species(HOH)% atom(:)% MMSymbol )
    if( MMSymbol(1) == "OX" ) then
        HOH_diss_parms% H_ptr(:) = [2,3]
    elseif ( MMSymbol(2) == "OX" ) then
        HOH_diss_parms% H_ptr(:) = [1,3]
    elseif ( MMSymbol(3) == "OX" ) then
        HOH_diss_parms% H_ptr(:) = [1,2]
    endif
end associate

end  subroutine set_H_pointers
!
!
!
!===============================
 subroutine  get_DWFF_Parameters
!===============================
implicit none

! atomic pairs
allocate( HOH_diss_parms% PairSymbols(3,2) )
HOH_diss_parms% PairSymbols(1,1:2) = ["OX","HX"]
HOH_diss_parms% PairSymbols(2,1:2) = ["OX","OX"]
HOH_diss_parms% PairSymbols(3,1:2) = ["HX","HX"]

! short range 2-body paramters
allocate( HOH_diss_parms% SR(2,3) )
HOH_diss_parms% SR(1,1) = OH_A_repulsive
HOH_diss_parms% SR(2,1) = OO_A_repulsive

HOH_diss_parms% SR(1,2) = HALF / OH_short_range
HOH_diss_parms% SR(2,2) = HALF / OO_short_range

HOH_diss_parms% SR(1,3) = OH_C6_dsprsion  ! 0.d0
HOH_diss_parms% SR(2,3) = OO_C6_dsprsion

! Coulomb 2-body paramters
allocate( HOH_diss_parms% Coul(3,4) )
HOH_diss_parms% Coul(1,1) = PointCharge_O * PointCharge_H    !........................ qiqj(O-H)
HOH_diss_parms% Coul(2,1) = PointCharge_O**2                 !........................ qiqj(O-O)
HOH_diss_parms% Coul(3,1) = PointCharge_H**2                 !........................ qiqj(H-H)

HOH_diss_parms% Coul(1,2) = DiffCharge_O * DiffCharge_H      !........................ qdiqdj(O-H)
HOH_diss_parms% Coul(2,2) = DiffCharge_O**2                  !........................ qdiqdj(O-O)
HOH_diss_parms% Coul(3,2) = DiffCharge_H**2                  !........................ qdiqdj(H-H)

HOH_diss_parms% Coul(1,3) = PointCharge_O*DiffCharge_H + DiffCharge_O*PointCharge_H  ! qdij_mix(O-H)
HOH_diss_parms% Coul(2,3) = two*PointCharge_O*DiffCharge_O                           ! qdij_mix(O-O)
HOH_diss_parms% Coul(3,3) = two*PointCharge_H*DiffCharge_H                           ! qdij_mix(H-H)

HOH_diss_parms% Coul(1,4) = half / OH_q_dsprsion             !........................ chrg_decay(O-H)
HOH_diss_parms% Coul(2,4) = half / OO_q_dsprsion             !........................ chrg_decay(O-O)
HOH_diss_parms% Coul(3,4) = half / HH_q_dsprsion             !........................ chrg_decay(H-H)

! Coulomb 3-body paramters
allocate( HOH_diss_parms% Angle(1,4) )
HOH_diss_parms% Angle(1,1) = lambda
HOH_diss_parms% Angle(1,2) = r0
HOH_diss_parms% Angle(1,3) = gama
HOH_diss_parms% Angle(1,4) = cos(theta_0)

end subroutine get_DWFF_Parameters
!
!
end module Build_DWFF

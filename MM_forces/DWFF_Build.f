module Build_DWFF

use iso_fortran_env
use constants_m
use MM_parms_module , only: DWFF_type
use type_m          , only: dynemoldir, warning
use util_m          , only: TO_UPPER_CASE, split_line
use MM_types        , only: Dissociative, MM_molecular, MM_atomic

    public :: include_DWFF_parameters, dump_DWFF_parameters, set_H_pointers
    public :: HOH_diss_parms

    private

    ! module variables ...
    type(dissociative) :: HOH_diss_parms

    real*8  :: PointCharge_O, PointCharge_H, DiffCharge_O, DiffCharge_H  
    real*8  :: OH_A_repulsive, OH_short_range, OH_q_dsprsion, OH_C6_dsprsion
    real*8  :: OO_A_repulsive, OO_short_range, OO_q_dsprsion, OO_C6_dsprsion
    real*8  :: HH_A_repulsive, HH_short_range, HH_q_dsprsion, HH_C6_dsprsion
    real*8  :: lambda, r0, gama, theta_0
    logical :: from_outside(8)

    ! module parameters ...
    real*8 :: temperature = 300.0d0
    real*8 :: pressure    = 1.0d0

    ! module variables ...
    real*8 :: mtxA(0:3,0:5) 

contains
!
!
!==================================
 subroutine include_DWFF_parameters
!==================================
implicit none

!local variables ...

select case (DWFF_type) 
    case("DIFFUSE")
        call Diffuse_DWFF_Parameters

    case("SPC_LIKE")
        call SPC_like_DWFF_Parameters

    case("QMMM")
        call QMMM_DWFF_Parameters
end select

call export_DWFF_Parameters

end subroutine include_DWFF_parameters
!
!
!
!
!==================================
 subroutine Diffuse_DWFF_Parameters 
!==================================
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
lambda   = 3.0d-17 ! <== Joule
theta_0  = 100.d0  ! <== degrees
r0       = 1.6d0   ! <== Angs
gama     = 1.3d0   ! <== Angs

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

end subroutine Diffuse_DWFF_Parameters 
!
!
!
!
!===================================
 subroutine SPC_like_DWFF_Parameters 
!===================================
implicit none
!-------------------------------------------------------------
!                PARAMETERS FROM THE ARTICLE
!
!       Towards a dissociative SPC-like water model: 
! probing the impact of intramolecular Coulombic contributions
!           Martin J. Wiedemair and Thomas S. Hofer
!           Phys. Chem. Chem. Phys., 2017, 19, 31910
!-------------------------------------------------------------

call define_mtxA 

!  Charges on Species
PointCharge_O = -0.898d0
PointCharge_H = +0.449d0

DiffCharge_O  = D_zero
DiffCharge_H  = D_zero

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
lambda   = 2.0d-17 ! <== Joule
theta_0  = 98.09d0 ! <== degrees
r0       = 1.55d0  ! <== Angs
gama     = 1.19d0  ! <== Angs

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

end subroutine SPC_like_DWFF_Parameters 
!
!
!
!===============================
 subroutine QMMM_DWFF_Parameters 
!===============================
implicit none
!-------------------------------------------------------------
!                PARAMETERS FROM THE ARTICLE
!
!       BASED ON THE PAPER: 
!       How to build a better pair potential for water
!            Bertrand Guillot and Yves Guissani
!           J. Chem. Phys. 2001, 114, 6720-6733
!-------------------------------------------------------------

call define_mtxA 

!  Charges on Species
PointCharge_O = -0.898d0
PointCharge_H = +0.449d0

DiffCharge_O  = D_zero
DiffCharge_H  = D_zero

!  Parameters of the Two-Body Potential
OH_A_repulsive = 2.283d-16  ! <== Joule
OH_short_range = PolynomialFunction( )       ! <== Angs
OH_q_dsprsion  = 1.5d0      ! <== Angs, from Guillot
OH_C6_dsprsion = D_zero     ! <== J*Angs^6

OO_A_repulsive = 4.250d-17  ! <== Joule
OO_short_range = 0.610d0    ! <== Angs
OO_q_dsprsion  = 1.5d0      ! <== Angs, from Guillot
OO_C6_dsprsion = 4.226d-18  ! <== J*Angs^6

HH_A_repulsive = D_zero     ! <== Joule
HH_short_range = D_zero     ! <== Angs
HH_q_dsprsion  = 1.5d0      ! <== Angs, from Guillot
HH_C6_dsprsion = D_zero     ! <== J*Angs^6

! Three-Body Parameters
lambda   = 2.0d-17 ! <== Joule
theta_0  = 98.09d0 ! <== degrees
r0       = 1.55d0  ! <== Angs
gama     = 1.19d0  ! <== Angs

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

end subroutine QMMM_DWFF_Parameters 
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
        HOH_diss_parms% O_ptr    = 1
        HOH_diss_parms% H_ptr(:) = [2,3]
    elseif ( MMSymbol(2) == "OX" ) then
        HOH_diss_parms% O_ptr    = 2
        HOH_diss_parms% H_ptr(:) = [1,3]
    elseif ( MMSymbol(3) == "OX" ) then
        HOH_diss_parms% O_ptr    = 3
        HOH_diss_parms% H_ptr(:) = [1,2]
    endif
end associate

end  subroutine set_H_pointers
!
!
!
!==================================
 subroutine  export_DWFF_Parameters
!==================================
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

HOH_diss_parms%contain_diffuse = .not. all(abs(HOH_diss_parms%Coul(:,2:3)) < mid_prec)
print*, "diffuse charge = ", HOH_diss_parms%contain_diffuse
print*,""

! Coulomb 3-body paramters
allocate( HOH_diss_parms% Angle(1,4) )
HOH_diss_parms% Angle(1,1) = lambda
HOH_diss_parms% Angle(1,2) = r0
HOH_diss_parms% Angle(1,3) = gama
HOH_diss_parms% Angle(1,4) = cos(theta_0)

! export charges
HOH_diss_parms% PointCharge_O = PointCharge_O
HOH_diss_parms% PointCharge_H = PointCharge_H

end subroutine export_DWFF_Parameters
!
!
!
!=======================================
 subroutine dump_DWFF_parameters(f_unit)
!=======================================
implicit none
integer, intent(in) :: f_unit

write(f_unit,"(A18,F8.4)")   "PointCharge_O  = ", PointCharge_O                                                                                        
write(f_unit,"(A18,F8.4)")   "PointCharge_H  = ", PointCharge_H
write(f_unit,"(A18,F8.4)")   "DiffCharge_O   = ", DiffCharge_O
write(f_unit,"(A18,F8.4)")   "DiffCharge_H   = ", DiffCharge_H
write(f_unit,"(A20,ES10.3)") "lambda  (Joule)  = ", lambda * factor3
write(f_unit,"(A20,F8.4)")   "theta_0 (degree) = ", theta_0 / deg_2_rad
write(f_unit,"(A20,F8.4)")   "r0 (Angstron)    = ", r0
write(f_unit,"(A20,F8.4)")   "gamma (Angstron) = ", gama

end subroutine dump_DWFF_parameters
!
!
!
end module Build_DWFF

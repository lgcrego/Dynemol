module Build_DWFF

use iso_fortran_env
use constants_m
use type_m      , only: dynemoldir, warning
use util_m      , only: TO_UPPER_CASE, split_line
use MM_types    , only: Dissociative, MM_molecular, MM_atomic

    public :: include_DWFF_parameters , HOH_diss_parms, read_DWFF_parameters, dump_DWFF_parameters, check_DWFF_parameters

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
if(.not. from_outside(1)) PointCharge_O = -0.904d0
if(.not. from_outside(2)) PointCharge_H = +0.452d0

if(.not. from_outside(3)) DiffCharge_O  = +0.226d0
if(.not. from_outside(4)) DiffCharge_H  = -0.113d0

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
if(.not. from_outside(5)) lambda   = 3.0d-17 ! <== Joule
if(.not. from_outside(6)) theta_0  = 100.d0  ! <== degrees
if(.not. from_outside(7)) r0       = 1.6d0   ! <== Angs
if(.not. from_outside(8)) gama     = 1.3d0   ! <== Angs

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

end subroutine get_DWFF_Parameters
!
!
!
!===============================
 subroutine read_DWFF_parameters
!===============================
implicit none

!local variables
character(len=120) :: line
character(len=80)  :: string
integer :: pos, fim, ioerr

from_outside = .false.

do 
    read(33,'(A)',iostat=ioerr) line

    ! End-of-file
    if (ioerr < 0) exit

    ! Read error
    if (ioerr > 0) then
        print*, "Error reading DWFF parameter file"
        stop
    end if

    ! Skip blank lines
    if (len_trim(line) == 0) cycle

    ! Skip comment lines
    if (line(1:1) == "!") cycle

    pos = index(line, "=")
    if (pos == 0) then
        print*, "Error parsing DWFF parameter file"
        stop
    end if

    fim = index(line, "!")
    if (fim == 0) fim = len_trim(line) + 1
    
    string = TO_UPPER_CASE( trim(adjustL(line(:pos-1))) )

    ioerr = 0
    select case(string)

         case("POINTCHARGE_O")
               read(line(pos+1:fim-1),*,iostat=ioerr)  Pointcharge_O
               from_outside(1) = .true.

         case("POINTCHARGE_H")
               read(line(pos+1:fim-1),*,iostat=ioerr)  Pointcharge_H
               from_outside(2) = .true.

         case("DIFFCHARGE_O")
               read(line(pos+1:fim-1),*,iostat=ioerr)  Diffcharge_O
               from_outside(3) = .true.

         case("DIFFCHARGE_H")
               read(line(pos+1:fim-1),*,iostat=ioerr)  Diffcharge_H
               from_outside(4) = .true.

         case("LAMBDA")
               read(line(pos+1:fim-1),*,iostat=ioerr)  lambda
               from_outside(5) = .true.

         case("THETA_0","THETA-0")
               read(line(pos+1:fim-1),*,iostat=ioerr)  theta_0
               from_outside(6) = .true.

         case("R0")
               read(line(pos+1:fim-1),*,iostat=ioerr)  R0
               from_outside(7) = .true.

         case("GAMA","GAMMA")
               read(line(pos+1:fim-1),*,iostat=ioerr)  gama
               from_outside(8) = .true.

        end select 

        if (ioerr /= 0) then
            print*, "Error parsing line:"
            print*, trim(line)
            stop
        end if
end do

end subroutine read_DWFF_parameters
!
!
!
!===============================
 subroutine dump_DWFF_parameters
!===============================
implicit none

!local variables
character(len=:) , allocatable :: keyword_list(:), string(:)
character(120)                 :: line
integer                        :: f_unit, ioerr

allocate( character(len=14) :: keyword_list(12) )

keyword_list = [ "OH_A_repulsive", "OH_short_range", "OH_q_dsprsion", "OH_C6_dsprsion", &
                 "OO_A_repulsive", "OO_short_range", "OO_q_dsprsion", "OO_C6_dsprsion", &
                 "HH_A_repulsive", "HH_short_range", "HH_q_dsprsion", "HH_C6_dsprsion"  ]

write(51,"(A18,F8.4)") "PointCharge_O  = ", PointCharge_O                                                                                        
write(51,"(A18,F8.4)") "PointCharge_H  = ", PointCharge_H
write(51,"(A18,F8.4)") "DiffCharge_O   = ", DiffCharge_O
write(51,"(A18,F8.4)") "DiffCharge_H   = ", DiffCharge_H

OPEN(newunit=f_unit, file=dynemoldir//'MM_forces/DWFF_Build.f', status='old')
do
    read(f_unit,'(A)',iostat=ioerr) line

    ! end-of-file
    if (ioerr < 0) exit

    ! read error
    if (ioerr > 0) then
        print*, "Error reading file: DWFF_Build.f"
        close(f_unit)
    end if
    
    allocate( string , source=split_line( line , token_length=14 ) )

    if (.not. allocated(string)) cycle

    if (size(string) > 0) then
        if (any(keyword_list == string(1))) then
            write(51,*) trim(line)
        end if

        if (trim(string(1)) == "HH_C6_dsprsion") then
            deallocate(string)
            exit
        end if
    end if

    deallocate(string)

end do
close(f_unit)

write(51,"(A20,ES10.3)") "lambda  (Joule)  = ", lambda * factor3
write(51,"(A20,F8.4)") "theta_0 (degree) = ", theta_0 / deg_2_rad
write(51,"(A20,F8.4)") "r0 (Angstron)    = ", r0
write(51,"(A20,F8.4)") "gamma (Angstron) = ", gama

end subroutine dump_DWFF_parameters
!
!
!
!===================================
subroutine check_DWFF_parameters(a)
!===================================
implicit none
type(MM_atomic), intent(in) :: a(:)

! Local variables
logical :: flag(size(a))
real(8), parameter :: tol = 1.d-10

flag = .false.

!-------------------------------------------------
! Check MM charges for DWFF atom types
!-------------------------------------------------
where (a%MMSymbol == "OX")
    flag = abs(a%MM_charge - PointCharge_O) > tol
elsewhere (a%MMSymbol == "HX")
    flag = abs(a%MM_charge - PointCharge_H) > tol
end where

!-------------------------------------------------
! Abort if inconsistencies are detected
!-------------------------------------------------
if (any(flag)) then
    call warning("halting: HOH.psf file has conflicting MM charges")
    stop
end if

end subroutine check_DWFF_parameters
!
!
!
end module Build_DWFF

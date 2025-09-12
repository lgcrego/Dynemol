module Reax_F_intra

    use constants_m
    use omp_lib
    use parameters_m , only: PBC 
    use for_force    , only: DWFF_intra, rcutsq
    use MD_read_m    , only: atom , molecule , MM 
    use BuildReaxWAT , only: HOH => HOH_diss_parms
                                

    public :: DW_f_intra 

    private

    ! module variables
    logical :: done = .false.
    real*8  :: A, B, C
    real*8 , allocatable :: f_ang(:,:)

    ! module parameters ...
    real, parameter :: c1 = 1.1283791671d0  ! <== 2/sqrt(PI)

contains
!
!
!========================
 subroutine DW_f_intra ()
!========================
implicit none

! local variables
real*8 , allocatable  :: tmp_force(:,:,:)
real*8  :: rkl(3)
real*8  :: rkl2, force, erg, DWFF_intra
integer :: i, j, atk, atl, pair_of_kind, ithr, numthr
character(len=2) :: type1, type2

if( .not. done ) call set_local_parameters



call  evaluate_3body_DWFF




do j = 1 , MM % N_of_atoms
    atom(j)% f_DWFF(:) = D_zero  
end do
DWFF_intra = D_zero

numthr = OMP_get_max_threads()
allocate( tmp_force (MM%N_of_atoms,3,numthr) , source = D_zero )

do i = 1 , MM % N_of_molecules

     if ( .not. molecule(i)% DWFF ) cycle

     do j = 1 , molecule(i) % NintraIJ

        ithr = OMP_get_thread_num() + 1

        atk  = molecule(i) % IntraIJ(j,1) 
        atl  = molecule(i) % IntraIJ(j,2) 

        if ( any([atom(atk)%flex, atom(atl)%flex]) ) then

            rkl(:) = atom(atk) % xyz(:) - atom(atl) % xyz(:)
            rkl(:) = rkl(:) - MM % box(:) * DNINT( rkl(:) * MM % ibox(:) ) * PBC(:)

            rkl2 = sum( rkl(:)**2 )

            if ( rkl2 < rcutsq ) then

                type1 = atom(atk)% MMSymbol
                type2 = atom(atl)% MMSymbol

                select case (trim(type1)//'-'//trim(type2))
                case ('HX-HX')
                    pair_of_kind = 3
                case ('OX-OX')
                    pair_of_kind = 2
                case ('HX-OX' , 'OX-HX')
                    pair_of_kind = 1
                end select

                call evaluate_2body_DWFF( pair_of_kind , rkl2 , force , erg )

                tmp_force(atk,1:3,ithr) = tmp_force(atk,1:3,ithr) + force * rkl(:)
                tmp_force(atl,1:3,ithr) = tmp_force(atl,1:3,ithr) - force * rkl(:)

                DWFF_intra = DWFF_intra + erg*factor3 

            end if
         end if
     end do
end do

! force units = J/mts = Newtons ...
! manual reduction (+: fnonbd , fnonch) ...
do i = 1, MM % N_of_atoms
    do j=1,numthr
        atom(i)% f_DWFF(1:3) = atom(i)% f_DWFF(1:3) + tmp_force(i,1:3,j)
    enddo
    atom(i)% f_DWFF(:) = atom(i)% f_DWFF(:) + f_ang(i,:)
end do

deallocate( tmp_force , f_ang )

end subroutine DW_F_intra
!
!
!
!==============================================================
 subroutine evaluate_2body_DWFF( k , rkl2 , DW_force , DW_erg )
!==============================================================
implicit none
integer , intent(in)  :: k
real*8  , intent(in)  :: rkl2 
real*8  , intent(out) :: DW_force
real*8  , intent(out) :: DW_erg

! local variables ...
real*8 :: irkl , ir2 , ir6 , ir8 , rkl
real*8 :: zeta , erfc_zeta , arg , exp_arg2
real*8 :: c2, c3, U0 , Ecoul , Fcoul , f_sr , E_sr

rkl  = SQRT(rkl2)
irkl = D_one / rkl

!----------------------------------------------------
! SR (short-range) applies only to O-H pair ...
!----------------------------------------------------
zeta = rkl * B
erfc_zeta = erfc(zeta) / zeta

ir2 = D_one / rkl2
ir6 = ir2 * ir2 * ir2
ir8 = ir6 * ir2

! SR Energy
E_sr = A*erfc_zeta - C*ir6
! SR Force
f_sr = A*( erfc_zeta + c1*exp(-zeta**2) )*ir2 - SIX*C*ir8

!----------------------------------------------------
! Coulomb applies to O-H and H-H pairs
! O-H ==> k = 1
! H-H ==> k = 3
!----------------------------------------------------

c2  = irsqPI * HOH% Coul(k,4)   
arg = rkl * HOH%Coul(k,4)
exp_arg2 = EXP(-arg**2)

! Energy
U0 = HOH%Coul(k,1) + HOH%Coul(k,2)*erf(arg) + HOH%Coul(k,3)* erf(arg*sqrt2)
Ecoul = coulomb * U0 * irkl

! Force
! Fcoul (not damped)
c3 = c2 * ( HOH%Coul(k,2) + sqrt2*HOH%Coul(k,3) * exp_arg2 )
Fcoul = coulomb * (U0*ir2*irkl - c3*exp_arg2*ir2)

!----------------------------------------------------
! total 
!----------------------------------------------------
DW_erg   = E_sr + Ecoul
DW_force = f_sr + Fcoul

end subroutine evaluate_2body_DWFF
!
!
!
!==============================================================
 subroutine evaluate_3body_DWFF( )
!==============================================================
implicit none

! local parameter ...
real*8 :: r0 = 1.6

! local_variables ...
real*8 , dimension (3):: rij, rik 
real*8  :: rijq, rijsq, rikq, riksq, fxyz
real*8  :: cos_theta, U3, U03, exp_arg, exponential
real*8  :: c1, c2, c3, f_ij, f_ik, inv_delta_0ij, inv_delta_0ik
integer :: i, j, k, l, n, ati, atj, atk 

!================================
! Angle - bending potential ...
!             J     K
!              \   /
!               \ / 
!                I 
!================================

allocate( f_ang(MM% N_of_atoms,3) , source = D_zero )

do i = 1 , MM % N_of_molecules

     if ( .not. molecule(i)% DWFF ) cycle

        ! MIND: the atomic sequence is JIK, with ATOM I AT THE VERTEX
        call get_bond_triplet(i, ati, atj, atk)

        if(.not. any([atom(atj)%flex, atom(ati)%flex, atom(atk)%flex])) cycle

            ! MIND: the atomic sequence is JIK, with ATOM I AT THE VERTEX
            ! rij = r_j - r_i
            rij(:) = atom(atj) % xyz(:) - atom(ati) % xyz(:)
            rij(:) = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
            rijq   = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            rijsq  = SQRT(rijq)
            if ( (rijsq+milli) > r0 ) cycle
            ! rik = r_k - r_i 
            rik(:) = atom(atk) % xyz(:) - atom(ati) % xyz(:)
            rik(:) = rik(:) - MM % box(:)*DNINT( rik(:) * MM % ibox(:) ) * PBC(:)
            rikq   = rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3)
            riksq  = SQRT(rikq)
            if ( (riksq+milli) > r0) cycle

            cos_theta = (rij(1)*rik(1) + rij(2)*rik(2) + rij(3)*rik(3)) / ( rijsq * riksq )

            inv_delta_0ij = 1.d0/(r0-rijsq)
            inv_delta_0ik = 1.d0/(r0-riksq)

            exp_arg = HOH%Angle(1,3)*( inv_delta_0ij + inv_delta_0ik )
            exponential = exp(-exp_arg)

            U03 = HOH%Angle(1,1) * (cos_theta - HOH%Angle(1,4)) * exponential
            U3  = U03 * (cos_theta - HOH%Angle(1,4))

            c1 = U3*HOH%Angle(1,3)    
            c2 = two*U03
            c3 = c2*cos_theta
             
            f_ij = c1*inv_delta_0ij**2 + c3/rijsq
            f_ik = c2/rijsq
            f_ang(atj,:) = f_ij*rij(:)/rijsq - f_ik*rik(:)/riksq

            f_ik = c1*inv_delta_0ik**2 + c3/riksq
            f_ij = c2/riksq
            f_ang(atk,:) = f_ik*rik(:)/riksq - f_ij*rij(:)/rijsq

            f_ang(ati,:) = -f_ang(atj,:) -f_ang(atk,:)

end do


end subroutine evaluate_3body_DWFF
!
!
!
!======================================================
 subroutine get_bond_triplet(i, ati, atj, atk)
! Given intraIJ array, determine the OX/HX atom indices
! and return them in ati, atj, atk.
!======================================================
 implicit none
 integer, intent(in)  :: i
 integer, intent(out) :: ati, atj, atk ! selected indices

 ! local variables
 integer :: idx1, idx2, idx3, idx4

 associate( intraIJ => molecule(i)% intraIJ )
       ! Unpack once
       idx1 = intraIJ(1,1)
       idx2 = intraIJ(1,2)
       idx3 = intraIJ(2,1)
       idx4 = intraIJ(2,2)

       if (atom(idx1)%MMSymbol == "OX") then
           ati = idx1
           atj = idx2
           atk = merge(idx3, idx4, atom(idx3)%MMSymbol == "HX")

       elseif (atom(idx2)%MMSymbol == "OX") then
           ati = idx2
           atj = idx1
           atk = merge(idx3, idx4, atom(idx3)%MMSymbol == "HX")

       else
           atj = idx1
           atk = idx2
           ati = merge(idx3, idx4, atom(idx3)%MMSymbol == "OX")
       end if
 end associate

end subroutine get_bond_triplet








!
!
!
!
!
!===================
 function ERFC ( X )
!===================
! ERFC = (1 - ERF), complementary error function ...
 implicit none
 real*8 :: ERFC
 real*8 :: A1, A2, A3, A4, A5, P, T, X, XSQ, TP
 parameter ( A1 = 0.254829592d0, A2 = -0.284496736d0 ) 
 parameter ( A3 = 1.421413741d0, A4 = -1.453122027d0 ) 
 parameter ( A5 = 1.061405429d0, P  =  0.3275911d0   ) 

 T    = d_one / ( d_one + P * X )
 XSQ  = X * X
 TP   = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
 ERFC = TP * EXP ( -XSQ )

end function ERFC
!
!
!
!================================
subroutine set_local_parameters()
!================================
implicit none

!local variables
real*8 :: chrg_decay

A = HOH% SR(1,1)
B = HOH% SR(1,2)
C = HOH% SR(1,3)

done = .true.

end subroutine set_local_parameters
!
!
!
end module Reax_F_intra

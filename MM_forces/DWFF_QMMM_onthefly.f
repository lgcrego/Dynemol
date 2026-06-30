module DWFF_QMMM

use iso_fortran_env
use constants_m
use type_m      , only : warning
use md_read_m   , only : atom, MM, special_pair_mtx
use Build_DWFF  , only : HOH => HOH_diss_parms

    public :: net_charge_prods
    public :: qd_qd, mix_q_qd

    private

    ! module variables ...
    real*8, allocatable :: qd_qd(:,:)
    real*8, allocatable :: mix_q_qd(:,:)

contains
!
!
!
!======================================
subroutine net_charge_prods(net_charge)
!======================================
implicit none
real*8, intent(in) :: net_charge(:)

!local variables
integer              :: k, l
character(len=2)     :: type1, type2

if (size(net_charge) /= MM%n_of_atoms) then
    error stop "net_charge size mismatch in net_charge_prods"
end if

if( .not. allocated(qd_qd) ) then
    allocate( qd_qd   ( MM%n_of_atoms , MM%n_of_atoms ) )
    allocate( mix_q_qd( MM%n_of_atoms , MM%n_of_atoms ) )
end if 

qd_qd    = D_zero
mix_q_qd = D_zero

!========================================================
! qd_qd and mix_q_qd are stored as full matrices.
!
! qd_qd is symmetric.
!
! mix_q_qd is symmetric for HX-HX and OX-OX pairs,
! but not for HX-OX pairs.
!========================================================
associate( qd => net_charge)

    do k = 1, MM%N_of_atoms-1
        do l = k+1, MM%N_of_atoms
    
            ! only for DWFF special pairs ...
            if ( special_pair_mtx(k,l) /= 3 ) cycle
    
               qd_qd(k,l) = qd(k) * qd(l)
               qd_qd(l,k) = qd_qd(k,l)
    
               type1 = atom(k)% MMSymbol
               type2 = atom(l)% MMSymbol
        
               select case (trim(type1)//'-'//trim(type2))
                    case ('HX-HX')
                        mix_q_qd(k,l) = HOH% PointCharge_H * (qd(k) + qd(l))
                        mix_q_qd(l,k) = mix_q_qd(k,l)
    
                    case ('OX-OX')
                        mix_q_qd(k,l) = HOH% PointCharge_O * (qd(k) + qd(l))
                        mix_q_qd(l,k) = mix_q_qd(k,l)
    
                    case ('HX-OX' , 'OX-HX')
    
                        if (type1 == "HX") &
                        then
                            mix_q_qd(k,l) = HOH% PointCharge_O * qd(k) &
                                          + HOH% PointCharge_H * qd(l)
    
                            mix_q_qd(l,k) = HOH% PointCharge_O * qd(l) &
                                          + HOH% PointCharge_H * qd(k)
                        else
                            mix_q_qd(k,l) = HOH% PointCharge_O * qd(l) &
                                          + HOH% PointCharge_H * qd(k)
    
                            mix_q_qd(l,k) = HOH% PointCharge_O * qd(k) &
                                          + HOH% PointCharge_H * qd(l)
                        end if
    
                    case default
                        CALL warning(" unexpected DWFF pair in subroutine net_charge_prods ")
                        write(*,*) k,l,type1,type2
                        error stop
    
               end select
        end do
    end do

end associate

end subroutine net_charge_prods

end module DWFF_QMMM


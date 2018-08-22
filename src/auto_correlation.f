#include "GPU.h"

!
! Calculates the auto correlation function: AC(t) = <PSI(0)|PSI(t)>,
! and it's Fourier transform:              FAC(w) = F[AC(t)],
! using FFTW.
!
! Calculates the occupation of the adiabatic (MO) states, and
! sums over the elements of the density matrix.
!
! Alberto Torres
!

! Notes:
!  - Works only in the molecular representations
!  - Devs.: - Call Auto_Correlation_init after calculating bra(t=0),
!             then Auto_Correlation after calc. ket,
!             and finally, after the dynamics completes, Auto_Correlation_end.
!           - Call MO_Occupation after MO orbitals have been calculated.

! To do: 
!  - moving basis
!  - work with restart = true
!     - MO_Occupation, done
!     - WF_AutoCorrelation, partially

module Auto_Correlation_m

    use blas95
    use Matrix_Math

    use type_m

    use constants_m  , only : twopi, h_bar, pico_2_sec
    use parameters_m , only : n_part, n_t, t_i, t_f, DOS_range, restart, nuclear_matter
    
    public :: Auto_Correlation_init, Auto_Correlation_restart, Auto_Correlation, Auto_Correlation_end, &
              MO_Occupation
    
    private

    logical                 :: atoms_move

    integer, parameter      :: out_AC  = 70, &
                               out_Occ = 71, &
                               out_Coh = 72

    complex*16, allocatable :: PSI_0(:,:), &   ! initial WaveFunction (<bra|)
                               AC(:,:)         ! WF Auto-Correlation


contains


!=============================================
subroutine Auto_Correlation_init( basis, bra )
!=============================================
implicit none
type(STO_basis), intent(in) :: basis(:)
complex*16,      intent(in) :: bra(:,:)

! local variables
integer :: n

! dimension
n = size(basis)


!------------------
! do atoms move? Or: does the basis change over time?
if( nuclear_matter == "MDynamics" ) then 
    atoms_move = .true.
    write(*,'(a)') "Auto_Correlation_init: not implemented for 'moving' basis set. Results will be spurious."
else
    atoms_move = .false.
end if


allocate( PSI_0(n, n_part), source = bra )

!------------------
! save PSI_0 for restart (binary format)
open(out_AC, file='AC.restart', status='unknown', form='unformatted')
write(out_AC) PSI_0
close(out_AC)

!------------------
! open file
open(out_AC, file='Auto_Correlation.dat', status='unknown', form='formatted')

! header
if(n_part == 1 ) then
    write(out_AC,'(a)') "# time(ps)     mod       phase     Re       Im       of <e(0)|e(t)>"
else
    write(out_AC,'(a)') "# time(ps)    mod_e     mod_h   phase_e   phase_h     Re_e      Re_h      Im_e      Im_h     of <q(0)|q(t)>, q=e,h"
end if


!------------------
! Calculate auto correlation at t=t_i:
allocate( AC(n_t, n_part) )
call Auto_Correlation( basis, bra, t_i, 1 )

end subroutine Auto_Correlation_init
!
!
!
!=====================================================================
subroutine Auto_Correlation_restart( n )
!=====================================================================
implicit none
integer, intent(in) :: n

! read initial saved PSI_0
allocate( PSI_0(n, n_part) )

open(out_AC, file='AC.restart', status='old', form='unformatted', access='sequential', action='read')
read(out_AC) PSI_0
close(out_AC)

! read Autocorrelation
allocate( AC(n_t, n_part) )

! open file
open(out_AC, file='Auto_Correlation.dat', status='old', access='append')

end subroutine Auto_Correlation_restart
!
!
!
!===================================================
subroutine Auto_Correlation( basis_t, ket_t, t, it )
!===================================================
implicit none
type(STO_basis), intent(in) :: basis_t(:)
complex*16,      intent(in) :: ket_t(:,:)
real*8,          intent(in) :: t
integer,         intent(in) :: it

! local variables
integer                 :: n, i

! dimension
n = size(basis_t)

! calculate auto correlation: AC = <bra_0|ket_t>
do i = 1, n_part
    AC(it,i) = dotu( PSI_0(:,i), ket_t(:,i) )
end do

! write results
write(out_AC,'(7f10.6)') t, abs(AC(it,:)), datan2(dimag(AC(it,:)),dreal(AC(it,:))), AC(it,:)
call flush(out_AC)


end subroutine Auto_Correlation
!
!
!
!==============================
subroutine Auto_Correlation_end
!==============================
implicit none
include 'fftw3.f'

! local variables
integer   :: i, n
integer*8 :: plan
real*8    :: T, w
real*8, allocatable :: FAC(:,:)         ! holds the Fourier transform of the Auto Correlation

n = 2*n_t

! deallocate and allocate
deallocate( PSI_0 )
allocate( FAC(n,n_part) )

! calculate fourier transform with FFTW
do i = 1, n_part

    call dfftw_plan_dft_c2r_1d( plan, n, AC(:,i), FAC(:,i), FFTW_ESTIMATE )  ! this is FFTW_BACKWARD

    call dfftw_execute_dft_c2r( plan, AC(:,i), FAC(:,i) )

    call dfftw_destroy_plan( plan )

end do

FAC = FAC/n  ! normalize

close(out_AC)
open(out_AC, file='Auto_Correlation_Fourier.dat', status='unknown')

T = 2*(t_f - t_i)*pico_2_sec       ! in seconds

! available frequency output units: index, Hz, rad/s, s, and eV
do i = n/2+1, n
    w = i-n-1                                    ! index
    w = w/T                                      ! f = index/T - 1/s (Hz)
    w = twopi*w                                  ! w = 2.pi.f - rad/s
!     w = merge( 0.d0, 1.d0/w/pico_2_sec, w==0.d0 )   ! T = 1/f - period units
    w = (pico_2_sec*h_bar)*w                     ! E = h_bar.w - eV (photon's energy for freq. w)
    write(out_AC, '(3es18.8)') w, FAC(i,:)
end do

write(out_AC,*) ''
do i = 1, n/2
    w = i-1                                      ! index
    w = w/T                                      ! f = index/T - 1/s (Hz)
    w = twopi*w                                  ! w = 2.pi.f - rad/s
!     w = merge( 0.d0, 1.d0/w/pico_2_sec, w==0.d0 )   ! T = 1/f - period (T) units
    w = (pico_2_sec*h_bar)*w                     ! E = h_bar.w - eV (photon's energy for freq. w)
    write(out_AC, '(3es18.8)') w, FAC(i,:)
end do

close(out_AC)
deallocate( AC, FAC)
call system('rm -f AC.restart')

end subroutine Auto_Correlation_end
!
!
!
!==================================================================
subroutine MO_Occupation( t, bra, ket, UNI_el, UNI_hl )
! Calculates the occupation of each level (adiabatic basis set is assumed)
!==================================================================
implicit none
real*8,        intent(in) :: t
complex*16,    intent(in) :: bra(:,:), ket(:,:)
type(R_eigen), intent(in) :: UNI_el
type(R_eigen), intent(in), optional :: UNI_hl

! local variables
logical, save       :: first_call = .true., All_states = .false.
integer             :: i, j, n, nE, p
integer, save       :: iEmin, iEmax
real*8              :: Emin, Emax
real*8, allocatable :: erg(:,:), occup_erg(:,:)
real*8, allocatable :: coh_diag(:), coh_off(:)
complex*16          :: rho

! check
if( n_part==2 .and. .not.present(UNI_hl)) then
    stop 'ERROR: Auto_Correlation.f: MO_Occupation: n_part = 2, but UNI_hl is not present'
end if

! dimension
n  = size(bra,1)


!------------------
! find min and max indexes inside limits

if (All_states) then
    iEmin = 1
    iEmax = n
else
    Emin = DOS_range%inicio    ! only orbitals within [Emin, Emax] window will be considered
    Emax = DOS_range%fim

    if (first_call) then

        iEmin = max( 1, maxloc(UNI_el%erg, 1, UNI_el%erg< Emin) + 1 )    ! index of the first orbital inside the energy limit
        iEmax = min( n, maxloc(UNI_el%erg, 1, UNI_el%erg<=Emax) )        ! and the last one

        if (n_part == 2) then
            iEmin = min( iEmin, maxloc(UNI_hl%erg, 1, UNI_hl%erg< Emin) + 1 )
            iEmax = max( iEmax, maxloc(UNI_hl%erg, 1, UNI_hl%erg<=Emax) )
        end if
    end if
end if

! nr. of orbitals inside energy window
nE = iEmax - iEmin + 1

if (first_call) then
    write(*,*)
    write(*,'(a,i6,a,i6,a,i6)') "Orbitals considered to calc. occupation: [",iEmin,",",iEmax,"] ->", nE
    write(*,*)

    if(.not. restart) then
        open(out_Occ, file='Occupation.bin', status='unknown', form='unformatted')
        write(out_Occ) n_t, nE, n_part, n, iEmin

        open(out_Coh, file='Coherences.dat', status='unknown')
        write(out_Coh,'(a)') "# time (ps);  Sum |rho_ij|^2 for i,j:  all;  diag;  off-diag   for el      (hl)"
    else
        open(out_Occ, file='Occupation.bin', status='old', access='append', form='unformatted')
        open(out_Coh, file='Coherences.dat', status='old', access='append')
    end if
end if


!------------------
! get energies

allocate(erg( nE, n_part ))

erg(:,1) = UNI_el%erg(iEmin:iEmax)     ! el energies

if(n_part==2) &
erg(:,2) = UNI_hl%erg(iEmin:iEmax)     ! hl energies


!------------------
! calculate occupation per level...

allocate( occup_erg(nE,n_part) )

!$omp parallel private(p,i,j) shared(n_part,iEmin,iEmax)
!$omp do collapse(2)
do p = 1, n_part
do i = iEmin, iEmax
    j = i - iEmin + 1
    occup_erg(j,p) = dreal( bra(i,p)*ket(i,p) )    ! imag. part is zero by construction (checked)
end do
end do
!$omp end do
!$omp end parallel

! write to file
write(out_Occ) t, erg(:,:), occup_erg(:,:)


!------------------
! Coherences

allocate( coh_diag(n_part), coh_off(n_part) )

coh_diag = 0.d0
coh_off  = 0.d0
!$omp  parallel do private(i,j,p,rho) default(shared) &
!$omp& collapse(2) reduction(+:coh_diag,coh_off)
do p = 1, n_part
do i = 1, n

    rho = bra(i,p)*ket(i,p)
    coh_diag(p) = coh_diag(p) + rho*conjg(rho)      ! diagonal

    do j = i+1, n
        rho =  bra(i,p)*ket(j,p)
        coh_off(p) = coh_off(p) + rho*conjg(rho)    ! off-diagonal
    end do

end do
end do
!$omp end parallel do

coh_off = 2.d0*coh_off  ! We used the Hermitian symmetry of rho above, so the factor 2

if (n_part==1) then
    write(out_Coh,'(f11.6,3f20.16)') t, coh_diag+coh_off, coh_diag, coh_off
else
    write(out_Coh,'(f11.6,3f20.16,a,3f20.16)') t  , coh_diag(1)+coh_off(1), coh_diag(1), coh_off(1), &
                                               ' ', coh_diag(2)+coh_off(2), coh_diag(2), coh_off(2)
end if

deallocate( occup_erg, erg )
first_call = .false.

end subroutine MO_Occupation
!
!
!
end module Auto_Correlation_m

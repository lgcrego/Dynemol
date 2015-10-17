#include "GPU.h"

!
! Calculate the auto correlation function: AC(t) = <PSI(0)|PSI(t)>,
! and it's Fourier transform:              AC(w) = F[AC(t)],
! using FFTW
!
! Alberto Torres
!

! Notes:
!  - Works only in the molecular representations
!  - Devs.: call Auto_Correlation_init after calculating bra(t=0),
!           then Auto_Correlation after calc. ket,
!           and finally, after the dynamics completes, Auto_Correlation_end.

! To do: 
!  - moving basis
!  - work with restart = true

module Auto_Correlation_m

    use blas95
    use Matrix_Math

    use type_m
    
    use constants_m  , only : twopi, h_bar, pico_2_sec
    use parameters_m , only : n_part, n_t, t_i, t_f, DOS_range, restart, nuclear_matter

    logical                      :: atoms_move

    integer, parameter           :: AC_unit = 7, spct_unit = 8

    complex*16, allocatable      :: PSI_0(:,:), AC(:,:)

    type(structure)              :: system_0
    type(STO_basis), allocatable :: basis_0(:)


contains


!==================================================================
subroutine Auto_Correlation_init( system, basis, bra )
!==================================================================
implicit none
type(structure), intent(in) :: system
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
else
    atoms_move = .false.
end if

! if( atoms_move ) then
!     system_0 = system
!     allocate( basis_0(n), source = basis )
! end if

allocate( PSI_0(n, n_part), source = bra )


!------------------
! save PSI_0 for restart (binary format)
open(AC_unit, file='AC.restart', status='unknown', form='unformatted')
write(AC_unit) PSI_0
close(AC_unit)

!------------------
! open file
open(AC_unit, file='Auto_Correlation.dat', status='unknown', form='formatted')

! header
if(n_part == 1 ) then
    write(AC_unit,'(a)') "# time(ps)     mod       phase     Re       Im       of <e(0)|e(t)>"
else
    write(AC_unit,'(a)') "# time(ps)    mod_e     mod_h   phase_e   phase_h     Re_e      Re_h      Im_e      Im_h     of <q(0)|q(t)>, q=e,h"
end if


!------------------
! Calculate auto correlation at t=t_i:
allocate( AC(n_t, n_part) )
call Auto_Correlation( system, basis, bra, t_i, 1 )

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

open(AC_unit, file='AC.restart', status='old', form='unformatted', access='sequential', action='read')
read(AC_unit) PSI_0
close(AC_unit)

! read Autocorrelation
allocate( AC(n_t, n_part) )

! open file
open(AC_unit, file='Auto_Correlation.dat', status='old', access='append')

end subroutine Auto_Correlation_restart
!
!
!
!=============================================================
subroutine Auto_Correlation( system_t, basis_t, ket_t, t, it )
!=============================================================
implicit none
type(structure), intent(in) :: system_t
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
write(AC_unit,'(7f10.6)') t, abs(AC(it,:)), datan2(dimag(AC(it,:)),dreal(AC(it,:))), AC(it,:)
call flush(AC_unit)

! print*,"Auto_Correlation done"
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
integer   :: i, j, n
integer*8 :: plan
real*8    :: T, w
real*8, allocatable :: FAC(:,:)         ! holds the Fourier transform of the Auto Correlation
! complex*16, allocatable :: test(:)

n = 2*n_t

! deallocate and allocate
deallocate( PSI_0 )
allocate( FAC(n,n_part) )

do i = 1, n_part

    ! calculate fourier transform with FFTW
    call dfftw_plan_dft_c2r_1d( plan, n, AC(:,i), FAC(:,i), FFTW_ESTIMATE )  ! this is FFTW_BACKWARD
    call dfftw_execute_dft_c2r( plan, AC(:,i), FAC(:,i) )
    call dfftw_destroy_plan( plan )

!!----------------------------------------------------
!!   test if back-transform of the transform returns original data
!     allocate( test(n/2) )
!     call dfftw_plan_dft_r2c_1d( plan, n, FAC(:,i), test, FFTW_PRESERVE_INPUT )
!     call dfftw_execute_dft_r2c(plan, FAC(:,i), test )
!     call dfftw_destroy_plan(plan)
!     test = test/n
!     do j = 1, n/2
!         if( abs(test(j)-AC(j,i)) > 1e-12 ) then
!             print*, j, abs(test(j)-AC(j,i)), test(j), AC(j,i)
!             stop 'ERROR: Auto_Correlation_end: F^(-1)[ F[ f(t) ](w) ](t) != f(t) --> back-transform of transform differs from original.'
!         end if
!     end do
!     deallocate( test )
!!   ok: tested forwards (t_f>t_i) and backwards (t_f<t_i) in time! Keeping code here
!!----------------------------------------------------

end do

FAC = FAC/n  ! normalize

close(AC_unit)
open(AC_unit, file='Auto_Correlation_Fourier.dat', status='unknown')

T = 2*(t_f - t_i)*pico_2_sec       ! in seconds

! available frequency output units: index, Hz, rad/s, s, and eV
do i = n/2+1, n
    w = i-n-1                                    ! index
    w = w/T                                      ! f = index/T - 1/s (Hz)
    w = twopi*w                                  ! w = 2.pi.f - rad/s
!     w = merge( 0.d0, 1.d0/w/pico_2_sec, w==0.d0 )   ! T = 1/f - period units
    w = (pico_2_sec*h_bar)*w                     ! E = h_bar.w - eV (photon's energy for freq. w)
    write(AC_unit, '(3es18.8)') w, FAC(i,:)
end do

write(AC_unit,*) ''
do i = 1, n/2
    w = i-1                                      ! index
    w = w/T                                      ! f = index/T - 1/s (Hz)
    w = twopi*w                                  ! w = 2.pi.f - rad/s
!     w = merge( 0.d0, 1.d0/w/pico_2_sec, w==0.d0 )   ! T = 1/f - period (T) units
    w = (pico_2_sec*h_bar)*w                     ! E = h_bar.w - eV (photon's energy for freq. w)
    write(AC_unit, '(3es18.8)') w, FAC(i,:)
end do

close(AC_unit)
deallocate( AC, FAC)

end subroutine Auto_Correlation_end
!
!
!
!==================================================================
subroutine MO_Occupation( t, bra, ket, UNI_el, UNI_hl )
! Calculates the occupation of each level (adiabatic basis set is assumed)
!==================================================================
implicit none
real*8,          intent(in) :: t
complex*16,      intent(in) :: bra(:,:), ket(:,:)
type(R_eigen),   intent(in) :: UNI_el
type(R_eigen),   intent(in), optional :: UNI_hl

! local variables
logical, save             :: first_call = .true.
integer                   :: i, j, n, m, p
integer, save             :: iEmin, iEmax
real*8                    :: Emin, Emax
! real*8                    :: norm
real*8, allocatable       :: erg(:,:), occup_erg(:,:)

! check
if( n_part==2 .and. .not.present(UNI_hl)) then
        stop 'ERROR: Auto_Correlation.f: Occupation: n_part = 2, but UNI_hl is not present'
end if

n  = size(bra,1)        ! dimension


!------------------
! find min and max indexes inside limits

! only orbitals with E inside these limits will be considered here
Emin = DOS_range%inicio
Emax = DOS_range%fim

if (first_call) then

    iEmin = max( 1, maxloc(UNI_el%erg,1, UNI_el%erg< Emin) + 1 )    ! index of the first orbital inside the energy limit
    iEmax = min( n, maxloc(UNI_el%erg,1, UNI_el%erg<=Emax) )        ! and the last one

    if( n_part == 2) then
        iEmin = min(iEmin, maxloc(UNI_hl%erg,1, UNI_hl%erg< Emin) + 1 )
        iEmax = max(iEmax, maxloc(UNI_hl%erg,1, UNI_hl%erg<=Emax) )
    end if
end if

! nr. of orbitals inside energy limits
m = iEmax - iEmin + 1

if (first_call) then
    write(*,*)
    write(*,'(a,i,a,i,a)') "Orbitals considered to calc. occupation: [",iEmin,",",iEmax,"] ->", iEmax - iEmin + 1
end if

!------------------
! get energies

allocate(erg( m, n_part ))

              erg(:,1) = UNI_el%erg(iEmin:iEmax)     ! el energies
if(n_part==2) erg(:,2) = UNI_hl%erg(iEmin:iEmax)     ! hl energies


!------------------
! calculate occupation per level...

allocate( occup_erg(m,n_part) )

!$omp parallel private(p,i,j) shared(n_part,m)
!$omp do collapse(2)
do p = 1, n_part
do i = 1, m
    j = i + iEmin - 1
    occup_erg(i,p) = dreal( bra(j,p)*ket(j,p) )    ! imag. part is zero by construction (checked)
end do
end do
!$omp end do
!$omp end parallel


!------------------
! write to file

if( first_call .and. (.not.restart) ) then
    open(spct_unit, file='Occupation.bin', status='unknown', form='unformatted')
    write(spct_unit) n_t, m, n_part
else
    open(spct_unit, file='Occupation.bin', status='old', access='append', form='unformatted')
end if

write(spct_unit) t, erg, occup_erg

close(spct_unit)


!------------------
! test

! do p = 1, n_part
!     norm=0.d0
!     do i = 1, n
!         norm = norm + occup_erg(i,p)
!     end do
!     if( t == t_i ) write(*,'(a,i1,a,2es15.5)') "Occupation: norm(",p,")=",norm,1.d0-norm
!     if( abs(1.d0 - norm) > 1.d-12 ) then
!         write(*,'(a,i1,a,es15.5)') "Occupation: not normalized! norm(",p,")=",1.d0-norm
!         stop
!     end if
! end do


deallocate( occup_erg, erg )
first_call = .false.

end subroutine MO_Occupation
!
!
!
end module Auto_Correlation_m

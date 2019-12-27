!
! Calculates Vibrational spectral density and Density of states
! via the Fourier transform of the velocity autocorrelation function
!
! Alberto Torres
!


! input module: use and set in parameters_MM.f
module VDOS_input

!   Sample length: Nr. of steps in each VACF sample: Nr. of v(0)'s
    integer :: VDOS_Nsteps_per_sample = 10000   ! 0 -> turn off

end module VDOS_input


module VDOS_m

    use constants_m  , only : low_prec
    use MD_read_m    , only : MM, atom
    use MM_input     , only : species
    use constants_m  , only : twopi, h_bar, pico_2_sec
    use parameters_m , only : n_t, t_i, t_f, restart, DRIVER

    use VDOS_input   , Lsample => VDOS_Nsteps_per_sample
    
    public :: VDOS_init, VDOS_Correlation, VDOS_end
    
    private

    ! module types ... 
    type velocity
        real*8 :: vel(3)
    end type velocity
   
    ! module parameters ... 
    integer, parameter          :: in = 80, out = 81, out2 = 82
   
    ! module variables ... 
    integer                     :: Nres,      &    ! Number of residues (:= MM % N_of_species)
                                   Natoms,    &
                                   Nsamples,  &
                                   isample,   &    ! sample index = [1,Nsamples] (this is the number of different 'v0')
                                   frame_in_s      ! frame index inside isample-th sample

    type(velocity), allocatable :: v0(:), mw_v0(:)

    real*8, allocatable         :: v0_norm(:), v0_norm_mw(:)
    real*8, allocatable         :: VACF(:,:), VACF_mw(:,:)   ! Velocity Auto Correlation Function

contains


! macro to expand dot product
#define VEC3_DOT_PROD(v,u) ( v(1)*u(1) + v(2)*u(2) + v(3)*u(3) )

! Parallelize if Natoms > SERIAL_THRESHOLD
#define SERIAL_THRESHOLD 500



!================
subroutine Get_v0
! Get initial velocity
!================
implicit none

! local variables
integer :: i, nr

v0_norm    = 0.d0
v0_norm_mw = 0.d0

!$omp parallel do private(i,nr) default(shared) reduction(+:v0_norm,v0_norm_mw) if(Natoms > SERIAL_THRESHOLD)
do i = 1, Natoms

    nr = atom(i) % nr

    v0(i) % vel = atom(i) % vel + low_prec
    v0_norm(nr) = v0_norm(nr) + VEC3_DOT_PROD( v0(i)%vel, v0(i)%vel )
    
    mw_v0(i) % vel = ( v0(i) % vel )*( atom(i) % mass )
    v0_norm_mw(nr) = v0_norm_mw(nr) + VEC3_DOT_PROD( mw_v0(i)%vel, v0(i)%vel )

end do
!$omp end parallel do

v0_norm(0)    = sum( v0_norm   (1:Nres) )
v0_norm_mw(0) = sum( v0_norm_mw(1:Nres) )

end subroutine Get_v0
!
!
!
!===================
subroutine VDOS_init
!===================
implicit none

! local variables
integer :: i

if (DRIVER == "q_dynamics"  .or. & 
    DRIVER == "avrg_confgs" .or. &
    DRIVER == "Genetic_Alg" .or. &
    DRIVER == "diagnostic") Lsample = 0   ! turn off

! check Lsample's values
if (Lsample == 0) then

    return

elseif (Lsample > n_t) then          ! can't be greater than the number of steps

    Lsample = n_t

elseif (mod(n_t,Lsample)/=0) then    ! should be a perfect divisor of n_t: no fractional Nsamples

    do while (mod(n_t,Lsample)/=0)
        Lsample = Lsample + 1
    end do

end if

Nres = MM % N_of_species

Natoms = MM % N_of_atoms

Nsamples = n_t / Lsample             ! VACF will be averaged over "Nsamples" samples

write(*,*)
write(*,'(a)')    " Velocity Autocorrelation: VDOS_init:"
write(*,'(a,i6)') "   Total nr. of steps         =", n_t
write(*,'(a,i6)') "   Sample length (time steps) =", Lsample
write(*,'(a,i3)') "   Nr. samples                =", Nsamples
write(*,'(a,i3)') "   Nr. residues               =", Nres
write(*,*)

allocate(    v0(Natoms), v0_norm   (0:Nres) )
allocate( mw_v0(Natoms), v0_norm_mw(0:Nres) )

allocate(    VACF( Lsample, 0:Nres ), source = 0.d0 )
allocate( VACF_mw( Lsample, 0:Nres ), source = 0.d0 )

if (restart) then
    call VDOS_restart
    return
end if

call Get_v0

isample    = 1   ! sample index = [1,Nsamples] (this is the number of different 'v0's)
frame_in_s = 1   ! frame index inside isample-th sample

!------------------
! save v0 for restart (binary format)
open(out, file='VDOS.restart', status='unknown', form='unformatted')
write(out) v0_norm, v0_norm_mw, v0
close(out)

open(out, file='DOS.trunk/VACF.dat', status='unknown')
write(out,'(a)',advance='no') "# frame   VACF: total    "
do i = 1, Nres
    write(out,'(a)',advance='no') species(i)%residue // '            '
end do
write(out,'(a)',advance='no') "MW-VACF: tot   "
do i = 1, Nres
    write(out,'(a)',advance='no') species(i)%residue // '            '
end do
write(out,*)

end subroutine VDOS_init
!
!
!
!======================
subroutine VDOS_restart
!======================
implicit none

open(in, file='VDOS.restart', status='old', form='unformatted', access='sequential', action='read')
read(in) v0_norm, v0_norm_mw, v0
close(in)

open(out, file='DOS.trunk/VACF.dat', status='old', access='append')

end subroutine VDOS_restart
!
!
!
!===================================
subroutine VDOS_Correlation( frame )
!===================================
implicit none

! inputs
integer, intent(in) :: frame          ! global time step index

! local variables
integer :: i, j
real*8, allocatable :: auto_corr(:), auto_corr_mw(:)    ! mw -> mass-weighted


if (Lsample == 0) return

allocate( auto_corr   (0:Nres), source = 0.d0)
allocate( auto_corr_mw(0:Nres), source = 0.d0)

! Calculate <v(0).v(t)>
!$omp parallel do private(i,j) default(shared) reduction(+:auto_corr,auto_corr_mw) if(Natoms > SERIAL_THRESHOLD)
do i = 1, Natoms
    j = atom(i) % nr
    auto_corr(j)    = auto_corr(j)    + VEC3_DOT_PROD(    v0(i)%vel, atom(i)%vel )  ! auto-correlalation
    auto_corr_mw(j) = auto_corr_mw(j) + VEC3_DOT_PROD( mw_v0(i)%vel, atom(i)%vel )  ! mass-weighted auto-correlalation
end do
!$omp end parallel do

auto_corr   (0) = sum( auto_corr   (1:Nres) )
auto_corr_mw(0) = sum( auto_corr_mw(1:Nres) )

VACF   (frame_in_s,:) = VACF   (frame_in_s,:) + auto_corr   (:) / v0_norm   (:)     ! accumulate normalized VACF
VACF_mw(frame_in_s,:) = VACF_mw(frame_in_s,:) + auto_corr_mw(:) / v0_norm_mw(:)     ! accumulate normalized MW-VACF

write(out,'(i7)',advance='no') frame_in_s
do i = 0, Nres
    write(out,'(es15.6)',advance='no') VACF(frame_in_s,i)/isample
end do
do i = 0, Nres
    write(out,'(es15.6)',advance='no') VACF_mw(frame_in_s,i)/isample
end do
write(out,*)

if( mod(frame,Lsample) == 0 ) then
    call Get_v0
    isample = isample + 1
    frame_in_s = 1
    write(out,*)
    return
end if

frame_in_s = frame_in_s + 1

end subroutine VDOS_Correlation
!
!
!
!==============================
subroutine VDOS_end
!==============================
implicit none
include 'fftw3.f'

! local variables
integer   :: i, j
integer*8 :: plan
real*8    :: T, norm

real*8, parameter :: Hz_to_eV = twopi*pico_2_sec*h_bar

real*8, allocatable :: Fvacf(:,:)         ! Fourier transform of the Auto Correlation
real*8, allocatable :: VDOS(:,:)          ! Vibrational DOS
real*8, allocatable :: N_flex_atoms(:)


if (Lsample == 0) return

! release resources not needed anymore
deallocate(v0, mw_v0)
close(out)  ! close VACF.dat

! count number of flexible atoms in each residue (for proprer normalization)
allocate( N_flex_atoms(0:Nres), source = 0.d0 )

do i = 1, Natoms
    j = atom(i)%nr
    if (atom(i)%flex)  N_flex_atoms(j) = N_flex_atoms(j) + 1
end do

N_flex_atoms(0) = sum(N_flex_atoms(1:Nres))  ! total


!======= Vibrational Spectral Density ========

allocate( Fvacf(Lsample,0:Nres) )

VACF = VACF/Nsamples  ! normalize average over samples

! using Real-Even symmetry: cosine fourier transform (refer to FFTW manual)
call dfftw_plan_r2r_1d( plan, Lsample, VACF(:,0), Fvacf(:,0), FFTW_REDFT00, FFTW_ESTIMATE )

do i = 0, Nres
    call dfftw_execute_r2r( plan, VACF(:,i), Fvacf(:,i) )
end do


allocate( VDOS(Lsample,0:Nres) )  ! here, this is really the VSD (spectral density)

VDOS(:,:) = Fvacf(:,:)**2

! normalize
do i = 0, Nres
    norm = (3*N_flex_atoms(i))/sum(VDOS(:,i))
    VDOS(:,i) = VDOS(:,i)*norm
end do

deallocate( VACF ) ! not needed anymore

! Sample's Period (simulation time of each sample)
T = 2*(t_f - t_i)*pico_2_sec/Nsamples  ! in seconds (factor 2 is to account for the real-even symmetry used above)

!------------------
! write to VSD-*.dat
do j = 0, Nres
    if (j==0) then
        open(out, file='DOS.trunk/VSD-total.dat', status='unknown')
    else
        open(out, file='DOS.trunk/VSD-'//trim(species(j)%residue)//'.dat', status='unknown')
    end if
    write(out,'(a)') "# E (eV)   spectral density   "
    do i = 1, Lsample
        write(out, '(2es)') ((i-1)/T)*Hz_to_eV, VDOS(i,j)
    end do
    close(out)
end do


!======== Vibrational Density of States ========

VACF_mw = VACF_mw/Nsamples  ! normalize average over samples

do i = 0, Nres
    call dfftw_execute_r2r( plan, VACF_mw(:,i), Fvacf(:,i) )
end do

call dfftw_destroy_plan( plan )

VDOS(:,:) = Fvacf(:,:)**2

! normalize
do i = 0, Nres
    norm = (3*N_flex_atoms(i))/sum(VDOS(:,i))
    VDOS(:,i) = VDOS(:,i)*norm
end do

!------------------
! write to VDOS-*.dat
do j = 0, Nres
    if (j==0) then
        open(out, file='DOS.trunk/VDOS-total.dat', status='unknown')
    else
        open(out, file='DOS.trunk/VDOS-'//trim(species(j)%residue)//'.dat', status='unknown')
    end if
    write(out,'(a)') "# E (eV)   VDOS"
    do i = 1, Lsample
        write(out, '(2es)') ((i-1)/T)*Hz_to_eV, VDOS(i,j)
    end do
    close(out)
end do

! clean up
deallocate( VDOS, Fvacf, N_flex_atoms )
call system('rm -f VDOS.restart')  ! if the calculation finishes, there's no more need of this file

end subroutine VDOS_end
!
!
!
end module VDOS_m

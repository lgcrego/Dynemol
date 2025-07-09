!
! Calculates Vibrational spectral density and Density of states
! via the Fourier transform of the velocity autocorrelation function
!
! Alberto Torres
! checked and modified by LGC Rego 09/July/2025


module VDOS_m

    use VDOS_tuning
    use constants_m  , only : low_prec
    use MD_read_m    , only : MM, atom
    use MM_input     , only : species, MM_atomic, MM_molecular
    use constants_m  , only : twopi, h_bar, pico_2_sec
    use parameters_m , only : n_t, t_i, t_f, restart, DRIVER
    
    public :: VDOS_init, VDOS_Correlation, VDOS_end
    
    private

    ! module types ... 
    type velocity
        real*8 :: vel(3)
    end type velocity
   
    ! module parameters ... 
    integer, parameter :: in = 80, out = 81, out2 = 82
   
    ! module variables ... 
    logical                      :: starting_a_new_sample 
    integer                      :: Nspc,      &    ! Number of species
                                    Natoms,    &
                                    Lsample,   &
                                    Nsamples,  &
                                    isample,   &    ! sample index = [1,Nsamples] (this is the number of different 'v0')
                                    frame_in_s      ! frame index inside isample-th sample
    type(MM_atomic), allocatable :: local_atomic(:)
    type(MM_molecular), allocatable :: local_species(:)
    type(velocity) , allocatable :: v0(:), mw_v0(:)
    real*8, allocatable          :: v0_norm(:), v0_norm_mw(:)
    real*8, allocatable          :: VACF(:,:), VACF_mw(:,:)   ! Velocity Auto Correlation Function

contains

! macro to expand dot product
#define VEC3_DOT_PROD(v,u) ( v(1)*u(1) + v(2)*u(2) + v(3)*u(3) )

! Parallelize if Natoms > SERIAL_THRESHOLD
#define SERIAL_THRESHOLD 500



!================================================
 subroutine render_projection_rules
! set ensuing subroutines for atomic projection 
!================================================
implicit none

! local variables 
integer :: i

allocate( local_atomic(MM%N_of_atoms) )

select case (projection_rule)
       case(1) ! <== Symbol
             do i = 1 , size(my_symbols)
                  where( atom% Symbol == my_symbols(i) ) local_atomic% my_species = i
             end do
             Nspc = size(my_symbols)
             allocate( local_species(Nspc) )
             local_species% residue = my_symbols

       case(2) ! <== residue
             do i = 1 , size(my_residues)
                  where( atom% residue == my_residues(i) ) local_atomic% my_species = i
             end do
             Nspc = size(my_residues)
             allocate( local_species(Nspc) )
             local_species% residue = my_residues

       case(3) ! <== MMSymbol
             do i = 1 , size(my_MMSymbols)
                  where( atom% MMSymbol == my_MMSymbols(i) ) local_atomic% my_species = i
             end do
             Nspc = size(my_MMSymbols)
             allocate( local_species(Nspc) )
             local_species% residue = my_MMSymbols
end select

end subroutine render_projection_rules
!
!
!
!================
subroutine Get_v0
! Get initial velocity
!================
implicit none

! local variables
integer :: i, indx

v0_norm    = 0.d0
v0_norm_mw = 0.d0

!$omp parallel do private(i,indx) default(shared) reduction(+:v0_norm,v0_norm_mw) if(Natoms > SERIAL_THRESHOLD)
do i = 1, Natoms

    indx = local_atomic(i) % my_species

    v0(i) % vel = atom(i) % vel + low_prec
    v0_norm(indx) = v0_norm(indx) + VEC3_DOT_PROD( v0(i)%vel, v0(i)%vel )
    
    mw_v0(i) % vel = ( v0(i) % vel )*( atom(i) % mass )
    v0_norm_mw(indx) = v0_norm_mw(indx) + VEC3_DOT_PROD( mw_v0(i)%vel, v0(i)%vel )

end do
!$omp end parallel do

v0_norm(0)    = sum( v0_norm   (1:Nspc) )
v0_norm_mw(0) = sum( v0_norm_mw(1:Nspc) )

starting_a_new_sample = .false.

end subroutine Get_v0
!
!
!===================
subroutine VDOS_init
!===================
implicit none

! local variables
integer :: i

Lsample = Nsteps_per_sample

if (DRIVER == "q_dynamics"  .or. & 
    DRIVER == "avrg_confgs" .or. &
    DRIVER == "Genetic_Alg" .or. &
    DRIVER == "diagnostic") Lsample = 0   ! turn off

! define Lsample's length
if (Lsample == 0) then
    return
elseif (Lsample > n_t) then          ! can't be greater than the number of steps
    Lsample = n_t
elseif (mod(n_t,Lsample)/=0) then    ! should be a perfect divisor of n_t: no fractional Nsamples
    do while (mod(n_t,Lsample)/=0)
        Lsample = Lsample + 1
    end do
end if

Natoms   = MM % N_of_atoms
Nsamples = n_t / Lsample             ! VACF will be averaged over "Nsamples" samples

call render_projection_rules

! print checklist
write(*,*)
write(*,'(a)')    " Velocity Autocorrelation: VDOS_init:"
write(*,'(a,i8)') "   Total nr. of steps         =", n_t
write(*,'(a,i8)') "   Sample length (time steps) =", Lsample
write(*,'(a,i8)') "   Nr. samples                =", Nsamples
write(*,'(a,i8)') "   Nr. species                =", Nspc
write(*,*)

allocate(    v0(Natoms), v0_norm   (0:Nspc) )
allocate( mw_v0(Natoms), v0_norm_mw(0:Nspc) )

allocate(    VACF( Lsample, 0:Nspc ), source = 0.d0 )
allocate( VACF_mw( Lsample, 0:Nspc ), source = 0.d0 )

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

open(out, file='dos.trunk/VACF.dat', status='unknown')
write(out,'(a)',advance='no') "# frame   VACF: total    "
do i = 1, Nspc
    write(out,'(a)',advance='no') local_species(i)%residue // '            '
end do
write(out,'(a)',advance='no') "MW-VACF: tot   "
do i = 1, Nspc
    write(out,'(a)',advance='no') local_species(i)%residue // '            '
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

open(out, file='dos.trunk/VACF.dat', status='old', access='append')

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

if (starting_a_new_sample) call Get_v0

allocate( auto_corr   (0:Nspc), source = 0.d0)
allocate( auto_corr_mw(0:Nspc), source = 0.d0)

! Calculate <v(0).v(t)>
!$omp parallel do private(i,j) default(shared) reduction(+:auto_corr,auto_corr_mw) if(Natoms > SERIAL_THRESHOLD)
do i = 1, Natoms
    j = local_atomic(i) % my_species
    auto_corr(j)    = auto_corr(j)    + VEC3_DOT_PROD(    v0(i)%vel, atom(i)%vel )  ! auto-correlalation
    auto_corr_mw(j) = auto_corr_mw(j) + VEC3_DOT_PROD( mw_v0(i)%vel, atom(i)%vel )  ! mass-weighted auto-correlalation
end do
!$omp end parallel do

auto_corr   (0) = sum( auto_corr   (1:Nspc) )
auto_corr_mw(0) = sum( auto_corr_mw(1:Nspc) )

VACF   (frame_in_s,:) = VACF   (frame_in_s,:) + auto_corr   (:) / v0_norm   (0)     ! accumulate normalized VACF
VACF_mw(frame_in_s,:) = VACF_mw(frame_in_s,:) + auto_corr_mw(:) / v0_norm_mw(0)     ! accumulate normalized MW-VACF

write(out,'(i8)',advance='no') frame_in_s
do i = 0, Nspc
    write(out,'(es15.6)',advance='no') VACF(frame_in_s,i)/isample
end do
do i = 0, Nspc
    write(out,'(es15.6)',advance='no') VACF_mw(frame_in_s,i)/isample
end do
write(out,*)

if( mod(frame,Lsample) == 0 ) then
    starting_a_new_sample = .true.
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

! count number of flexible atoms in each species (for proprer normalization)
allocate( N_flex_atoms(0:Nspc), source = 0.d0 )

do i = 1, Natoms
    j = local_atomic(i)%my_species
    if (atom(i)%flex)  N_flex_atoms(j) = N_flex_atoms(j) + 1
end do

N_flex_atoms(0) = sum(N_flex_atoms(1:Nspc))  ! total

!======= Vibrational Spectral Density ========

allocate( Fvacf(Lsample,0:Nspc) )

VACF = VACF/Nsamples  ! normalize average over samples

! using Real-Even symmetry: cosine fourier transform (refer to FFTW manual)
call dfftw_plan_r2r_1d( plan, Lsample, VACF(:,0), Fvacf(:,0), FFTW_REDFT00, FFTW_ESTIMATE )
do i = 0, Nspc
    call dfftw_execute_r2r( plan, VACF(:,i), Fvacf(:,i) )
end do

allocate( VDOS(Lsample,0:Nspc) )  ! here, this is really the VSD (spectral density)

VDOS(:,:) = Fvacf(:,:)

! normalize
norm = (3*N_flex_atoms(0))/sum(VDOS(:,0))
do i = 0, Nspc
    VDOS(:,i) = VDOS(:,i)*norm
end do

deallocate( VACF ) ! not needed anymore

! Sample's Period (simulation time of each sample)
T = 2*(t_f - t_i)*pico_2_sec/Nsamples  ! in seconds (factor 2 is to account for the real-even symmetry used above)

!------------------
! write to VelPS-*.dat
do j = 0, Nspc
    if (j==0) then
        open(out, file='dos.trunk/VelPS-total.dat', status='unknown')
    else
        open(out, file='dos.trunk/VelPS-'//trim(local_species(j)%residue)//'.dat', status='unknown')
    end if
    write(out,'(a)') "# E (eV)   VelPS (spectral density of velocities)   "
    do i = 1, Lsample
        write(out, '(2es)') ((i-1)/T)*Hz_to_eV, VDOS(i,j)
    end do
    close(out)
end do

!======== Vibrational Density of States ========
! mass-weighted auto-correlalation

VACF_mw = VACF_mw/Nsamples  ! normalize average over samples

do i = 0, Nspc
    call dfftw_execute_r2r( plan, VACF_mw(:,i), Fvacf(:,i) )
end do

! Destroy the plan to free up resources
call dfftw_destroy_plan( plan )

VDOS(:,:) = Fvacf(:,:)

! normalize
norm = (3*N_flex_atoms(0))/sum(VDOS(:,0))
do i = 0, Nspc
    VDOS(:,i) = VDOS(:,i)*norm
end do

!------------------
! write to VDOS-*.dat
do j = 0, Nspc
    if (j==0) then
        open(out, file='dos.trunk/VDOS-total.dat', status='unknown')
    else
        open(out, file='dos.trunk/VDOS-'//trim(local_species(j)%residue)//'.dat', status='unknown')
    end if
    write(out,'(a)') "# E (eV)   VDOS (density of vibrational states)"
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

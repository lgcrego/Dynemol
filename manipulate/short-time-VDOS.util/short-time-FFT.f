module gauss
contains

    elemental real*8 function gaussian( x, s )
        ! inputs
        real*8, intent(in) :: x, s

        ! local
!         real*8, parameter :: inv_sqrt_2pi = 0.398942280401432677939946059934381868475858631164934657665925  ! 1/sqrt(2*pi)
!         real*8 :: inv_s

!         inv_s = 1.d0/s

        gaussian = dexp( -0.5d0*((x/s)**2) ) ! * inv_s*inv_sqrt_2pi
    end function gaussian

    elemental real*8 function igaussian( n, s )
        ! inputs
        integer, intent(in) :: n, s

!         ! local
!         real*8, parameter :: inv_sqrt_2pi = 0.398942280401432677939946059934381868475858631164934657665925  ! 1/sqrt(2*pi)
!         real*8 :: inv_s

!         inv_s = 1.d0/s

        igaussian = dexp( -0.5d0*dfloat(n**2)/dfloat(s)**2 )
    end function igaussian

end module gauss


program post_process

use gauss
implicit none
include 'fftw3.f'

character(len=64) :: option, val

character(len=1)    :: skip
character(len=128)  :: file_in, file_out

logical :: file_exists, normalize, file_is_given

integer, parameter  :: in = 7, out = 8
integer             :: n, m, i, j, k, k_thres, p, mi, mf

integer*8           :: plan

real*8              :: dT, avg, time_step, t0, fmax
real*8, allocatable :: t(:), s(:), ws(:), norm(:), f(:,:), fd(:,:)

complex*16, allocatable :: Fs(:)

real*8, parameter   :: Hz_2_eV = 4.135667516e-15

! options:
logical :: padding, find_maxima
integer :: jump                       ! number of time-steps to jump when outputting
real*8  :: Energy_threshold, sigma

! defaults:
Energy_threshold = 0.5d0  ! eV
sigma            = 75.d-3 ! ps
jump = 1
padding     = .false.
find_maxima = .false.

#define freq(i)  ( 1.d12*((i)-1)/dT )
#define erg(i)   ( freq(i)*Hz_2_eV )
#define time(i)  ( t0 + (i)*time_step )
#define log2(x)  ( dlog(x)/dlog(2.d0) )

file_is_given = .false.
file_in  = 'vacf.dat'
file_out = 'vsd.dat'


do i = 1, iargc()
    call getarg(i, option)

    select case (option)

        case ("-h","--help","-help")
            write(*,*)
            write(*,'(a)') "Options:"
            write(*,'(a)') "  -h                    this help"
            write(*,'(a)') "  -normalize            normalize all components of the spectrogram at t=0 to 1"
            write(*,'(a)') "  -in  <file1>          read input data from file (default is 'vacf.dat')"
            write(*,'(a)') "  -out <file2>          write output data to file (default is 'vsd.dat')"
            write(*,'(a)') "  -sigma <float>        the time width of the ST-FFT window (default= 75.d-3)"
            write(*,'(a)') "  -threshold <float>    freqs. with energy (eV, if time is in ps) above this are discarded (default= 0.5)"
            write(*,'(a)') "  -jump <int>           number of time-steps to jump when outputting (default= 1)"
            write(*,'(a)') "  -padding              make number of points = 2^n (best for FFT) by padding zeros"
            write(*,'(a)') "  -find_maxima          find the maxima (saved to E_x_t_max.dat)"
            write(*,*)
            stop

        case ("-normalize")
            normalize = .true.

        case ("-in")
            call getarg( i+1, file_in )

        case ("-out")
            call getarg( i+1, file_out )

        case ("-threshold")
            call getarg( i+1, val )
            read(val,*) Energy_threshold

        case ("-sigma")
            call getarg( i+1, val )
            read(val,*) sigma

        case ("-jump")
            call getarg( i+1, val )
            read(val,*) jump

        case ("-padding")
            padding = .true.

        case ("-find_maxima")
            find_maxima = .true.

    end select
end do

! check
inquire(file=file_in, exist=file_exists)
if (.not. file_exists) then
    write(*,'(3a)') "File ",file_in," doesn't exists. Skiping."
    stop
end if

! open and count
open(in, file=file_in, status='old', readonly)
read(in,*) !skip first line
n = 0
do
    read(in,*, end=666)
    n = n + 1
end do
666 rewind(in)
write(*,'(a,i6)') 'n = ',n

if (padding) then
    p = 2*n                             ! extra terms at the begining and at the end
    j = ceiling(log2(dfloat(n+2*p)))  ! power of 2 is better for FFT
    m = 2**j
    p = (m - n)/2                     ! readjusting p
else
    p = 0
    m = n
end if
k = m/2 + 1                       ! size of the Fourier r2c transform (refer to FFTW documentation)
mi = p + 1
mf = mi + n - 1                   ! = p+n

write(*,'(a,i6,a,i2)') 'm = ',m,' = 2^',j
write(*,'(a,i6)') 'p = ',p
write(*,'(a,i6)') 'k = ',k
write(*,'(a,2i6)') 'mf = ',mf,p+n

allocate( t(n), s(m) )
s = 0.d0

! read signal
read(in,*) !skip first line
do i = 1, n
    read(in,*) t(i), s(mi + i - 1)
end do
close(in)

t0 = t(1)
time_step = (t(n) - t0)/dfloat(n-1)
dT = m*time_step
deallocate(t)

! duplicate: s(-t) = s(t) for the terms at the begining
do i = 1, p
    s(i) = s(mi + p - i + 1)
end do

! duplicate: s(T+t) = s(T-t) for the terms at the end
do i = 1, p
    s(mf + i) = s(mf - i)
end do
    
!     open(out, file='test.dat', status='unknown', buffered='yes')
!     do i = 1, m
!         write(out,*) time(i-mi), s(i)
!     end do
!     close(out)

! spectogram: using short-time Fourier transform

open(out, file=file_out, status='unknown', buffered='yes')

allocate( ws(m), Fs(k) )
ws = 0.d0
    
k_thres = 1
do while ( erg(k_thres) <= Energy_threshold )
    k_thres = k_thres + 1
end do
write(*,*) "k,E threshold =", k_thres, erg(k_thres)

do i = mi, mf ! time loop

    do j = 1, m ! convolution loop
        ws(j) = s(j)*gaussian( time(j)-time(i), sigma )
    end do
        
    ! align signal with padded constant line
!     ws(1:n) = ws(1:n) - ws(n)

    ! Remove average, i.e., zero freq.
    avg = 0.d0
    do j = 1, m
        avg = avg + ws(j)
    end do
    avg = avg/m
    ws(:) = ws(:) - avg

    if (i==mi) call dfftw_plan_dft_r2c_1d( plan, m, ws, Fs, FFTW_MEASURE) !FFTW_ESTIMATE )
    call dfftw_execute_dft_r2c( plan, ws, Fs )
    if (i==mf) call dfftw_destroy_plan( plan )
        
    if (normalize) then
        if (i==mi) then
            write(*,*) "yes", i, mi
            allocate( norm(k_thres) )
            do j = 1, k_thres
                norm(j) = abs(Fs(j))**2
            end do
        end if
    else
        Fs = Fs/m
    end if

    if (mod(i,jump)==0) then
        do j = 1, k_thres
            if (normalize) then
                write(out,'(3es)') time(i-mi), erg(j), abs(Fs(j))**2/norm(j)
            else
                write(out,'(3es)') time(i-mi), erg(j), abs(Fs(j))**2
            end if
        end do
        write(out,*)
    end if

end do ! time loop
    
close(out)

! clean up
deallocate( s, Fs, ws )
if (normalize) deallocate( norm )


if (find_maxima) then

    m = (mf - mi + 1)/jump

    allocate( f(k_thres,m), source=0.d0 )    ! transposed!
    allocate( t(m), s(k_thres), fd(2,k_thres) )

    open(in, file=file_out, status='old', readonly)
    
    fmax=0.d0
    do i = 1, m
        do k = 1, k_thres
            read(in,*) t(i), s(k), f(k,i)
            fmax = max(fmax,f(k,i))
        end do
        read(in,*)
    end do
    
    close(in)

    time_step = (t(m) - t(1))/dfloat(m-1)
    dt = 12.d0*time_step
    
    write(*,*)
    write(*,*) "Find maxima:"
    write(*,*) "m,k=",m,k
    write(*,*) "time step=",time_step
    write(*,*) "t=",t(m-3:m)
    write(*,*) "f=",s(k_thres-3:k_thres)
    write(*,*) "fmax=", fmax
!     write(*,*) f(1:3,1:3)
    

    p = 5     ! order of the derivatives we will calculate
    p = p/2

    open(out, file='E_x_t_max.dat', status='unknown', buffered='yes')
    write(out,'(a)') "# E       t        VDOS"

    do i = 1, m
    
        ! first and second derivatives (in energy/freq.)
        !dir$ ivdep
        do k = 1+p, k_thres-p
            fd(1,k) = ( f(k-2,i) -  8.d0*f(k-1,i)                +  8.d0*f(k+1,i) - f(k+2,i) ) / dt
            fd(2,k) = (-f(k-2,i) + 16.d0*f(k-1,i) - 30.d0*f(k,i) + 16.d0*f(k+1,i) - f(k+2,i) ) / (dt*time_step)
        end do
        
        do k = 1+p+1, k_thres-p-1
            if ( fd(1,k-1)*fd(1,k+1)<0.d0 .and. (fd(2,k)<0.d0) ) then
                write(out,'(2f14.8,*(es16.8))') s(k), t(i), f(k,i) !, fd(:,k)
            end if
        end do
        write(out,*)
        
    end do
    
    close(out)
    
    deallocate( t, s, f )

end if

end program

! timing macros
#define TIME_SET    integer :: c1,c2,cr; call system_clock(count_rate=cr)
#define TIME_INIT   call system_clock(c1)
#define TIME_END(S) call system_clock(c2); print*, S,real(c2-c1)/real(cr)," s"; call flush(6)

! print macros
#define PRINT_INIT(S) write(*,*); write(*,*) '========================================================'; write(*,'(a)') S; write(*,*)
#define PRINT_END     write(*,*); write(*,*) '--------------------------------------------------------'; write(*,*)

! Max. for ifort compiler
#define BLK_SIZE 67108864



module Occupation

implicit none

public :: Post_proccess_Occupation

! module variables ...
real*8 :: it_step 

 contains

subroutine Post_proccess_Occupation
implicit none

integer, parameter   :: in = 7, out = 8
integer              :: n_t, n_grid1, n_grid2, nE, Ntotal, n_part, it, iE, iEmin, ip, isigma, i, j, k, ii, jj
integer              :: indexE, indexEi(2), indexEf(2), Nwindow, N_coh_list
integer, allocatable :: E_coh_ij_list(:)

real*8               :: sigma, dE_window, dt, dE
real*8               :: Ei, Ef, E_resolution, Ewindow, dEg1, dEg2, Ewmin, Ewmax
real*8, allocatable  :: t(:), ergs(:,:,:), occ(:,:,:), rho(:,:,:), occ_der(:,:,:), &
                        f(:), occ_grid(:,:,:)
           
character :: option,             &
             calc_occ_pl,        &
             calc_occ_der_pl,    &
             calc_coherences_pl, &
             calc_occ_smear,     &
             calc_rho_abxdE_ab,  &
             calc_rho_ab_smear

logical :: done, need_grid(2)

TIME_SET

need_grid = .false.
done      = .false.

 calc_occ_pl        = '1'
 calc_occ_der_pl    = '2'
 calc_coherences_pl = '3'
 calc_occ_smear     = '4'
 calc_rho_abxdE_ab  = '5'
 calc_rho_ab_smear  = '6'


do while (.not. done)

    call system( "clear" )

    write(*,'(/a)') '>>>    Choose what you want to calculate (multiple options allowed) and hit 0 to continue    <<<'

    write(*,'(/a,a)') calc_occ_pl,       '  : Occupation per MO                    -->  t x E_i x rho_ii'
    write(*,'(/a,a)') calc_occ_der_pl,   '  : Occupation change in time per MO     -->  t x E_i x d(rho_ii)/dt'
    write(*,'(/a,a)') calc_coherences_pl,'  : Coherences between selected MO`s     -->  t x i x j x |rho_ij| x {E_i(t=0) - E_j(t=0)}'
    write(*,'(/a,a)') calc_occ_smear,    '  : Nice occupations (option 1 smeared on a grid, will need grid_1)'
    write(*,'(/a,a)') calc_rho_abxdE_ab, '  : Coherences vs. E difference (grid_2) -->  t x (E_a - E_b) x |rho(E_a,E_b)|'
    write(*,'(/a,a)') calc_rho_ab_smear, '  : Nice Coherences vs. E difference (grid_2)'

    write(*,'(/a,a)') '0  : In the end, choose this option to carry out the selected calculations'

    write(*,'(/a,a)',advance='no') '>>>   '
    read (*,'(a,a)') option

    select case ( option )

        case ('0')
            done = .true.

        case ('1')
            calc_occ_pl = 'X'

        case ('2')
            calc_occ_der_pl = 'X'
            calc_occ_pl     = 'X'

        case ('3')
            calc_coherences_pl = 'X'
            calc_occ_pl        = 'X'

        case ('4')
            calc_occ_smear = 'X'
            calc_occ_pl    = 'X'
            need_grid(1)   = .true.

        case ('5')
            calc_rho_abxdE_ab = 'X'
            calc_occ_pl       = 'X'
            need_grid(:)      = .true.

        case ('6')
            calc_rho_ab_smear = 'X'
            calc_rho_abxdE_ab = 'X'
            calc_occ_pl       = 'X'
            need_grid(:)      = .true.

    end select
end do



!------------------------------------------------------------------------
TIME_INIT
PRINT_INIT(" Read data")

open(in, file='Occupation.bin', status='old', form='unformatted', access='sequential',  buffered='yes', blocksize=BLK_SIZE, readonly)
read(in) n_t, nE, n_part, Ntotal, iEmin

write(*,*) "nr. time steps = ", n_t
write(*,*) "nr. particles  = ", n_part
write(*,*) "nr. levels     = ", nE
write(*,*) "nr. levels tot = ", Ntotal
write(*,*) "first level    = ", iEmin
write(*,*)

allocate( t(n_t), ergs(nE,n_part,n_t), occ(nE,n_part,n_t) )   ! mind the indexing order of z
allocate( f(n_part) )

write(*,*) "Reading file..."
do it = 1, n_t
    read(in,end=666) t(it), ergs(:,:,it), occ(:,:,it)
end do

666 continue
if( it < n_t+1 ) then
    n_t = it-1
    write(*,*) "Data incomplete: new bounds:"
    write(*,*) "nr. time steps = ", n_t
    write(*,*) "nr. levels     = ", nE
    write(*,*) "nr. particles  = ", n_part
end if

close(in)
deallocate(f)
TIME_END("...done in ")
write(*,*)
!------------------------------------------------------------------------
! time step 
write(*,'(a)',advance='no') " >> time step for treating the data (it_step = integer value): "
read(*,*) it_step
write(*,*)

!------------------------------------------------------------------------
! Energy window
write(*,'(a)') " Define the energy window (only energies inside will be considered):"
write(*,'(a)',advance='no') " >> Energy window min. (in eV): "
read(*,*) Ewmin
write(*,'(a)',advance='no') " >> Energy window max. (in eV): "
read(*,*) Ewmax

write(*,*)
write(*,'(a)') "Finding MO level indices for Emin and Emax ..."

Nwindow = 0.d0

! find min and max indexes inside Ewindow
do ip = 1, n_part
    indexEi(ip) = maxloc(ergs(:,ip,1), 1, ergs(:,ip,1)< Ewmin)     ! index of the first orbital inside the energy window
    indexEf(ip) = maxloc(ergs(:,ip,1), 1, ergs(:,ip,1)<=Ewmax)     ! and the last one
    Nwindow = max(Nwindow, indexEf(ip) - indexEi(ip) + 1)
    if (ip==1) then
        write(*,*) "El:"
    else
        write(*,*) "Hl:"
    end if
    write(*,*) " Lower  state = ", indexEi(ip) + iEmin-1, ergs(indexEi(ip),ip,1)
    write(*,*) " Higher state = ", indexEf(ip) + iEmin-1, ergs(indexEf(ip),ip,1)
end do
write(*,*) "Nr. of states inside energy window = ", Nwindow

PRINT_END
!------------------------------------------------------------------------



!------------------------------------------------------------------------
if ( calc_occ_pl == 'X' ) then
    TIME_INIT
    PRINT_INIT(" Occupation (per level)")

    write(*,*) "Writing to Occup_el(hl).dat ..."

    open(out, file='Occup_el.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)

    do it = 1, n_t, it_step
    do iE = indexEi(1), indexEf(1)
        write(out,'(f11.6,i6,f11.6,2es13.3e3)') t(it), iE+iEmin-1, ergs(iE,1,it), occ(iE,1,it)
    end do;  write(out,*)
    end do

    close(out)

    if (n_part==2) then
        open(out, file='Occup_hl.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)

        do it = 1, n_t, it_step
        do iE = indexEi(2), indexEf(2)
            write(out,'(f11.6,i6,f11.6,2es13.3e3)') t(it), iE+iEmin-1, ergs(iE,2,it), occ(iE,2,it)
        end do;  write(out,*)
        end do

        close(out)
    end if

    TIME_END("...done in ")
    PRINT_END
end if
!------------------------------------------------------------------------



!------------------------------------------------------------------------
if ( calc_occ_der_pl == 'X' ) then

    TIME_INIT
    PRINT_INIT(" Derivative (per level): d(occ)/dt")

    dt = (t(n_t)-t(1))/n_t
    write(*,'(a,f11.8)') "  dt=",    dt

    allocate( occ_der(Nwindow,n_part,n_t) )

    write(*,*) "Calculating..."

    !$omp parallel do private(it,iE,ip) shared(occ,occ_der,dt)
    do it = 2, n_t-1
    do ip = 1, n_part
    do iE = indexEi(ip), indexEf(ip)
        j = iE - indexEi(ip) + 1
        occ_der(j,ip,it) = (occ(iE,ip,it+1) - occ(iE,ip,it-1))/(2.d0*dt)     ! centered difference
    end do
    end do
    end do
    !$omp end parallel do

    TIME_END("...done in ")

    TIME_INIT
    write(*,*) "Writing to Occup_der_el(hl).dat ... "

    open(out, file='Occup_der_el.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)
    write(out,'(a)') "# time (ps)    E_level      d(occ)/dt"

    do it = 2, n_t-1, it_step
        do iE = indexEi(1), indexEf(1)
            j = iE - indexEi(1) + 1
            write(out,'(2f11.6,1es18.8e3)') t(it), ergs(iE,1,it), occ_der(j,1,it)
        end do
        write(out,*)
    end do
    close(out)

    if (n_part==2) then
        open(out, file='Occup_der_hl.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)
        write(out,'(a)') "# time (ps)    E_level      d(occ)/dt"
        do it = 2, n_t-1, it_step
            do iE = indexEi(2), indexEf(2)
                j = iE - indexEi(2) + 1
                write(out,'(3f11.6,2es18.8e3)') t(it), ergs(iE,2,it), occ_der(j,2,it)
            end do
            write(out,*)
        end do
        close(out)
    end if

    deallocate( occ_der )

    TIME_END("...done in ")
    PRINT_END

end if
!------------------------------------------------------------------------



!------------------------------------------------------------------------
if ( calc_coherences_pl == 'X' ) then

    TIME_INIT
    PRINT_INIT(" Coherences (per level): |rho_ij|")

    ! get number of selected leves
    write(*,'(a)',advance='no') " >> Enter the number of levels: "
    read(*,'(i)') N_coh_list
    
    allocate( E_coh_ij_list(N_coh_list), f(n_part) )
    
    ! get their indices
    write(*,'(a)') " >> Enter the levels indices (from system-ergs.dat): "

    do i = 1, N_coh_list
        write(*,'(a,i2,a)',advance='no') '      # ',i,': '
        read(*,'(i)') E_coh_ij_list(i)
    end do

    write(*,*) "Calculating and writing to Sel_Coh.dat ... "

    open(out, file='Sel_Coh.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)
    write(out,'(a)') "# time (ps);  index i;   index j;    |rho_ij| for el(hl);   Ei - Ej (eV)"

    i = indexEi(ip)

    do it = 1, n_t, it_step
    do i  = 1, size(E_coh_ij_list)

        ii = E_coh_ij_list(i) - iEmin + 1
    
        do j  = i+1, size(E_coh_ij_list)
    
            jj = E_coh_ij_list(j) - iEmin + 1

            dE = ergs(jj,1,it) - ergs(ii,1,it)

            f(:) = occ(ii,:,it)*occ(jj,:,it)

            write(out,'(f11.6,2i6,2es13.3e3,f11.6)') t(it), E_coh_ij_list(i), E_coh_ij_list(j), f(:), dE

        end do
    end do;  write(out,*)
    end do

    TIME_END("...done in ")
    close(out)
    deallocate( E_coh_ij_list, f ) 
    PRINT_END

end if
!------------------------------------------------------------------------



!------------------------------------------------------------------------
if ( need_grid(1) ) then

    TIME_INIT
    PRINT_INIT(" Convert occuppation from MO-levels 'non-uniform grid' to uniform energy grid")

    ! get grid_1
    write(*,'(a)') "  Grid_1 is needed for energies:"
    write(*,'(a)',advance='no') " >> Enter the begining of the energy grid (in eV): "
    read(*,*) Ei
    write(*,'(a)',advance='no') " >> Enter the end of the energy grid (in eV): "
    read(*,*) Ef
    write(*,'(a)',advance='no') " >> Enter the resolution of the energy grid (in eV): "
    read(*,*) E_resolution

    Ewindow = Ef - Ei                      ! size of the energy grid
    n_grid1 = Ewindow/E_resolution         ! number of bins
    dEg1    = Ewindow/n_grid1              ! bin size (grid)

    write(*,*) "Calculating..."

    ! the grid
    allocate( occ_grid(n_grid1,n_part,n_t), source = 0.d0 )

    !$omp parallel do schedule(dynamic) private(it,iE,ip,indexE) shared(n_t,n_part,n_grid1,indexEi,indexEf,ergs,Ewindow,occ,occ_grid)
    do it = 1, n_t
    do ip = 1, n_part
    do iE = indexEi(ip), indexEf(ip)

        ! index conversion: [levels] erg(i) -> occ_grid(indexE) [grid]
        indexE = int( ((ergs(iE,ip,it) - Ei)/Ewindow)*n_grid1 )
    
        ! just an histogram (levels can coincide into the same bin)
        occ_grid(indexE,ip,it) = occ_grid(indexE,ip,it) + occ(iE,ip,it)

    end do
    end do
    end do
    !$omp end parallel do

    deallocate( ergs, occ )

    TIME_END("...done in ")
    PRINT_END

end if
!------------------------------------------------------------------------
#define E_GRID1(i) ( Ei + (i-1)*dEg1 )



!------------------------------------------------------------------------
if ( calc_occ_smear == 'X' ) then

    TIME_INIT
    PRINT_INIT(" Smear occupation (grid_1)")
    
    write(*,'(a)',advance='no') "  >> Enter the broadening width (sigma) in eV: "
    read(*,*) sigma

    isigma = int(10*sigma/dEg1 + 0.5)

    write(*,*) "Calculating and writing to Occup_smear.dat ..."

    open(out, file='Occup_smear.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)

#define SOFT(expr) dsqrt(expr)
! #define SOFT(expr) (expr)**(1.d0/3.d0)
! #define SOFT(expr) (log10(expr) + 1.d0)
! #define SOFT(expr) (expr)

    allocate( f(n_part) )

    do it = 1, n_t, it_step
    do iE = 1, n_grid1
    do ip = 1, n_part
        f(ip) = 0.d0
        do k = max(1, iE-isigma), min(n_grid1, iE+isigma)
            f(ip) = f(ip) + occ_grid(k,ip,it) * SOFT( dexp( -((iE-k)*dEg1)**2 / (2*sigma**2) ) )
        end do
    end do
    write(out,'(2f11.6,2es13.3e3)') t(it), E_GRID1(iE), f(:)
    end do
    write(out,*)
    end do

    close(out)
    deallocate(f)
    TIME_END("...done in ")
    PRINT_END

end if
!------------------------------------------------------------------------



!------------------------------------------------------------------------
if ( need_grid(2) ) then

    ! get grid_2
    write(*,'(a)') "  Grid_2 is needed for energy differences:"
    write(*,'(a)') "    The begining is at 0 eV."
    write(*,'(a)',advance='no') " >> Enter the max. Delta_E to consider (in eV): "
    read(*,*) dE_window                    ! maximum Delta_E to consider

    n_grid2 = dE_window/E_resolution       ! number of bins
    dEg2    = dE_window/n_grid2            ! bin size (grid)

end if
!------------------------------------------------------------------------
#define E_GRID2(i) ( (i-1)*dEg2 )



!------------------------------------------------------------------------
if ( calc_rho_abxdE_ab == 'X' ) then

    TIME_INIT
    PRINT_INIT(" Coherences (grid 2): |rho_ab| x Delta_E_ab")

    write(*,*) "Calculating..."

    allocate( rho(n_grid2,n_part,n_t), source = 0.d0 )

    !$omp parallel do schedule(dynamic) private(it,iE,ip,k,dE,j) shared(occ_grid,rho,dE_window)
    do it = 1, n_t-1
    do iE = 1, n_grid1
        do k = iE+1, n_grid1

            dE = E_GRID1(k) - E_GRID1(iE)

            if( dE>dE_window ) cycle

            j = int( (dE/dE_window)*n_grid2 + 1 )  ! index corresponding to dE on grid 2

            do ip = 1, n_part
                rho(j,ip,it) = rho(j,ip,it) + occ_grid(iE,ip,it)*occ_grid(k,ip,it)
            end do

        end do
    end do
    end do
    !$omp end parallel do

    TIME_END("...done in ")

    TIME_INIT
    write(*,*) "Writing to Occup_diff.dat ..."

    open(out, file='Occup_diff.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)

    do it = 1, n_t, it_step
    do iE = 1, n_grid2
        write(out,'(2f11.6,2es18.8e3)') t(it), E_GRID2(iE), rho(iE,:,it)
    end do;  write(out,*)
    end do

    close(out)
    deallocate( occ_grid )
    TIME_END("...done in ")
    PRINT_END

end if
!------------------------------------------------------------------------



!------------------------------------------------------------------------
if ( calc_rho_ab_smear == 'X' ) then

    TIME_INIT
    PRINT_INIT(" Smear Coherences vs. E difference (grid_2)")

    allocate(f(n_part))

    write(*,*) "Calculating..."

    isigma = int(20*sigma/dEg2 + 0.5)

    ! reusing array occ
    if(allocated(occ)) deallocate(occ)
    allocate( occ(n_grid2,n_part,n_t), source = 0.d0 )

#undef SOFT
#define SOFT(expr) dsqrt(expr)
! #define SOFT(expr) (expr)**(1.d0/3.d0)
! #define SOFT(expr) (log10(expr) + 1.d0)
! #define SOFT(expr) (expr)

    !$omp parallel do private(it,iE,ip,k) shared(isigma,occ,rho,sigma)
    do it = 1, n_t
    do iE = 1, n_grid2
    do ip = 1, n_part
        do k = max(1, iE-isigma), min(n_grid2, iE+isigma)
            occ(iE,ip,it) = occ(iE,ip,it) + rho(k,ip,it) * SOFT(dexp( -((iE-k)*dEg2)**2 / (2*sigma**2) ))
        end do
    end do
    end do
    end do
    !$omp end parallel do

    TIME_END("...done in ")
    write(*,*) "Writing to Occup_diff_smear.dat ..."

    open(out, file='Occup_diff_smear.dat', status='unknown', buffered='yes', blocksize=BLK_SIZE)

    do it = 1, n_t, it_step
    do iE = 1, nE
        write(out,'(2f11.6,2es13.3e3)') t(it), E_GRID2(iE), occ(iE,:,it)
    end do;  write(out,*)
    end do

    close(out)
    deallocate(f)
    TIME_END("...done in ")
    PRINT_END

end if
!------------------------------------------------------------------------

end subroutine Post_proccess_Occupation

end module Occupation

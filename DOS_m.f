module DOS_m

    use type_m
    use constants_m
    use Structure_Builder
    use Overlap_Builder
    use EHT_parameters , the_chemical_atom => atom

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

 contains
!
!
!
subroutine  TDOS(erg)

real*8  , ALLOCATABLE , intent(in)  :: erg(:)

real*8  , allocatable :: e_grid(:) , DOS_grid(:) , DOS(:)
real*8                :: erg_MO(size(erg))
real*8                :: gauss_norm , two_sigma2
integer               :: n_of_DOS_states
integer , parameter   :: npoints = 1500

gauss_norm = 1.d0 / (sigma*sqrt2PI)
two_sigma2 = 2.d0 * sigma*sigma
step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

! find the energies in the range [DOS_range%inicio,DOS_range%fim] 
j = 1
do i = 1 , size(erg)
    if( (erg(i) >= DOS_range%inicio) == (erg(i) <= DOS_range%fim) ) then
        erg_MO(j) = erg(i)
        j = j + 1
    end if
end do

! number of states in the range [DOS_range%inicio,DOS_range%fim]
n_of_DOS_states = j-1


ALLOCATE(e_grid(npoints), DOS_grid(npoints), DOS(npoints))

forall(k=1:npoints) e_grid(k) = (k-1)*step + DOS_range%inicio

! the total density of states
DOS(:) = 0.d0
do j = 1 , n_of_DOS_states

    forall(k=1:npoints) DOS_grid(k) = gauss_norm*dexp(-(e_grid(k)-erg_MO(j))*(e_grid(k)-erg_MO(j))/two_sigma2)

    DOS(:) = DOS(:) + DOS_grid(:)

end do

OPEN(unit=3,file='TDOS.dat',status='unknown')
    do i = 1 , npoints
        write(3,10) e_grid(i) , DOS(i)
    end do
CLOSE(3)

DEALLOCATE(e_grid,DOS_grid,DOS)

10   FORMAT(2F12.5)

print*, '>> TDOS done <<'

end subroutine TDOS
!
!
!
subroutine  PDOS(system,zL,zR,erg)

type(structure) , intent(in)               :: system
complex*16      , ALLOCATABLE , intent(in) :: zL(:,:) , zR(:,:)
real*8          , ALLOCATABLE , intent(in) :: erg(:)

real*8  , allocatable :: e_grid(:) , DOS_grid(:) , DOS(:)
real*8                :: gauss_norm , two_sigma2 , projection
integer               :: list_of_DOS_states(size(erg))
integer               :: i , j , n_of_atoms , ioerr
integer               :: atom(system%atoms) 
integer , parameter   :: npoints = 1500

gauss_norm = 1.d0 / (sigma*sqrt2PI)
two_sigma2 = 2.d0 * sigma*sigma
step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

ALLOCATE(e_grid(npoints), DOS_grid(npoints), DOS(npoints))

! find the energies in the range [DOS_range%inicio,DOS_range%fim]
j = 1
do i = 1 , size(erg)
    if( (erg(i) >= DOS_range%inicio) == (erg(i) <= DOS_range%fim) ) then
        list_of_DOS_states(j) = i
        j = j + 1
    end if
end do

! number of states in the range [DOS_range%inicio,DOS_range%fim]
n_of_DOS_states = j-1

! reads the list of atoms
atom(:) = 0
OPEN(unit=9,file='PDOS_list-of-atoms',status='unknown')

j=1
do i = 1 , system%atoms
    if( system%fragment(i) == 'M' ) then
        write(9,*) i , system%symbol(i)
        atom(j) = i
        j = j + 1
    end if
end do

CLOSE(9)

! number of atoms of species ....
n_of_atoms = j-1

forall(k=1:npoints) e_grid(k) = (k-1)*step + DOS_range%inicio

DOS(:) = 0.d0
do l = 1 , n_of_atoms

    i1   = system%BasisPointer(atom(l)) + 1
    i2   = system%BasisPointer(atom(l)) + the_chemical_atom(system%AtNo(atom(l)))%DOS 

    do n = 1 , n_of_DOS_states

        j = list_of_DOS_states(n)
        projection = 0.d0
        do i = i1 , i2
            projection = projection + zL(j,i)*zR(i,j)
        end do
        
        erg_MO = erg(list_of_DOS_states(n)) 
        forall(k=1:npoints) DOS_grid(k) = gauss_norm*dexp(-(e_grid(k)-erg_MO)*(e_grid(k)-erg_MO)/two_sigma2)

        DOS(:) = DOS(:) + projection*DOS_grid(:)

    end do

end do

OPEN(unit=3,file='PDOS.dat',status='unknown')
    do i = 1 , npoints
        write(3,10) e_grid(i) , DOS(i)
    end do
CLOSE(3)

DEALLOCATE(e_grid,DOS_grid,DOS)

10   FORMAT(2F10.5)

print*, '>> PDOS done <<'

end subroutine PDOS
!
!
!
end module DOS_m

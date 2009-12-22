module DOS_m

    use type_m
    use constants_m
    use Structure_Builder
    use Overlap_Builder
    use Semi_Empirical_Parms , the_chemical_atom => atom

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

 contains
!
!
!
!===================================================
subroutine  Total_DOS( erg , TDOS , internal_sigma )
!===================================================
real*8        , ALLOCATABLE , intent(in)    :: erg(:)
type(f_grid)                , intent(inout) :: TDOS
real*8        , OPTIONAL    , intent(in)    :: internal_sigma

! . local variables
real*8  , allocatable :: DOS_partial(:) , erg_MO(:)
real*8                :: gauss_norm , sgm , two_sigma2
integer               :: i1 , i2 , n_of_DOS_states

if( present(internal_sigma) ) then 
    sgm = internal_sigma
else
    sgm = sigma
end if

npoints = size(TDOS%grid)

gauss_norm = 1.d0  !  1.d0 / (sgm*sqrt2PI)  <== for gauss_norm = 1 the gaussians are not normalized ...
two_sigma2 = 2.d0 * sgm*sgm
step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

! states in the range [DOS_range%inicio,DOS_range%fim]
i1 = maxloc(erg , 1 , erg <= DOS_range%inicio) 
i2 = maxloc(erg , 1 , erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1

allocate( erg_MO(n_of_DOS_states) )

! find the energies in the range [DOS_range%inicio,DOS_range%fim]
erg_MO = erg( i1 : i2 )

allocate( DOS_partial(npoints) )

forall(k=1:npoints) TDOS%grid(k) = (k-1)*step + DOS_range%inicio

! the total density of states
TDOS%func(:) = 0.d0
do j = 1 , n_of_DOS_states

    forall(k=1:npoints) DOS_partial(k) = gauss_norm*dexp(-(TDOS%grid(k)-erg_MO(j))*(TDOS%grid(k)-erg_MO(j))/two_sigma2)

    TDOS%func(:) = TDOS%func(:) + DOS_partial(:)

end do

TDOS%average = TDOS%average + TDOS%func

DEALLOCATE(DOS_partial)

print*, '>> TDOS done <<'

end subroutine Total_DOS
!
!
!
!===========================================================================
subroutine  Partial_DOS( system, QM , PDOS , residue , nr , internal_sigma )
!===========================================================================
type(structure)             , intent(in)    :: system
type(eigen)                 , intent(in)    :: QM
type(f_grid)  , allocatable , intent(inout) :: PDOS(:)
character(3)                , intent(in)    :: residue
integer       , OPTIONAL    , intent(in)    :: nr
real*8        , OPTIONAL    , intent(in)    :: internal_sigma

! . local variables
real*8                :: gauss_norm , sgm , two_sigma2 , projection
real*8  , allocatable :: DOS_partial(:) 
integer , allocatable :: list_of_DOS_states(:)
integer               :: i , j , i1 , i2 , n_of_atoms , ioerr
integer               :: atom(system%atoms) 

if( present(internal_sigma) ) then 
    sgm = internal_sigma
else
    sgm = sigma
end if

npoints = size( PDOS(nr)%grid )

gauss_norm = 1.d0 !  / (sgm*sqrt2PI)    <== for gauss_norm = 1 the gaussians are not normalized ...
two_sigma2 = 2.d0 * sgm*sgm
step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

allocate( DOS_partial(npoints) )

! number of states in the range [DOS_range%inicio,DOS_range%fim]
i1 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%inicio) 
i2 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1

allocate( list_of_DOS_states(n_of_DOS_states) )

! states in the range [DOS_range%inicio,DOS_range%fim]
forall( i=1:n_of_DOS_states ) list_of_DOS_states(i) = i1 + i

! reads the list of atoms
atom(:) = 0
j=1
do i = 1 , system%atoms
    if( (system%residue(i) == RESIDUE) ) then
        atom(j) = i
        j = j + 1
    end if
end do

! number of atoms of species ....
n_of_atoms = j-1

forall(k=1:npoints) PDOS(nr)%grid(k) = (k-1)*step + DOS_range%inicio

PDOS(nr)%func(:) = 0.d0
do l = 1 , n_of_atoms

    i1   = system%BasisPointer(atom(l)) + 1
    i2   = system%BasisPointer(atom(l)) + the_chemical_atom(system%AtNo(atom(l)))%DOS 

    do n = 1 , n_of_DOS_states

        j = list_of_DOS_states(n)

        projection = 0.d0
        do i = i1 , i2
            projection = projection + QM%L(j,i)*QM%R(i,j)
        end do
        
        erg_MO = QM%erg(list_of_DOS_states(n)) 
        forall(k=1:npoints) DOS_partial(k) = gauss_norm*dexp(-(PDOS(nr)%grid(k)-erg_MO)*(PDOS(nr)%grid(k)-erg_MO)/two_sigma2)

        PDOS(nr)%func(:) = PDOS(nr)%func(:) + projection*DOS_partial(:)

    end do

end do

PDOS(nr)%average = PDOS(nr)%average + PDOS(nr)%func

DEALLOCATE(DOS_partial,list_of_DOS_states)

print*, '>> ',residue,' PDOS done <<'

end subroutine Partial_DOS
!
!
!
end module DOS_m

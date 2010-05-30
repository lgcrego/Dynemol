module DOS_m

    use type_m
    use constants_m
    use Structure_Builder
    use Overlap_Builder
    use Semi_Empirical_Parms , the_chemical_atom => atom

 contains
!
!
!
!===================================================
subroutine  Total_DOS( erg , TDOS , internal_sigma )
!===================================================
implicit none
real*8        , ALLOCATABLE , intent(in)    :: erg(:)
type(f_grid)                , intent(inout) :: TDOS
real*8        , OPTIONAL    , intent(in)    :: internal_sigma

! local variables ...
real*8  , allocatable :: DOS_partial(:) , erg_MO(:) , peaks(:)
real*8                :: gauss_norm , sgm , two_sigma2 , step , sub_occupation
integer               :: i1 , i2 , n_of_DOS_states , npoints, k ,j 

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

n_of_DOS_states = i2 - i1 + 1

allocate( erg_MO(n_of_DOS_states) )

! find the energies in the range [DOS_range%inicio,DOS_range%fim]
erg_MO = erg( i1 : i2 )

allocate( peaks      (npoints) )
allocate( DOS_partial(npoints) )

forall(k=1:npoints) TDOS%grid(k) = (k-1)*step + DOS_range%inicio

! the total density of states
TDOS%peaks(:) = 0.d0
TDOS%func (:) = 0.d0
do j = 1 , n_of_DOS_states

    peaks = 0.d0
    where( dabs(TDOS%grid-erg_MO(j)) < (step/two) ) peaks = D_one

    TDOS%peaks = TDOS%peaks + peaks

    DOS_partial = 0.d0
    where( ((TDOS%grid-erg_MO(j))**2/two_sigma2) < 25.d0 ) DOS_partial = gauss_norm*exp( -(TDOS%grid-erg_MO(j))**2 / two_sigma2 )

    TDOS%func(:) = TDOS%func(:) + DOS_partial(:)

end do

TDOS%average = TDOS%average + TDOS%func

! occupation of PDOS(nr) ...
TDOS%occupation(1) = two * TDOS%peaks(1)
do k = 2 , npoints 
    TDOS%occupation(k) = TDOS%occupation(k-1) + two*TDOS%peaks(k) 
end do
sub_occupation  = two * (i1 - 1)
TDOS%occupation = sub_occupation + TDOS%occupation

DEALLOCATE( peaks , DOS_partial )

print*, '>> TDOS done <<'

end subroutine Total_DOS
!
!
!
!===========================================================================
subroutine  Partial_DOS( system, QM , PDOS , nr , internal_sigma )
!===========================================================================
implicit none
type(structure)             , intent(in)    :: system
type(C_eigen)               , intent(in)    :: QM
type(f_grid)  , allocatable , intent(inout) :: PDOS(:)
integer       , OPTIONAL    , intent(in)    :: nr
real*8        , OPTIONAL    , intent(in)    :: internal_sigma

! local variables ...
real*8                :: gauss_norm , sgm , two_sigma2 , projection , step , erg_MO , sub_occupation
real*8  , allocatable :: peaks(:) , DOS_partial(:) 
integer , allocatable :: list_of_DOS_states(:)
integer               :: i , j , i1 , i2 , n_of_atoms , ioerr , npoints , k , l , n_of_DOS_states , n
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

allocate( peaks      (npoints) )
allocate( DOS_partial(npoints) )

! number of states in the range [DOS_range%inicio,DOS_range%fim] ...
i1 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%inicio) 
i2 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1 + 1

allocate( list_of_DOS_states(n_of_DOS_states) )

! states in the range [DOS_range%inicio,DOS_range%fim] ...
forall( i=1:n_of_DOS_states ) list_of_DOS_states(i) = i1 + (i-1)

! reads the list of atoms ...
atom(:) = 0
j=1
do i = 1 , system%atoms
    if( (system%residue(i) == PDOS(nr)%residue) ) then
        atom(j) = i
        j = j + 1
    end if
end do

! number of atoms of species residue ...
n_of_atoms = j-1

forall(k=1:npoints) PDOS(nr)%grid(k) = (k-1)*step + DOS_range%inicio

PDOS(nr)%peaks(:) = 0.d0
PDOS(nr)%func(:)  = 0.d0
do l = 1 , n_of_atoms

    i1   = system%BasisPointer(atom(l)) + 1
    i2   = system%BasisPointer(atom(l)) + the_chemical_atom(system%AtNo(atom(l)))%DOS 

    do n = 1 , n_of_DOS_states

        j = list_of_DOS_states(n)

        projection = 0.d0
        do i = i1 , i2
            projection = projection + QM%L(j,i)*QM%R(i,j)
        end do
        
        erg_MO = QM%erg(j) 

        peaks = 0.d0
        where( dabs(PDOS(nr)%grid-erg_MO) < (step/two) ) peaks = projection

        PDOS(nr)%peaks = PDOS(nr)%peaks + peaks

        DOS_partial = 0.d0
        where( ((PDOS(nr)%grid-erg_MO)**2/two_sigma2) < 25.d0 ) DOS_partial = gauss_norm*exp( -(PDOS(nr)%grid-erg_MO)**2 / two_sigma2 )

        PDOS(nr)%func(:) = PDOS(nr)%func(:) + projection*DOS_partial(:)

    end do

end do

PDOS(nr)%average = PDOS(nr)%average + PDOS(nr)%func

! occupation of PDOS(nr) ...
PDOS(nr)%occupation(1) = two * PDOS(nr)%peaks(1)
do k = 2 , npoints 
    PDOS(nr)%occupation(k) = PDOS(nr)%occupation(k-1) + two*PDOS(nr)%peaks(k) 
end do
sub_occupation = undernith_occupation(system,QM,atom,list_of_DOS_states(1))
PDOS(nr)%occupation = sub_occupation + PDOS(nr)%occupation

DEALLOCATE( peaks , DOS_partial , list_of_DOS_states )

print*, '>> ',PDOS(nr)%residue,' PDOS done <<'

end subroutine Partial_DOS
!
!
!
!=========================================================
 function  undernith_occupation( system, QM , atom , top )
!=========================================================
implicit none
type(structure)             , intent(in)    :: system
type(C_eigen)               , intent(in)    :: QM
integer       , OPTIONAL    , intent(in)    :: atom(:)
integer                     , intent(in)    :: top

! local variables ...
real*8                :: undernith_occupation
real*8  , allocatable :: state_projection(:) 
integer               :: j , l , i1 , i2 , n_of_atoms 

! state projection up to highest undernith state ...
allocate( state_projection(top-1) , source = 0.d0 )

! number of atoms of species = residue ...
n_of_atoms = size(atom)

do j = 1 , top-1

    do l = 1 , n_of_atoms

        i1 = system%BasisPointer(atom(l)) + 1
        i2 = system%BasisPointer(atom(l)) + the_chemical_atom(system%AtNo(atom(l)))%DOS 

        state_projection(j) = state_projection(j) + sum( QM%L(j,i1:i2)*QM%R(i1:i2,j) )
        
    end do

end do

! occupation of states up to state (top-1) ...
undernith_occupation = two * sum(state_projection)

deallocate( state_projection )

end function undernith_occupation
!
!
!
end module DOS_m

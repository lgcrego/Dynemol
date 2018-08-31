module PDOS_tool_m

    use type_m
    use omp_lib    
    use constants_m
    use parameters_m            , only : sigma , DOS_range        
    use Solvated_M              , only : DeAllocate_PDOS 
    use Semi_Empirical_Parms    , only : the_chemical_atom => atom

    public  :: Partial_DOS

    private

    ! modulo variables ...
    real*8                :: gauss_norm , two_sigma2 , step
    integer , allocatable :: list_of_DOS_states(:)
    integer , allocatable :: atom(:) 
    integer               :: n_of_DOS_states 

 contains
!
!
!
!=====================================================
 subroutine  Partial_DOS( system, QM , nr , instance )
!=====================================================
implicit none
type(structure) , intent(in)  :: system
type(R_eigen)   , intent(in)  :: QM
integer         , intent(in)  :: nr
character(*)    , intent(in)  :: instance

! local variables ...
real*8       , allocatable :: tmp_PDOS_peaks(:) , tmp_PDOS_func(:) 
real*8                     :: sgm , N_of_residues
integer                    :: i , i1 , i2 , j , n_of_atoms , npoints , k , l 
character(15)              :: string
type(f_grid) , allocatable :: PDOS(:)

CALL DeAllocate_PDOS( PDOS , flag="alloc" )

npoints = size( PDOS(nr)%grid )

sgm         = sigma
gauss_norm  = 1.d0 !  / (sgm*sqrt2PI)    <== for gauss_norm = 1 the gaussians are not normalized ...
two_sigma2  = 2.d0 * sgm*sgm

step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

! number of states in the range [DOS_range%inicio,DOS_range%fim] ...
i1 = maxloc(QM%erg , 1 , QM%erg <  DOS_range%inicio) + 1
i2 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1 + 1

allocate( list_of_DOS_states(n_of_DOS_states) )

! states in the range [DOS_range%inicio,DOS_range%fim] ...
forall( i=1:n_of_DOS_states ) list_of_DOS_states(i) = i1 + (i-1)

! reads the list of atoms ...
allocate( atom(system%atoms) , source=I_zero ) 
j=1
do i = 1 , system%atoms
    ! only quantum species contribute to PDOS ...
    if( (system%residue(i) == PDOS(nr)%residue) .AND. (system%QMMM(i) == "QM") ) then
        atom(j) = i
        j = j + 1
    end if
end do

! number of atoms of species residue ...
n_of_atoms = j-1

forall(k=1:npoints) PDOS(nr)%grid(k) = (k-1)*step + DOS_range%inicio

allocate( tmp_PDOS_peaks (npoints) , source = D_zero )
allocate( tmp_PDOS_func  (npoints) , source = D_zero )

!$OMP parallel 
    !$OMP DO reduction(+ : tmp_PDOS_peaks , tmp_PDOS_func )
    do l = 1 , n_of_atoms

        CALL tmp_PDOS( system , QM , l , PDOS(nr)%grid , tmp_PDOS_peaks , tmp_PDOS_func )

    end do
    !$OMP END DO
!$OMP end parallel

PDOS(nr)%peaks = tmp_PDOS_peaks
PDOS(nr)%func  = tmp_PDOS_func

deallocate( tmp_PDOS_peaks , tmp_PDOS_func )

PDOS(nr)%average = PDOS(nr)%average + PDOS(nr)%func

DEALLOCATE( atom , list_of_DOS_states )

! save PDOS ...
N_of_residues = size( PDOS )
string = "PDOS-"//instance//"-"//PDOS(nr)%residue//".dat" 
OPEN( unit=3 , file=string , status='unknown' )
    do i = 1 , size(PDOS(nr)%func)
        write(3,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i) , PDOS(nr)%peaks(i) 
    end do
CLOSE(3)

print*, '>> ',PDOS(nr)%residue,' PDOS done <<'

CALL DeAllocate_PDOS( PDOS , flag="dealloc" )

10   FORMAT(3F12.5)

end subroutine Partial_DOS
!
!
!
!=====================================================================================
subroutine tmp_PDOS( system , QM , l , tmp_PDOS_grid , tmp_PDOS_peaks, tmp_PDOS_func ) 
!=====================================================================================
implicit none
type(structure)  , intent(in)   :: system
type(R_eigen)    , intent(in)   :: QM
integer          , intent(in)   :: l
real*8           , intent(in)   :: tmp_PDOS_grid  (:)
real*8           , intent(out)  :: tmp_PDOS_peaks (:)
real*8           , intent(out)  :: tmp_PDOS_func  (:)

! local variables ...
real*8           :: projection , erg_MO 
integer          :: i , j , n , i1 , i2 , grid_size

grid_size = size( tmp_PDOS_grid )

i1 = system%BasisPointer(atom(l)) + 1
i2 = system%BasisPointer(atom(l)) + the_chemical_atom(system%AtNo(atom(l)))%DOS 

do n = 1 , n_of_DOS_states

    j = list_of_DOS_states(n)

    projection = 0.d0
    do i = i1 , i2
        projection = projection + QM%L(j,i)*QM%R(i,j)
    end do

    erg_MO = QM%erg(j) 

    do i = 1 , grid_size

        if(dabs(tmp_PDOS_grid(i)-erg_MO) < (step/two) ) tmp_PDOS_peaks(i) = tmp_PDOS_peaks(i) + projection 

        if( ((tmp_PDOS_grid(i)-erg_MO)**2/two_sigma2) < 25.d0 ) &
        tmp_PDOS_func(i) = tmp_PDOS_func(i) + projection*gauss_norm*exp( -(tmp_PDOS_grid(i)-erg_MO)**2 / two_sigma2 )

    end do    

end do

end subroutine tmp_PDOS
!
!
!
end module PDOS_tool_m

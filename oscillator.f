module Oscillator_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Structure_Builder
    use Overlap_Builder
    use EHT_parameters , chemical_atom => atom

 contains
!
!
!
subroutine  Optical_Transitions(system, basis, M_matrix, zL, zR, erg)

type(structure) , intent(in)  :: system
type(STO_basis) , intent(in)  :: basis(:)
type(dipole)    , intent(in)  :: M_matrix(:,:)
complex*16      , allocatable , intent(in)  :: zL(:,:) , zR(:,:)
real*8          , allocatable , intent(in)  :: erg(:)

integer :: i, j, xyz, j_bra, j_ket, dim_basis, dim_bra, dim_ket

integer      , dimension(:)   , allocatable :: bra_POINTER, ket_POINTER
real*8       , dimension(:,:) , allocatable :: a, bra, ket, matrix, R_vector, Transition_Strength
type(dipole) , dimension(:,:) , allocatable :: origin_Dependent, origin_Independent, Transition

real*8               :: resonance
real*8               :: gauss_norm , two_sigma2 , step
real*8 , allocatable :: e_grid(:) , peak_ij(:) , spec_peaks(:) , spec_broad(:)

integer , parameter  :: npoints = 1500
real*8  , parameter  :: one = 1.d0 , zero = 0.d0
real*8  , parameter  :: osc_const = 1.65338d-4  ! <== (2/3)*(m_e/hbar*hbar) ; unit = 1/( eV * (a_B)^2 ) 

!=====================================================================================
! define bra-&-ket states

dim_basis = size(zR(:,1))
dim_bra = count( (erg >= empty%inicio)    == (erg <= empty%fim)    )
dim_ket = count( (erg >= occupied%inicio) == (erg <= occupied%fim) )

allocate(bra_POINTER(dim_bra))
allocate(ket_POINTER(dim_ket)) 
j_bra = 1
j_ket = 1
do i = 1 , size(erg) 
    if( (erg(i) >= empty%inicio) == (erg(i) <= empty%fim) ) then
        bra_POINTER(j_bra) = i 
        j_bra = j_bra + 1
    else if( (erg(i) >= occupied%inicio) == (erg(i) <= occupied%fim) ) then
        ket_POINTER(j_ket) = i
        j_ket = j_ket + 1
    end if
end do

! atomic positions measured from the Center of Charge
allocate(R_vector(system%atoms,3))
forall(xyz=1:3) R_vector(:,xyz) = system%coord(:,xyz) - system%Center_of_Charge(xyz)

allocate( bra                 (dim_bra,dim_basis)   )
allocate( ket                 (dim_basis,dim_ket)   )
allocate( a                   (dim_bra,dim_basis)   )
allocate( matrix              (dim_basis,dim_basis) )
allocate( origin_Dependent    (dim_bra,dim_ket)     )
allocate( origin_Independent  (dim_bra,dim_ket)     )
allocate( Transition          (dim_bra,dim_ket)     )
allocate( Transition_Strength (dim_bra,dim_ket)     )

! Origin dependent DP = sum{C_dagger * vec{R} * S_ij * C}

ket = real( zR(:,ket_POINTER) )
do xyz = 1 , 3

    forall(i=1:dim_basis) bra(:,i) = real( zL(bra_POINTER,i) ) * R_vector(basis(i)%atom,xyz)

    forall(i=1:dim_bra,j=1:dim_ket) origin_Dependent(i,j)%DP(xyz) = sum( bra(i,:) * ket(:,j) )

end do

! Origin independent DP = sum{C_dagger * vec{M_matrix(i,j)} * C}

a   =  real( zL(bra_POINTER,:) )
!ket =  real( transpose(zL(ket_POINTER,:)) )
forall( i=1:dim_ket , j=1:dim_basis )  ket(j,i) =  real( zL(ket_POINTER(i),j) )

do xyz = 1 , 3  

    matrix = M_matrix%DP(xyz)

    CALL gemm(a,matrix,bra,'N','N',one,zero)    

    forall(i=1:dim_bra,j=1:dim_ket) origin_Independent(i,j)%DP(xyz) = sum( bra(i,:) * ket(:,j) )

end do

! Dipole Transition Matrix <bra|r|ket>

forall(i=1:dim_bra,j=1:dim_ket) 

    Transition(i,j)%DP = origin_Dependent(i,j)%DP + origin_Independent(i,j)%DP

    Transition_Strength(i,j) = sum(Transition(i,j)%DP**2)

end forall

deallocate(a,bra,ket,matrix,R_vector)
deallocate(origin_Dependent,origin_Independent,Transition)

!=====================================================================================
! the gaussians are not normalized

two_sigma2 = 2.d0 * sigma*sigma

step = (erg(bra_POINTER(dim_bra))-erg(ket_POINTER(1))) / float(npoints-1)

ALLOCATE(e_grid(npoints), peak_ij(npoints), spec_peaks(npoints), spec_broad(npoints))

forall(k=1:npoints) e_grid(k) = (k-1)*step 

! the optical spectrum : peaks and broadened lines
spec_peaks = 0.d0
spec_broad = 0.d0
do i=1,dim_bra
    do j=1,dim_ket 

        resonance = erg(bra_POINTER(i))-erg(ket_POINTER(j))
        Transition_Strength(i,j) = osc_const * resonance * Transition_Strength(i,j)

        peak_ij = 0.d0
        where( dabs(e_grid-resonance) < step ) peak_ij = Transition_Strength(i,j)
        spec_peaks = spec_peaks + peak_ij    

        peak_ij = 0.d0
        where( ((e_grid-resonance)**2/two_sigma2) < 25.d0 ) &
        peak_ij = dexp( -(e_grid-resonance)**2 / two_sigma2 )

        spec_broad = spec_broad + Transition_Strength(i,j) * peak_ij

    end do
end do

! save the peak and broadened specs
OPEN(unit=3,file='spectrum.dat',status='unknown')
    do i = 1 , npoints
        write(3,10) e_grid(i) , spec_broad(i) , spec_peaks(i)
    end do
10   FORMAT(3F13.9)
CLOSE(3)

deallocate(e_grid,peak_ij,spec_peaks,spec_broad)
deallocate(bra_POINTER,ket_POINTER,Transition_Strength)
!=====================================================================================

print*, '>> Optical Spectrum done <<'

end subroutine Optical_Transitions
!
!
!
end module Oscillator_m

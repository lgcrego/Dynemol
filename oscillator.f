module Oscillator_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Multipole_Core

    public :: Optical_Transitions , Transition_Dipole_Builder

    private

 contains
!
!
!
!-------------------------------------------------
subroutine  Optical_Transitions(system, basis, QM)
!-------------------------------------------------
type(structure) , intent(in) :: system
type(STO_basis) , intent(in) :: basis(:)
type(eigen)     , intent(in) :: QM

! . local variables: transition dipole
integer                        :: i , j , dim_bra , dim_ket
real*8           , allocatable :: Transition_Strength(:,:)

! . local variables: resonance spectrum
type(transition)               :: Trans_DP
real*8           , allocatable :: e_grid(:) , peak_ij(:) , spec_peaks(:) , spec_broad(:)
real*8                         :: gauss_norm , two_sigma2 , step , resonance
real*8           , parameter   :: one = 1.d0 , zero = 0.d0
real*8           , parameter   :: osc_const = 1.65338d-4  ! <== (2/3)*(m_e/hbar*hbar) ; unit = 1/( eV * (a_B)^2 ) 
integer          , parameter   :: npoints = 1500

!-------------------------------------------------------------
! . Dipole Transition Matrix <bra|r|ket> = <empty|r|occupied>
!-------------------------------------------------------------
trans_DP%bra_range = empty
trans_DP%ket_range = occupied

CALL Transition_Dipole_Builder(system, basis, QM, Trans_DP)

dim_bra = size(trans_DP%bra_POINTER)
dim_ket = size(trans_DP%ket_POINTER)

allocate( Transition_Strength (dim_bra,dim_ket) )

forall(i=1:dim_bra,j=1:dim_ket)  Transition_Strength(i,j) = sum(Trans_DP%matrix(i,j)%DP**2)

! . the gaussians are not normalized
two_sigma2 = 2.d0 * sigma*sigma

step = (QM%erg(trans_DP%bra_POINTER(dim_bra))-QM%erg(trans_DP%ket_POINTER(1))) / float(npoints-1)

ALLOCATE(e_grid(npoints), peak_ij(npoints), spec_peaks(npoints), spec_broad(npoints))

forall(k=1:npoints) e_grid(k) = (k-1)*step 

! . the optical spectrum : peaks and broadened lines
spec_peaks = 0.d0
spec_broad = 0.d0
do i=1,dim_bra
    do j=1,dim_ket 

        resonance = QM%erg(trans_DP%bra_POINTER(i)) - QM%erg(trans_DP%ket_POINTER(j))
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

! . save the peak and broadened specs
OPEN(unit=3,file='spectrum.dat',status='unknown')
    do i = 1 , npoints
        write(3,10) e_grid(i) , spec_broad(i) , spec_peaks(i)
    end do
10   FORMAT(3F13.9)
CLOSE(3)

deallocate(e_grid,peak_ij,spec_peaks,spec_broad)
deallocate(trans_DP%bra_POINTER,trans_DP%ket_POINTER,trans_DP%matrix,Transition_Strength)

print*, '>> Optical Spectrum done <<'

end subroutine Optical_Transitions
!
!
!
!
!------------------------------------------------------------
subroutine  Transition_Dipole_Builder(system, basis, QM, DP)
!------------------------------------------------------------
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
type(eigen)     , intent(in)    :: QM
type(transition), intent(inout) :: DP

! . local variables
integer                                        :: i, j, xyz, j_bra, j_ket, dim_basis, dim_bra, dim_ket, dim_bigger
real*8          , parameter                    :: one = 1.d0 , zero = 0.d0
real*8          , dimension(:,:) , allocatable :: a, left, right, matrix, R_vector
type(R3_vector) , dimension(:,:) , allocatable :: origin_Dependent, origin_Independent

! . define bra-&-ket states
select case (DP%flag)

    case('Redfield')

        DP%bra_indx_range%inicio = max( electrons%inicio , 1           )
        DP%ket_indx_range%inicio = max( holes    %inicio , 1           )
        DP%bra_indx_range%fim    = min( electrons%fim    , size(basis) )
        DP%ket_indx_range%fim    = min( holes    %fim    , size(basis) )

        allocate(DP%bra_POINTER(DP%bra_indx_range%fim - DP%bra_indx_range%inicio +1))
        allocate(DP%ket_POINTER(DP%ket_indx_range%fim - DP%ket_indx_range%inicio +1))

        DP%bra_POINTER = (/(i , i = DP%bra_indx_range%inicio,DP%bra_indx_range%fim)/) 
        DP%ket_POINTER = (/(i , i = DP%ket_indx_range%inicio,DP%ket_indx_range%fim)/) 
  
        dim_bra = size(DP%bra_POINTER)
        dim_ket = size(DP%ket_POINTER)

    case default

        dim_bra = count( (QM%erg >= DP%bra_range%inicio) == (QM%erg <= DP%bra_range%fim) )
        dim_ket = count( (QM%erg >= DP%ket_range%inicio) == (QM%erg <= DP%ket_range%fim) )

        allocate(DP%bra_POINTER(dim_bra))
        allocate(DP%ket_POINTER(dim_ket)) 
        j_bra = 1
        j_ket = 1
        do i  = 1 , size(QM%erg) 
            if( (QM%erg(i) >= DP%bra_range%inicio) == (QM%erg(i) <= DP%bra_range%fim) ) then
                DP%bra_POINTER(j_bra) = i 
                j_bra = j_bra + 1
            else if( (QM%erg(i) >= DP%ket_range%inicio) == (QM%erg(i) <= DP%ket_range%fim) ) then
                DP%ket_POINTER(j_ket) = i
                j_ket = j_ket + 1
            end if
        end do

end select

! . atomic positions measured from the Center of Charge
allocate(R_vector(system%atoms,3))
forall(xyz=1:3) R_vector(:,xyz) = system%coord(:,xyz) - system%Center_of_Charge(xyz)

dim_basis  = size(QM%R(:,1))
 
allocate( matrix              ( dim_basis  , dim_basis)  )
allocate( origin_Dependent    ( dim_bra    , dim_ket)    )
allocate( origin_Independent  ( dim_bra    , dim_ket)    )
allocate( DP%matrix           ( dim_bra    , dim_ket)    )

!...........................................................................
! . Origin dependent DP     =   sum{ C[ai] * {R_i + R_j}/2 * S_ij * C[jb] }  ....

! . Origin independent DP   =   sum{ C[ai] * {DP_matrix_AO[ij](i) + DP_matrix_AO[ji](j)}/2 * C[jb] }  ....
!...........................................................................

allocate( a     ( dim_bra , dim_basis)  )
allocate( left  ( dim_bra , dim_basis)  )
allocate( right ( dim_basis  , dim_ket) )

forall(i=1:dim_ket) right(:,i) = real( QM%R(:,DP%ket_POINTER(i)) )
do xyz = 1 , 3
    forall( j=1:dim_basis , i=1:dim_bra )  left(i,j) = real( QM%L(DP%bra_POINTER(i),j) ) * R_vector(basis(j)%atom,xyz) / two
    forall( j=1:dim_ket   , i=1:dim_bra )  origin_Dependent(i,j)%DP(xyz) = sum( left(i,:) * right(:,j) )
end do

forall( i=1:dim_bra )                  a(i,:)     =  real( QM%L(DP%bra_POINTER(i),:) )
forall( i=1:dim_ket , j=1:dim_basis )  right(j,i) =  real( QM%L(DP%ket_POINTER(i),j) )
do xyz = 1 , 3  
    matrix = DP_matrix_AO%DP(xyz)
    CALL gemm(a,matrix,left,'N','N',one,zero)    
    forall( j=1:dim_ket , i=1:dim_bra ) origin_Independent(i,j)%DP(xyz) = sum( left(i,:) * right(:,j) ) / two
end do

deallocate( a , left , right )
!...........................................................................

allocate( a     ( dim_ket , dim_basis)  )
allocate( left  ( dim_ket , dim_basis)  )
allocate( right ( dim_basis  , dim_bra) )

forall(i=1:dim_bra) right(:,i) = real( QM%R(:,DP%bra_POINTER(i)) )
do xyz = 1 , 3
    forall( j=1:dim_basis , i=1:dim_ket )  left(i,j) = real( QM%L(DP%ket_POINTER(i),j) ) * R_vector(basis(j)%atom,xyz) / two
    forall( j=1:dim_ket   , i=1:dim_bra )  origin_Dependent(i,j)%DP(xyz) = origin_Dependent(i,j)%DP(xyz) + sum( left(j,:) * right(:,i) )
end do

forall( i=1:dim_ket )                  a(i,:)     =  real( QM%L(DP%ket_POINTER(i),:) )
forall( i=1:dim_bra , j=1:dim_basis )  right(j,i) =  real( QM%L(DP%bra_POINTER(i),j) )
do xyz = 1 , 3  
    matrix = DP_matrix_AO%DP(xyz)
    CALL gemm(a,matrix,left,'N','N',one,zero)    
    forall( j=1:dim_ket , i=1:dim_bra ) origin_Independent(i,j)%DP(xyz) = origin_Independent(i,j)%DP(xyz) + sum( left(j,:) * right(:,i) ) / two
end do

deallocate( a , left , right )
!...........................................................................

! . Dipole Transition Matrix <bra|r|ket>
forall( i=1:dim_bra , j=1:dim_ket )  DP%matrix(i,j)%DP = origin_Dependent(i,j)%DP + origin_Independent(i,j)%DP

deallocate(matrix,R_vector)
deallocate(origin_Dependent,origin_Independent)

end subroutine Transition_Dipole_Builder
!
!
!
end module Oscillator_m

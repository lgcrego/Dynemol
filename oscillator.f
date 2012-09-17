module Oscillator_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use DP_main_m                   , only  : DP_matrix_AO
    use parameters_m                , only  : empty ,               &
                                              occupied ,            &
                                              rho_range ,           &
                                              sigma

    public :: Optical_Transitions   ,  Transition_Dipole_Builder

    private

 contains
!
!
!
!=============================================================================
subroutine  Optical_Transitions( system , basis , QM , SPEC , internal_sigma )
!=============================================================================
type(structure)                , intent(in)     :: system
type(STO_basis)                , intent(in)     :: basis(:)
type(R_eigen)                  , intent(in)     :: QM
type(f_grid)                   , intent(inout)  :: SPEC
real*8          , OPTIONAL     , intent(in)     :: internal_sigma

! . local variables: transition dipole
integer                        :: i , j , dim_bra , dim_ket
real*8           , allocatable :: Transition_Strength(:,:)

! . local variables: resonance spectrum
type(transition)               :: Trans_DP
real*8           , allocatable :: SPEC_peaks(:) , SPEC_func(:)
real*8                         :: gauss_norm , sgm , two_sigma2 , step , resonance , osc_const
real*8           , parameter   :: one = 1.d0 , zero = 0.d0
real*8           , parameter   :: osc_const_parameter = 1.65338d-4  ! <== (2/3)*(m_e/h_bar*h_bar) ; unit = 1/( eV * (a_B)^2 ) 

!-------------------------------------------------------------
! . Dipole Transition Matrix <bra|r|ket> = <empty|r|occupied>
!-------------------------------------------------------------


osc_const = osc_const_parameter

npoints = size( SPEC%grid )

trans_DP%bra_range = empty
trans_DP%ket_range = occupied

CALL Transition_Dipole_Builder(system, basis, QM, Trans_DP)

dim_bra = size(trans_DP%bra_PTR)
dim_ket = size(trans_DP%ket_PTR)

allocate( Transition_Strength (dim_bra,dim_ket) )

forall(i=1:dim_bra,j=1:dim_ket)  Transition_Strength(i,j) = sum(Trans_DP%matrix(i,j)%DP**2)


! . the gaussians are not normalized ...
if( present(internal_sigma) ) then 
    sgm = internal_sigma
else
    sgm = sigma
end if
two_sigma2 = 2.d0 * sgm*sgm

step = (QM%erg(trans_DP%bra_PTR(dim_bra))-QM%erg(trans_DP%ket_PTR(1))) / float(npoints-1)

forall(k=1:npoints) SPEC%grid(k) = (k-1)*step 

! . the optical spectrum : peaks and broadened lines ...
allocate( SPEC_peaks(npoints) , source = 0.d0 )
allocate( SPEC_func (npoints) , source = 0.d0 )

!$OMP parallel do private( resonance ) reduction( + : SPEC_peaks , SPEC_func ) 
do j=1,dim_ket 
    do i=1,dim_bra

        resonance = QM%erg(trans_DP%bra_PTR(i)) - QM%erg(trans_DP%ket_PTR(j))
        Transition_Strength(i,j) = osc_const * resonance * Transition_Strength(i,j)

        do k = 1 , npoints

            if( dabs(SPEC%grid(k)-resonance) < step ) SPEC_peaks(k) = SPEC_peaks(k) + Transition_Strength(i,j)

            if( ((SPEC%grid(k)-resonance)**2/two_sigma2) < 25.d0 ) &
            SPEC_func(k) = SPEC_func(k) + Transition_Strength(i,j) * dexp( -(SPEC%grid(k)-resonance)**2 / two_sigma2 )

        end do

    end do
end do
!$OMP end parallel do

SPEC%peaks = SPEC_peaks
SPEC%func  = SPEC_func

SPEC%average = SPEC%average + SPEC%func


deallocate( SPEC_peaks , SPEC_func )



deallocate( trans_DP%bra_PTR , trans_DP%ket_PTR , trans_DP%matrix , Transition_Strength )




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
type(R_eigen)   , intent(in)    :: QM
type(transition), intent(inout) :: DP

! . local variables
integer                                        :: i, j, xyz, j_bra, j_ket, dim_basis, dim_bra, dim_ket
real*8          , parameter                    :: one = 1.d0 , zero = 0.d0
real*8          , dimension(:,:) , allocatable :: a, left, right, matrix, R_vector
type(R3_vector) , dimension(:,:) , allocatable :: origin_Dependent, origin_Independent

! . define bra-&-ket states
select case (DP%flag)

    case('Redfield')

        DP % bra_indx_range % inicio = max( rho_range % inicio , 1           )
        DP % bra_indx_range % fim    = min( rho_range % fim    , size(basis) )
        dim_bra                      = DP % bra_indx_range % fim - DP % bra_indx_range % inicio + 1
   
        DP % ket_indx_range % inicio = DP % bra_indx_range % inicio 
        DP % ket_indx_range % fim    = DP % bra_indx_range % fim    
        dim_ket                      = dim_bra

        allocate( DP % bra_PTR (dim_bra) )
        allocate( DP % ket_PTR (dim_ket) )

        DP % bra_PTR = [ (i , i = DP % bra_indx_range % inicio , DP % bra_indx_range % fim) ] 
        DP % ket_PTR = [ (i , i = DP % ket_indx_range % inicio , DP % ket_indx_range % fim) ] 
  
    case default

        dim_bra = count( (QM%erg >= DP%bra_range%inicio) .AND. (QM%erg <= DP%bra_range%fim) )
        dim_ket = count( (QM%erg >= DP%ket_range%inicio) .AND. (QM%erg <= DP%ket_range%fim) )

        allocate(DP%bra_PTR(dim_bra))
        allocate(DP%ket_PTR(dim_ket)) 
        j_bra = 1
        j_ket = 1
        do i  = 1 , size(QM%erg) 
            if( (QM%erg(i) >= DP%bra_range%inicio) .AND. (QM%erg(i) <= DP%bra_range%fim) ) then
                DP%bra_PTR(j_bra) = i 
                j_bra = j_bra + 1
            else if( (QM%erg(i) >= DP%ket_range%inicio) .AND. (QM%erg(i) <= DP%ket_range%fim) ) then
                DP%ket_PTR(j_ket) = i
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

!..........................................................................................................
! . Origin dependent DP     =   sum{ C[ai] * {R_i + R_j}/2 * S_ij * C[jb] }  ....

! . Origin independent DP   =   sum{ C[ai] * {DP_matrix_AO[ij](i) + DP_matrix_AO[ji](j)}/2 * C[jb] }  ....
!..........................................................................................................

allocate( a     ( dim_bra , dim_basis)  )
allocate( left  ( dim_bra , dim_basis)  )
allocate( right ( dim_basis  , dim_ket) )

forall(i=1:dim_ket) right(:,i) = QM%R(:,DP%ket_PTR(i)) 

!$OMP parallel 
    do xyz = 1 , 3

        !$OMP single
        do j = 1 , dim_basis 
            !$OMP task shared(left,QM,R_vector)
            do i = 1 , dim_bra   
                left(i,j) = QM%L(DP%bra_PTR(i),j) * R_vector(basis(j)%atom,xyz) / two
            end do
            !$OMP end task
        end do
        !$OMP end single

        !$OMP single
        do j = 1 , dim_ket   
            !$OMP task shared(origin_Dependent,left,right)
            do i = 1 , dim_bra 
                origin_Dependent(i,j)%DP(xyz) = sum( left(i,:) * right(:,j) )
            end do
            !$OMP end task
        end do
        !$OMP end single

    end do
!$OMP end parallel

forall(i=1:dim_bra) a(i,:) = QM%L(DP%bra_PTR(i),:)

!$OMP parallel 
    !$OMP single
    do i = 1 , dim_ket 
        !$OMP task shared(right)
        do j = 1 , dim_basis   
            right(j,i) = QM%L(DP%ket_PTR(i),j)
        end do
        !$OMP end task
    end do
    !$OMP end single
!$OMP end parallel

do xyz = 1 , 3  
    matrix = DP_matrix_AO(:,:,xyz)
    CALL gemm(a,matrix,left,'N','N',one,zero)    
    !$OMP parallel 
        !$OMP single
        do j = 1 , dim_ket 
            !$OMP task shared(origin_Independent,left,right)
            do i = 1 , dim_bra  
                origin_Independent(i,j)%DP(xyz) = sum( left(i,:) * right(:,j) ) / two
            end do
            !$OMP end task
        end do
        !$OMP end single
    !$OMP end parallel
end do

deallocate( a , left , right )
!...........................................................................

allocate( a     ( dim_ket , dim_basis)  )
allocate( left  ( dim_ket , dim_basis)  )
allocate( right ( dim_basis  , dim_bra) )

forall(i=1:dim_bra) right(:,i) = QM%R(:,DP%bra_PTR(i))

!$OMP parallel
    do xyz = 1 , 3

        !$OMP single
        do j = 1 , dim_basis  
            !$OMP task shared(left,QM,R_vector)
            do i = 1 , dim_ket   
                left(i,j) = QM%L(DP%ket_PTR(i),j) * R_vector(basis(j)%atom,xyz) / two
            end do
            !$OMP end task
        end do
        !$OMP end single

        !$OMP single
        do j = 1 , dim_ket   
            !$OMP task shared(origin_Dependent,left,right)
            do i = 1 , dim_bra   
                origin_Dependent(i,j)%DP(xyz) = origin_Dependent(i,j)%DP(xyz) + sum( left(j,:) * right(:,i) )
            end do
            !$OMP end task
        end do
        !$OMP end single

    end do
!$OMP end parallel    

forall( i=1:dim_ket ) a(i,:) = QM%L(DP%ket_PTR(i),:)

!$OMP parallel    
    !$OMP single
    do i = 1 , dim_bra 
        !$OMP task shared(right)
        do j = 1 , dim_basis   
            right(j,i) =  QM%L(DP%bra_PTR(i),j)
        end do
        !$OMP end task
    end do
    !$OMP end single
!$OMP end parallel    

do xyz = 1 , 3  
    matrix = DP_matrix_AO(:,:,xyz)
    CALL gemm(a,matrix,left,'N','N',one,zero)    
    !$OMP parallel    
        !$OMP single
        do j = 1 , dim_ket  
            !$OMP task
            do i = 1 , dim_bra  
                origin_Independent(i,j)%DP(xyz) = origin_Independent(i,j)%DP(xyz) + sum( left(j,:) * right(:,i) ) / two
            end do
            !$OMP end task
        end do
        !$OMP end single
    !$OMP end parallel    
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

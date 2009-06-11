 module FMO_m

    use type_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Allocation_m
    use QCModel_Huckel
    use projectors
    use EHT_parameters
    use Overlap_Builder
    use Structure_Builder

    integer , public , protected :: spin(0:n_part) , orbital(0:n_part)

    public :: FMO_analysis

    private

 contains
!
!
!
!-------------------------------------------------------------------
 subroutine FMO_analysis( system, basis, zR, FMO_L, FMO_R, erg_FMO )
!-------------------------------------------------------------------
 type(structure)               , intent(in)  :: system
 type(STO_basis)               , intent(in)  :: basis(:)
 complex*16      , allocatable , intent(in)  :: zR(:,:)
 complex*16      , allocatable , intent(out) :: FMO_L(:,:) , FMO_R(:,:)  
 real*8          , allocatable , intent(out) :: erg_FMO(:)

! . local variables
 type(structure)               :: FMO_system
 type(STO_basis) , allocatable :: FMO_basis(:)
 real*8          , allocatable :: wv_FMO(:,:) 
 real*8                        :: entropy            
 character(len=1)              :: fragment

 fragment = 'D'

! orbitals to be propagated
 orbital(0) = HOMO_state    ; spin(0) = +1 
 orbital(1) = initial_state ; spin(1) = +1 

! FMO_system = fragment

 FMO_system%atoms = count(system%fragment == fragment)

 CALL Allocate_Structures(FMO_system%atoms,FMO_system)

 forall(i=1:3)
 FMO_system%coord(:,i) =  pack(system%coord(:,i) , system%fragment == fragment ) 
 end forall
 FMO_system%AtNo       =  pack( system%AtNo      , system%fragment == fragment ) 
 FMO_system%k_WH       =  pack( system%k_WH      , system%fragment == fragment )
 FMO_system%symbol     =  pack( system%symbol    , system%fragment == fragment )
 FMO_system%fragment   =  pack( system%fragment  , system%fragment == fragment )
 FMO_system%copy_No    =  0

 CALL Basis_Builder( FMO_system, FMO_basis )
 
 CALL eigen_FMO( FMO_system, FMO_basis, wv_FMO, erg_FMO )

 CALL projector( FMO_L, FMO_R, zR, basis%fragment, fragment, wv_FMO )

! "entropy" of the FMO states with respect to the system 
 OPEN(unit=9,file='entropy.dat',status='unknown')
 do i = 1 , size(FMO_basis)
     entropy = - sum( cdabs(FMO_L(:,i))*dlog(cdabs(FMO_L(:,i))) ) 
     write(9,*) erg_FMO(i) , entropy
 end do
 CLOSE(9)   

 DeAllocate( FMO_basis , wv_FMO)

 do i = 0 , n_part
    Print 59, orbital(i) , erg_FMO(orbital(i))
 end do   

 print*, ''
 print*, '>> FMO analysis done <<'

 include 'formats.h'

 end subroutine FMO_analysis
! 
!
!
!
 subroutine  eigen_FMO( system, basis, wv_FMO, erg_FMO )

 type(structure)               , intent(in)  :: system
 type(STO_basis)               , intent(in)  :: basis(:)
 real*8          , ALLOCATABLE , intent(out) :: wv_FMO(:,:)
 real*8          , ALLOCATABLE , intent(out) :: erg_FMO(:)

 integer               :: N_of_molecule_electrons, i, j
 real*8  , ALLOCATABLE :: s_FMO(:,:) , h_FMO(:,:) 

 ALLOCATE( s_FMO(size(basis),size(basis)), h_FMO(size(basis),size(basis)),  erg_FMO(size(basis)) )

!-----------------------------------------------------------------------

 CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

 DO j = 1 , size(basis)
   DO i = 1 , j 

      h_FMO(i,j) = huckel( i, j, S_FMO(i,j), basis )     !! <== define h_FMO
 
   END DO
 END DO

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

 CALL SYGVD(h_FMO,s_FMO,erg_FMO,1,'V','U',info)

 If (info /= 0) write(*,*) 'info = ',info,' in SYGVD/eigen_FMO '

!---------------------------------------------------------------------

 ALLOCATE( wv_FMO(size(basis),size(basis)) )

 wv_FMO = transpose(h_FMO)

 DeAllocate( s_FMO , h_FMO )

! save energies of the FMO system 
 OPEN(unit=9,file='molecule-ergs.dat',status='unknown')
    N_of_molecule_electrons = sum(atom(system%AtNo)%Nvalen)
    write(9,*) float(N_of_molecule_electrons) / 2.0
    do i = 1 , size(basis)
        write(9,*) i , erg_FMO(i)
    end do
 CLOSE(9)   

 print*, '>> eigen_FMO done <<'

 end subroutine
!
!
!
 end module FMO_m

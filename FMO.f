 module FMO_m

    use type_m
    use parameters_m                , only : n_part ,                   &
                                             Survival ,                 &
                                             initial_state ,            &
                                             hole_state
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Allocation_m                , only : Allocate_Structures ,      &
                                             Deallocate_Structures
    use QCModel_Huckel              , only : Huckel
    use Semi_Empirical_Parms        , only : atom
    use Overlap_Builder             , only : Overlap_Matrix
    use Structure_Builder           , only : Basis_Builder

    integer      , allocatable , public , protected :: orbital(:)
    character(2) , allocatable , public , protected :: eh_tag(:) 

    public :: FMO_analysis

    private

    ! module variables ...
    logical , save  :: done = .false.

 contains
!
!
!
!=================================================================
 subroutine FMO_analysis( system, basis, CR, FMO , MO , fragment )
!=================================================================
 implicit none
 type(structure)                                , intent(in)    :: system
 type(STO_basis)                                , intent(in)    :: basis(:)
 complex*16         , optional  , allocatable   , intent(in)    :: CR(:,:)
 type(C_eigen)                                  , intent(out)   :: FMO
 real*8             , optional  , allocatable   , intent(inout) :: MO(:)
 character(*)       , optional                  , intent(in)    :: fragment

! local variables ...
 type(structure)               :: FMO_system
 type(STO_basis) , allocatable :: FMO_basis(:)
 real*8          , allocatable :: wv_FMO(:,:) 
 real*8                        :: entropy            
 integer                       :: i

! orbitals to be propagated ...
If( .NOT. done ) then

    allocate( orbital(n_part) , eh_tag(n_part) )

    orbital(1) = initial_state ; eh_tag(1) = "el"
    orbital(2) = hole_state    ; eh_tag(2) = "hl"  

    done = .true.

end If

! FMO_system = fragment

 FMO_system%atoms = count(system%fragment == fragment)

 CALL Allocate_Structures(FMO_system%atoms,FMO_system)

 forall(i=1:3)
 FMO_system%coord(:,i) =  pack(system%coord(:,i) , system%fragment == fragment ) 
 end forall
 FMO_system%AtNo       =  pack( system%AtNo      , system%fragment == fragment ) 
 FMO_system%Nvalen     =  pack( system%Nvalen    , system%fragment == fragment ) 
 FMO_system%k_WH       =  pack( system%k_WH      , system%fragment == fragment )
 FMO_system%symbol     =  pack( system%symbol    , system%fragment == fragment )
 FMO_system%fragment   =  pack( system%fragment  , system%fragment == fragment )
 FMO_system%MMSymbol   =  pack( system%MMSymbol  , system%fragment == fragment )
 FMO_system%residue    =  pack( system%residue   , system%fragment == fragment )
 FMO_system%nr         =  pack( system%nr        , system%fragment == fragment )
 FMO_system%copy_No    =  0

 CALL Basis_Builder( FMO_system , FMO_basis )
 
 CALL eigen_FMO( FMO_system , FMO_basis , wv_FMO , FMO , fragment )

! get wv_FMO orbital in local representation and leave subroutine ... 
 if( present(MO) ) then

    ! MO vector used at Chebyshev propagator ...
    allocate( MO(size(FMO_basis)) )
    MO(:) = wv_FMO(orbital(1),:)

    print*, ''
    print*, '>> FMO analysis done <<'

    return

 end if

! wv_FMO needed only for time-propagation of wv_FMO state ...
 If( Survival ) then
     
     CALL projector( FMO , CR , basis%fragment , fragment , wv_FMO )

    ! "entropy" of the FMO states with respect to the system 
    OPEN(unit=9,file='entropy.dat',status='unknown')
    do i = 1 , size(FMO_basis)
        entropy = - sum( cdabs(FMO%L(:,i))*dlog(cdabs(FMO%L(:,i))) ) 
        write(9,*) FMO%erg(i) , entropy
    end do
    CLOSE(9)   

    DeAllocate( wv_FMO )

 end IF

 DeAllocate( FMO_basis )
 CALL DeAllocate_Structures( FMO_system )

 print*, ''
 print*, '>> FMO analysis done <<'

 include 'formats.h'

 end subroutine FMO_analysis
!
!
!
!----------------------------------------------------------------
 subroutine projector( FMO, CR, basis_fragment, fragment, wv_FMO)
!----------------------------------------------------------------
 implicit none
 type(C_eigen)                           , intent(inout) :: FMO
 complex*16       , ALLOCATABLE , target , intent(in)    :: CR(:,:)
 character(len=1)                        , intent(in)    :: basis_fragment(:)
 character(len=1)                        , intent(in)    :: fragment
 real*8           , ALLOCATABLE          , intent(in)    :: wv_FMO(:,:)

! local variables ...
 complex*16 , pointer  :: CR_FMO(:,:) => null()
 integer               :: ALL_size , FMO_size , i , j , p1 , p2
 real*8                :: check

 ALL_size = size( CR(:,1) )                     ! <== basis size of the entire system
 FMO_size = size( wv_FMO(1,:) )                 ! <== basis size of the FMO system

 Allocate( FMO%L (ALL_size,FMO_size) )
 Allocate( FMO%R (ALL_size,FMO_size) )
 Allocate( CR_FMO(FMO_size,ALL_size) )

 p1 =  minloc( [(i,i=1,ALL_size)] , 1,basis_fragment == fragment )
 p2 =  maxloc( [(i,i=1,ALL_size)] , 1,basis_fragment == fragment )

! . the fragment basis MUST correspond to a contiguous window ... 
 CR_FMO => CR(p1:p2,:)

!--------------------------------------------------------------------------------------
!             writes the isolated FMO eigenfunctions in the MO basis 
! the isolated orbitals are stored in the "ROWS of wv_FMO" and in the "COLUMNS of FMO"

 FMO%L = ( 0.d0 , 0.d0 )
 FMO%R = ( 0.d0 , 0.d0 )

 forall( j=1:FMO_size, i=1:ALL_size )

    FMO%L(i,j) = sum( wv_FMO(j,:) * CR_FMO(:,i) )

 end forall    

 FMO%R = FMO%L

 check = dreal( sum( FMO%L(1:ALL_size,:)*FMO%R(1:ALL_size,:) ) )

 if( dabs(check-FMO_size) < low_prec ) then
     print*, '>> projection done <<'
 else
     Print 58 , check 
     print*, '---> problem in projector <---'
 end if

!-----------------------------------------------------------------------------------------

 nullify( CR_FMO )

 include 'formats.h'

 end subroutine projector
! 
!
!
!--------------------------------------------------------------
 subroutine  eigen_FMO( system, basis, wv_FMO, FMO , fragment )
!--------------------------------------------------------------
 implicit none
 type(structure)               , intent(in)  :: system
 type(STO_basis)               , intent(in)  :: basis(:)
 real*8          , ALLOCATABLE , intent(out) :: wv_FMO(:,:)
 type(C_eigen)                 , intent(out) :: FMO       
 character(*)    , optional    , intent(in)  :: fragment

! local variables ... 
 integer               :: N_of_FMO_electrons, i, j , info
 real*8  , ALLOCATABLE :: s_FMO(:,:) , h_FMO(:,:)

 ALLOCATE( s_FMO(size(basis),size(basis)), h_FMO(size(basis),size(basis)),  FMO%erg(size(basis)) )

!-----------------------------------------------------------------------

 CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

 DO j = 1 , size(basis)
   DO i = 1 , j 

      h_FMO(i,j) = huckel( i, j, S_FMO(i,j), basis )     !! <== define h_FMO
 
   END DO
 END DO

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

 CALL SYGVD(h_FMO,s_FMO,FMO%erg,1,'V','U',info)

 If (info /= 0) write(*,*) 'info = ',info,' in SYGVD/eigen_FMO '

 FMO % Fermi_State = sum( system%Nvalen ) / two
!---------------------------------------------------------------------

 ALLOCATE( wv_FMO(size(basis),size(basis)) )

 wv_FMO = transpose(h_FMO)

 DeAllocate( s_FMO , h_FMO )

! save energies of the FMO system 
 If( present(fragment) .AND. fragment=="H" ) then
    OPEN(unit=9,file='hl_FMO-ergs.dat',status='unknown')
 else
    OPEN(unit=9,file='el_FMO-ergs.dat',status='unknown')
 end IF

 N_of_FMO_electrons = sum( system%Nvalen )
 write(9,*) float(N_of_FMO_electrons) / 2.0
 do i = 1 , size(basis)
    write(9,*) i , FMO%erg(i)
 end do
 CLOSE(9)   

 print*, '>> eigen_FMO done <<'

 end subroutine
!
!
!
 end module FMO_m

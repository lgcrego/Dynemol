 module FMO_m

    use type_m
    use parameters_m                , only : n_part ,                   &
                                             Survival ,                 &
                                             initial_state ,            &
                                             hole_state ,               &
                                             LCMO
    use f95_precision
    use blas95
    use lapack95
    use Allocation_m                , only : Allocate_Structures ,      &
                                             Deallocate_Structures
    use Semi_Empirical_Parms        , only : atom
    use Overlap_Builder             , only : Overlap_Matrix
    use Structure_Builder           , only : Basis_Builder
    use LCMO_m                      , only : LCMO_Builder

    integer      , allocatable , public :: orbital(:)
    character(2) , allocatable , public :: eh_tag(:) 

    public :: FMO_analysis

    private

    ! module variables ...
    logical , save  :: done = .false.

 contains
!
!
!
!=================================================================
 subroutine FMO_analysis( system, basis, CR, FMO , MO , instance )
!=================================================================
 implicit none
 type(structure)                                , intent(inout) :: system
 type(STO_basis)                                , intent(inout) :: basis(:)
 real*8             , optional  , allocatable   , intent(in)    :: CR(:,:)
 type(R_eigen)                                  , intent(out)   :: FMO
 real*8             , optional  , allocatable   , intent(inout) :: MO(:)
 character(*)       , optional                  , intent(in)    :: instance

! local variables ...
 type(structure)               :: FMO_system
 type(STO_basis) , allocatable :: FMO_basis(:)
 real*8          , allocatable :: wv_FMO(:,:) 
 integer                       :: i
 character(1)                  :: fragment
 character(1)    , allocatable :: system_fragment(:) , basis_fragment(:)

 CALL preprocess( system, basis, system_fragment , basis_fragment , fragment , instance )

!FMO_system = fragment ...
 FMO_system%atoms = count(system%fragment == fragment)

 CALL Allocate_Structures(FMO_system%atoms,FMO_system)

! Notice: not everything needs to be cloned into FMO ... 
 forall(i=1:3)
 FMO_system%coord(:,i) =  pack(system%coord(:,i) , system%fragment == fragment ) 
 end forall
 FMO_system%AtNo       =  pack( system%AtNo      , system%fragment == fragment ) 
 FMO_system%Nvalen     =  pack( system%Nvalen    , system%fragment == fragment ) 
 FMO_system%k_WH       =  pack( system%k_WH      , system%fragment == fragment )
 FMO_system%symbol     =  pack( system%symbol    , system%fragment == fragment )
 FMO_system%fragment   =  pack( system%fragment  , system%fragment == fragment )
 FMO_system%MMSymbol   =  pack( system%MMSymbol  , system%fragment == fragment )
 FMO_system%QMMM       =  pack( system%QMMM      , system%fragment == fragment )
 FMO_system%residue    =  pack( system%residue   , system%fragment == fragment )
 FMO_system%nr         =  pack( system%nr        , system%fragment == fragment )
 FMO_system%V_shift    =  pack( system%V_shift   , system%fragment == fragment )
 FMO_system%copy_No    =  0

! check point ...
 If( any(FMO_system%QMMM /= "QM") ) stop ">> FMO fragment contains MM atoms <<"

 CALL Basis_Builder( FMO_system , FMO_basis )

 CALL eigen_FMO( FMO_system , FMO_basis , wv_FMO , FMO , fragment )

 If ( LCMO ) CALL LCMO_Builder( wv_FMO , FMO%erg , instance )
! the following subroutine can be used to check the LCMO packets ... 
! call check_casida_builder( FMO_system , FMO_basis , wv_FMO , FMO )

 If( present(MO) ) then

    ! get wv_FMO orbital in local representation and leave subroutine ... 
    ! MO vector used at Chebyshev propagator ...
    allocate( MO(size(FMO_basis)) )

    select case(instance)
        case ("E","D")
        MO(:) = wv_FMO(orbital(1),:)
        case ("H")
        MO(:) = wv_FMO(orbital(2),:)
    end select

 else If( Survival ) then

    ! wv_FMO needed only for time-propagation of wv_FMO state ...
    CALL projector( FMO , CR , basis%fragment , fragment , wv_FMO )

 end IF

 DeAllocate( FMO_basis , wv_FMO )
 CALL DeAllocate_Structures( FMO_system )

 system % fragment = system_fragment
 basis  % fragment = basis_fragment

 deallocate( system_fragment , basis_fragment )

 print*, '>> FMO analysis done <<'

 include 'formats.h'

 end subroutine FMO_analysis
!
!
!
!==============================================================================================
 subroutine preprocess( system, basis, system_fragment , basis_fragment , fragment , instance )
!==============================================================================================
implicit none
 type(structure)                            , intent(inout) :: system
 type(STO_basis)                            , intent(inout) :: basis(:)
 character(1)    , allocatable              , intent(out)   :: system_fragment(:) 
 character(1)    , allocatable              , intent(out)   :: basis_fragment(:)
 character(1)                               , intent(out)   :: fragment
 character(*)                   , optional  , intent(in)    :: instance

! orbitals to be propagated ...
 If( (.NOT. done) .AND. Survival ) then

    If( any(system%Hl)          .AND. n_part == 1 )  stop ">> n_part = 1 , but El/Hl is .true.   <<"
    If( (.NOT. any(system%Hl))  .AND. n_part == 2 )  stop ">> n_part = 2 , but El/Hl is .false.  <<"

    allocate( orbital(n_part) , eh_tag(n_part) )

    orbital(1) = initial_state ; eh_tag(1) = "el"
    If( any(system%Hl) ) then
    orbital(2) = hole_state    ; eh_tag(2) = "hl"  
    End If

    done = .true.

 end If

! setting the fragment ...
 allocate(system_fragment (system%atoms) , source = system % fragment )
 allocate( basis_fragment (size(basis) ) , source =  basis % fragment )

 If( .not. present(instance) ) then

        ! entire system ...
        fragment = "#"       
        system % fragment  = fragment
        basis  % fragment  = fragment

    else

    fragment = instance

    select case (instance)

        case( "D" )
            where( system%El ) system % fragment  = instance
            where( basis%El  ) basis  % fragment  = instance

        case( "E" )
            where( system%El ) system % fragment  = instance
            where( basis%El  ) basis  % fragment  = instance

        case( "H" )
            where( system%Hl ) system % fragment  = instance
            where( basis%Hl  ) basis  % fragment  = instance
    end select 

 end if

 end subroutine preprocess
!
!
!
!----------------------------------------------------------------
 subroutine projector( FMO, CR, basis_fragment, fragment, wv_FMO)
!----------------------------------------------------------------
 implicit none
 type(R_eigen)                           , intent(inout) :: FMO
 real*8           , ALLOCATABLE , target , intent(in)    :: CR(:,:)
 character(len=1)                        , intent(in)    :: basis_fragment(:)
 character(len=1)                        , intent(in)    :: fragment
 real*8           , ALLOCATABLE          , intent(in)    :: wv_FMO(:,:)

! local variables ...
 real*8     , pointer  :: CR_FMO(:,:) => null()
 integer               :: ALL_size , FMO_size , i , j , p1 , p2
 real*8                :: check

 ALL_size = size( CR(:,1) )                     ! <== basis size of the entire system
 FMO_size = size( wv_FMO(1,:) )                 ! <== basis size of the FMO system

 Allocate( FMO%L (ALL_size,FMO_size) , source=D_zero )
 Allocate( FMO%R (ALL_size,FMO_size) , source=D_zero )
 Allocate( CR_FMO(FMO_size,ALL_size) )

 p1 =  minloc( [(i,i=1,ALL_size)] , 1,basis_fragment == fragment )
 p2 =  maxloc( [(i,i=1,ALL_size)] , 1,basis_fragment == fragment )

!the fragment basis MUST correspond to a contiguous window ... 
 CR_FMO => CR(p1:p2,:)

!--------------------------------------------------------------------------------------
!             writes the isolated FMO eigenfunctions in the MO basis 
! the isolated orbitals are stored in the "ROWS of wv_FMO" and in the "COLUMNS of FMO"

 forall( j=1:FMO_size, i=1:ALL_size )

    FMO%L(i,j) = sum( wv_FMO(j,:) * CR_FMO(:,i) )

 end forall    

 FMO%R = FMO%L

 check = sum( FMO%L(1:ALL_size,:)*FMO%R(1:ALL_size,:) ) 

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
 type(R_eigen)                 , intent(out) :: FMO       
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

 FMO % Fermi_State = sum(system%Nvalen)/two + mod( sum(system%Nvalen) , 2 )
!---------------------------------------------------------------------

 ALLOCATE( wv_FMO(size(basis),size(basis)) )

 wv_FMO = transpose(h_FMO)

 DeAllocate( s_FMO , h_FMO )

! save energies of the FMO system 
 If( present(fragment) .AND. (fragment=="H") ) then
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

 end subroutine eigen_FMO
!
!
!
!====================================
 pure function Huckel(i,j,S_ij,basis)
!====================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel
 real*8  :: k_eff , k_WH , c1 , c2 , c3

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

 k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

 huckel = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / two

 IF(i == j) huckel = basis(i)%IP

 end function Huckel
!
!
!
!--------------------------------------------------------------
 subroutine  check_casida_builder( system, basis, wv_FMO, FMO )
!--------------------------------------------------------------
 implicit none
 type(structure) , intent(in)  :: system
 type(STO_basis) , intent(in)  :: basis(:)
 real*8          , intent(in)  :: wv_FMO(:,:)
 type(R_eigen)   , intent(in)  :: FMO       

! local variables ... 
 integer               :: i, j , nn
 real*8                :: erg , pop
 real*8  , ALLOCATABLE :: s_FMO(:,:) , h_FMO(:,:) , tmp_S(:) , tmp_E(:)

nn = size(basis)

ALLOCATE( s_FMO(nn,nn) , h_FMO(nn,nn) , tmp_S(nn) , tmp_E(nn) ) 

!-----------------------------------------------------------------------

CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

DO j = 1 , nn
   DO i = 1 , j 

      h_FMO(i,j) = huckel( i, j, S_FMO(i,j), basis )     !! <== define h_FMO
      h_FMO(j,i) = h_FMO(i,j)
 
   END DO
END DO

pop = D_zero
do i = 1 , nn
   
    tmp_S = matmul( S_FMO , wv_FMO(i,:) )
    tmp_E = matmul( h_FMO , wv_FMO(i,:) )

    pop = dot_product( wv_FMO(i,:) , tmp_S ) + pop
    erg = dot_product( wv_FMO(i,:) , tmp_E )

    print*, i, erg , FMO % erg(i)
end do

print*, pop

deallocate( s_FMO , h_FMO , tmp_S , tmp_E )

end subroutine  check_casida_builder
! 
!
! 
end module FMO_m

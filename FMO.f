 module FMO_m

    use type_m
    use f95_precision
    use blas95
    use lapack95
    use parameters_m                , only : driver ,                   &
                                             n_part ,                   &
                                             Survival ,                 &
                                             EnvField_ ,                &
                                             Environ_type ,             &
                                             Induced_ ,                 & 
                                             electron_state ,           &
                                             hole_state ,               &
                                             LCMO ,                     &
                                             verbose
    use Allocation_m                , only : Allocate_Structures ,      &
                                             Deallocate_Structures
    use Dielectric_Potential        , only : Q_phi
    use DP_potential_m              , only : DP_phi
    use DP_main_m                   , only : DP_matrix_AO
    use Semi_Empirical_Parms        , only : atom
    use tuning_m                    , only : eh_tag , orbital 
    use Overlap_Builder             , only : Overlap_Matrix
    use Structure_Builder           , only : Basis_Builder
    use Hamiltonians                , only : X_ij , Huckel_with_Fields 
    use LCMO_m                      , only : LCMO_Builder

    public :: FMO_analysis , eh_tag , orbital

    private

    type(R_eigen) :: Dual , tmp


 contains
!
!
!
!=================================================================
 subroutine FMO_analysis( system, basis, UNI, FMO , MO , instance )
!=================================================================
 implicit none
 type(structure)                          , intent(inout) :: system
 type(STO_basis)                          , intent(inout) :: basis(:)
 type(R_eigen)  , optional                , intent(in)    :: UNI
 type(R_eigen)                            , intent(out)   :: FMO
 real*8         , optional  , allocatable , intent(inout) :: MO(:)
 character(*)   , optional                , intent(in)    :: instance

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
    CALL projector( FMO , UNI , basis%fragment , fragment )

 end IF

 DeAllocate( FMO_basis , wv_FMO )
 CALL DeAllocate_Structures( FMO_system )

 system % fragment = system_fragment
 basis  % fragment = basis_fragment

 deallocate( system_fragment , basis_fragment )

 Print*, '>> FMO analysis done <<'

 include 'formats.h'

 end subroutine FMO_analysis
!
!
!
!==============================================================================================
 subroutine preprocess( system, basis, system_fragment , basis_fragment , fragment , instance )
!==============================================================================================
implicit none
 type(structure)                           , intent(inout) :: system
 type(STO_basis)                           , intent(inout) :: basis(:)
 character(1)    , allocatable             , intent(out)   :: system_fragment(:) 
 character(1)    , allocatable             , intent(out)   :: basis_fragment(:)
 character(1)                              , intent(out)   :: fragment
 character(*)                  , optional  , intent(in)    :: instance

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
!----------------------------------------------------------
 subroutine projector( FMO, UNI, basis_fragment, fragment )
!----------------------------------------------------------
 implicit none
 type(R_eigen)                  , intent(inout) :: FMO
 type(R_eigen)                  , intent(in)    :: UNI
 character(len=1)               , intent(in)    :: basis_fragment(:)
 character(len=1)               , intent(in)    :: fragment

! local variables ...
 integer :: UNI_size , FMO_size , i , j , p1 , p2 , k
 real*8  :: check

 UNI_size = size( UNI%R  (:,1) )   ! <== basis size of the entire system
 FMO_size = size( Dual%R (:,1) )   ! <== basis size of the FMO system

!-----------------------------------------------------------------------------
! cast the FMO eigenvectors in UNI eigen-space

 allocate( tmp%L(UNI_size,FMO_size), source=D_zero )
 allocate( tmp%R(UNI_size,FMO_size), source=D_zero )
 k = 0 
 do i = 1 , UNI_size
    if( basis_fragment(i) == fragment ) then

        k = k + 1
        tmp%R(i,:) = Dual%R(k,:)
        tmp%L(i,:) = Dual%L(:,k) 

    end if
 end do
 ! tmp is deallocated ...
 CALL move_alloc( tmp%R , Dual%R )
 CALL move_alloc( Tmp%L , Dual%L )

!---------------------------------------------------------------------------
!        writes the isolated FMO eigenfunctions in the MO basis 
! orbitals are stored in the "ROWS of FMO%L" and in the "COLUMNS of FMO%R"

 Allocate( FMO%L (FMO_size,UNI_size) , source=D_zero )
 Allocate( FMO%R (UNI_size,FMO_size) , source=D_zero )

 CALL gemm( UNI%L , Dual%R , FMO%R , 'N' , 'N' )
 CALL gemm( Dual%L , UNI%R , FMO%L , 'T' , 'N' )

 check = 0.d0
 do i = 1 , FMO_size
    ! %L*%R = A^T.S.C.C^T.S.A = 1
    check = check + sum( FMO%L(i,:)*FMO%R(:,i) ) 
 end do

 if( dabs(check-FMO_size) < low_prec ) then
     Print*, '>> projection done <<'
 else
     Print 58 , check 
     Print*, '---> problem in projector <---'
 end if

!--------------------------------------------------------------------------------------------------

 deallocate( Dual%L , Dual%R )

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
 integer               :: N_of_FMO_electrons, i, N , info
 real*8  , ALLOCATABLE :: s_FMO(:,:) , h_FMO(:,:) , dumb_S(:,:)

 N = size(basis)

 ALLOCATE( s_FMO(N,N)  , h_FMO(N,N) ,  FMO%erg(N) )

 CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

 If( EnvField_ .OR. Induced_ ) then
     h_FMO = even_more_extended_Huckel( system , basis , S_FMO )
 else
     h_FMO = Build_Huckel( basis , S_FMO )
 end If

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

 ALLOCATE( dumb_S(N,N) , source = s_FMO )

 CALL SYGVD(h_FMO , dumb_S , FMO%erg , 1 , 'V' , 'U' , info)

 If (info /= 0) write(*,*) 'info = ',info,' in SYGVD/eigen_FMO '

 FMO % Fermi_State = sum(system%Nvalen)/two + mod( sum(system%Nvalen) , 2 )

 ALLOCATE( wv_FMO(N,N) )

 wv_FMO = transpose(h_FMO)

 If( driver == "slice_AO" ) then 
     ALLOCATE(Dual%L(N,N) , source = wv_FMO) 
     ALLOCATE(Dual%R(N,N)) 
     CALL symm( s_FMO , h_FMO , Dual%R )
 end IF

 DeAllocate( s_FMO , h_FMO , dumb_S )

! save energies of the FMO system 
 If( present(fragment) .AND. (fragment=="H") ) then
    OPEN(unit=9,file='hl_FMO-ergs.dat',status='unknown')
 else
    OPEN(unit=9,file='el_FMO-ergs.dat',status='unknown')
 end IF

 N_of_FMO_electrons = sum( system%Nvalen )
 write(9,*) float(N_of_FMO_electrons) / 2.0
 do i = 1 , N
    write(9,*) i , FMO%erg(i)
 end do
 CLOSE(9)   

 Print*, '>> eigen_FMO done <<'

 end subroutine eigen_FMO
!
!
!
!===================================================
 function Build_Huckel( basis , S_matrix ) result(h)
!===================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

! local variables ... 
integer :: i , j , N
real*8  , allocatable   :: h(:,:)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

N = size(basis)
ALLOCATE( h(N,N) , source = D_zero )

do j = 1 , N
  do i = 1 , j 

        h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j)

        h(j,i) = h(i,j)

    end do
end do

end function Build_Huckel
!
!
!
!=========================================================================
 function even_more_extended_huckel( system , basis , S_matrix ) result(h)
!=========================================================================
implicit none
type(structure) , intent(in) :: system
type(STO_basis) , intent(in) :: basis(:)
real*8          , intent(in) :: S_matrix(:,:)

! local variables ...
integer               :: i , j , ia , ib , ja , jb , N
real*8                :: Rab , DP_4_vector(4)
real*8  , ALLOCATABLE :: h(:,:) 

N = size(basis)
Allocate( h(N,N) , source = D_zero )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!$OMP parallel do &
!$OMP   default(shared) &
!$OMP   schedule(dynamic, 1) &
!$OMP   private(ib, ia, Rab, jb, ja, j, i, DP_4_vector)
do ib = 1, system%atoms
    do ia = ib+1, system%atoms

        Rab = GET_RAB(system%coord(ib,:), system%coord(ia,:))
        if (Rab > cutoff_Angs) then
           cycle
        end if

        select case (Environ_Type)
           case('DP_MM','DP_QM')
               DP_4_vector = DP_phi( system , ia , ib )
           case default
               DP_4_vector =  Q_phi( system , ia , ib )
        end select

        do jb = 1, atom(system%AtNo(ib))% DOS
            do ja = 1, atom(system%AtNo(ia))% DOS

               j = system% BasisPointer(ib) + jb
               i = system% BasisPointer(ia) + ja

               h(i,j) = huckel_with_FIELDS(i , j , S_matrix(i,j) , basis , DP_4_vector )

               h(j,i) = h(i,j)

            end do
        end do

    end do
end do  
!$OMP END PARALLEL DO

forall( i=1:N ) h(i,i) = X_ij( i , i , basis ) 

end function even_more_extended_huckel
!
!
!
!==================================================
pure function GET_RAB(a_coord, b_coord) result(rab)
!==================================================
 implicit none

 ! args
 real*8, intent(in) :: a_coord(:)
 real*8, intent(in) :: b_coord(:)

 ! result
 real*8 :: rab

 rab = SUM((a_coord - b_coord) ** 2)
 rab = SQRT(rab)
end function GET_RAB
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
 integer               :: i, nn
 real*8                :: erg , pop
 real*8  , ALLOCATABLE :: s_FMO(:,:) , h_FMO(:,:) , tmp_S(:) , tmp_E(:)

nn = size(basis)

ALLOCATE( s_FMO(nn,nn) , h_FMO(nn,nn) , tmp_S(nn) , tmp_E(nn) ) 

!-----------------------------------------------------------------------

CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

If( EnvField_ .OR. Induced_ ) then
    h_FMO = even_more_extended_Huckel( system , basis , S_FMO )
else
    h_FMO = Build_Huckel( basis , S_FMO )
end If

pop = D_zero
do i = 1 , nn
   
    tmp_S = matmul( S_FMO , wv_FMO(i,:) )
    tmp_E = matmul( h_FMO , wv_FMO(i,:) )

    pop = dot_product( wv_FMO(i,:) , tmp_S ) + pop
    erg = dot_product( wv_FMO(i,:) , tmp_E )

    Print*, i, erg , FMO % erg(i)
end do

Print*, pop

deallocate( s_FMO , h_FMO , tmp_S , tmp_E )

end subroutine  check_casida_builder
!
!
!
end module FMO_m

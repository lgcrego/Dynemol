 module FMO_m

    use IFPORT
    use type_m
    use f95_precision
    use blas95
    use lapack95
    use parameters_m                , only : driver ,               &
                                             n_part ,               &
                                             Survival ,             &
                                             EnvField_ ,            &
                                             Environ_type ,         &
                                             Induced_ ,             & 
                                             electron_state ,       &
                                             hole_state ,           &
                                             LCMO ,                 &
                                             SOC ,                  &
                                             verbose
    use Allocation_m                , only : Allocate_Structures ,  &
                                             Deallocate_Structures
    use Dielectric_Potential        , only : Q_phi
    use DP_potential_m              , only : DP_phi
    use DP_main_m                   , only : DP_matrix_AO
    use Semi_Empirical_Parms        , only : atom
    use tuning_m                    , only : eh_tag , orbital 
    use Overlap_Builder             , only : Overlap_Matrix
    use Structure_Builder           , only : Basis_Builder
    use Hamiltonians                , only : X_ij , Huckel_with_Fields , spin_orbit_h
    use LCMO_m                      , only : LCMO_Builder
    use Hamiltonians                , only : spin_orbit_h

    public :: FMO_analysis , eh_tag , orbital

    private

    integer       :: UNI_size , FMO_size 
    type(C_eigen) :: Dual 
    type(C_eigen) :: tmp

 contains
!
!
!
!=====================================================================
 subroutine FMO_analysis( system , basis , UNI , FMO , AO , instance )
!=====================================================================
 implicit none
 type(structure)           , intent(inout) :: system
 type(STO_basis)           , intent(inout) :: basis(:)
 type(C_eigen)  , optional , intent(in)    :: UNI
 type(C_eigen)  , optional , intent(out)   :: FMO
 type(C_eigen)  , optional , intent(inout) :: AO
 character(*)   , optional , intent(in)    :: instance

! local variables ...
type(structure)               :: FMO_system
type(STO_basis) , allocatable :: FMO_basis(:)
integer                       :: i
character(1)                  :: fragment
character(1)    , allocatable :: system_fragment(:) , basis_fragment(:)
logical                       :: TorF

CALL preprocess( system , basis , system_fragment , basis_fragment , fragment , instance )

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
FMO_system%T_xyz      =  system%T_xyz
FMO_system%copy_No    =  0

! check point ...
If( any(FMO_system%QMMM /= "QM") ) then
    TorF = systemQQ("sed '11i >>> FMO fragment contains MM atoms <<<' warning.signal |cat")                                  
    stop     
end If

CALL Basis_Builder( FMO_system , FMO_basis )

CALL eigen_FMO( FMO_system , FMO_basis , fragment )

! If ( LCMO ) CALL LCMO_Builder( Dual%L, Dual%erg , instance )
! the following subroutine can be used to check the LCMO states ... 
! call check_casida_builder( FMO_system , FMO_basis , Dual%L, Dual%erg )

If( Survival .AND. (.not. present(AO)) ) then

    ! Psi_0 = FMO used at AO/MO propagator ...
    CALL projector( FMO , UNI , basis%fragment , fragment , instance = 'MO' )

elseIf( present(AO) ) then

    if( SOC ) then
        print*, "present(AO) not implemented for SOC=.true. ==> Problem in routine FMO_analysis (FMO.f)"
        stop
    end if
    
    ! Psi_0 in local representation ... 
    ! used at Chebyshev propagator ...
    CALL projector( basis_fragment = basis%fragment , fragment = fragment , instance = 'AO' )

    allocate( AO%L(UNI_size,2) , AO%R(UNI_size,2) )

    select case(instance)
        case ("E","D")
            AO%L(:,1) = tmp%L(orbital(1),:)
            AO%R(:,1) = tmp%R(:,orbital(1))
        case ("H")
            AO%L(:,2) = tmp%L(orbital(2),:)
            AO%R(:,2) = tmp%R(:,orbital(2))
    end select

    deallocate( tmp%L , tmp%R )

end IF

DeAllocate( FMO_basis )
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
!================================================================================================
 subroutine preprocess( system , basis , system_fragment , basis_fragment , fragment , instance )
!================================================================================================
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
    system % fragment = fragment
    basis  % fragment = fragment

else

    fragment = instance

    select case (instance)

        case( "D" )
            where( system%El ) system % fragment = instance
            where( basis%El  ) basis  % fragment = instance

        case( "E" )
            where( system%El ) system % fragment = instance
            where( basis%El  ) basis  % fragment = instance

        case( "H" )
            where( system%Hl ) system % fragment = instance
            where( basis%Hl  ) basis  % fragment = instance
                
    end select

end if

end subroutine preprocess
!
!
!
!========================================================================
 subroutine projector( FMO , UNI , basis_fragment , fragment , instance )
!========================================================================
implicit none
type(C_eigen)    , optional , intent(out) :: FMO
type(C_eigen)    , optional , intent(in)  :: UNI
character(len=1)            , intent(in)  :: basis_fragment(:)
character(len=1)            , intent(in)  :: fragment
character(len=2)            , intent(in)  :: instance

! local variables ...
complex*16 , allocatable :: aux(:,:)
real*8  :: check
integer :: i , j , k

UNI_size = size( basis_fragment )   ! <== basis size of the entire system
FMO_size = size( Dual%R (:,1)   )   ! <== basis size of the FMO system

select case ( instance ) 

    case('MO') 

        !--------------------------------------------------------------------------
        ! cast the FMO eigenvectors in UNI eigen-space
        !--------------------------------------------------------------------------

        allocate( aux(FMO_size,UNI_size), source=C_zero )

        k = 0 
        do i = 1 , UNI_size
            if( basis_fragment(i) == fragment ) then

                k = k + 1
                aux(k,:) = UNI%R(i,:)

            end if
        end do

        !-------------------------------------------------------------------------
        ! isolated FMO eigenfunctions in MO basis for use in AO/MO propagator
        ! orbitals are stored in the "ROWS of FMO%L" and in the "COLUMNS of FMO%R"
        !-------------------------------------------------------------------------
        Allocate( FMO%erg(FMO_size)          , source=D_zero )
        Allocate( FMO%L  (FMO_size,UNI_size) , source=C_zero )
        Allocate( FMO%R  (UNI_size,FMO_size) , source=C_zero )

        CALL gemm( Dual%L , aux , FMO%L )

        aux = dconjg(FMO%L)
        FMO%R = transpose(aux)
        
        deallocate( aux )

        forall(k=1:FMO_size) FMO%erg(k) = dreal(sum( FMO%L(k,:)*UNI%erg(:)*FMO%R(:,k) ))

        FMO% Fermi_state = Dual% Fermi_state

        check = D_zero
        do i = 1 , FMO_size
           ! %L*%R = A^T.S.C.C^T.S.A = 1
           check = check + dreal(sum( FMO%L(i,:) * FMO%R(:,i) ))
        end do

        if( dabs(check-FMO_size) < low_prec ) then
            Print*, '>> projection done <<'
        else
            Print 58 , check 
            Print*, '---> problem in projector <---'
        end if

    case('AO') 

        !-------------------------------------------------------------------------
        ! isolated FMO eigenfunctions in AO basis for use in Taylor propagator
        ! orbitals are stored in the "ROWS of AO%L" and in the "COLUMNS of AO%R"
        !-------------------------------------------------------------------------

        allocate( aux(UNI_size,FMO_size), source=C_zero )
        k = 0
        do i = 1 , UNI_size
           if( basis_fragment(i) == fragment ) then

               k = k + 1
               aux(i,:) = Dual%L(:,k) 

           end if
        end do
        ! aux is deallocated ...
        CALL move_alloc( aux , Dual%L )

        allocate( tmp%L(FMO_size,UNI_size) , source=transpose(Dual%L) )
        allocate( tmp%R(UNI_size,FMO_size) , source=Dual%L            )

end select

deallocate( Dual%L , Dual%R , Dual%erg )

include 'formats.h'

end subroutine projector
! 
!
!
!==================================================
 subroutine  eigen_FMO( system , basis , fragment )
!==================================================
implicit none
type(structure)               , intent(in)  :: system
type(STO_basis)               , intent(in)  :: basis(:)
character(*)    , optional    , intent(in)  :: fragment

! local variables ... 
complex*16 , ALLOCATABLE :: h_FMO(:,:) , h_spin(:,:) , S_complex(:,:) , dumb_S(:,:)
real*8     , ALLOCATABLE :: s_FMO(:,:) , h_Huckel(:,:)
real*8                   :: Fermi_level
integer                  :: N_of_FMO_electrons , i , j , N , info

N = size(basis)

ALLOCATE( s_FMO(N,N)  , h_Huckel(N,N) ,  Dual%erg(N) )

CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

If( EnvField_ .OR. Induced_ ) then
    h_Huckel = even_more_extended_Huckel( system , basis , S_FMO )
else
    h_Huckel = Build_Huckel( basis , S_FMO )
end If

if( SOC ) then

!    CALL spin_orbit_h( basis , h_spin , s_FMO )
!    allocate( h_FMO(N,N) , source = dcmplx( h_Huckel , D_zero ) + h_spin )
!    deallocate( h_spin )
    allocate( h_FMO(N,N) , source = dcmplx( h_Huckel , D_zero ) )

else

    allocate( h_FMO(N,N) , source = dcmplx( h_Huckel , D_zero ) )

end if

deallocate(h_Huckel)

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

ALLOCATE( S_complex(N,N) , source = dcmplx( s_FMO , D_zero ) )
ALLOCATE( dumb_S(N,N) , source = S_complex )

CALL HEGVD( h_FMO , dumb_S , Dual%erg , 1 , 'V' , 'L' , info )
if ( info /= 0 ) write(*,*) 'info = ',info,' in HEGVD/eigen/FMO '

ALLOCATE(Dual%R(N,N)) 
ALLOCATE(Dual%L(N,N)) 
Dual%L = transpose(conjg(h_FMO))

CALL hemm( S_complex , h_FMO , Dual%R )

DeAllocate( S_complex , s_FMO , h_FMO , dumb_S )

Dual% Fermi_state = sum(system%Nvalen)/two + mod( sum(system%Nvalen) , 2 )

! save energies of the FMO system 
If( present(fragment) .AND. (fragment=="H") ) then
    OPEN(unit=9,file='hl_FMO-ergs.dat',status='unknown')
else
    OPEN(unit=9,file='el_FMO-ergs.dat',status='unknown')
end IF

N_of_FMO_electrons = sum( system%Nvalen )
Fermi_level = merge(dfloat(N_of_FMO_electrons),dfloat(N_of_FMO_electrons)/TWO,SOC)
write(9,*) "# Fermi level = " , Fermi_level
if( SOC ) then
    write(9,100) "#" , "Level" , "Energy" , "Sz"
    do i = 1 , N
        write(9,*) i , Dual%erg(i) , dreal( sum( Dual%L(i,:) * dcmplx( basis( : ) % s , D_zero ) * Dual%R(:,i) ) )
    end do
else
    write(9,100) "#" , "Level" , "Energy"
    do i = 1 , N
        write(9,*) i , Dual%erg(i)
    end do
end if
close(9)

if( SOC ) then
    If( present(fragment) .AND. (fragment=="H") ) then
        OPEN(unit=9,file='hl_FMO-complet.dat',status='unknown')
    else
        OPEN(unit=9,file='el_FMO-complet.dat',status='unknown')
    end IF
    write(9,101) "Orb" , "Energy" , "Sz" , "s up" , "py up" , "pz up" , "px up" , "s down" , "py down" , "pz down" , "px down"
    do i = 1 , N
        write(9,102) i , Dual%erg(i) , dreal( sum( Dual%L(i,:) * Dual%R(:,i) * basis(:) % s ) ) , &
                     ( dreal( Dual%L(i,j) * Dual%R(j,i) ) , j = 1 , N )
    end do
    close(9)
end if

Print*, '>> eigen_FMO done <<'

100 format(a1,a11,a19,a24)
101 format(a3,a11,9a13)
102 format(i3,f11.6,9f13.9)

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
real*8  , allocatable :: h(:,:)
integer               :: i , j , N , N2

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
!----------------------------------------------------------

N  = size(basis)

ALLOCATE( h(N,N) , source = D_zero )

if( SOC ) then

    N2 = N/2

    do j = 1 , N2
        do i = j , N2

            ! spin up orbital block
            h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j) 
            ! spin down orbital block
            h(i+N2,j+N2) = h(i,j)

        end do
    end do

else

    do j = 1 , N
        do i = j , N

            h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j) 

        end do
    end do

end if

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
!==================================================================
 subroutine  check_casida_builder( system, basis, wv_FMO, erg_FMO )
!==================================================================
 implicit none
 type(structure) , intent(in)  :: system
 type(STO_basis) , intent(in)  :: basis(:)
 real*8          , intent(in)  :: wv_FMO(:,:)
 real*8          , intent(in)  :: erg_FMO(:)

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

    tmp_E = matmul( h_FMO , wv_FMO(i,:) )

    pop = dot_product( wv_FMO(i,:) , tmp_S ) + pop
    erg = dot_product( wv_FMO(i,:) , tmp_E )

    Print*, i, erg , erg_FMO(i)
end do

Print*, pop

deallocate( s_FMO , h_FMO , tmp_S , tmp_E )

end subroutine  check_casida_builder
!
!
!
end module FMO_m

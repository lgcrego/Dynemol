module qdynamics_m
! MPI-free module; includes MPI-free version of QCModel_Huckel.f

 use f95_precision
 use blas95
 use lapack95
 use omp_lib
 use type_m
 use constants_m
 use parameters_m           , only : spectrum , DP_Moment , &
                                     survival , DP_Field_ , &
                                     NetCharge , Induced_ , &
                                     verbose
 use Solvated_M             , only : DeAllocate_TDOS ,      &
                                     DeAllocate_PDOS ,      &
                                     DeAllocate_SPEC 
 use FMO_m                  , only : FMO_analysis ,         &
                                     eh_tag 
 use DOS_m                  , only : Total_DOS ,            &
                                     Partial_DOS
 use Structure_Builder      , only : Generate_Structure ,   &
                                     Basis_Builder ,        &
                                     Unit_Cell ,            &
                                     Extended_Cell ,        &
                                     ExCell_Basis
 use Overlap_Builder        , only : Overlap_Matrix
 use polarizability_m       , only : Induced_DP_phi
 use DP_main_m              , only : Dipole_Matrix , DP_matrix_AO
 use DP_potential_m         , only : Molecular_DPs , DP_phi
 use Oscillator_m           , only : Optical_Transitions
 use Schroedinger_m         , only : Huckel_dynamics ,      &
                                     ElHl_dynamics ,        &
                                     DeAllocate_QDyn
 use Data_Output            , only : Dump_stuff ,           &
                                     Net_Charge

 public :: qdynamics

 private

 contains
!
!
!
!====================
subroutine qdynamics
!====================
implicit none 

! local variables ...
 integer                        :: nr , N_of_residues
 logical                        :: DIPOLE_ , el_hl_
 type(R_eigen)                  :: UNI
 type(R_eigen)                  :: el_FMO , hl_FMO
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(f_time)                   :: QDyn

 
! preprocessing stuff ...................................

el_hl_  = any( Unit_Cell%Hl)
DIPOLE_ = ( spectrum .OR. DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )
CALL DeAllocate_QDyn( QDyn , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

!.........................................................

! Quantum Dynamics ...

 CALL Generate_Structure(1)

 If( NetCharge ) allocate( Net_Charge(Extended_Cell%atoms) )

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 If( DP_field_ )CALL Molecular_DPs( Extended_Cell )

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( DIPOLE_  ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R )  

 If( spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( survival ) then
    
    select case( el_hl_ )

        case( .false. )

            CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, el_FMO , instance="E")

            CALL Huckel_dynamics( Extended_Cell, ExCell_basis, UNI, el_FMO , QDyn=QDyn )

        case( .true. )

            CALL FMO_analysis( Extended_Cell, ExCell_basis, UNI%R, el_FMO , instance="E")

            CALL FMO_analysis( Extended_Cell , ExCell_basis , UNI%R, hl_FMO , instance="H")

            CALL ElHl_dynamics( Extended_Cell , ExCell_basis , UNI , el_FMO , hl_FMO , QDyn )

        end select

 end If

 CALL Dump_stuff( TDOS , PDOS , SPEC , QDyn )

 CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
 CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
 CALL DeAllocate_SPEC( SPEC , flag="dealloc" )
 CALL DeAllocate_QDyn( QDyn , flag="dealloc" )

 include 'formats.h'

end subroutine qdynamics
!
!
!=============================================
! MPI-free version of QCModel_Huckel.f
 subroutine EigenSystem( system , basis , QM )
!=============================================
use Matrix_math
implicit none
type(structure)                             , intent(in)    :: system
type(STO_basis)                             , intent(in)    :: basis(:)
type(R_eigen)                               , intent(inout) :: QM

! local variables ...
real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:)  
real*8  , ALLOCATABLE :: h(:,:) , S_matrix(:,:)
real*8  , ALLOCATABLE :: dumb_S(:,:) , tool(:,:) , S_eigen(:) 
integer               :: i , j , N , info 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

N = size(basis)

CALL Overlap_Matrix( system , basis , S_matrix )

If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(N))

Allocate(      h(N,N) )
Allocate( dumb_S(N,N) )

! clone S_matrix because SYGVD will destroy it ...
dumb_s = S_matrix

If( DP_field_ .OR. Induced_ ) then

!$OMP PARALLEL DO schedule( GUIDED , 10 )
    do j = 1 , N
        do i = j , N
    
            h(i,j) = huckel_with_FIELDS(i,j,S_matrix(i,j),basis)

        end do
    end do  
!$OMP END PARALLEL DO

else
    do j = 1 , N
        do i = j , N

            h(i,j) = huckel(i,j,S_matrix(i,j),basis) 

        end do
    end do
end If

CALL SYGVD( h , dumb_S , QM%erg , 1 , 'V' , 'L' , info )

If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '
!---------------------------------------------------
!   ROTATES THE HAMILTONIAN:  H --> H*S_inv 
!
!   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S.|C> 
!
!   Rv = <AO|MO> coefficients
!---------------------------------------------------

Allocate( Lv(size(basis),size(basis)) )

Lv = h

Deallocate(h)

Allocate( Rv(size(basis), size(basis)) )

!CALL gemm(S_matrix,Lv,Rv,'N','N',D_one,D_zero)
call Multiply( S_matrix, Lv, Rv )

DEALLOCATE( S_matrix )
!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(size(basis),size(basis))) 
! eigenvectors in the rows of QM%L
QM%L = transpose(Lv) 
Deallocate( Lv )

If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(size(basis),size(basis)))
! eigenvectors in the columns of QM%R
QM%R = Rv
Deallocate( Rv )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! save energies of the TOTAL system ...
OPEN(unit=9,file='system-ergs.dat',status='unknown')
    do i = 1 , N
        write(9,*) i , QM%erg(i)
    end do
CLOSE(9)  

If( verbose ) Print*, '>> EigenSystem done <<'

end subroutine EigenSystem
!
!
!
!============================================
 pure function Huckel( i , j , S_ij , basis )
!============================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel
 real*8  :: k_eff , k_WH , c1 , c2 , c3 , c4

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 if (i == j) then
    huckel = basis(i)%IP + basis(i)%V_shift
 else
    c1 = basis(i)%IP - basis(j)%IP
    c2 = basis(i)%IP + basis(j)%IP

    c3 = (c1/c2)*(c1/c2)

    c4 = (basis(i)%V_shift + basis(j)%V_shift)*HALF

    k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

    k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

    huckel = k_eff*S_ij*c2/two + c4*S_ij
 endif

 end function Huckel
!
!
!
!========================================================
 pure function Huckel_with_FIELDS( i , j , S_ij , basis )
!========================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: DP(4)
 real*8  :: r0(3) , Ri(3) , vector(3)
 real*8  :: Huckel_with_FIELDS
 real*8  :: k_eff , k_WH , c1 , c2 , c3
 logical :: flag

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
   
 DP = D_zero

 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

 k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

 huckel_with_FIELDS = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / two

 IF( i == j ) huckel_with_FIELDS = basis(i)%IP 

 flag = ( abs(S_ij) > mid_prec ) 

 IF( flag ) then
    
     r0(1) = ( basis(i)%x + basis(j)%x ) / two
     r0(2) = ( basis(i)%y + basis(j)%y ) / two
     r0(3) = ( basis(i)%z + basis(j)%z ) / two

     Ri(1) = basis(i)%x
     Ri(2) = basis(i)%y
     Ri(3) = basis(i)%z

     vector = DEBYE_inv * DP_matrix_AO(i,j,:) + S_ij * ( Ri - r0 )    ! <== in Angs

     If( Induced_  ) DP = Induced_DP_phi( i , j , basis )

     If( DP_field_ ) DP = DP + DP_phi( i , j , basis )

     huckel_with_FIELDS = huckel_with_FIELDS + S_ij*DP(1) + dot_product( vector(1:3) , DP(2:4) )

 end if

end function Huckel_with_FIELDS
!
!
!
end module qdynamics_m

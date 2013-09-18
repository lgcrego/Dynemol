 module Embedded_FF_Alpha

    use type_m
    use omp_lib
    use constants_m
    use parameters_m                , only : DP_Field_,       &
                                             Induced_ ,       &
                                             verbose
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Overlap_Builder             , only : Overlap_Matrix
    use QCModel_Huckel              , only : Huckel , Huckel_with_Fields
    use DP_main_m                   , only : DP_matrix_AO , Dipole_Moment

    public :: AlphaPolar

    private


    ! module variables ...
    real*8  , ALLOCATABLE :: H0(:,:) , S(:,:)


 contains

!=======================================
 subroutine AlphaPolar( system , basis ) 
!=======================================
implicit none
type(structure)             , intent(inout) :: system
type(STO_basis)             , intent(in)    :: basis(:)

! local variables ...
integer                         :: mm , i , j , xyz
real*8                          :: alpha_ii(3) 
real*8          , ALLOCATABLE   :: H(:,:) , DP_AO(:,:) 
type(R_eigen)                   :: UNI
type(R3_vector)                 :: Induced(-2:2)

! local parameters ...
real*8 , parameter :: base_Field = 5.0d-4           ! <== in Volts/Angs ...
real*8 , parameter :: Debye_unit = 4.803204d0       ! <== 1e*d[Angs]*4.803204 = p[Debye]
real*8 , parameter :: factor     = 14.39965173d0    ! <== 1e/Volt = 14.399 Angs

mm = size(basis)

! build the field independent H and S matrices ...
CALL Build_H0_and_S( system , basis )

ALLOCATE( H(mm,mm) , source=D_zero )

! field dependent hamiltonian and Induced DP moments (DP is calculated in Debyes) ...

ALLOCATE( DP_AO(mm,mm) , source=D_zero )

! for each molecular axis F_xyz ...
do xyz = 1 , 3 

    select case ( xyz )

        case (1)
        forall( i=1:mm , j=1:mm ) DP_AO(i,j) = DP_matrix_AO(i,j,xyz) + basis(i)%x*S(i,j)

        case (2)
        forall( i=1:mm , j=1:mm ) DP_AO(i,j) = DP_matrix_AO(i,j,xyz) + basis(i)%y*S(i,j)

        case (3)
        forall( i=1:mm , j=1:mm ) DP_AO(i,j) = DP_matrix_AO(i,j,xyz) + basis(i)%z*S(i,j)

    end select

    do i = -2 , 2

        If( i /= 0 ) then

            H(:,:) = H0(:,:) +  DP_AO(:,:)*base_Field*float(i)  
            CALL Eigenstates( H , UNI )

            ! Dipole moment in Debye ...
            CALL Dipole_Moment( system , basis , UNI%L , UNI%R , DP_total=Induced(i)%DP )

        end if

    end do

    ! diagonal elements of the Alpha tensor , JCP 109, 7756 (1998) ...
    Alpha_ii(xyz) = two/three * (Induced(1)%DP(xyz) - Induced(-1)%DP(xyz)) - D_one/twelve * (Induced(2)%DP(xyz) - Induced(-2)%DP(xyz))

    ! diagonal elements of the polarizability tensor in Angs^{3} ...
    Alpha_ii(xyz) = ( (Alpha_ii(xyz)/Debye_unit) / base_Field ) * factor    

end do

Print 188 , Alpha_ii , sum( Alpha_ii ) / three
Print 189 , Alpha_ii / (a_Bohr*a_Bohr*a_Bohr) , sum( Alpha_ii ) / (three*a_Bohr*a_Bohr*a_Bohr)

DEALLOCATE( H , H0 , S , DP_AO )

include 'formats.h'

end subroutine AlphaPolar
!
!
!
!
!===========================================
 subroutine Build_H0_and_S( system , basis )
!===========================================
implicit none
type(structure)  , intent(in)    :: system
type(STO_basis)  , intent(in)    :: basis(:)

! local variables ...
integer :: i , j 

verbose = .NOT. verbose

CALL Overlap_Matrix( system , basis , S )

ALLOCATE( H0(size(basis),size(basis)) , source=D_zero)

If( DP_field_ .OR. Induced_ ) then

    !$OMP PARALLEL DO schedule( GUIDED , 10 )
    do j = 1 , size(basis)
        do i = 1 , j
     
            H0(i,j) = huckel_with_FIELDS( i , j , S(i,j) , basis )

        end do
    end do  
    !$OMP END PARALLEL DO

else
 
    do j = 1 , size(basis)
        do i = 1 , j
     
            H0(i,j) = Huckel( i , j , S(i,j) , basis )

        end do
    end do  

end if

end subroutine Build_H0_and_S
!
!
!
!================================
 subroutine Eigenstates( H , QM )
!================================
implicit none
real*8          , allocatable   , intent(inout) :: H(:,:) 
type(R_eigen)                   , intent(inout) :: QM

! local variables ...
integer               :: mm , info
real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) 
real*8  , ALLOCATABLE :: dumb_s(:,:) 

 mm = size( H(:,1) )

 ALLOCATE( dumb_S(mm,mm) , source=S )

 If( .NOT. allocated(QM%erg) ) ALLOCATE( QM%erg(mm) )

 CALL SYGVD( H , dumb_S , QM%erg , 1 , 'V' , 'U' , info )

 DEALLOCATE(dumb_S)

 ALLOCATE( Lv(mm,mm) )

 Lv = H

 ALLOCATE( Rv(mm,mm) )

 CALL gemm( S , Lv , Rv , 'N' , 'N' , D_one , D_zero )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

 If( .NOT. allocated(QM%L) ) ALLOCATE( QM%L(mm,mm) ) 
! eigenvectors in the rows of QM%L
 QM%L = transpose(Lv) 
 DEALLOCATE( Lv )

 If( .NOT. ALLOCATED(QM%R) ) ALLOCATE( QM%R(mm,mm) )
! eigenvectors in the columns of QM%R
 QM%R = Rv
 DEALLOCATE( Rv )

!  the order of storage is the ascending order of eigenvalues
!----------------------------------------------------------

end subroutine Eigenstates
!
!
!
end module Embedded_FF_Alpha

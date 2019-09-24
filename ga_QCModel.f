module GA_QCModel_m

    use type_m
    use constants_m
    use f95_precision
    use blas95
    use lapack95
    use parameters_m            , only : Alpha_Tensor , EnvField_ , Induced_
    use Semi_Empirical_Parms    , only : element => atom  
    use Structure_Builder       , only : Extended_Cell 
    use Overlap_Builder         , only : Overlap_Matrix
    use Hamiltonians            , only : X_ij , even_more_extended_Huckel
    use Multipole_Routines_m    , only : rotationmultipoles ,   &
                                         multipole_messages ,   &
                                         multipoles1c ,         &
                                         multipoles2c 

    public :: MO_erg_diff, GA_eigen, GA_DP_Analysis, Mulliken, AlphaPolar, Bond_Type, MO_character, Localize, Exclude, i_

    private 

    interface Mulliken
        module procedure R_Mulliken
        module procedure C_Mulliken
    end interface

    ! module variables ...
    Real*8  , allocatable :: DP_matrix_AO(:,:,:)
    Real*8  , allocatable :: H0(:,:) , S(:,:)        ! <== to be used by AlphaPolar ...
    integer , allocatable :: occupancy(:)
    integer               :: i_= 0

contains
!
!
!
!====================================================================
 function MO_erg_diff( GA , up , down , dE_ref , weight ) result(cost)
!====================================================================
implicit none
type(R_eigen)            , intent(in) :: GA
integer                  , intent(in) :: up
integer                  , intent(in) :: down
real                     , intent(in) :: dE_ref
real          , optional , intent(in) :: weight

!local variables ...
real   :: w
real*8 :: cost , delta_E

delta_E = GA%erg(up) - GA%erg(down)

w = merge( weight , 1.0 , present(weight) ) 
cost = (delta_E - dE_ref) * w

i_ = i_ + 1

end function MO_erg_diff
!
!
!
!===============================================================
 function Exclude( GA , basis , MO , atom , residue , threshold)
!===============================================================
implicit none
type(R_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
integer         , optional  , intent(in) :: atom(:)
character(len=*), optional  , intent(in) :: residue
real            , optional  , intent(in) :: threshold

! local variables ...
integer               :: i
real*8                :: exclude , population
logical , allocatable :: mask(:) , mask_1(:) , mask_2(:)

allocate( mask  (size(basis)) , source=.false. )
allocate( mask_1(size(basis)) , source=.false. )
allocate( mask_2(size(basis)) , source=.false. )

!====================================================
IF( .NOT. present(atom) ) then
    mask_1 = .true.
else
    do i = 1 , size(atom) 
        where( basis%atom == atom(i) ) mask_1 = .true.
    end do
end IF
!====================================================
IF( .NOT. present(residue) ) then
    mask_2 = .true.
else
    where( basis%residue == residue ) mask_2 = .true.
end IF
!====================================================

! the total mask ...
mask = ( mask_1 .AND. mask_2 )

population = sqrt( sum( GA%L(MO,:) * GA%R(:,MO) , mask ) )

If( present(threshold) ) then
    exclude = merge( D_zero , large , population < threshold )
else
    ! default value is assumed ...
    exclude = merge( D_zero , large , population < 1.0d-3 )
end if

deallocate( mask )

i_ = i_ + 1

end function exclude
!
!
!
!
!===========================================================================
 function Localize( GA , basis , MO , atom , EHSymbol , residue , threshold)
!===========================================================================
implicit none
type(R_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
integer         , optional  , intent(in) :: atom(:)
character(len=*), optional  , intent(in) :: EHSymbol
character(len=*), optional  , intent(in) :: residue
real            , optional  , intent(in) :: threshold

! local variables ...
integer               :: i
real*8                :: Localize , population
logical , allocatable :: mask(:) , mask_1(:) , mask_2(:) , mask_3(:)

allocate( mask  (size(basis)) , source=.false. )
allocate( mask_1(size(basis)) , source=.false. )
allocate( mask_2(size(basis)) , source=.false. )
allocate( mask_3(size(basis)) , source=.false. )

!====================================================
IF( .NOT. present(atom) ) then
    mask_1 = .true.
else
    do i = 1 , size(atom) 
        where( basis%atom == atom(i) ) mask_1 = .true.
    end do
end IF
!====================================================
IF( .NOT. present(residue) ) then
    mask_2 = .true.
else
    where( basis%residue == residue ) mask_2 = .true.
end IF
!====================================================
IF( .NOT. present(EHSymbol) ) then
    mask_3 = .true.
else
    where( basis%EHSymbol == EHSymbol ) mask_3 = .true.
end IF
!====================================================

! the total mask ...
mask = ( mask_1 .AND. mask_2 .AND. mask_3)

population = sqrt( sum( GA%L(MO,:) * GA%R(:,MO) , mask ) )

If( present(threshold) ) then
    localize = merge( D_zero , large , population > threshold )
else
    ! default value is assumed, 85% of localization ...
    localize = merge( D_zero , large , population > 0.85 )
end if

deallocate( mask )

i_ = i_ + 1

end function Localize
!
!
!
!=============================================
 function MO_character( GA , basis , MO , AO )
!=============================================
implicit none
type(R_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
character(len=*), optional  , intent(in) :: AO

! local variables ...
integer               :: l , m
real*8                :: MO_character , population
logical , allocatable :: mask(:) 

 allocate( mask  (size(basis)) , source=.false. )

 select case( AO ) 
      
    case( 's', 'S' )  

        l = 0 ; m = 0

    case( 'py', 'Py' , 'PY' )

        l = 1 ; m = -1

    case( 'pz', 'Pz' , 'PZ' )

        l = 1 ; m = 0
        
    case( 'px', 'Px' , 'PX' )

        l = 1 ; m = +1
        
    case( 'dxy', 'Dxy' , 'DXY' )

        l = 2 ; m = -2
        
    case( 'dyz', 'Dyz' , 'DYZ' )

        l = 2 ; m = -1
        
    case( 'dz2', 'Dz2' , 'DZ2' )

        l = 2 ; m = 0 
        
    case( 'dxz', 'Dxz' , 'DXZ' )

        l = 2 ; m = +1

    case( 'dx2y2', 'Dx2y2' , 'DX2Y2' )

        l = 2 ; m = +2 
  
    case default

        stop " >> error in [MO_character] subroutine check input arguments <<"

 end select

 where( (basis%l == l) .AND. (basis%m == m) ) mask = .true.

 population = sqrt( sum( GA%L(MO,:) * GA%R(:,MO) , mask ) )

 MO_character = merge( D_zero , large , population> HALF )

deallocate( mask )

i_ = i_ + 1

end function MO_character
!
!
!
!======================================================================
 function Bond_Type( system , GA , MO , atom1 , atom2 , AO , instance )
!======================================================================
implicit none
type(structure)  , intent(in) :: system
type(R_eigen)    , intent(in) :: GA
integer          , intent(in) :: MO
integer          , intent(in) :: atom1
integer          , intent(in) :: atom2
character(*)     , intent(in) :: AO
character(len=1) , intent(in) :: instance

real*8 :: bond_type 

! local variables ...
integer :: indx1 , indx2
real*8  :: bond_signal

select case( AO ) 

    case( 's', 'S' )

        indx1 = system% BasisPointer(atom1) + 1
        indx2 = system% BasisPointer(atom2) + 1

    case( 'py', 'Py' , 'PY' )

        indx1 = system% BasisPointer(atom1) + 2
        indx2 = system% BasisPointer(atom2) + 2

    case( 'pz', 'Pz' , 'PZ' )

        indx1 = system% BasisPointer(atom1) + 3
        indx2 = system% BasisPointer(atom2) + 3

    case( 'px', 'Px' , 'PX' )

        indx1 = system% BasisPointer(atom1) + 4
        indx2 = system% BasisPointer(atom2) + 4

    case( 'dxy', 'Dxy' , 'DXY' )

        indx1 = system% BasisPointer(atom1) + 5
        indx2 = system% BasisPointer(atom2) + 5

    case( 'dyz', 'Dyz' , 'DYZ' )

        indx1 = system% BasisPointer(atom1) + 6
        indx2 = system% BasisPointer(atom2) + 6

    case( 'dz2', 'Dz2' , 'DZ2' )

        indx1 = system% BasisPointer(atom1) + 7
        indx2 = system% BasisPointer(atom2) + 7

    case( 'dxz', 'Dxz' , 'DXZ' )

        indx1 = system% BasisPointer(atom1) + 8
        indx2 = system% BasisPointer(atom2) + 8

    case( 'dx2y2', 'Dx2y2' , 'DX2Y2' )

        indx1 = system% BasisPointer(atom1) + 9
        indx2 = system% BasisPointer(atom2) + 9

    case default

        stop " >> error in [bond] subroutine check input arguments <<"

end select

! pre-processing ...
bond_type = D_zero
If( dabs(GA%R(indx1,MO)) < mid_prec .OR. dabs(GA%R(indx2,MO)) < mid_prec ) return

! actual calculations start here ...
bond_signal = sign( 1.d0 , GA%R(indx1,MO) * GA%R(indx2,MO) )

select case ( instance )

    case( '+' )  ! <== Bonding ...

        bond_type = merge( D_zero , large , bond_signal > D_zero )

    case( '-' )  ! <== Anti-Bonding ...

        bond_type = merge( D_zero , large , bond_signal < D_zero )

    case default

        stop " >> error in [Bond_Type] subroutine check input arguments <<"

end select

i_ = i_ + 1

end function
!
!
!
!================================================================================
 function R_Mulliken( GA , basis , MO , atom , AO , EHSymbol , residue , weight )
!================================================================================
implicit none
type(R_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
integer         , optional  , intent(in) :: atom(:)
character(len=*), optional  , intent(in) :: AO
character(len=*), optional  , intent(in) :: EHSymbol
character(len=*), optional  , intent(in) :: residue
real            , optional  , intent(in) :: weight

! local variables ...
integer               :: i , l , m
real*8                :: R_Mulliken
logical , allocatable :: mask(:) , mask_1(:) , mask_2(:) , mask_3(:)  , mask_4(:)  

allocate( mask  (size(basis)) , source=.false. )
allocate( mask_1(size(basis)) , source=.false. )
allocate( mask_2(size(basis)) , source=.false. )
allocate( mask_3(size(basis)) , source=.false. )
allocate( mask_4(size(basis)) , source=.false. )

!====================================================
IF( .NOT. present(atom) ) then
    mask_1 = .true.
else
    do i = 1 , size(atom) 
        where( basis%atom == atom(i) ) mask_1 = .true.
    end do
end IF
!====================================================
IF( .NOT. present(AO) ) then
    mask_2 = .true.
else
    select case( AO ) 
     
       case( 's', 'S' )
     
           l = 0 ; m = 0
     
       case( 'py', 'Py' , 'PY' )
     
           l = 1 ; m = -1

       case( 'pz', 'Pz' , 'PZ' )
     
           l = 1 ; m = 0
           
       case( 'px', 'Px' , 'PX' )
     
           l = 1 ; m = +1
           
       case( 'dxy', 'Dxy' , 'DXY' )
     
           l = 2 ; m = -2
           
       case( 'dyz', 'Dyz' , 'DYZ' )
     
           l = 2 ; m = -1
           
       case( 'dz2', 'Dz2' , 'DZ2' )
     
           l = 2 ; m = 0 
           
       case( 'dxz', 'Dxz' , 'DXZ' )
     
           l = 2 ; m = +1
     
       case( 'dx2y2', 'Dx2y2' , 'DX2Y2' )
     
           l = 2 ; m = +2 
     
       case default
     
           stop " >> error in [Mulliken] subroutine check input arguments <<"

    end select

    where( (basis%l == l) .AND. (basis%m == m) ) mask_2 = .true.

end IF
!====================================================
IF( .NOT. present(EHSymbol) ) then
    mask_3 = .true.
else
    where( basis%EHSymbol == EHSymbol ) mask_3 = .true.
end IF
!====================================================
IF( .NOT. present(residue) ) then
    mask_4 = .true.
else
    where( basis%residue == residue ) mask_4 = .true.
end IF
!====================================================

! the total mask ...
mask = ( mask_1 .AND. mask_2 .AND. mask_3 .AND. mask_4 )

! perform the population analysis ...
R_Mulliken = real( sum( GA%L(MO,:) * GA%R(:,MO) , mask ) )


If( .not. present(weight)) then

    i_ = i_ + 1                     ! <= updat me

else If( weight > 0 ) then 

    ! apply weight ...
    R_Mulliken = R_Mulliken * weight

    i_ = i_ + 1                     ! <= also updat me

end If                              ! <= otherwise, dont update me and dont apply weight

deallocate( mask , mask_1 , mask_2 , mask_3 , mask_4 )

end function R_Mulliken
!
!
!
!===========================================================================
 function C_Mulliken( GA , basis , MO , atom , AO_ang , EHSymbol , residue )
!===========================================================================
implicit none
type(C_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
integer         , optional  , intent(in) :: atom
integer         , optional  , intent(in) :: AO_ang
character(len=*), optional  , intent(in) :: EHSymbol
character(len=*), optional  , intent(in) :: residue

! local variables ...
complex*16            :: C_Mulliken
logical , allocatable :: mask(:) , mask_1(:) , mask_2(:) , mask_3(:) , mask_4(:)

allocate( mask  (size(basis)) , source=.false. )
allocate( mask_1(size(basis)) , source=.false. )
allocate( mask_2(size(basis)) , source=.false. )
allocate( mask_3(size(basis)) , source=.false. )
allocate( mask_4(size(basis)) , source=.false. )

!====================================================
IF( .NOT. present(atom) ) then
    mask_1 = .true.
else
    where( basis%atom == atom ) mask_1 = .true.
end IF
!====================================================
IF( .NOT. present(AO_ang) ) then
    mask_2 = .true.
else
    where( basis%l == AO_ang ) mask_2 = .true.
end IF
!====================================================
IF( .NOT. present(EHSymbol) ) then
    mask_3 = .true.
else
    where( basis%EHSymbol == EHSymbol ) mask_3 = .true.
end IF
!====================================================
IF( .NOT. present(residue) ) then
    mask_4 = .true.
else
    where( basis%residue == residue ) mask_4 = .true.
end IF
!====================================================


mask = ( mask_1 .AND. mask_2 .AND. mask_3 .AND. mask_4)

! perform the population analysis ...
C_Mulliken = sum( GA%L(MO,:) * GA%R(:,MO) , mask )

deallocate( mask , mask_1 , mask_2 , mask_3 , mask_4 )

i_ = i_ + 1

end function C_Mulliken
!
!
!
!===================================================
 subroutine  GA_eigen( system , basis , FMO , flag )
!===================================================
 implicit none
 type(structure)              , intent(in)    :: system
 type(STO_basis)              , intent(in)    :: basis(:)
 type(R_eigen)                , intent(out)   :: FMO       
 integer         , optional   , intent(inout) :: flag 

! local variables ... 
 integer               :: i , j , N , info
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) , s_FMO(:,:) , h_FMO(:,:) , dumb_S(:,:) 

 real*8  , ALLOCATABLE :: S_eigen(:) , tool(:,:)


! local parameters ... 
 real*8  , parameter   :: one = 1.d0 , zero = 0.d0

 N = size(basis)

 ALLOCATE( s_FMO   (size(basis),size(basis)) )
 ALLOCATE( h_FMO   (size(basis),size(basis)) )
 ALLOCATE( dumb_S  (size(basis),size(basis)) )
 ALLOCATE( FMO%erg (size(basis)            ) )

!-----------------------------------------------------------------------

 CALL Overlap_Matrix( system, basis, S_FMO, "GA-CG" )

! clone S_matrix because SYGVD will destroy it ... 
 dumb_S = S_FMO

 If( EnvField_ .OR. Induced_ ) then
     h_FMO = even_more_extended_Huckel( system , basis , S_FMO )
 else
     h_FMO = Build_Huckel( basis , S_FMO )
 end If

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

 CALL SYGVD( h_FMO , dumb_S , FMO%erg , 1 , 'V' , 'U' , info )

 If (info /= 0) write(*,*) 'info = ',info,' in GA_Eigen '
 If ( present(flag) ) flag = info

 !--------------------------------------------------------
 ! Overlap Matrix Factorization: S^(1/2) ...
 Allocate( S_eigen(N) )

 dumb_s = S_FMO

 CALL SYEVD(dumb_S , S_eigen , 'V' , 'L' , info)

 Allocate( tool(N,N) , source = transpose(dumb_S) )

 forall( i=1:N ) tool(:,i) = sqrt(S_eigen) * tool(:,i)

 !now S_matrix = S^(1/2) Lowdin Orthogonalization matrix ...
 CALL gemm(dumb_S , tool , S_FMO , 'N' , 'N')

 DEALLOCATE( S_eigen , dumb_S , tool )

 !---------------------------------------------------
 !RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S^(1/2).|C> 
 !
 !normalizes the L&R eigenvectors as < L(i) | R(i) > = 1
 !---------------------------------------------------

 Allocate( Lv(N,N) )
 Allocate( Rv(N,N) )

 Lv = h_FMO
 Deallocate( h_FMO )

 ! Rv = S^(1/2) * Lv ...
 CALL symm( S_FMO , Lv , Rv )

 ALLOCATE(FMO%R(N,N))
 ! eigenvectors in the columns of QM%R
 FMO%R = Rv

 ALLOCATE(FMO%L(N,N))
 ! eigenvectors in the rows of QM%L
 FMO%L = transpose(FMO%R)

 Deallocate( Lv , Rv , S_FMO )

 end subroutine GA_eigen
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
!======================================================================
 subroutine GA_DP_Analysis( system , basis , L_vec , R_vec , Total_DP )
!======================================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: L_vec(:,:) , R_vec(:,:)
real*8          , intent(out)   :: Total_DP(3)

CALL Build_Dipole_Matrix( system , basis )

CALL Dipole_Moment( system , basis , L_vec , R_vec , Total_DP )

end subroutine GA_DP_Analysis
!
!
!
!==================================================
 subroutine AlphaPolar( system , basis , Alpha_ii ) 
!==================================================
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
real*8           , intent(out)   :: Alpha_ii(3) 

! local variables ...
integer                         :: mm , i , j , xyz
real*8          , ALLOCATABLE   :: H(:,:) , DP_AO(:,:) 
type(R_eigen)                   :: UNI
type(R3_vector)                 :: Induced(-2:2)

! local parameters ...
real*8 , parameter :: conversion_factor = 20.23100999d0  ! <== see AlphaTensor.f
real*8 , parameter :: base_Field        = 5.0d-4

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
            CALL Alpha_eigen( H , UNI )

            ! Dipole moment in Debye ...
            CALL Dipole_Moment( system , basis , UNI%L , UNI%R , Total_DP=Induced(i)%DP )

        end if

    end do

    ! diagonal elements of the Alpha tensor , JCP 109, 7756 (1998) ...
    Alpha_ii(xyz) = two/three * (Induced(1)%DP(xyz) - Induced(-1)%DP(xyz)) - D_one/twelve * (Induced(2)%DP(xyz) - Induced(-2)%DP(xyz))

    ! diagonal elements of the polarizability tensor in a_B^{3} ...
    Alpha_ii(xyz) = ( Alpha_ii(xyz)/base_Field ) * conversion_factor    

end do

DEALLOCATE( H , H0 , S , DP_AO , DP_matrix_AO )

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

CALL Overlap_Matrix( system , basis , S , "GA-CG" )

ALLOCATE( H0(size(basis),size(basis)) , source=D_zero)

H0 = Build_Huckel( basis , S )

end subroutine Build_H0_and_S
!
!
!
!================================
 subroutine Alpha_eigen( H , QM )
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

end subroutine Alpha_eigen
!
!
!
!================================================
 subroutine Build_DIPOLE_Matrix( system , basis )
!================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)

! local variables
real*8  :: expa, expb, Rab 
integer :: a , b , ia , ib , ja , jb 
integer :: na , la , ma 
integer :: nb , lb , mb
integer :: lmult , i , j

real*8  , parameter :: tol = 1.d-10 
integer , parameter :: mxl = 5 , mxmult = 3 , mxlsup = max(mxl,mxmult)
real*8  , parameter :: cutoff_Angs = 10.d0

real*8 , dimension((mxmult+1)**2,-mxl:mxl,-mxl:mxl)        :: qlm
real*8 , dimension(-mxlsup:mxlsup,-mxlsup:mxlsup,0:mxlsup) :: rl , rl2

lmult = 1 ! <== DIPOLE MOMENT

allocate( DP_matrix_AO(size(basis),size(basis),3) , source=D_zero )

do ib = 1 , system%atoms
do ia = 1 , system%atoms  

! calculate rotation matrix for the highest l

    call RotationMultipoles( system , ia , ib , Rab , lmult , rl , rl2 )

    If(Rab > cutoff_Angs) goto 10

    do jb = 1 , element(system%AtNo(ib))%DOS  ;  b = system%BasisPointer(ib) + jb
    do ja = 1 , element(system%AtNo(ia))%DOS  ;  a = system%BasisPointer(ia) + ja

        na = basis(a)%n ;   la = basis(a)%l ;   ma = basis(a)%m
        nb = basis(b)%n ;   lb = basis(b)%l ;   mb = basis(b)%m

        CALL Multipole_Messages(na,nb,la,lb)

!---------------------------------------------------------------------------------------------------- 
!       sum over zeta coefficients
        do i = 1 , basis(a)%Nzeta
        do j = 1 , basis(b)%Nzeta
   
            expa = basis(a)%zeta(i)
            expb = basis(b)%zeta(j)

            if( ia==ib ) then

!               CALLS THE SUBROUTINE FOR THE MULTIPOLES OF ONE-CENTER DISTRIBUTIONS

                qlm = 0.d0   ! check this !!!!

                call multipoles1c(na, la, expa, nb, lb, expb, lmult, qlm)

            else 

!               CALLS THE SUBROUTINE FOR THE MULTIPOLES OF TWO-CENTER DISTRIBUTIONS

                qlm = 0.d0   

                call multipoles2c(na, la, expa, nb, lb, expb, Rab, lmult, rl, qlm)

            end if

!           p_x(a,b) 
            DP_matrix_AO(a,b,1) = DP_matrix_AO(a,b,1) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(4,ma,mb)
!           p_y(a,b)
            DP_matrix_AO(a,b,2) = DP_matrix_AO(a,b,2) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(2,ma,mb)
!           p_z(a,b)
            DP_matrix_AO(a,b,3) = DP_matrix_AO(a,b,3) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(3,ma,mb)

        end do
        end do
!---------------------------------------------------------------------------------------------------- 
    enddo
    enddo
10 end do
end do

end subroutine Build_DIPOLE_Matrix
!
!
!
!=====================================================================
 subroutine Dipole_Moment( system , basis , L_vec , R_vec , Total_DP )
!=====================================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: L_vec(:,:) , R_vec(:,:)
real*8          , intent(out)   :: Total_DP(3)

! local variables ...
integer                       :: i, states, xyz, n_basis, Fermi_state
real*8                        :: Nuclear_DP(3), Electronic_DP(3), Center_of_Charge(3)
real*8          , allocatable :: R_vector(:,:)
real*8          , allocatable :: a(:,:), b(:,:)
type(R3_vector) , allocatable :: origin_Dependent(:), origin_Independent(:)

! local parameters ...
real*8          , parameter   :: Debye_unit = 4.803204d0
real*8          , parameter   :: one = 1.d0 , zero = 0.d0

Center_of_Charge = C_of_C(system)

! atomic positions measured from the Center of Charge ...
allocate(R_vector(system%atoms,3))
forall(xyz=1:3) R_vector(:,xyz) = system%coord(:,xyz) - Center_of_Charge(xyz)

! Nuclear dipole ; if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0)
Nuclear_DP = D_zero

! Electronic dipole 
n_basis     = size(basis)
Fermi_state = system% N_of_electrons/2 + mod( system% N_of_electrons , 2 ) 
If( .not. allocated(occupancy)) then
    allocate(occupancy(Fermi_state), source = 2)
    occupancy(Fermi_state) =  2 - mod( sum( system%Nvalen ) , 2 )
end If

 
allocate( a(n_basis,n_basis) )
allocate( b(n_basis,n_basis) )
allocate( origin_Dependent(Fermi_state) )
allocate( origin_Independent(Fermi_state) )

do xyz = 1 , 3

!   origin dependent DP = sum{C_dagger * vec{R} * S_ij * C}

    forall(states=1:Fermi_state)

        forall(i=1:n_basis) a(states,i) = L_vec(states,i) * R_vector(basis(i)%atom,xyz)

        origin_Dependent(states)%DP(xyz) = occupancy(states) * sum( a(states,:) * R_vec(:,states) )

    end forall    
 
!   origin independent DP = sum{C_dagger * vec{DP_matrix_AO(i,j)} * C}

    b = DP_matrix_AO(:,:,xyz)
       
    CALL gemm(L_vec,b,a,'N','N',one,zero)    

    forall(states=1:Fermi_state) origin_Independent(states)%DP(xyz) = occupancy(states) * sum(a(states,:)*L_vec(states,:))

end do

forall(xyz=1:3) Electronic_DP(xyz) = sum( origin_Dependent%DP(xyz) + origin_Independent%DP(xyz) )
 
Total_DP = ( Nuclear_DP - Electronic_DP ) * Debye_unit

deallocate( R_vector , a , b   )
deallocate( origin_Dependent   )
deallocate( origin_Independent )

! AlphaPolar will deallocate it ...
If( .not. Alpha_Tensor ) deallocate( DP_matrix_AO )

end subroutine Dipole_Moment
!
!
!
!================================================
 function C_of_C(a)    result( Center_of_Charge )
!================================================
implicit none
type(structure) , intent(in) :: a

! local variables
real*8 , allocatable :: Qi_Ri(:,:) 
real*8               :: total_valence , Center_of_Charge(3)
integer              :: i , j

! sum_i = (q_i * vec{r}_i) / sum_i q_i ...

 allocate(Qi_Ri(a%atoms,3))

 forall(j=1:3,i=1:a%atoms) Qi_Ri(i,j) = element(a%AtNo(i))% Nvalen * a%coord(i,j)

 total_valence = sum( element(a%AtNo(:))% Nvalen )

 forall(j=1:3) Center_of_Charge(j) = sum(Qi_Ri(:,j)) / total_valence

 deallocate(Qi_Ri)

end function C_of_C
!
!
!
end module GA_QCModel_m

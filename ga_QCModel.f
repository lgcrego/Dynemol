module GA_QCModel_m

    use type_m
    use constants_m
    use f95_precision
    use blas95
    use lapack95
    use parameters_m            , only : Alpha_Tensor
    use Semi_Empirical_Parms    , only : element => atom  
    use Structure_Builder       , only : Extended_Cell 
    use Overlap_Builder         , only : Overlap_Matrix
    use Multipole_Routines_m    , only : rotationmultipoles ,   &
                                         multipole_messages ,   &
                                         multipoles1c ,         &
                                         multipoles2c 

    public ::  GA_eigen , GA_DP_Analysis , Mulliken , AlphaPolar

    private 

    interface Mulliken
        module procedure R_Mulliken
        module procedure C_Mulliken
    end interface

    ! module variables ...
    Real*8  , allocatable :: DP_matrix_AO(:,:,:)
    Real*8  , allocatable :: H0(:,:) , S(:,:)        ! <== to be used by AlphaPolar ...
    integer , allocatable :: occupancy(:)

contains
!
!
!================================================================================
 pure function R_Mulliken( GA , basis , MO , atom , AO_ang , EHSymbol , residue )
!================================================================================
implicit none
type(R_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
integer         , optional  , intent(in) :: atom
integer         , optional  , intent(in) :: AO_ang
character(len=*), optional  , intent(in) :: EHSymbol
character(len=*), optional  , intent(in) :: residue

! local variables ...
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

! the total mask ...
mask = ( mask_1 .AND. mask_2 .AND. mask_3 .AND. mask_4 )

! perform the population analysis ...
R_Mulliken = sum( GA%L(MO,:) * GA%R(:,MO) , mask )

deallocate( mask , mask_1 , mask_2 , mask_3 , mask_4 )

end function R_Mulliken
!
!
!
!================================================================================
 pure function C_Mulliken( GA , basis , MO , atom , AO_ang , EHSymbol , residue )
!================================================================================
implicit none
type(C_eigen)               , intent(in) :: GA
type(STO_basis)             , intent(in) :: basis(:)
integer                     , intent(in) :: MO
integer         , optional  , intent(in) :: atom
integer         , optional  , intent(in) :: AO_ang
character(len=*), optional  , intent(in) :: EHSymbol
character(len=*), optional  , intent(in) :: residue

! local variables ...
real*8                :: C_Mulliken
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
 integer               :: i , j , info
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) , s_FMO(:,:) , h_FMO(:,:) , dumb_S(:,:) 

! local parameters ... 
 real*8  , parameter   :: one = 1.d0 , zero = 0.d0

 ALLOCATE( s_FMO   (size(basis),size(basis)) )
 ALLOCATE( h_FMO   (size(basis),size(basis)) )
 ALLOCATE( dumb_S  (size(basis),size(basis)) )
 ALLOCATE( FMO%erg (size(basis)            ) )

!-----------------------------------------------------------------------

 CALL Overlap_Matrix( system, basis, S_FMO, "GA-CG" )

! clone S_matrix because SYGVD will destroy it ... 
 dumb_S = S_FMO

 DO j = 1 , size(basis)
   DO i = 1 , j 

      h_FMO(i,j) = Huckel_bare( i, j, S_FMO(i,j), basis )     !! <== define h_FMO
 
   END DO
 END DO

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

 CALL SYGVD(h_FMO,dumb_S,FMO%erg,1,'V','U',info)

 If (info /= 0) write(*,*) 'info = ',info,' in SYGVD/eigen_FMO '
 If ( present(flag) ) flag = info

 DEALLOCATE(dumb_S)

!     ---------------------------------------------------
!   ROTATES THE HAMILTONIAN:  H --> H*S_inv 
!
!   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S.|C> 
!
!   Rv = <AO|MO> coefficients
!     ---------------------------------------------------

 ALLOCATE(Lv(size(basis),size(basis)))

    Lv = h_FMO

 DEALLOCATE(h_FMO)

 ALLOCATE(Rv(size(basis),size(basis)))

    CALL gemm(S_FMO,Lv,Rv,'N','N',one,zero)  

 DEALLOCATE( S_FMO )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

 ALLOCATE(FMO%L(size(basis),size(basis))) 
! eigenvectors in the rows of QM%L
    FMO%L = transpose(Lv)                 
 DEALLOCATE( Lv )

 ALLOCATE(FMO%R(size(basis),size(basis)))
! eigenvectors in the columns of QM%R
    FMO%R = Rv             
 DEALLOCATE( Rv )

!  the order of storage is the ascending order of eigenvalues
!----------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 end subroutine GA_eigen
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

do j = 1 , size(basis)
    do i = 1 , j
     
        H0(i,j) = Huckel_bare( i , j , S(i,j) , basis )

    end do
end do  

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
Nuclear_DP(xyz) = D_zero

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
!=================================================
 pure function Huckel_bare( i , j , S_ij , basis )
!=================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: k_eff , k_WH , Huckel_bare , c1 , c2 , c3 , c4

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 if (i == j) then
    huckel_bare = basis(i)%IP + basis(i)%V_shift
 else
    c1 = basis(i)%IP - basis(j)%IP
    c2 = basis(i)%IP + basis(j)%IP

    c3 = (c1/c2)*(c1/c2)

    c4 = (basis(i)%V_shift + basis(j)%V_shift)*HALF

    k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

    k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

    huckel_bare = k_eff*S_ij*c2/two + c4*S_ij
 endif

 end function Huckel_bare
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

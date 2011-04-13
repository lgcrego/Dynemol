module DP_excited_m

    use type_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use parameters_m            , only : OPT_basis ,                    &
                                         hole_state ,                   &
                                         excited_state => initial_state         
    use Allocation_m            , only : Allocate_Structures
    use Semi_Empirical_Parms    , only : atom ,                         &
                                         Include_OPT_parameters
    use Overlap_Builder         , only : Overlap_Matrix
    use Structure_Builder       , only : Basis_Builder
    use Multipole_Routines_m    , only : rotationmultipoles ,           &
                                         multipole_messages ,           &
                                         multipoles1c ,                 &
                                         multipoles2c

    public :: el_hl_StaticDPs 

    private

    type(R3_vector) , allocatable :: moiety_DP_matrix_AO(:,:)

 contains
!
!
!
!===============================================
 function el_hl_StaticDPs( system , instance )
!===============================================
 implicit none
 type(structure) , intent(in)  :: system
 character(*)    , intent(in)  :: instance

! local variables ...
 type(structure)                 :: FMO_system
 type(STO_basis) , allocatable   :: FMO_basis(:)
 type(R_eigen)                   :: FMO
 integer                         :: i
 real*8                          :: el_hl_StaticDPs(3) , DP_FMO(3)

! specialFMO_system <= molecule with FMO = .true. ...

 FMO_system%atoms = count( system%FMO )

 CALL Allocate_Structures( FMO_system%atoms , FMO_system )

 forall(i=1:3)
 FMO_system%coord(:,i) =  pack( system%coord(:,i) , system%FMO ) 
 end forall
 FMO_system%AtNo       =  pack( system%AtNo       , system%FMO ) 
 FMO_system%Nvalen     =  pack( system%Nvalen     , system%FMO ) 
 FMO_system%k_WH       =  pack( system%k_WH       , system%FMO )
 FMO_system%symbol     =  pack( system%Symbol     , system%FMO )
 FMO_system%fragment   =  pack( system%fragment   , system%FMO )
 FMO_system%residue    =  pack( system%residue    , system%FMO )
 FMO_system%nr         =  pack( system%nr         , system%FMO )
 FMO_system%MMSymbol   =  pack( system%MMSymbol   , system%FMO )
 FMO_system%FMO        =  pack( system%FMO        , system%FMO )
 FMO_system%copy_No    =  0

 FMO_system%N_of_electrons = sum( FMO_system%Nvalen )

 CALL Basis_Builder( FMO_system , FMO_basis )

 If( OPT_basis ) CALL Include_OPT_parameters( FMO_basis )

 CALL eigen_FMO( FMO_system , FMO_basis , FMO )

 CALL Build_DIPOLE_Matrix( FMO_system , FMO_basis )
 
 CALL DP_moments( FMO_system , FMO_basis , FMO%L , FMO%R , instance , DP_FMO )

 el_hl_StaticDPs = DP_FMO 

 DeAllocate( FMO_basis )
 DeAllocate( moiety_DP_matrix_AO )
 DeAllocate( FMO%L , FMO%R , FMO%erg )

 end function el_hl_StaticDPs
!
!
!
!=============================================================================
 subroutine DP_moments( system , basis , L_vec , R_vec , instance , Total_DP )
!=============================================================================
implicit none
type(structure) , intent(in)  :: system
type(STO_basis) , intent(in)  :: basis(:)
real*8          , intent(in)  :: L_vec(:,:) 
real*8          , intent(in)  :: R_vec(:,:)
character(*)    , intent(in)  :: instance
real*8          , intent(out) :: Total_DP(3)

! local variables ...
integer                       :: i , j , xyz , n_basis 
real*8                        :: total_valence , target_state , Q_center(3)
real*8          , allocatable :: a(:) , b(:,:) , R_vector(:,:) , Qi_Ri(:,:)
type(R3_vector)               :: origin_Dependent , origin_Independent

!-----------------------------------------------------------------------------------------------
! calculate Center_of_Charge for special_FMO

allocate(Qi_Ri(system%atoms,3) , source = D_zero )     

forall( j=1:3 , i=1:system%atoms ) Qi_Ri(i,j) = system%Nvalen(i) * system%coord(i,j)

total_valence = sum( system%Nvalen )

forall(j=1:3) Q_center(j) = sum( Qi_Ri(:,j) ) / total_valence

deallocate( Qi_Ri )
!-----------------------------------------------------------------------------------------------

! atomic positions measured from the Center of Charge ...
allocate(R_vector(system%atoms,3))
forall(xyz=1:3) R_vector(:,xyz) = system%coord(:,xyz) - Q_center(xyz)

! Nuclear dipole ; if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0)

n_basis = size(basis)

allocate( a(n_basis)         , source = D_zero )
allocate( b(n_basis,n_basis) , source = D_zero )

select case ( instance )

    case( "hole" )

        target_state = hole_state

    case( "electron" )

        target_state = excited_state

end select


do xyz = 1 , 3

    ! origin dependent DP = sum{bra * vec{R} * S_ij * ket}
    
    forall(i=1:n_basis) a(i) = L_vec(target_state,i) * R_vector(basis(i)%atom,xyz)

    origin_Dependent%DP(xyz) =  sum( a(:)*R_vec(:,target_state) ) 

    ! origin independent DP = sum{C_dagger * vec{DP_matrix_AO(i,j)} * C}

    b = moiety_DP_matrix_AO%DP(xyz)

    forall(i=1:n_basis) a(i) = sum( L_vec(target_state,:)*b(:,i) )

    origin_Independent%DP(xyz) = sum( a(:)*L_vec(target_state,:) ) 
       
end do

Total_DP = origin_Dependent%DP + origin_Independent%DP 

deallocate( a , b , R_vector )

end subroutine DP_moments
!
!
!
!=============================================
 subroutine  eigen_FMO( system , basis , FMO )
!=============================================
 implicit none
 type(structure) , intent(in)  :: system
 type(STO_basis) , intent(in)  :: basis(:)
 type(R_eigen)   , intent(out) :: FMO       

! local variables ... 
 integer               :: i , j , info
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) , s_FMO(:,:) , h_FMO(:,:) , dumb_S(:,:) 


 ALLOCATE( s_FMO   (size(basis),size(basis)) )
 ALLOCATE( h_FMO   (size(basis),size(basis)) )
 ALLOCATE( dumb_S  (size(basis),size(basis)) )
 ALLOCATE( FMO%erg (size(basis)            ) )

!-----------------------------------------------------------------------

 CALL Overlap_Matrix( system, basis, S_FMO, purpose='FMO' )

! clone S_matrix because SYGVD will destroy it ... 
 dumb_S = S_FMO

 DO j = 1 , size(basis)
   DO i = 1 , j 

      h_FMO(i,j) = huckel_Molecule( i, j, S_FMO(i,j), basis )     !! <== define h_FMO
 
   END DO
 END DO

!-------- solve generalized eH eigenvalue problem H*Q = E*S*Q

 CALL SYGVD(h_FMO,dumb_S,FMO%erg,1,'V','U',info)

 If (info /= 0) write(*,*) 'info = ',info,' in SYGVD/eigen_FMO '

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

    CALL gemm(S_FMO,Lv,Rv,'N','N',D_one,D_zero)  

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

 end subroutine eigen_FMO
!
!
!
!
!=====================================================
 pure function Huckel_Molecule( i , j , S_ij , basis )
!=====================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: k_eff , k_WH , Huckel_Molecule , c1 , c2 , c3

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / 2.d0

 k_eff = k_WH + c3 + c3 * c3 * (1.d0 - k_WH)

 huckel_Molecule = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / 2.d0

 IF(i == j) huckel_Molecule = basis(i)%IP

 end function Huckel_Molecule
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
real*8  :: expa, expb, xab , yab , zab , Rab 
integer :: AtNo_a , AtNo_b
integer :: a , b , ia , ib , ja , jb , i , j
integer :: na , la , ma 
integer :: nb , lb , mb
integer :: lmult

real*8  , parameter :: tol = 1.d-10 
integer , parameter :: mxl = 5 , mxmult = 3 , mxlsup = max(mxl,mxmult)
real*8  , parameter :: cutoff_Angs = 10.d0

real*8 , dimension((mxmult+1)**2,-mxl:mxl,-mxl:mxl)        :: qlm
real*8 , dimension(-mxlsup:mxlsup,-mxlsup:mxlsup,0:mxlsup) :: rl , rl2

lmult = 1 ! <== DIPOLE MOMENT

allocate( moiety_DP_matrix_AO(size(basis),size(basis)) )

forall(i=1:3) moiety_DP_matrix_AO(:,:)%dp(i) = 0.d0

do ib = 1 , system%atoms
do ia = 1 , system%atoms  

! calculate rotation matrix for the highest l

    call RotationMultipoles(system,ia,ib,xab,yab,zab,Rab,lmult,rl,rl2)

    If(Rab*a_Bohr > cutoff_Angs) goto 10

    do jb = 1 , atom(system%AtNo(ib))%DOS  ;  b = system%BasisPointer(ib) + jb
    do ja = 1 , atom(system%AtNo(ia))%DOS  ;  a = system%BasisPointer(ia) + ja

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

                qlm = 0.d0   

                call multipoles1c(na, la, expa, nb, lb, expb, lmult, qlm)

            else 

!               CALLS THE SUBROUTINE FOR THE MULTIPOLES OF TWO-CENTER DISTRIBUTIONS

                qlm = 0.d0   

                call multipoles2c(na, la, expa, nb, lb, expb, xab, yab, zab, Rab, lmult, rl, rl2, qlm)

            end if

!           p_x(a,b) 
            moiety_DP_matrix_AO(a,b)%dp(1) = moiety_DP_matrix_AO(a,b)%dp(1) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(4,ma,mb)
!           p_y(a,b)
            moiety_DP_matrix_AO(a,b)%dp(2) = moiety_DP_matrix_AO(a,b)%dp(2) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(2,ma,mb)
!           p_z(a,b)
            moiety_DP_matrix_AO(a,b)%dp(3) = moiety_DP_matrix_AO(a,b)%dp(3) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(3,ma,mb)

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
end module DP_excited_m

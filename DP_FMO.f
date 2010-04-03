 module DP_FMO_m

    use type_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Allocation_m            , only : Allocate_Structures
    use Semi_Empirical_Parms    , only : atom ,                 &
                                         Include_OPT_parameters
    use Overlap_Builder         , only : Overlap_Matrix
    use Structure_Builder       , only : Basis_Builder
    use Multipole_Core          , only : rotationmultipoles ,   &
                                         multipole_messages ,   &
                                         multipoles1c ,         &
                                         multipoles2c

    public :: DP_FMO_analysis

    private

    type(R3_vector) , allocatable :: DP_matrix_AO(:,:)

 contains
!
!
!
!=============================================================
 subroutine DP_FMO_analysis( system , Q_center , DP_FMO , nr )
!=============================================================
 type(structure) , intent(in)    :: system
 real*8          , intent(in)    :: Q_center(3)
 real*8          , intent(inout) :: DP_FMO(3)
 integer         , intent(in)    :: nr 

! local variables ...
 type(structure)                 :: FMO_system
 type(STO_basis) , allocatable   :: FMO_basis(:)
 type(R_eigen)                   :: FMO

! FMO_system = solvent molecule with residue # nr ...

 FMO_system%atoms = count( system%nr == nr )

 CALL Allocate_Structures( FMO_system%atoms , FMO_system )

 forall(i=1:3)
 FMO_system%coord(:,i) =  pack(system%coord(:,i) , system%nr == nr ) 
 end forall
 FMO_system%AtNo       =  pack( system%AtNo      , system%nr == nr ) 
 FMO_system%k_WH       =  pack( system%k_WH      , system%nr == nr )
 FMO_system%symbol     =  pack( system%Symbol    , system%nr == nr )
 FMO_system%fragment   =  pack( system%fragment  , system%nr == nr )
 FMO_system%nr         =  pack( system%nr        , system%nr == nr )
 FMO_system%MMSymbol   =  pack( system%MMSymbol  , system%nr == nr )
 FMO_system%copy_No    =  0

 FMO_system%N_of_electrons = sum( atom(FMO_system%AtNo(:))%Nvalen )

 CALL Basis_Builder( FMO_system , FMO_basis )

 If( OPT_basis ) CALL Include_OPT_parameters( FMO_basis )

 CALL DP_eigen_FMO( FMO_system , FMO_basis , FMO )

 CALL Build_DIPOLE_Matrix( FMO_system , FMO_basis )

 CALL Dipole_Moment( FMO_system , FMO_basis , FMO%L , FMO%R , Q_center , DP_FMO )

 DeAllocate( FMO_basis )
 DeAllocate( DP_matrix_AO )
 DeAllocate( FMO%L , FMO%R , FMO%erg )

 end subroutine DP_FMO_analysis
!
!
!
!==============================================
 subroutine  DP_eigen_FMO( system, basis, FMO )
!==============================================
 type(structure) , intent(in)  :: system
 type(STO_basis) , intent(in)  :: basis(:)
 type(R_eigen)   , intent(out) :: FMO       

! local variables ... 
 integer               :: i , j
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) , s_FMO(:,:) , h_FMO(:,:) , dumb_S(:,:) 

! local parameters ... 
 real*8  , parameter   :: one = 1.d0 , zero = 0.d0

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

 end subroutine DP_eigen_FMO
!
!
!
!
!================================================
 subroutine Build_DIPOLE_Matrix( system , basis )
!================================================
implicit real*8 (a-h,o-z)
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)

! local variables
real*8  :: expa, expb, xab , yab , zab , Rab 
integer :: AtNo_a , AtNo_b
integer :: a , b , ia , ib , ja , jb 
integer :: na , la , ma 
integer :: nb , lb , mb
integer :: lmult

real*8  , parameter :: tol = 1.d-10 
integer , parameter :: mxl = 5 , mxmult = 3 , mxlsup = max(mxl,mxmult)
real*8  , parameter :: cutoff_Angs = 10.d0

real*8 , dimension((mxmult+1)**2,-mxl:mxl,-mxl:mxl)        :: qlm
real*8 , dimension(-mxlsup:mxlsup,-mxlsup:mxlsup,0:mxlsup) :: rl , rl2

lmult = 1 ! <== DIPOLE MOMENT

allocate( DP_matrix_AO(size(basis),size(basis)) )

forall(i=1:3) DP_matrix_AO(:,:)%dp(i) = 0.d0

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

                qlm = 0.d0   ! check this !!!!

                call multipoles1c(na, la, expa, nb, lb, expb, lmult, qlm)

            else 

!               CALLS THE SUBROUTINE FOR THE MULTIPOLES OF TWO-CENTER DISTRIBUTIONS

                qlm = 0.d0   

                call multipoles2c(na, la, expa, nb, lb, expb, xab, yab, zab, Rab, lmult, rl, rl2, qlm)

            end if

!           p_x(a,b) 
            DP_matrix_AO(a,b)%dp(1) = DP_matrix_AO(a,b)%dp(1) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(4,ma,mb)
!           p_y(a,b)
            DP_matrix_AO(a,b)%dp(2) = DP_matrix_AO(a,b)%dp(2) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(2,ma,mb)
!           p_z(a,b)
            DP_matrix_AO(a,b)%dp(3) = DP_matrix_AO(a,b)%dp(3) + basis(a)%coef(i)*basis(b)%coef(j)*qlm(3,ma,mb)

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
!========================================================================================
 subroutine Dipole_Moment( system , basis , L_vec , R_vec , Center_of_Charge , Total_DP )
!========================================================================================

type(structure) , intent(in)  :: system
type(STO_basis) , intent(in)  :: basis(:)
real*8          , intent(in)  :: L_vec(:,:) , R_vec(:,:)
real*8          , intent(in)  :: Center_of_Charge(3)
real*8          , intent(out) :: Total_DP(3)

! local variables ...
integer                       :: i, j, states, xyz, n_basis, Fermi_state
real*8                        :: Nuclear_DP(3), Electronic_DP(3) 
real*8          , allocatable :: R_vector(:,:)
real*8          , allocatable :: a(:,:), b(:,:)
type(R3_vector) , allocatable :: origin_Dependent(:), origin_Independent(:)

! local parameters ...
real*8          , parameter   :: Debye_unit = 4.803204d0
real*8          , parameter   :: one = 1.d0 , zero = 0.d0


! atomic positions measured from the Center of Charge ...
allocate(R_vector(system%atoms,3))
forall(xyz=1:3) R_vector(:,xyz) = system%coord(:,xyz) - Center_of_Charge(xyz)

! Nuclear dipole ; if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0)
forall(xyz=1:3) Nuclear_DP(xyz) = sum( atom( system%AtNo(:) )%Nvalen * R_vector(:,xyz) )

! Electronic dipole 
n_basis      =  size(basis)
Fermi_state  =  system%N_of_electrons / 2
 
allocate( a(n_basis,n_basis) )
allocate( b(n_basis,n_basis) )
allocate( origin_Dependent(Fermi_state) )
allocate( origin_Independent(Fermi_state) )

do xyz = 1 , 3

!   origin dependent DP = sum{C_dagger * vec{R} * S_ij * C}

    forall(states=1:Fermi_state)

        forall(i=1:n_basis) a(states,i) = L_vec(states,i) * R_vector(basis(i)%atom,xyz)

        origin_Dependent(states)%DP(xyz) = 2.d0 * sum( a(states,:) * R_vec(:,states) )

    end forall    
 
!   origin independent DP = sum{C_dagger * vec{DP_matrix_AO(i,j)} * C}

    b = DP_matrix_AO%DP(xyz)
       
    CALL gemm(L_vec,b,a,'N','N',one,zero)    

    forall(states=1:Fermi_state) origin_Independent(states)%DP(xyz) = 2.d0 * sum(a(states,:)*L_vec(states,:))

end do

forall(xyz=1:3) Electronic_DP(xyz) = sum( origin_Dependent%DP(xyz) + origin_Independent%DP(xyz) )
 
Total_DP = ( Nuclear_DP - Electronic_DP ) * Debye_unit

deallocate(R_vector,a,b)
deallocate(origin_Dependent)
deallocate(origin_Independent)

end subroutine Dipole_Moment
!
!
!
!=====================================================
 pure function Huckel_Molecule( i , j , S_ij , basis )
!=====================================================
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: k_eff , k_WH , Huckel_Molecule

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
end module DP_FMO_m

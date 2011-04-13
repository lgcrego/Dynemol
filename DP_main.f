module DP_main_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Semi_Empirical_Parms    , only  : atom
    use parameters_m            , only  : verbose ,                             &
                                          DP_Moment ,                           &
                                          static ,                              &
                                          hole_state ,                          &
                                          excited_state => initial_state        ! for static calculations initial state = excited state
    use Multipole_Routines_m    , only  : rotationmultipoles ,                  &
                                          multipole_messages ,                  &
                                          multipoles1c ,                        &
                                          multipoles2c ,                        &
                                          util_multipoles
    use DP_excited_m            , only  : el_hl_StaticDPs
    use FMO_m                   , only  : eh_tag


    type(R3_vector) , allocatable , public , protected :: DP_matrix_AO(:,:)

    public :: Dipole_Matrix 
    public :: Dipole_Moment

    private

contains
!
!
!
!=====================================================================
 subroutine Dipole_Matrix( system , basis , L_vec , R_vec , Total_DP )
!=====================================================================
implicit none
type(structure)             , intent(inout) :: system
type(STO_basis)             , intent(in)    :: basis(:)
complex*16      , optional  , intent(in)    :: L_vec(:,:) , R_vec(:,:)
real*8          , optional  , intent(out)   :: Total_DP(3) 

! local variables ...
real*8  :: Sparsity(3)
integer :: NonZero(3) , M_size , i

If( verbose ) Print 153
!----------------------------------------------------------
!       initialize DIPOLE MATRIX M(i,j)[x,y,z]

 CALL Util_Multipoles

! size of M matrix ...
M_size = sum(atom(system%AtNo)%DOS)

If( allocated(DP_matrix_AO) ) deallocate( DP_matrix_AO )
allocate(DP_matrix_AO(M_size,M_size))

CALL Build_DIPOLE_Matrix(system,basis)

forall(i=1:3) NonZero(i) = count(DP_matrix_AO(:,:)%dp(i) /= 0.d0)

Sparsity(:) = dfloat(NonZero(:))/dfloat((M_size**2))

If( verbose ) Print 73, Sparsity  

! execute only for static calculations ...
if ( DP_Moment .AND. static ) CALL Dipole_Moment( system , basis , L_vec , R_vec , DP_total=Total_DP ) 

!----------------------------------------------------------
If( verbose ) then
    Print*, '>> Dipole Moment done <<'
    Print 155
end If

 include 'formats.h'

end subroutine Dipole_Matrix
!
!
!
!==================================================================================================
 subroutine Dipole_Moment( system , basis , L_vec , R_vec , AO_bra , AO_ket , Dual_ket , DP_total )
!==================================================================================================
implicit none
type(structure)             , intent(inout) :: system
type(STO_basis)             , intent(in)    :: basis(:)
complex*16                  , intent(in)    :: L_vec    (:,:) 
complex*16                  , intent(in)    :: R_vec    (:,:)
complex*16      , optional  , intent(in)    :: AO_bra   (:,:) 
complex*16      , optional  , intent(in)    :: AO_ket   (:,:) 
complex*16      , optional  , intent(in)    :: Dual_ket (:,:)
real*8          , optional  , intent(out)   :: DP_total (3) 

! local variables ...
integer                       :: i, j, states, xyz, n_basis, Fermi_state
real*8                        :: Nuclear_DP(3), Electronic_DP(3), hole_DP(3), excited_DP(3), Total_DP(3)
real*8          , allocatable :: R_vector(:,:)
complex*16      , allocatable :: a(:,:), b(:,:)
logical         , allocatable :: AO_mask(:)
type(R3_vector) , allocatable :: origin_Dependent(:), origin_Independent(:)

! local parameters ...
real*8          , parameter   :: Debye_unit = 4.803204d0

! define system for DP_Moment calculation ...
allocate( AO_mask(size(basis)) , source = .true. )
AO_mask = merge( basis%FMO , AO_mask , any(basis%FMO) )

! if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0)
CALL Center_of_Charge( system , R_vector )

Nuclear_DP = D_zero 

! Electronic dipole 
n_basis      =  size(basis)
Fermi_state  =  sum( system%Nvalen ) / two
 
allocate( a(n_basis,n_basis)              , source = C_zero )
allocate( b(n_basis,n_basis)              , source = C_zero )
allocate( origin_Dependent  (Fermi_state) )
allocate( origin_Independent(Fermi_state) )

do xyz = 1 , 3

!   origin dependent DP = sum{C_dagger * vec{R} * S_ij * C}

    forall( states=1:Fermi_state )

        forall( i=1:n_basis ) a(states,i) = L_vec(states,i) * R_vector(basis(i)%atom,xyz)

        origin_Dependent(states)%DP(xyz)  = two * real( sum( a(states,:)*R_vec(:,states) , AO_mask ) )

    end forall    

!   origin independent DP = sum{C_dagger * vec{DP_matrix_AO(i,j)} * C}

    b = DP_matrix_AO%DP(xyz)

    CALL gemm( L_vec , b , a , 'N' , 'N' , C_one , C_zero )    

    forall( states=1:Fermi_state ) origin_Independent(states)%DP(xyz) = two * real( sum( a(states,:)*L_vec(states,:) , AO_mask ) )

end do

deallocate(a,b)

!--------------------------------------------------------------------------------------
! Build DP_Moment ...
!--------------------------------------------------------------------------------------

! contribution from the valence states ...
forall(xyz=1:3) Electronic_DP(xyz) = sum( origin_Dependent%DP(xyz) + origin_Independent%DP(xyz) )

! contribution from the hole and electronic-wavepackets ... 
! excited-state case: hole_state /= 0 ...
If( hole_state /= I_zero ) then

    If( static ) then

        hole_DP    = el_hl_StaticDPs( system , instance="hole" )

        excited_DP = el_hl_StaticDPs( system , instance="electron" )

    else

        If( (eh_tag(1) /= "el") .OR. (eh_tag(2) /= "hl") ) pause ">>> check call to wavepacket_DP <<<"

        hole_DP    =  wavepacket_DP( basis , AO_mask , R_vector , AO_bra(:,2) , AO_ket(:,2) , Dual_ket(:,2) )

        excited_DP =  wavepacket_DP( basis , AO_mask , R_vector , AO_bra(:,1) , AO_ket(:,1) , Dual_ket(:,1) )

    end If

    Electronic_DP = Electronic_DP - hole_DP + excited_DP

end If

Total_DP = ( Nuclear_DP - Electronic_DP ) * Debye_unit

!--------------------------------------------------------------------------------------

Print 154, Total_DP, sqrt(sum(Total_DP*Total_DP))

If( present(DP_total) ) DP_total = Total_DP

deallocate(R_vector , AO_mask)
deallocate(origin_Dependent)
deallocate(origin_Independent)

include 'formats.h'

end subroutine Dipole_Moment
!
!
!
!=============================================================================
 pure function wavepacket_DP( basis , mask , R_vector , bra , ket , Dual_ket )
!=============================================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
logical         , intent(in)    :: mask(:)
real*8          , intent(in)    :: R_vector(:,:)
complex*16      , intent(in)    :: bra(:) 
complex*16      , intent(in)    :: ket(:) 
complex*16      , intent(in)    :: Dual_ket(:)
real*8                          :: wavepacket_DP(3)

! local variables ...
integer                         :: i, j, xyz, n_basis
complex*16      , allocatable   :: a(:), b(:,:)
type(R3_vector)                 :: origin_Dependent, origin_Independent

n_basis = size(basis)
 
allocate( a(n_basis)         , source = C_zero )
allocate( b(n_basis,n_basis) , source = C_zero )

do xyz = 1 , 3

        ! origin dependent DP = sum{bra * vec{R} * S_ij * ket}

        forall( i=1:n_basis ) a(i) = bra(i) * R_vector(basis(i)%atom,xyz)

        origin_Dependent%DP(xyz) = real( sum( a(:)*Dual_ket(:) , mask ) )

        ! origin independent DP = sum{bra * vec{DP_matrix_AO(i,j)} * ket}

        b = DP_matrix_AO%DP(xyz)

        CALL gemv( b , ket , a , C_one , C_zero )    

        origin_Independent%DP(xyz) = real( sum( bra(:)*a(:) , mask ) )

end do

wavepacket_DP = origin_Dependent%DP + origin_Independent%DP 

deallocate( a , b )

end function wavepacket_DP
!
!
!
!===========================================
 subroutine Center_of_Charge( a , R_vector )
!===========================================
implicit none
type(structure)                 , intent(inout) :: a
real*8          , allocatable   , intent(out)   :: R_vector(:,:)

! local variables ...
integer               :: i , j
real*8                :: total_valence
real*8  , allocatable :: Qi_Ri(:,:) 
logical , allocatable :: mask(:)

! define system for DP_Moment calculation ...
allocate( mask(a%atoms) , source = .true. )
mask = merge( mask , a%FMO , count(a%FMO) == I_zero )

! sum_i = (q_i * vec{r}_i) / sum_i q_i ...

allocate( Qi_Ri(a%atoms,3) , source = D_zero )

forall( j=1:3 , i=1:a%atoms , mask(i) ) Qi_Ri(i,j) = a%Nvalen(i) * a%coord(i,j)

total_valence = sum( a%Nvalen , mask )

forall(j=1:3) a%Center_of_Charge(j) = sum( Qi_Ri(:,j) , mask ) / total_valence

! atomic positions measured from the Center of Charge
allocate( R_vector(a%atoms,3) , source = D_zero )
forall( j=1:3 , i=1:a%atoms , mask(i) ) R_vector(i,j) = a%coord(i,j) - a%Center_of_Charge(j)

deallocate( Qi_Ri , mask )

end subroutine Center_of_Charge
!
!
!
!============================================
subroutine Build_DIPOLE_Matrix(system, basis)
!============================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)

! local variables
real*8  :: expa, expb, xab , yab , zab , Rab 
integer :: AtNo_a , AtNo_b
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

forall(i=1:3) DP_matrix_AO(:,:)%dp(i) = 0.d0

do ib = 1  , system%atoms
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

                qlm = 0.d0   ! check this !!!!!

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
end module DP_main_m

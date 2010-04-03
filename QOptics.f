module QOptics_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use Oscillator_m
    use FMO_m           , only : orbital
    use Multipole_Core

    public :: Evolution , RK_setup

    private

    type(R3_vector)    , allocatable  , save :: D(:,:) 
    type(transition)                  , save :: P
    complex*16         , allocatable  , save :: D_P_scalar(:,:)
    real*8             , allocatable  , save :: MO_erg(:)

 contains
!
!
!
!-------------------------------------------------------
subroutine Evolution( time , f , dfdt )
!-------------------------------------------------------
 implicit none
 real*8     , intent(in)    :: time
 complex*16 , intent(in)    :: f(:,:)
 complex*16 , intent(inout) :: dfdt(:,:)


 CALL Redfield_Equations( time , f , dfdt )



 end subroutine Evolution
!
!
!
!-----------------------------------------------
subroutine Redfield_Equations( time , f , dfdt )
!-----------------------------------------------
 implicit none
 real*8             , intent(in)    :: time
 complex*16         , intent(in)    :: f(:,:)
 complex*16         , intent(inout) :: dfdt(:,:)

! . local parameters ...
 complex*16         , parameter     :: one = (1.d0,0.d0) , zero = (0.d0,0.d0)

! . local variables ...
 integer                            :: i , j , xyz , dim_f
 complex*16         , allocatable   :: dfdt_tmp(:,:) , tmp(:,:)

 dim_f = size( f(:,1) )

! terms 1 & 3 ... OK

 CALL gemm( D_P_scalar , f , dfdt , 'N' , 'N' , one , zero )  

 CALL gemm( f , D_P_scalar , dfdt , 'N' , 'T' , one , one )  

 dfdt = -dfdt

! term 2 ... OK

 allocate( dfdt_tmp (dim_f,dim_f) )
 allocate( tmp      (dim_f,dim_f) )

 do xyz = 1 , 3

    dfdt_tmp = f

    tmp = dcmplx(D%DP(xyz))
    CALL trmm( tmp , dfdt_tmp , 'L' , 'L' , 'T' )

    tmp = dcmplx(P%matrix%DP(xyz))
    CALL trmm( tmp , dfdt_tmp , 'R' , 'L' )

    dfdt = dfdt + dfdt_tmp

 end do

! term 4 ... OK

 do xyz = 1 , 3

    dfdt_tmp = f

    tmp = dcmplx(P%matrix%DP(xyz))
    CALL trmm( tmp , dfdt_tmp , 'L' , 'L' , 'T' )

    tmp = dcmplx(D%DP(xyz))
    CALL trmm( tmp , dfdt_tmp , 'R' , 'L' )

    dfdt = dfdt + dfdt_tmp

 end do

 deallocate( dfdt_tmp , tmp )

 
 forall( j=1:dim_f , i=1:dim_f ) dfdt(i,j) = dfdt(i,j) * cdexp( zi * ( MO_erg(i)-MO_erg(j) ) * time / h_bar )

 forall( j=1:dim_f , i=1:dim_f ) dfdt(i,j) = dfdt(i,j) * 1.d1
 
include 'formats.h'

end subroutine Redfield_Equations
!
!
!
!---------------------------------------------------------------
subroutine RK_setup( system , basis , QM , FMO , rho , dim_rho )
!---------------------------------------------------------------
 type(structure)                    , intent(in)    :: system
 type(STO_basis)                    , intent(in)    :: basis(:)
 type(C_eigen)                      , intent(in)    :: QM
 type(C_eigen)                      , intent(in)    :: FMO
 complex*16         , allocatable   , intent(out)   :: rho(:,:)
 integer                            , intent(out)   :: dim_rho

! . local parameter ...
 real*8             , parameter     :: Gamma_ij = 1.89808d-6       ! <== (erg_ij^3 e^2)/( 6 pi epsilon_0 c^3 h_bar^4) ; unit = 1 / (eV^3 * ps * Angs^2) 

! . local variables ...
 integer                            :: i , j , xyz
 type(C3_vector)    , allocatable   :: D_P(:,:)

!-------------------------------------------------------- 
! . transition dipole matrices ... 

 P%flag = 'Redfield'

 CALL Transition_Dipole_Builder( system , basis , QM , P )

 dim_rho = size( P%matrix(:,1) )

 forall( i=1:dim_rho ) P%matrix(i,i)%DP = 0.d0

!-------------------------------------------------------- 
! . D[i,j] = E[i,j]^3 * Gamma[i,j] * <i|r|j> is lower trinagular 

 allocate( D(dim_rho,dim_rho) )

 forall( j=1:dim_rho , i=1:dim_rho ) D(i,j)%DP = 0.d0

 forall( j=1:dim_rho , i=1:dim_rho , i>j )

    D(i,j)%DP = (   (QM%erg(P%bra_PTR(i)) - QM%erg(P%ket_PTR(j)))   *   &
                    (QM%erg(P%bra_PTR(i)) - QM%erg(P%ket_PTR(j)))   *   &
                    (QM%erg(P%bra_PTR(i)) - QM%erg(P%ket_PTR(j)))   )   *   Gamma_ij    *   P % matrix(i,j) % DP
 end forall
!-------------------------------------------------------- 
! . D_P[i,j] = P*D^T

 allocate( D_P        (dim_rho,dim_rho) )
 allocate( D_P_scalar (dim_rho,dim_rho) )

 forall( j=1:dim_rho , i=1:dim_rho ) D_P(i,j)%DP = 0.d0

 forall( xyz=1:3 , j=1:dim_rho , i=2:dim_rho )

    D_P(i,j)%DP(xyz) = sum( P%matrix(i,1:i-1)%DP(xyz) * D(j,1:i-1)%DP(xyz) )

 end forall

 forall( j=1:dim_rho , i=1:dim_rho ) D_P_scalar(i,j) = sum( D_P(i,j)%DP(:) )

 deallocate(D_P)

!-------------------------------------------------------- 
! .  define the initial density matrix rho ...

 allocate( rho(dim_rho,dim_rho) )

! . f = rho(0)[i,j] = <i|L> <L|j> = ket[i] X bra[j] 

 forall( j=1:dim_rho , i=1:dim_rho )
    rho(i,j)    =   FMO%L( P%bra_PTR(j),orbital(1) ) * cdexp( + zi * QM%erg(P%bra_PTR(j)) * t_i / h_bar ) &
                *   FMO%R( P%ket_PTR(i),orbital(1) ) * cdexp( - zi * QM%erg(P%ket_PTR(i)) * t_i / h_bar )
 end forall

 trace = real( sum( (/( rho(i,i) , i=1,dim_rho )/) ) )

 Print 185 , 1.d0 - trace
!-------------------------------------------------------- 
! . save MO energies for future use in evolution ...

 allocate( MO_erg(dim_rho) )

 MO_erg = QM%erg(P%ket_PTR)


 include 'formats.h'

end subroutine RK_setup
!
!
!
!
end module QOptics_m

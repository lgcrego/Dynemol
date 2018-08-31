module QMMM_m
   
    use constants_m
    use parameters_m         , only : PBC
    use for_force            , only : pot_total 
    use MD_read_m            , only : atom , molecule , MM
    use polarizability_m     , only : Induced_DP
    use Semi_Empirical_Parms , only : chemical_element => atom

    public :: QMMM_FORCE

    private

    ! module parameters ...
    real*8  :: D_2_eAngs = 0.20819434d0
    real*8  :: beta_QQ   = 1.0
    real*8  :: mu_QQ     = 6.0
    real*8  :: beta_QP   = 0.8
    real*8  :: mu_QP     = 10.0

    ! module variables ...
    real*8 , allocatable :: dq(:)

contains
!
!
!
!==================================
 subroutine QMMM_FORCE( NetCharge )
!==================================
implicit none
real*8  , intent(in) :: NetCharge(:)

! local variables ...
real*8  :: FourVector(4) , U_Coul
integer :: i , j , ati , atj 

return

print*,"QMMM_FORCE"

if( .NOT. allocated(dq) ) allocate( dq(size(atom)) )
dq = NetCharge

forall( i=1:size(atom) ) atom(i) % fcoupling(:) = D_zero
U_Coul = D_zero

!================================================================================================
do ati = 1 , MM % N_of_atoms 
    do atj = ati+1 , MM % N_of_atoms

            if( atom(atj)% flex .OR. atom(ati)% flex ) then
        
            FourVector(1:4) = CoulFourVector( ati , atj )

            atom(ati)% fcoupling(1:3) = atom(ati)% fcoupling(1:3) + FourVector(1:3)
            atom(atj)% fcoupling(1:3) = atom(atj)% fcoupling(1:3) - FourVector(1:3)

            U_Coul = U_Coul + FourVector(4)

        end if
    end do
end do
!================================================================================================

! Append total force with Excited State Coulombic terms; force units = J/mts = Newtons ...
forall( i=1:MM % N_of_atoms ) atom(i)% ftotal(:) = atom(i)% ftotal(:) + atom(i)% fcoupling(:) * Angs_2_mts

pot_total = pot_total + U_Coul * mol*micro*factor3/MM%N_of_molecules 

end subroutine QMMM_FORCE
!
!
!
!=========================================================
 pure  function CoulFourVector( i , j ) result(FourVector)
!=========================================================
implicit none
integer , intent(in) :: i
integer , intent(in) :: j

! local variables ...
real*8  , dimension (4) :: FourVector
real*8  , dimension (3) :: rij , a , QiPj , QjPi , F_QQ , F_QP
real*8                  :: rijq , rijsq , Q_i, Q_j , QQ_ij , U_QQ , U_QP
real*8                  :: g , step , barrier

rij(:) = atom(i) % xyz(:) - atom(j) % xyz(:)
rij(:) = rij(:) - MM % box(:) * DINT( rij(:) * MM % ibox(:) )
rijq   = sum( rij(:) * rij(:) )
rijsq  = sqrt( rijq )

! charge/charge interaction ...
!================================================================================================
QQ_ij = dq(i)*dq(j) + atom(i)%charge*dq(j) + atom(j)%charge*dq(i)

! force ...
F_QQ = QQ_ij * rij(1:3) / (rijq * rijsq) 

! energy ...
U_QQ = QQ_ij / rijsq
!================================================================================================

! charge/induced-dipole interaction ...
!================================================================================================
Q_i = atom(i)%charge + dq(i)
! a = (p.r)r ...
a(1:3) = dot_product(Induced_DP(j,1:3),rij(1:3)) * rij(1:3)
! a = 3*(p.r)*r / ( |r|^2 ) - p ...
QiPj(1:3) = Q_i * (THREE * a(1:3) * ( D_ONE / rijq ) - Induced_DP(j,1:3) )

Q_j = atom(j)%charge + dq(j)
! a = (p.r)r ...
a(1:3) = dot_product(Induced_DP(i,1:3),rij(1:3)) * rij(1:3)
! a = 3*(p.r)*r / ( |r|^2 ) - p ...
QjPi(1:3) = Q_j * (THREE * a(1:3) * ( D_ONE / rijq ) - Induced_DP(i,1:3) )

! force ...
F_QP = (QiPj - QjPi) / (rijq * rijsq)  

! energy ...
a(:) = Q_i*Induced_DP(j,:) - Q_j*Induced_DP(i,:)
U_QP = dot_product( a(:) , rij(:) ) / ( rijq * rijsq ) 
!================================================================================================


! applying smooth cutoffs ...
!================================================================================================
! charge-charge ...
g       = exp(beta_QQ*(rijsq - mu_QQ))
step    = g / (g + D_one)
barrier = beta_QQ*step - beta_QQ*step*step  

U_QQ = U_QQ * step  
F_QQ = F_QQ * step - U_QQ * barrier * rij(1:3)/rijsq 

! charge-dipole ...
g       = exp(beta_QP*(rijsq - mu_QP))
step    = g / (g + D_one)
barrier = beta_QP*step - beta_QP*step*step  

U_QP = U_QP * step  
F_QP = F_QP * step - U_QP * barrier * rij(1:3)/rijsq 
!================================================================================================

FourVector(1:3) = Coulomb*F_QQ + Coulomb*D_2_eAngs*F_QP
FourVector(4)   = Coulomb*U_QQ + Coulomb*D_2_eAngs*U_QP

end function CoulFourVector
!
!
!
end module QMMM_m

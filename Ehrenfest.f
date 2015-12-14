! Program for computing Ehrenfest forces from Huckel Hamiltonian
module Ehrenfest_Builder

    use type_m
    use constants_m
    use parameters_m            , only  : verbose , n_part
    use MD_read_m               , only  : atom
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: EhrenfestForce 

    private

    !module variables ...
    real*8                :: delta = 1.d-8
    real*8  , allocatable :: grad_S(:,:)

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]

contains
!
!
!
!=============================================================================
 subroutine EhrenfestForce( system , basis , MO_bra , MO_ket , QM_el , QM_hl )
!=============================================================================
 implicit none
 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 complex*16                 , intent(in)    :: MO_bra(:,:)
 complex*16                 , intent(in)    :: MO_ket(:,:)
 type(R_eigen)              , intent(in)    :: QM_el
 type(R_eigen)   , optional , intent(in)    :: QM_hl

! local variables ... 
 integer            :: i 
 real*8             :: wp_occ(size(basis)) 

! local parameters ...
 real*8 , parameter :: eVAngs_2_Newton = 1.602176565d-9

forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero

do i = 1 , system% atoms

    ! electron contribution ...
    wp_occ(:) = MO_bra(:,1) * MO_ket(:,1)

    atom(i)% Ehrenfest = real(Ehrenfest( system, basis, wp_occ , QM_el , i )) * eVAngs_2_Newton 

    ! hole conntribution ...
    If( n_part > 1 ) then
        
        wp_occ(:) = MO_bra(:,2) * MO_ket(:,2)
        
        ! single or separate QM's ...
        If( present(QM_hl) ) then
            atom(i)% Ehrenfest = atom(i)% Ehrenfest - real(Ehrenfest( system, basis, wp_occ , QM_hl , i )) * eVAngs_2_Newton 
        else
            atom(i)% Ehrenfest = atom(i)% Ehrenfest - real(Ehrenfest( system, basis, wp_occ , QM_el , i )) * eVAngs_2_Newton 
        end If
      
    end If    
       
end do

deallocate( grad_S )

include 'formats.h'

end subroutine EhrenfestForce
!
!
!
!=====================================================================
 function Ehrenfest( system, basis, wp_occ , QM , site ) result(Force)
!=====================================================================
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
real*8           , intent(in)    :: wp_occ(:)
type(R_eigen)    , intent(in)    :: QM
integer          , intent(in)    :: site 

! local variables ...
real*8  :: xb , yb , zb
real*8  :: delta_b(3) 
integer :: ib , xyz
integer :: i , j , n

! local arrays ...
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) 
real*8                :: Force(3)

verbose = .false.
If( .NOT. allocated(grad_S) ) allocate( grad_S(size(basis),size(basis)) )

!force on atom site ...
ib = site 

! save coordinate ...
xb = system% coord (ib,1) 
yb = system% coord (ib,2) 
zb = system% coord (ib,3) 

do xyz = 1 , 3

       delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )

       system% coord (ib,1) = xb + delta_b(1)
       system% coord (ib,2) = yb + delta_b(2)
       system% coord (ib,3) = zb + delta_b(3)

       CALL Overlap_Matrix( system , basis , S_fwd )

       system% coord (ib,1) = xb - delta_b(1)
       system% coord (ib,2) = yb - delta_b(2)
       system% coord (ib,3) = zb - delta_b(3)

       CALL Overlap_Matrix( system , basis , S_bck )

       grad_S = (S_fwd - S_bck) / (TWO*delta) 

       Force(xyz) = C_zero

       do i = 1 , size(basis)
       do j = 1 , size(basis)

            do n = 1 , size(basis) 

                Force(xyz) = Force(xyz) - wp_occ(n) * ( Huckel_stuff(i,j,basis) - QM%erg(n) ) * grad_S(i,j) * QM%L(n,i) * QM%L(n,j)
       
            end do

       end do
       end do
            
end do 

! recover original system ...
system% coord (ib,1) = xb 
system% coord (ib,2) = yb 
system% coord (ib,3) = zb 

end function Ehrenfest
!
!
!
!
!======================================
 function Huckel_stuff( i , j , basis )
!======================================
implicit none
integer         , intent(in) :: i , j
type(STO_basis) , intent(in) :: basis(:)

!local variables ...
real*8 :: Huckel_stuff
real*8 :: k_eff , k_WH , c1 , c2 , c3

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

if (i == j) then

    huckel_stuff = basis(i)%IP

else

    c1 = basis(i)%IP - basis(j)%IP
    c2 = basis(i)%IP + basis(j)%IP

    c3 = (c1/c2)*(c1/c2)

    k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

    k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

    Huckel_stuff = k_eff * c2 * HALF

end if 

end function Huckel_stuff
!
!
end module Ehrenfest_Builder

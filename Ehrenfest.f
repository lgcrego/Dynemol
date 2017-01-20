! Program for computing Ehrenfest forces from Huckel Hamiltonian
module Ehrenfest_Builder

    use type_m
    use constants_m
    use parameters_m            , only  : driver , verbose , n_part , QMMM
    use Structure_Builder       , only  : Unit_Cell 
    use Overlap_Builder         , only  : Overlap_Matrix
    use Allocation_m            , only  : DeAllocate_Structures    

    public :: EhrenfestForce 

    private

    !module variables ...
    real*8                      :: delta = 1.d-8
    real*8      , allocatable   :: grad_S(:,:) , S(:,:)
    complex*16  , allocatable   :: rho_eh(:,:)

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]

contains
!
!
!
!=============================================================================
 subroutine EhrenfestForce( system , basis , MO_bra , MO_ket , QM_el , QM_hl )
!=============================================================================
 use MD_read_m              , only  : atom
 implicit none
 type(structure)            , intent(inout) :: system
 type(STO_basis)            , intent(in)    :: basis(:)
 complex*16                 , intent(in)    :: MO_bra(:,:)
 complex*16                 , intent(in)    :: MO_ket(:,:)
 type(R_eigen)              , intent(in)    :: QM_el
 type(R_eigen)   , optional , intent(in)    :: QM_hl

! local variables ... 
 integer :: i , j
 logical :: mask

! local parameters ...
 real*8 , parameter :: eVAngs_2_Newton = 1.602176565d-9

forall( i=1:system% atoms ) atom(i)% Ehrenfest(:) = D_zero

If( .NOT. allocated(rho_eh) ) allocate( rho_eh(size(basis),size(basis)) )

! preprocess overlap matrix for Pulay calculations ...
CALL Overlap_Matrix( system , basis )
CALL Overlap_Matrix( system , basis , S )

! build up electron-hole density matrix ...
forall( i=1:size(basis) , j=1:size(basis) ) rho_eh(i,j) = MO_ket(j,1)*MO_bra(i,1) - MO_ket(j,2)*MO_bra(i,2)

select case( driver )

    case( "slice_AO" )

        do i = 1 , system% atoms
            atom(i)% Ehrenfest = Ehrenfest_AO( system, basis, QM_el, i ) * eVAngs_2_Newton 
        end do

    case( "slice_ElHl" )

        ! electron and hole contributions ...
        do i = 1 , system% atoms
            atom(i)% Ehrenfest = Ehrenfest_ElHl( system, basis, MO_bra, MO_ket, QM_el, QM_hl, i ) * eVAngs_2_Newton 
        end do

end select

deallocate( S )
If( allocated(grad_S) ) deallocate( grad_S )

include 'formats.h'

end subroutine EhrenfestForce
!
!
!
!===============================================================
 function Ehrenfest_AO( system, basis, QM , site ) result(Force)
!===============================================================
use Semi_empirical_parms , only: atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
type(R_eigen)    , intent(in)    :: QM
integer          , intent(in)    :: site 

! local variables ...
integer :: i , j , m , n , xyz
integer :: k , ik , DOS_atom_k , BasisPointer_k 
real*8  :: Force_AD , Force_NAD , rho_nm

! local arrays ...
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) 
real*8                :: Force(3) , tmp_coord(3) , delta_b(3) 
real*8                :: X_ij

verbose = .false.
If( .NOT. allocated(grad_S) ) allocate( grad_S(size(basis),size(basis)) )

!force on atom site ...
k = site 
DOS_atom_k     =  atom( system% AtNo(k) )% DOS
BasisPointer_k =  system% BasisPointer(k) 

! save coordinate ...
tmp_coord = system% coord(k,:)

do xyz = 1 , 3

       delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )
 
       system% coord (k,:) = tmp_coord + delta_b
       CALL Overlap_Matrix( system , basis , S_fwd )

       system% coord (k,:) = tmp_coord - delta_b
       CALL Overlap_Matrix( system , basis , S_bck )
       grad_S = (S_fwd - S_bck) / (TWO*delta) 

!       grad_S = (S_fwd - S) / delta 

       Force_AD = D_zero

       ! adiabatic component of the Force ...
       !==============================================================================================
       !$omp parallel do schedule(dynamic,10) private(n,ik,i,j) default(shared) reduction(+:Force_AD)
       do n = 1 , size(basis) 
         do ik = 1 , DOS_atom_k
            i = BasisPointer_k + ik
            do j = 1 , size(basis)
                Force_AD = Force_AD - rho_eh(n,n) * ( Huckel_stuff(i,j,basis) - QM%erg(n) ) * grad_S(i,j) * QM%L(n,i) * QM%L(n,j)
            end do
         end do
       end do
       !$end parallel do
       !==============================================================================================

       Force_NAD = D_zero

       ! non-adiabatic component of the Force ...
       !==============================================================================================
       !$omp parallel do schedule(dynamic,10) private(n,m,rho_nm,ik,i,j,X_ij) default(shared) reduction(+:Force_NAD)
       do n = 1   , size(basis)
       do m = n+1 , size(basis)

         rho_nm = real( rho_eh(m,n) ) 

         do ik = 1 , DOS_atom_k
            i = BasisPointer_k + ik
            do j = 1 , size(basis)

                X_ij = Huckel_stuff(i,j,basis) 

                Force_NAD = Force_NAD - rho_nm * grad_S(j,i)                                                &
                          * ((X_ij-QM%erg(m))*QM%L(m,j)*QM%L(n,i) + (X_ij-QM%erg(n))*QM%L(m,i)*QM%L(n,j))

            end do
         end do

       end do
       end do
       !$end parallel do
       !==============================================================================================

       Force(xyz) = two * (Force_AD + Force_NAD)

end do 

! recover original system ...
system% coord (k,:) = tmp_coord

end function Ehrenfest_AO
!
!
!
!=========================================================================================
 function Ehrenfest_ElHl( system, basis, MO_bra, MO_ket,QM_el, QM_hl, site ) result(Force)
!=========================================================================================
use Semi_empirical_parms , only: atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
complex*16       , intent(in)    :: MO_bra(:,:)
complex*16       , intent(in)    :: MO_ket(:,:)
type(R_eigen)    , intent(in)    :: QM_el
type(R_eigen)    , intent(in)    :: QM_hl
integer          , intent(in)    :: site 

! local variables ...
integer :: i , j , m , n , xyz
integer :: k , ik , DOS_atom_k , BasisPointer_k 
real*8  :: Force_AD , Force_NAD , rho_nm_el , rho_nm_hl

! local arrays ...
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) 
real*8                :: Force(3) , tmp_coord(3) , delta_b(3) , wp_occ(size(basis),2)

verbose = .false.
If( .NOT. allocated(grad_S) ) allocate( grad_S(size(basis),size(basis)) )

!force on atom site ...
k = site 
DOS_atom_k     =  atom( system% AtNo(k) )% DOS
BasisPointer_k =  system% BasisPointer(k) 

! some preprocessing ...
forall(j=1:2) wp_occ(:,j) = MO_bra(:,j)*MO_ket(:,j)

! save coordinate ...
tmp_coord = system% coord(k,:)

do xyz = 1 , 3

       delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )
    
       system% coord (k,:) = tmp_coord + delta_b
       CALL Overlap_Matrix( system , basis , S_fwd )

       system% coord (k,:) = tmp_coord - delta_b
       CALL Overlap_Matrix( system , basis , S_bck )

       grad_S = (S_fwd - S_bck) / (TWO*delta) 

       Force_AD = D_zero

       ! adiabatic component of the Force ...
       !==============================================================================================
       !$omp parallel do schedule(dynamic,10) private(n,ik,i,j) default(shared) reduction(+:Force_AD)
       do n = 1 , size(basis) 
         do ik = 1 , DOS_atom_k
            i = BasisPointer_k + ik
            do j = 1 , size(basis)
                Force_AD = Force_AD                                                                                             &
                         - wp_occ(n,1) * ( Huckel_stuff(i,j,basis) - QM_el%erg(n) ) * grad_S(i,j) * QM_el%L(n,i) * QM_el%L(n,j) &
                         + wp_occ(n,2) * ( Huckel_stuff(i,j,basis) - QM_hl%erg(n) ) * grad_S(i,j) * QM_hl%L(n,i) * QM_hl%L(n,j) 
            end do
         end do
       end do
       !$end parallel do
       !==============================================================================================

       Force_NAD = D_zero

       ! non-adiabatic component of the Force ...
       !==============================================================================================
       !$omp parallel do schedule(dynamic,10) private(n,m,rho_nm_el,rho_nm_hl,ik,i,j) default(shared) reduction(+:Force_NAD)
       do n = 1   , size(basis)
       do m = n+1 , size(basis)
         rho_nm_el = two * real(MO_bra(m,1)*MO_ket(n,1)) 
         rho_nm_hl = two * real(MO_bra(m,2)*MO_ket(n,2)) 

         do ik = 1 , DOS_atom_k
            i = BasisPointer_k + ik
            do j = 1 , size(basis)

                Force_NAD = Force_NAD                                                                                                     &
                          - rho_nm_el * ( QM_el%erg(n)*QM_el%L(m,j)*QM_el%L(n,i) + QM_el%erg(m)*QM_el%L(m,i)*QM_el%L(n,j) ) * grad_S(i,j) &
                          + rho_nm_hl * ( QM_hl%erg(n)*QM_hl%L(m,j)*QM_hl%L(n,i) + QM_hl%erg(m)*QM_hl%L(m,i)*QM_hl%L(n,j) ) * grad_S(i,j) 

            end do
         end do

       end do
       end do
       !$end parallel do
       !==============================================================================================

       Force(xyz) = two*Force_AD + Force_NAD

end do 

! recover original system ...
system% coord (k,:) = tmp_coord

end function Ehrenfest_ElHl
!
!
!
!
!===========================================
 pure function Huckel_stuff( i , j , basis )
!===========================================
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

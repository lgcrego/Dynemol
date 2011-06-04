module Coulomb_m

    use constants_m         
    use type_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack

    public  :: Build_Coulomb_potential , wormhole_to_Coulomb

    private

    ! module variables ...
    complex*16 , allocatable , save  :: AO_bra(:,:) , AO_ket(:,:)

contains
!
!
!
!=======================================================================================
 subroutine Build_Coulomb_potential( S_matrix , basis , V_coul , V_coul_El , V_coul_Hl )
!=======================================================================================
implicit none
real*8                          , intent(in)  :: S_matrix  (:,:)
type(STO_basis)                 , intent(in)  :: basis     (:)
complex*16      , allocatable   , intent(out) :: V_coul    (:,:) 
complex*16      , allocatable   , intent(out) :: V_coul_El (:) 
complex*16      , allocatable   , intent(out) :: V_coul_Hl (:)

! local variables ...
complex*16                  :: coeff_El_ij , coeff_El_kl , coeff_HL_ij , coeff_Hl_kl
real*8                      :: mid_xyz_ij(3) , mid_xyz_kl(3) , dif_xyz_kl(3) , delta_xyz(3) , distance_ij_kl
real*8                      :: factor
integer                     :: i , j , k , l , basis_size

basis_size = size( basis(:) )

! if there is no electron-hole pair save S_matrix and leave the subroutine ...
if( .NOT. allocated(AO_bra) ) then
    
    allocate( V_coul_El (basis_size)            , source = C_zero)
    allocate( V_coul_Hl (basis_size)            , source = C_zero)
    allocate( V_coul    (basis_size,basis_size) , source = C_zero)

    return

end if

allocate( V_coul_El (basis_size)            , source = C_zero )
allocate( V_coul_Hl (basis_size)            , source = C_zero )
allocate( V_coul    (basis_size,basis_size) , source = C_zero )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  two, three, and four-center Coulomb potential matrix elements for Electrons and Holes ...

do j = 1 , basis_size
do i = 1 , j

    if( abs(S_matrix(i,j)) > high_prec ) then

        ! coefficients for V_Electron potential ...
        coeff_El_ij = AO_bra(i,1) * AO_ket(j,1) * S_matrix(i,j)

        ! coefficients for V_Hole potential ...
        coeff_Hl_ij = AO_bra(i,2) * AO_ket(j,2) * S_matrix(i,j)
        
        mid_xyz_ij(1) = (basis(i)%x + basis(j)%x)/two
        mid_xyz_ij(2) = (basis(i)%y + basis(j)%y)/two
        mid_xyz_ij(3) = (basis(i)%z + basis(j)%z)/two

        do l = 1 , basis_size
        do k = 1 , l

            factor = merge( two , D_one , k/=l )

            ! condition not yet defined ...
            if( (basis(i)%atom==basis(j)%atom) .AND. (basis(k)%atom==basis(l)%atom) .AND. (basis(i)%atom==basis(k)%atom) ) exit

            if( abs(S_matrix(k,l)) > high_prec ) then

                ! coefficients for V_Electron potential ...
                coeff_Hl_kl = AO_bra(k,2) * AO_ket(l,2) * S_matrix(k,l)

                ! coefficients for V_Hole potential ...
                coeff_El_kl = AO_bra(k,1) * AO_ket(l,1) * S_matrix(k,l)

                mid_xyz_kl(1) = (basis(k)%x + basis(l)%x)/two
                mid_xyz_kl(2) = (basis(k)%y + basis(l)%y)/two
                mid_xyz_kl(3) = (basis(k)%z + basis(l)%z)/two

                delta_xyz(:)  = mid_xyz_ij(:) - mid_xyz_kl(:) 

                distance_ij_kl = sqrt( sum(delta_xyz(:)*delta_xyz(:)) ) / a_Bohr

                ! electron-hole "exchange" mechanism ...
                if( distance_ij_kl < mid_prec ) then
                    dif_xyz_kl(1) = basis(k)%x - basis(l)%x
                    dif_xyz_kl(2) = basis(k)%y - basis(l)%y
                    dif_xyz_kl(3) = basis(k)%z - basis(l)%z
                    
                    distance_ij_kl = sqrt( sum(dif_xyz_kl(:)*dif_xyz_kl(:)) ) / a_Bohr
                end if

                if( i == j ) then
                    ! diagonal elements of V_Electron ...
                    V_coul_El(i) = - coeff_El_ij * factor*Real(coeff_Hl_kl) / distance_ij_kl

                    ! diagonal elements of V_Hole ...
                    V_coul_Hl(i) = - coeff_Hl_ij * factor*Real(coeff_El_kl) / distance_ij_kl
                else
                    ! UPPER TRIANGLE: OFF-diagonal elements for V_Electron potential ...
                    V_coul(i,j) =  - coeff_El_ij * factor*Real(coeff_Hl_kl) / distance_ij_kl

                    ! LOWER TRIANGLE: OFF-diagonal elements for V_Hole potential ...
                    V_coul(j,i) =  - coeff_Hl_ij * factor*Real(coeff_El_kl) / distance_ij_kl
                end if

            end if      !<== [|S(k,l)| > high_prec]
        end do
        end do
    end if              !<== [|S(i,j)| > high_prec]
end do
end do

V_coul_El (:)   = V_coul_El (:)  * Hartree_2_eV * 1.d2
V_coul_Hl (:)   = V_coul_Hl (:)  * Hartree_2_eV * 1.d2
V_coul    (:,:) = V_coul    (:,:)* Hartree_2_eV * 1.d2

end subroutine Build_Coulomb_potential
!
!
!
!===========================================
 subroutine wormhole_to_Coulomb( bra , ket )
!===========================================
implicit none
complex*16  , intent(in) :: bra(:,:)
complex*16  , intent(in) :: ket(:,:)

! local variables ...
integer :: row , column

! save AO_bra and AO_ket in Coulomb_m ...

if( .NOT. allocated(AO_bra) ) then

    row     = size( bra(:,1) )
    column  = size( bra(1,:) )

    allocate( AO_bra(row,column) , source=bra )
    allocate( AO_ket(row,column) , source=ket )

else

    AO_bra = bra
    AO_ket = ket

end if

end subroutine wormhole_to_Coulomb
!
!
!
!
end module Coulomb_m

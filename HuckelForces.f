! Program for computing Hellman-Feynman-Pulay forces form Huckel Hamiltonian ...
module HuckelForces_m

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only  : PBC , verbose
    use Overlap_Builder         , only  : Overlap_Matrix
    use PBC_m                   , only  : Generate_Periodic_Structure
    use Semi_Empirical_Parms    , only  : ChemAtom => atom
    use Allocation_m            , only  : DeAllocate_Structures    
    use Ehrenfest_Builder       , only  : RotationOverlap, solap, rotar, util_overlap, dlmn

    public :: HuckelForces

    private

    !module variables ...
    real*8  :: delta = 1.d-8

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]

contains
!
!
!
!==============================================
 subroutine HuckelForces( system , basis , QM )
!==============================================
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(in)    :: basis(:)
 type(R_eigen)   , intent(in)    :: QM

! local variables ... 
 integer                         :: i , i1 , i2 , n , n_MO , Fermi_level , method
 real*8          , allocatable   :: bra(:), ket(:), Force(:,:), force_atom(:,:)
 type(structure)                 :: pbc_system
 type(STO_basis) , allocatable   :: pbc_basis(:)

 CALL util_overlap     

 n_MO = size(QM%erg)
 allocate( bra  ( size(basis)              ) )
 allocate( ket  ( size(basis)              ) )
 allocate( Force( 3*system% atoms , 0:n_MO ) , source = D_zero )

 ! if no PBC: pbc_system = system ...
 CALL Generate_Periodic_Structure( system, pbc_system, pbc_basis ) 

 write(*,'(/a)') ' Choose the method : '
 write(*,'(/a)') ' (1) = Hellman-Feynman-Pulay '
 write(*,'(/a)') ' (2) = Numerical derivative of PES '
 read (*,'(I)') method

 select case( method )

 case( 1 )
 !=========================================================================
 ! Hellman-Feynman-Pulay ...

 do n = 1 , n_MO

    ! bra = ket ...
    bra = QM%L(n,:)
    ket = bra 
    do i = 1 , system% atoms

        i1 = (i-1)*3 + 1
        i2 = (i-1)*3 + 3

        Force( i1:i2 ,n ) = Hellman_Feynman_Pulay( system, basis, bra, ket, QM%erg(n), i )

    end do
 end do

 case( 2 )
 !=========================================================================
 ! numerical derivative of the PES ...

 verbose = .false.
 allocate( force_atom( n_MO , 3 ) )

 do i = 1 , system% atoms

        force_atom = grad_E( system, basis, i )

        i1 = (i-1)*3 + 1
        i2 = (i-1)*3 + 3

        forall( n=1:n_MO ) Force( i1:i2 , n ) = force_atom( n , : )

 end do
 !=========================================================================

 end select

! center of mass force ...
 do n = 1 , n_MO
    Print 200, n , sum( Force(:,n) )
 end do
 
 ! Force(:,0) = total force at ground state ...
 Fermi_level = system% N_of_electrons / 2
 forall( i=1:size(Force(:,0)) ) Force(i,0) = sum( Force(i,1:Fermi_level) )

 do i = 1 , 3*system%atoms
    write(30,*) i , Force(i,0)
 end do

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
OPEN( unit=3 , file='EHM_forces.nmd' , status='unknown' )

write(3,*) "EHM Force Analysis"

write( 3 , '(A6 ,1000A3)'   ) "names "         , (system % Symbol(i)   , i = 1 , system% atoms)
write( 3 , '(A9 ,1000A4)'   ) "resnames "      , (system % residue(i)  , i = 1 , system% atoms)
write( 3 , '(A6 ,1000A2)'   ) "chids "         , [("A"                 , i = 1 , system% atoms)]             
write( 3 , '(A7 ,1000I4)'   ) "resids "        , (system % nr(i)       , i = 1 , system% atoms)
write( 3 , '(A6 ,1000A2)'   ) "betas "         , [("0"                 , i = 1 , system% atoms)]             
write( 3 , '(A12,3000F8.4)' ) "coordinates "   , (system % coord(i,:)  , i = 1 , system% atoms)

do n = 0 , n_MO
    write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , n , Force(:,n) 
end do

close(3)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 CALL Deallocate_Structures(pbc_system)
 deallocate(pbc_basis)

 include 'formats.h'

end subroutine HuckelForces
!
!
!
!=================================================================================
function Hellman_Feynman_Pulay( system, basis, bra, ket, erg, site ) result(Force)
!=================================================================================
use util_m , factorial => fact
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in) :: basis(:)
real*8           , intent(in) :: bra(:)
real*8           , intent(in) :: ket(:)
real*8           , intent(in) :: erg
integer          , intent(in) :: site 

! local variables ...
real*8  :: xb , yb , zb 
real*8  :: delta_b(3) 
integer :: ib , xyz
integer :: i , j 

! local arrays ...
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) , grad_S(:,:)
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

       Force(xyz) = D_zero
       do i = 1 , size(basis)
       do j = 1 , size(basis)

           Force(xyz) = Force(xyz) - ( Huckel_stuff(i,j,basis) - erg ) * grad_S(i,j) * bra(i) * ket(j)

       end do
       end do

end do 

! recover original system ...
system% coord (ib,1) = xb 
system% coord (ib,2) = yb 
system% coord (ib,3) = zb 

end function Hellman_Feynman_Pulay
!
!
!
!
!====================================================
 function grad_E( system, basis, site ) result(Force)
!====================================================
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
integer          , intent(in)    :: site 

! local variables ...
integer :: ib , xyz 
real*8  :: xb , yb , zb
real*8  :: delta_b(3) 
real*8  :: erg_fwd(size(basis)) , erg_bck(size(basis)) , Force(size(basis),3)

Force = D_zero

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

            CALL LocalEigenSystem( system , basis , erg_fwd )

            system% coord (ib,1) = xb - delta_b(1)
            system% coord (ib,2) = yb - delta_b(2)
            system% coord (ib,3) = zb - delta_b(3)

            CALL LocalEigenSystem( system , basis , erg_bck )

            Force(:,xyz) = - (erg_fwd(:) - erg_bck(:)) / (TWO*delta) 

    end do 

! recover original system ...
system% coord (ib,1) = xb 
system% coord (ib,2) = yb 
system% coord (ib,3) = zb 

end function grad_E
!
!
!
!
!===================================================
 subroutine LocalEigenSystem( system , basis , erg )
!===================================================
implicit none
type(structure)  , intent(in)  :: system
type(STO_basis)  , intent(in)  :: basis(:)
real*8           , intent(out) :: erg(:)

! local variables ...
real*8  , ALLOCATABLE :: h(:,:) 
real*8  , ALLOCATABLE :: S_matrix(:,:)
integer               :: i , j , info

CALL Overlap_Matrix( system , basis , S_matrix )

allocate( h(size(basis),size(basis)) )

do j = 1 , size(basis)
    do i = j, size(basis)

        h(i,j) = Huckel_stuff( i , j , basis ) * S_matrix(i,j)

    end do
end do    

CALL SYGVD( h , S_matrix , erg , 1 , 'V' , 'L' , info )

If ( info /= 0 ) write(*,*) 'info = ',info,' in LocalEigenSystem '

deallocate( h )

end subroutine LocalEigenSystem
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
!
!
end module HuckelForces_m

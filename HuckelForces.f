! Program for computing Hellman-Feynman-Pulay forces form Huckel Hamiltonian ...
module HuckelForces_m

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use constants_m
    use parameters_m            , only : verbose , EnvField_
    use Overlap_Builder         , only : Overlap_Matrix
    use Allocation_m            , only : DeAllocate_Structures    
    use Hamiltonians            , only : X_ij , even_more_extended_Huckel

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

 n_MO = size(QM%erg)
 allocate( bra  ( size(basis)              ) )
 allocate( ket  ( size(basis)              ) )
 allocate( Force( 3*system% atoms , 0:n_MO ) , source = D_zero )

 write(*,'(/a)') ' Choose the method : '
 write(*,'(/a)') ' (1) = Hellman-Feynman-Pulay '
 write(*,'(/a)') ' (2) = Numerical derivative of PES '
 read (*,'(I)') method

 select case( method )

 case( 1 )
 !=========================================================================
 ! Hellman-Feynman-Pulay ...

 CALL Overlap_Matrix( system , basis )

 do n = 1 , n_MO

    ! bra = ket ...
    bra = QM%L(n,:)
    ket = bra 
    do i = 1 , system% atoms

        If( system% QMMM(i) /= "QM" ) cycle

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

 do i = 1 , 3*system% atoms
    write(30,*) i , Force(i,0)
 end do

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
OPEN( unit=3 , file='HFP_forces.nmd' , status='unknown' )

write(3,*) "HFP Force Analysis"

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

include 'formats.h'

end subroutine HuckelForces
!
!
!
!=================================================================================
function Hellman_Feynman_Pulay( system, basis, bra, ket, erg, site ) result(Force)
!=================================================================================
use Semi_Empirical_parms , only : atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in) :: basis(:)
real*8           , intent(in) :: bra(:)
real*8           , intent(in) :: ket(:)
real*8           , intent(in) :: erg
integer          , intent(in) :: site 

! local variables ...
integer :: i , j , k , ik , xyz

! local arrays ...
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) , grad_S(:,:)
real*8                :: Force(3) , tmp_coord(3) , delta_b(3) 

verbose = .false.
If( .NOT. allocated(grad_S) ) allocate( grad_S(size(basis),size(basis)) )

!force on atom site ...
k = site 

! save coordinate ...
tmp_coord = system% coord(k,:)

do xyz = 1 , 3

       delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )

       system% coord (k,:) = tmp_coord + delta_b
       CALL Overlap_Matrix( system , basis , S_fwd , purpose='Pulay')

       system% coord (k,:) = tmp_coord - delta_b
       CALL Overlap_Matrix( system , basis , S_bck , purpose='Pulay')

       grad_S = (S_fwd - S_bck) / (TWO*delta) 

       Force(xyz) = D_zero

       do ik = 1 , atom( system% AtNo(k) )% DOS  
            i = system% BasisPointer(k) + ik
            do j = 1 , size(basis)

                Force(xyz) = Force(xyz) - ( X_ij(i,j,basis) - erg ) * grad_S(i,j) * bra(i) * ket(j)

            end do
       end do

end do 

Force = two * Force

! recover original system ...
system% coord (k,:) = tmp_coord

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
integer :: k , xyz 
real*8  :: delta_b(3) , tmp_coord(3)
real*8  :: erg_fwd(size(basis)) , erg_bck(size(basis)) , Force(size(basis),3)

Force = D_zero

!force on atom site ...
k = site 

! save coordinate ...
tmp_coord = system% coord(k,:)

    do xyz = 1 , 3

            delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )

            system% coord (k,:) = tmp_coord + delta_b

            CALL LocalEigenSystem( system , basis , erg_fwd )

            system% coord (k,:) = tmp_coord - delta_b

            CALL LocalEigenSystem( system , basis , erg_bck )

            Force(:,xyz) = - (erg_fwd(:) - erg_bck(:)) / (TWO*delta) 

    end do 

! recover original system ...
system% coord (k,:) = tmp_coord

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
integer               :: info

CALL Overlap_Matrix( system , basis , S_matrix )

allocate( h(size(basis),size(basis)) )

If( EnvField_ ) then
    h(:,:) = even_more_extended_Huckel( system , basis , S_matrix )
else
    h(:,:) = Build_Huckel( basis , S_matrix )
end If

CALL SYGVD( h , S_matrix , erg , 1 , 'V' , 'L' , info )

If ( info /= 0 ) write(*,*) 'info = ',info,' in LocalEigenSystem '

deallocate( h )

end subroutine LocalEigenSystem
!
!
!
!===================================================
 function Build_Huckel( basis , S_matrix ) result(h)
!===================================================
implicit none
type(STO_basis) , intent(in)    :: basis(:)
real*8          , intent(in)    :: S_matrix(:,:)

! local variables ... 
integer              :: i , j , N
real*8 , allocatable :: h(:,:)

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

N = size(basis)
ALLOCATE( h(N,N) , source = D_zero )

do j = 1 , N
  do i = j , N

        h(i,j) = X_ij( i , j , basis ) * S_matrix(i,j)

    end do
end do

end function Build_Huckel
!
!
!
!
end module HuckelForces_m

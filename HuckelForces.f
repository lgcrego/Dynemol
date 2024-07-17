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

    public :: HuckelForces , Force

    private

    !module variables ...
    real*8  :: delta = 1.d-8

    !module parameters ...
    integer , parameter   :: xyz_key(3) = [1,2,3]
    real*8  , allocatable :: Force(:,:) , aux(:,:)

contains
!
!
!
!==================================================================
 subroutine HuckelForces( system , basis , QM , approach , eh_PES )
!==================================================================
 implicit none
 type(structure) , intent(inout)         :: system
 type(STO_basis) , intent(in)            :: basis(:)
 type(R_eigen)   , intent(in)            :: QM
 character(*)    , optional , intent(in) :: approach
 integer         , optional , intent(in) :: eh_PES(2)

! local variables ... 
 integer              :: i , i1 , i2 , n , n_MO , Fermi_level , method , n_eh(2)
 real*8 , allocatable :: bra(:), ket(:), force_atom(:,:)

 allocate( aux( size(basis) , 9 ) )

 n_MO = size(QM%erg)
 allocate( bra ( size(basis) ) )
 allocate( ket ( size(basis) ) )

 if( present(approach) ) &
 then
     method = 4  ! <== XS geometry optimization
 else
     write(*,'(/a)') ' Choose the method : '
     write(*,'(/a)') ' (1) = Hellman-Feynman-Pulay '
     write(*,'(/a)') ' (2) = Numerical derivative of PES '
     write(*,'(/a)') ' (3) = el-hl pair force '
     write(*,'(/a)',advance='no') '>>>   '
     read (*,'(I)') method
 end if

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 select case( method )

        case( 1 )
                 !=========================================================================
                 ! Hellman-Feynman-Pulay ...

                 allocate( Force( 3*system% atoms , 0:n_MO ) , source = D_zero )

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
                 allocate( Force( 3*system% atoms , 0:n_MO ) , source = D_zero )
                 allocate( force_atom( n_MO , 3 ) )

                 do i = 1 , system% atoms

                        force_atom = grad_E( system, basis, i )

                        i1 = (i-1)*3 + 1
                        i2 = (i-1)*3 + 3

                        forall( n=1:n_MO ) Force( i1:i2 , n ) = force_atom( n , : )

                 end do

        case( 3 )
                 !=========================================================================
                 ! el-hl pair force ...
                 
                 write(*,'(/a)') '> enter MOs associated with the electron and hole:'
                 write(*,'(a)',advance='no') 'n_el = '
                 read (*,*) n_eh(1)
                 write(*,'(a)',advance='no') 'n_hl = '
                 read (*,*) n_eh(2)

                 ! resetting array dimensions ...
                 if( allocated(Force) ) deallocate( Force )
                 allocate( Force( 3*system% atoms , 0:2 ) , source = D_zero )

                 CALL Overlap_Matrix( system , basis )

                 do n = 1 , 2

                      bra = QM%L(n_eh(n),:)
                      ket = bra 

                      do i = 1 , system% atoms

                          If( system% QMMM(i) /= "QM" .OR. system% flex(i) == .false. ) cycle

                          i1 = (i-1)*3 + 1
                          i2 = (i-1)*3 + 3

                          Force( i1:i2 ,n ) = Hellman_Feynman_Pulay( system, basis, bra, ket, QM%erg(n_eh(n)), i )

                      end do
                 end do

        case( 4 )
                 !=========================================================================
                 ! el-hl pair force , used for XS geometry optimization ...
                 
                 If( present(eh_PES) ) n_eh = eh_PES
                 
                 if( allocated(Force) ) deallocate( Force )
                 allocate( Force( 3*system% atoms , 2 ) , source = D_zero )

                 select case( approach )

                 case( "HFP" , "hfp" )
                     ! HFP approach ...
                     CALL Overlap_Matrix( system , basis )
                     do n = 1 , 2
                          bra = QM%L(n_eh(n),:)
                          ket = bra 
                          do i = 1 , system% atoms
                              If( system% QMMM(i) /= "QM" .OR. system% flex(i) == .false. ) cycle
                              i1 = (i-1)*3 + 1
                              i2 = (i-1)*3 + 3
                              Force( i1:i2 ,n ) = Hellman_Feynman_Pulay( system, basis, bra, ket, QM%erg(n_eh(n)), i )
                          end do
                          end do
                     !---------------------------------------------------------------

                 case( "FDM" , "fdm" )
                     ! Finite Difference approach ...
                     allocate( force_atom( 2 , 3 ) )
                     do i = 1 , system% atoms
                            If( system% QMMM(i) /= "QM" .OR. system% flex(i) == .false. ) cycle
                            force_atom = grad_E( system, basis, i )
                            i1 = (i-1)*3 + 1
                            i2 = (i-1)*3 + 3
                            forall( n=1:2 ) Force( i1:i2 , n ) = force_atom( n , : )
                     end do
                     !---------------------------------------------------------------

                 case default
                     CALL warning("halting: force method must be HFP or FDM, check XS_ERG_class.f")
                     stop
                     !---------------------------------------------------------------

                 end select

                 ! mission accomplished ...
                 deallocate(aux)
                 RETURN

 end select
 
 deallocate(aux)

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
select case (method)
       case(1:2)
                ! center of mass force (must be zero) ...
                do n = 1 , n_MO
                   Print 200, n , sum( Force(:,n) )
                end do

                ! Force(:,0) = total force at ground state ...
                Fermi_level = system% N_of_electrons / 2
                forall( i=1:size(Force(:,0)) ) Force(i,0) = sum( Force(i,1:Fermi_level) )
 
       case(3)
                ! center of mass force (must be zero) ...
                do n = 1 , 2
                   Print 200, n_eh(n) , sum( Force(:,n) )
                end do

                ! Force(:,0) = el-hl pair force ...
                forall( i=1:size(Force(:,0)) ) Force(i,0) = Force(i,1) - Force(i,2)
end select
 

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
OPEN( unit=3 , file='ancillary.trunk/HFP_forces.nmd' , status='unknown' )

write(3,*) "HFP Force Analysis"

write( 3 , '(A6 ,1000A3)'   ) "names "         , (system % Symbol(i)   , i = 1 , system% atoms)
write( 3 , '(A9 ,1000A4)'   ) "resnames "      , (system % residue(i)  , i = 1 , system% atoms)
write( 3 , '(A6 ,1000A2)'   ) "chids "         , [("A"                 , i = 1 , system% atoms)]             
write( 3 , '(A7 ,1000I4)'   ) "resids "        , (system % nr(i)       , i = 1 , system% atoms)
write( 3 , '(A6 ,1000A2)'   ) "betas "         , [("0"                 , i = 1 , system% atoms)]             
write( 3 , '(A12,3000F8.4)' ) "coordinates "   , (system % coord(i,:)  , i = 1 , system% atoms)

select case (method)
       case(1:2)
               do n = 0 , n_MO
                   write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , n , Force(:,n) 
               end do
 
       case(3)
                write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , 0 , Force(:,0) 
                do n = 2 , 1 , -1
                   write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , n_eh(n) , Force(:,n) 
                end do
end select

close(3)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
print*, ""
print*, '>> saving HFP_forces.nmd in directory ancillary.trunk/ <<'

include 'formats.h'

end subroutine HuckelForces
!
!
!
!=====================================================================================
function Hellman_Feynman_Pulay( system, basis, bra, ket, erg, site ) result(Force_HFP)
!=====================================================================================
use Semi_Empirical_parms , only : atom
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in) :: basis(:)
real*8           , intent(in) :: bra(:)
real*8           , intent(in) :: ket(:)
real*8           , intent(in) :: erg
integer          , intent(in) :: site 

! local variables ...
integer :: i , j , k , ik , xyz , DOSk , BPk

! local arrays ...
real*8  , allocatable :: S_fwd(:,:) , S_bck(:,:) , grad_S(:,:)
real*8                :: Force_HFP(3) , tmp_coord(3) , delta_b(3) 

verbose = .false.
If( .NOT. allocated(grad_S) ) allocate( grad_S(size(basis),size(basis)) )
grad_S = D_zero
Force_HFP = D_zero

!force on atom site ...
k    = site 
DOSk = atom( system% AtNo(k) )% DOS
BPk  = system% BasisPointer(k)

! save coordinate ...
tmp_coord = system% coord(k,:)

do xyz = 1 , 3

       delta_b = delta * merge(D_one , D_zero , xyz_key == xyz )

       system% coord (k,:) = tmp_coord + delta_b
       CALL Overlap_Matrix( system , basis , S_fwd )

       system% coord (k,:) = tmp_coord - delta_b
       CALL Overlap_Matrix( system , basis , S_bck )

       grad_S = (S_fwd - S_bck) / (TWO*delta) 

       do concurrent ( ik=1:DOSk , j=1:size(basis) )
            i = BPk + ik
            aux(j,ik) =  - ( X_ij(i,j,basis) - erg ) * grad_S(i,j) * bra(i) * ket(j)
       end do

       Force_HFP(xyz) = sum( aux(:,1:DOSk) )

end do 

Force_HFP = two * Force_HFP

! recover original system ...
system% coord (k,:) = tmp_coord

end function Hellman_Feynman_Pulay
!
!
!
!
!========================================================
 function grad_E( system, basis, site ) result(Force_xyz)
!========================================================
implicit none
type(structure)  , intent(inout) :: system
type(STO_basis)  , intent(in)    :: basis(:)
integer          , intent(in)    :: site 

! local variables ...
integer :: k , xyz 
real*8  :: delta_b(3) , tmp_coord(3)
real*8  :: erg_fwd(size(basis)) , erg_bck(size(basis)) , Force_xyz(size(basis),3)

Force_xyz = D_zero

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

            Force_xyz(:,xyz) = - (erg_fwd(:) - erg_bck(:)) / (TWO*delta) 

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

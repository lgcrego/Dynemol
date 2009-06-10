 module projectors

    use type_m
    use constants_m
    use Structure_Builder
    use EHT_parameters , the_chemical_atom => atom

 contains
! 
!
!
!-----------------------------------------------
 subroutine projector( FMO_L, FMO_R, zR, wv_FMO)
!-----------------------------------------------
 complex*16 , ALLOCATABLE , intent(out) :: FMO_L(:,:) , FMO_R(:,:)
 complex*16 , ALLOCATABLE , intent(in)  :: zR(:,:)
 real*8     , ALLOCATABLE , intent(in)  :: wv_FMO(:,:)

! local variables 
 integer :: ALL_size , FMO_size , cluster_size , i , j
 real*8  :: check

 ALL_size = size( zR(:,1) )                     ! <== basis size of the entire system
 FMO_size = size( wv_FMO(1,:) )                 ! <== basis size of the FMO system
 cluster_size  = ALL_size - FMO_size            ! <== basis size of the cluster 

 Allocate( FMO_L (ALL_size,FMO_size) )
 Allocate( FMO_R (ALL_size,FMO_size) )

!-----------------------------------------------------------------------------
!    writes the isolated molecule wavefunction |k> in the MO basis 
!    the isolated orbitals are stored in the ROWS of wv_FMO

 FMO_L = ( 0.d0 , 0.d0 )
 FMO_R = ( 0.d0 , 0.d0 )

 forall( j=1:FMO_size, i=1:ALL_size )

    FMO_L(i,j) = sum( wv_FMO(j,:) * zR(cluster_size+1:ALL_size,i) )

 end forall    

 FMO_R = FMO_L

 check = dreal( sum( FMO_L(1:ALL_size,:)*FMO_R(1:ALL_size,:) ) )

 if( dabs(check-FMO_size) < low_prec ) then
     print*, '>> projection done <<'
 else
     Print 58 , check 
     pause '---> problem in projector <---'
 end if

!-----------------------------------------------------------------------------

 include 'formats.h'

 end subroutine
!
!
!
!-----------------------------------------------
 function pop_Slater(basis,za,zb,fragment)
!-----------------------------------------------
 type(STO_basis)  , intent(in) :: basis(:)
 complex*16       , intent(in) :: za(:,:) , zb(:,:)
 character(len=2) , optional   :: fragment

! local variables
 real*8                                :: pop_Slater
 complex*16 , dimension(size(za(1,:))) :: pop
 integer                               :: n , i 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!        projection of | k(t) >  onto the atom k_atom 

 pop(:) = (0.d0,0.d0)

 if( present(fragment) ) then
    forall(n=1:n_part , i=1:size(basis) , basis(i)%fragment == fragment)  pop(n) = pop(n) + za(i,n)*zb(i,n)
 else
    forall(n=1:n_part , i=1:size(basis))  pop(n) = pop(n) + za(i,n)*zb(i,n)
 end if

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 pop_Slater = real(sum(pop))

 end function

 end module projectors

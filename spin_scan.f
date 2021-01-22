module Spin_Align_m

    use type_m
    use constants_m
    use parameters_m       , only : verbose , SOC , B_field
    use QCModel_Huckel     , only : EigenSystem
    use Structure_Builder  , only : Generate_Structure ,   &
                                    Basis_Builder 

    public :: ScanSpinOrientation

    private 

    ! module variables ...
    logical              :: done = .false.
    real*8 , allocatable :: atomQ(:) 

contains
!
!
!
!================================================
 subroutine ScanSpinOrientation( system , basis )
!================================================
implicit none
type(structure)              , intent(inout) :: system
type(STO_basis), allocatable , intent(inout) :: basis(:)

!local parameters ...
integer              :: theta_steps = 360
character(len=1)     :: xyz(3) = ['x','y','z']

!local variable ...
integer              :: i , j , k , MO , ati 
real*8               :: degree , theta
real*8               :: CC(3) 
real*8 , allocatable :: Sz(:,:) , erg(:,:)
type(C_eigen)        :: UNI
logical              :: exist

verbose = NO

inquire(file="frames.pdb", EXIST=exist)
if( exist ) call systemQQ( "rm frames.pdb" )

MO = 11

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! ROTATE MOLECULE

allocate( Sz(theta_steps,3) , erg(theta_steps,3) )

degree = (PI/180.d0)       

do j = 1 , 3

       CALL preprocess( system , basis , MO )

       do k = 0 , theta_steps-1
       
              theta = k * degree 
       
              CALL rotate( system%coord , angle=degree , axis = xyz(j) )
              
              CALL dump_geometry( system , theta , k , axis=xyz(j) )
       
              CALL Basis_Builder( system , basis )
           
              CALL EigenSystem( system , basis , UNI )
              
              Sz (k+1,j) = real( sum( UNI%L(MO,:)*UNI%R(:,MO)*dcmplx(basis(:)%S) ) ) 
              erg(k+1,j) = UNI%erg(MO)  
              
              ! atomic net-charge ...
              atomQ = d_zero
              do ati = 1 , system%atoms
                 atomQ(ati) = atomQ(ati) + abs( sum( UNI%L(MO,:)*UNI%R(:,MO) , basis(:)%atom == ati ) )
                 end do
       
              forall( i=1:3 ) CC(i) = sum(atomQ*system%coord(:,i))
              
              forall( i=1:3 ) system%coord(:,i) = system%coord(:,i)-CC(i)
       
              end do
              end do

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
open (10, file='spin.trunk/Sz.dat' , status='unknown')
open (11, file='spin.trunk/erg.dat', status='unknown')
do k = 0 , theta_steps-1
     theta = k * degree 
     write(10,'(4F10.6)') theta*(180.0/PI) , (Sz (k+1,i),i=1,3)
     write(11,'(4F10.5)') theta*(180.0/PI) , (erg(k+1,i),i=1,3)  
     end do
close(10)
close(11)

end subroutine ScanSpinOrientation
!
!
!
!=====================================
 subroutine rotate( a , angle , axis )
!=====================================
implicit none
real*8           , intent(inout) :: a(:,:)
real*8           , intent(in)    :: angle
character(len=1) , intent(in)    :: axis

!local variables
integer              :: i , j , atoms
real*8               :: R(3,3)
real*8 , allocatable :: b(:,:)

allocate(b , source=a)

atoms = size(a(:,1))

select case (axis)
       case('x') 
           !------------------------
           R      =  0.d0
           R(1,1) =  1.d0
           R(2,2) =  dcos(angle)
           R(2,3) = -dsin(angle)
           R(3,2) =  dsin(angle)
           R(3,3) =  dcos(angle)
           !------------------------
       case('y') 
           !------------------------
           R      =  0.d0
           R(2,2) =  1.d0
           R(1,1) =  dcos(angle)
           R(1,3) =  dsin(angle)
           R(3,1) = -dsin(angle)
           R(3,3) =  dcos(angle)
           !------------------------
       case('z') 
           !------------------------
           R      =  0.d0 
           R(3,3) =  1.d0
           R(1,1) =  dcos(angle)
           R(1,2) = -dsin(angle)
           R(2,1) =  dsin(angle)
           R(2,2) =  dcos(angle)
           !------------------------
end select     

do i = 1,atoms 
do j = 1,3  
   a(i,j) = sum( R(j,:)*b(i,:) ) 
   end do
   end do

end subroutine rotate
!
!
!
!=================================================
 subroutine preprocess( system , basis , MO_indx )
! translates system to Center of Charge (CC) ...
!=================================================
implicit none
type(structure) , intent(inout) :: system
type(STO_basis) , intent(in)    :: basis(:)
integer         , intent(in)    :: MO_indx

!local variable ...
integer                     :: MO , xyz , ati
real*8                      :: CC(3)
real*8 , allocatable , save :: coord_in_store(:,:)
type(C_eigen)               :: UNI

if( .not. done ) then

    allocate( atomQ(system%atoms) )
    allocate( coord_in_store , source = system%coord )

    !----- SEED MAGNETIC FIELD -------
    B_field = [ 0.0d0 , 0.0d0 , 0.5d0 ] 

    done = yes

else

    system%coord = coord_in_store

end if

CALL EigenSystem( system , basis , UNI )

MO = MO_indx

! atomic net-charge ...
do ati = 1 , system%atoms
   atomQ(ati) = atomQ(ati) + abs( sum( UNI%L(MO,:)*UNI%R(:,MO) , basis(:)%atom == ati ) )
   end do

forall( xyz=1:3 ) CC(xyz) = sum(atomQ*system%coord(:,xyz))
forall( xyz=1:3 ) system%coord(:,xyz) = system%coord(:,xyz)-CC(xyz)

deallocate(UNI%L,UNI%R,UNI%erg)
                                      
end  subroutine 
!
!
!
!========================================
 subroutine dump_geometry(a,theta,k,axis)
!========================================
implicit none
type(structure)  , intent(in) :: a
real*8           , intent(in) :: theta
integer          , intent(in) :: k
character(len=1) , intent(in) :: axis

! local variables ...
integer :: i , j

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
OPEN (unit=4, file='frames.pdb', status='unknown', access='append')

if(k==0) write(4,6) 'rotating molecule'
write(4,7) 'rotating molecule , axis = ', axis, theta*(180.0/PI)
write(4,1) 'CRYST1' , a%T_xyz(1) , a%T_xyz(2) , a%T_xyz(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
write(4,993) 'MODEL' ,  k

do i = 1 , a%atoms

            write(4,2)  'HETATM'                        ,  &    ! <== non-standard atom
                        i                               ,  &    ! <== global number
                        a%Symbol(i)                     ,  &    ! <== atom type
                        ' '                             ,  &    ! <== alternate location indicator
                        a%residue(i)                    ,  &    ! <== residue name
                        ' '                             ,  &    ! <== chain identifier
                        a%nr(i)                         ,  &    ! <== residue sequence number
                        ' '                             ,  &    ! <== code for insertion of residues
                        ( a%coord(i,j) , j=1,3 )        ,  &    ! <== xyz coordinates
                        1.00                            ,  &    ! <== occupancy
                        0.00                            ,  &    ! <== temperature factor
                        ' '                             ,  &    ! <== segment identifier
                        ' '                             ,  &    ! <== here only for tabulation purposes
                        a%Symbol(i)                             ! <== chemical element symbol

end do

write(4,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , a%atoms , 0 , a%atoms , 0
write(4,8) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a20)
7 FORMAT(a27,a3,F9.3)
8 FORMAT(a3)
993 FORMAT(a5,i8)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

end subroutine dump_geometry
!
!
!
end module Spin_Align_m

module Topology_routines

use types_m
use constants_m
use Read_Parms  , only: Atomic_Mass

public :: connect , dump_topol , InputIntegers , get_topology

private

! module variables ...
integer, allocatable :: InputIntegers(:,:) 
integer, allocatable :: bond_list(:,:) , angle_list(:,:) , dihedral_list(:,:)
logical, allocatable :: bond_matrix(:,:) , angle_matrix(:,:,:) , dihedral_matrix(:,:,:,:)
logical              :: done = .false.

contains
!
!
!=========================
 subroutine connect( sys )
!=========================
implicit none
type(universe)                  , intent(inout)   :: sys 

! local varibles ...
character(len=1)                :: yn
character(len=2)                :: S1 , S2
integer                         :: i , j
real*8                          :: cutoff
logical                         :: flag

CALL system( "clear" )

write(*,'(/a)',advance='no') ">>> Connect Specific Bonds ? (y/n) "
read (*,'(a)') yn

if ( yn == "y" ) then

    allocate( sys % topol (sys%N_of_atoms,sys%N_of_atoms) , source = .false. )

    do

        write(*,'(1x,3/a)') "Choose chemical elements whose bonds are to be connected (@ to exit) : "
        read*, S1
        If( S1 == "@" ) exit
        read*, S2

        write(*,'(1x,a)') "cut-off distance for the bond: "
        read*, cutoff

        do j = 1 , sys%N_of_atoms 
            do i = j+1 , sys%N_of_atoms 

                flag = (sys%atom(i)%Symbol==S1 .AND. sys%atom(j)%Symbol==S2)       &
                    .OR.                                                           &  
                       (sys%atom(i)%Symbol==S2 .AND. sys%atom(j)%Symbol==S1) 

                If( flag .AND. sqrt(sum((sys%atom(i)%xyz - sys%atom(j)%xyz)**2)) < cutoff ) then
        
                    sys % topol(i,j) = .true. 
                    sys % topol(j,i) = .true.

                end If

            end do
        end do

    end do

end If

CALL system( "clear" )
    
end subroutine connect

!
!
!========================================
 subroutine dump_topol( sys , file_unit )
!========================================
implicit none
type(universe) , intent(in) :: sys 
integer        , intent(in) :: file_unit

! local varibles ...
integer :: i , j

do j = 1 , sys%N_of_atoms

    if( any(sys%topol(:,j) ) ) then

        write(file_unit,'(/A6,I5)', advance='no') "CONECT" , j 

        do i = 1 , sys % N_of_atoms
            If( sys % topol(i,j) ) write(file_unit,'(I5)', advance='no') i
        end do

    end if

end do
    
write(file_unit,'(/)')
    
end subroutine dump_topol
!
!
!
!=========================================
subroutine get_topology( sys , file_type )
!=========================================
implicit none 
type(universe) , intent(inout) :: sys
integer        , intent(in)    :: file_type

! local variables ...
logical :: TorF
integer :: i, j, n, total_bonds, total_angs, total_diheds 
integer :: mol_conect

!----------------------------------------------
!         generate topology files 
!----------------------------------------------

select case( file_type )
       case(1)
           OPEN(unit=10,file='seed.itp',status='unknown')
       case(2)
           OPEN(unit=10,file='seed.psf',status='unknown')
       end select

! get the Atomic_Masses ...
sys%atom%mass = Atomic_Mass(sys%atom%AtNo)

do n = 1 , maxval(sys%atom%nresid)
   
       !----------------------------------------
       ! heading 
       !----------------------------------------
       write(10,101) "[ moleculetype ]"
       write(10,*) , "SPECIES", 3
       write(10,*) 
       
       write(10,102) "[ atoms ]"
       do i = 1 , sys%N_of_atoms

           write(10,5) i                    ,  &  ! <== serial number within the residue
                       sys%atom(i)%Symbol   ,  &  ! <== force field descriptor of atom type
                       sys%atom(i)%nresid   ,  &  ! <== residue identifier
                       sys%atom(i)%resid    ,  &  ! <== residue name
                       sys%atom(i)%MMSymbol ,  &  ! <== atom type
                       sys%atom(i)%nrcg     ,  &  ! <== charge group
                       sys%atom(i)%charge   ,  &  ! <== charge of atom type       
                       sys%atom(i)%mass           ! <== mass of chemical element 
       end do
end do
write(10,*)
write(10,*)

!----------------------------------------
! start topology connections
!----------------------------------------
mol_conect = 0
do i = 1 , sys%total_conect
    j = size( pack( InputIntegers(i,:) , InputIntegers(i,:) /= 0 ) )
    mol_conect = j + mol_conect
end do

if( mol_conect == 0 ) then
    Write(*,*)
    Write(*,*) "======================= W A R N I N G ======================="
    Write(*,*) "No CONECT in the pdb file; cannot generate bonds, angles, etc"
    Write(*,*) "============================================================="
    Write(*,*)
endif

!----------------------------------------
! Assign CONECT to a logical bond matrix ...
!----------------------------------------
CALL get_bond_matrix( sys , total_bonds ) 

!--------------------------------------------------------
! Assign bond matrix to a logical i,k,j angle matrix ...
!--------------------------------------------------------
if( mol_conect > 3 ) then 
    CALL get_angle_matrix( sys , total_angs ) 
end if

!-----------------------------------------------------------
! Assign angle matrix to a logical i--l dihedral matrix ...
!-----------------------------------------------------------
if( mol_conect > 4 ) then 
    CALL get_dihedral_matrix( sys , total_diheds ) 
end if

CALL generate_topology_list( sys%atom )

select case( file_type )
    case(1)
        CALL write_seed_itp
    case(2)
        CALL write_seed_psf
    end select

TorF = Checking_Topology( bond_list , angle_list , dihedral_list )
If( TorF ) then
    Print*, "error detected in Topology , check Topology.log"
    stop
    End If

deallocate( bond_matrix, angle_matrix, dihedral_matrix )

5   FORMAT(i6,a8,i8,a7,a8,i8,F9.4,F9.4)
101 FORMAT(a16)
102 FORMAT(a9)

end subroutine get_topology
!
!
!
!
!
!
!==============================================
subroutine get_bond_matrix( sys , total_bonds )
!==============================================
implicit none 
type(universe) , intent(in)  :: sys
integer        , intent(out) :: total_bonds

! local variables ...
integer :: i , j , k , l , N_atoms 
logical :: flag1 , flag2

N_atoms = sys%N_of_atoms

allocate( bond_matrix(N_atoms,N_atoms), source=.false. )

do i = 1, N_atoms
  do j = i, N_atoms
    do k = 1, sys%Total_conect
      flag1 = ( i == InputIntegers(k,1) )
      if( flag1 .eqv. .true.) then
        do l = 2, 5
          flag2 = ( j == InputIntegers(k,l) )
          if( flag2 .eqv. .true.) bond_matrix(i,j) = .true.
        end do
      end if
    end do
    bond_matrix(j,i) = bond_matrix(i,j)
  end do
end do

total_bonds = (size( pack( bond_matrix(:,:), bond_matrix(:,:) .eqv. .true. )))/ 2

end subroutine get_bond_matrix
!
!
!
!==============================================
subroutine get_angle_matrix( sys , total_angs )
!==============================================
implicit none 
type(universe) , intent(in)  :: sys
integer        , intent(out) :: total_angs

! local variables ...
integer :: i , j , k , m , N_atoms 
logical :: flag1 , flag2

N_atoms = sys%N_of_atoms

allocate( angle_matrix(N_atoms,N_atoms,N_atoms), source=.false. )

m = 0
do i = 1, N_atoms
  do j = i, N_atoms
    do k = 1, N_atoms
      flag1 = ( bond_matrix(i,k) .or. bond_matrix(k,i) ) .eqv. .true.
      flag2 = ( bond_matrix(j,k) .or. bond_matrix(k,j) ) .eqv. .true.
      if( flag1 .and. flag2 ) then
        angle_matrix(i,k,j) = .true.
        angle_matrix(j,k,i) = angle_matrix(i,k,j)
        angle_matrix(i,k,i) = .false.
        m = m + 1
      end if
    end do
  end do
end do
total_angs = m 

end subroutine get_angle_matrix
!
!
!
!===================================================
subroutine get_dihedral_matrix( sys , total_diheds )
!===================================================
implicit none 
type(universe) , intent(in)  :: sys
integer        , intent(out) :: total_diheds

! local variables ...
integer :: i , j , k , l , N_atoms
logical :: flag1 , flag2

N_atoms = sys%N_of_atoms

allocate( dihedral_matrix(N_atoms,N_atoms,N_atoms,N_atoms),source=.false. )

total_diheds=0
do i = 1, N_atoms
  do j = 1, N_atoms
    do k = 1, N_atoms
      if( angle_matrix(i,j,k) .eqv. .true. ) then
        do l = 1, N_atoms
          flag1 = bond_matrix(l,k) .or. bond_matrix(k,l)
          flag2 = bond_matrix(l,i) .or. bond_matrix(i,l)
          if( flag1 ) then 
            dihedral_matrix(i,j,k,l) = .true.
            dihedral_matrix(i,j,k,j) = .false.
            dihedral_matrix(l,k,j,i) = dihedral_matrix(i,j,k,l) 
          end if
          if( flag2 ) then 
            dihedral_matrix(l,i,j,k) = .true.
            dihedral_matrix(j,i,j,k) = .false.
            dihedral_matrix(k,j,i,l) = dihedral_matrix(l,i,j,k)
          end if
          if(dihedral_matrix(i,j,k,l) .eqv. .true.) total_diheds=total_diheds+1
        end do
      end if
    end do
  end do
end do
total_diheds = total_diheds/2 

end subroutine get_dihedral_matrix
!
!
!
!=======================================
subroutine generate_topology_list( aux )
!=======================================
implicit none 
type(atomic) , intent(in) :: aux(:)

! local variables ...
integer :: i , j , k , l , m , N_atoms , ati , atj , atk , atl , max_bonds
integer , allocatable :: tmp_list(:,:)

N_atoms = size(aux)

max_bonds = N_atoms*(N_atoms-1)/2

!-----------------------------------------------------------
! Writing stuff ...  
!-----------------------------------------------------------

allocate(tmp_list( max_bonds , 2 ))
m = 0
do i = 1, N_atoms - 1
  do j = i, N_atoms
    ati = aux(i) % my_intra_id
    atj = aux(j) % my_intra_id
    if( bond_matrix(ati,atj)  .eqv. .true. ) then
             m = m + 1
             tmp_list(m,1) = i
             tmp_list(m,2) = j
             endif
  end do
end do
allocate(bond_list , source=tmp_list(1:m,:))
deallocate(tmp_list)

allocate(tmp_list( max_bonds , 3 ))
m = 0
do i = 1, N_atoms
  do j = i, N_atoms
    do k = 1, N_atoms
      ati = aux(i) % my_intra_id
      atj = aux(j) % my_intra_id
      atk = aux(k) % my_intra_id
      if( angle_matrix(ati,atk,atj) .eqv. .true. ) then
                m = m + 1
                tmp_list(m,1) = i
                tmp_list(m,2) = k
                tmp_list(m,3) = j
                endif
    end do
  end do
end do
allocate(angle_list , source=tmp_list(1:m,:))
deallocate(tmp_list)

allocate(tmp_list( max_bonds , 4 ))
m = 0
do i = 1 , N_atoms
  do j = 1, N_atoms
    do k = 1, N_atoms
      do l = 1, N_atoms
        if( i < l ) then
          ati = aux(i) % my_intra_id
          atj = aux(j) % my_intra_id
          atk = aux(k) % my_intra_id
          atl = aux(l) % my_intra_id
          if( dihedral_matrix(ati,atj,atk,atl) .eqv. .true. ) then
                       m = m + 1
                       tmp_list(m,1) = i
                       tmp_list(m,2) = j
                       tmp_list(m,3) = k
                       tmp_list(m,4) = l
                       endif
        end if
      end do
    end do
  end do
end do
allocate(dihedral_list , source=tmp_list(1:m,:))
deallocate(tmp_list)

end subroutine generate_topology_list
!
!
!
!========================
subroutine write_seed_itp
!========================
implicit none 

! local variables ...
integer :: j , k , Nbonds , Nangs , Ndiheds

!-----------------------------------------------------------
! Writing stuff ...  
!-----------------------------------------------------------

Nbonds = size(bond_list(:,1))
write(10,105) "[ bonds ]"
do k = 1 , Nbonds
     write(10,102) (bond_list(k,j) , j=1,2) , 1
end do

write(10,*) " "
write(10,*) " "
Nangs = size(angle_list(:,1))
write(10,106) "[ angles ]"
do k = 1 , Nangs
     write(10,103) (angle_list(k,j) , j=1,3) , 1
end do

write(10,*) " "
write(10,*) " "
Ndiheds = size(dihedral_list(:,1))
write(10,107) "[ dihedrals ]"
do k = 1 , Ndiheds
     write(10,104) (dihedral_list(k,j) , j=1,4) ,3
end do 

close(10)

102 format(3I4)
103 format(4I4)
104 format(5I4)
105 FORMAT(a9)
106 FORMAT(a10)
107 FORMAT(a13)

end subroutine write_seed_itp
!
!
!
!========================
subroutine write_seed_psf
!========================
implicit none 

! local variables ...
integer :: j , k , n , ioerr , Nbonds , Nangs , Ndiheds

!-----------------------------------------------------------
! Writing stuff ...  
!-----------------------------------------------------------

Nbonds = size(bond_list(:,1))
write(10,105) Nbonds , "   !NBOND: bonds"
do k = 1 , ceiling( Nbonds / four ) - 1
  write(10 ,100, iostat=ioerr )  ( ( bond_list((k-1)*4+n,j) , j=1,2 ) , n=1,4 )
end do
write(10 ,100, iostat=ioerr )  ( ( bond_list((k-1)*4+n,j) , j=1,2 ) , n=1,merge(4,mod(NBonds,4),mod(NBonds,4)==0) )

write(10,*) " "
write(10,*) " "
Nangs = size(angle_list(:,1))
write(10,106) Nangs , "   !NTHETA: angles"
do k = 1 , ceiling(Nangs/three)-1
  write(10 ,101, iostat=ioerr )  ( ( angle_list((k-1)*3+n,j) , j=1,3 ) , n=1,3 )
end do
write(10 ,101, iostat=ioerr )  ( ( angle_list((k-1)*3+n,j) , j=1,3 ) , n=1,merge(3,mod(NAngs,3),mod(NAngs,3)==0) )


write(10,*) " "
write(10,*) " "
Ndiheds = size(dihedral_list(:,1))
write(10,107) Ndiheds , "   !NPHI: dihedrals"
do k = 1 , ceiling(Ndiheds/two)-1
   write(10 ,102, iostat=ioerr )  ( ( dihedral_list((k-1)*2+n,j) , j=1,4 ) , n=1,2 )
end do 
write(10 ,102, iostat=ioerr )  ( ( dihedral_list((k-1)*2+n,j) , j=1,4 ) , n=1,merge(2,mod(Ndiheds,2),mod(Ndiheds,2)==0) )

close(10)

100 FORMAT(t10,I4,t16,I4,t30,I4,t36,I4,t50,I4,t56,I4,t70,I4,t76,I4)
101 FORMAT(t10,I4,t16,I4,t22,I4,t35,I4,t41,I4,t47,I4,t60,I4,t66,I4,t72,I4)
102 FORMAT(t10,I4,t16,I4,t22,I4,t28,I4,t41,I4,t47,I4,t53,I4,t58,I4)
105 FORMAT(t10,I4,a16)
106 FORMAT(t10,I4,a18)
107 FORMAT(t10,I4,a19)

end subroutine write_seed_psf
!
!
!
!================================================================
 function Checking_Topology( bonds , angs , diheds ) result(TorF)
!================================================================
implicit none
integer , intent(in) :: bonds (:,:)
integer , intent(in) :: angs  (:,:)
integer , intent(in) :: diheds(:,:)
logical              :: TorF
 
! local variables ... 
integer               :: i , x , y , z
integer               :: Nbonds , Nangs , Ndiheds , KeyLeft , KeyRight
integer , allocatable :: BondKeys(:) , AngKeys(:)
logical               :: flag

Nbonds  =  size(bonds (:,1)) 
Nangs   =  size(angs  (:,1))
Ndiheds =  size(diheds(:,1))

! checking bonds topology ...
allocate( BondKeys(Nbonds) )
do i = 1 , Nbonds

     x = bonds(i,1)  ;  y = bonds(i,2) 
     BondKeys(i) = PairingFunction( x , y , verbose = .true. ) 

end do

! checking angs topology ...
do i = 1 , Nangs

     flag = .false.

     x = angs(i,1)  ;  y = angs(i,2) 
     KeyLeft = PairingFunction(x,y) 
     If( .not. any(KeyLeft == BondKeys) ) call error_message(i,angs,flag,instance="ang")

     x = angs(i,2)  ;  y = angs(i,3) 
     KeyRight = PairingFunction(x,y) 
     If( .not. any(KeyRight == BondKeys) ) call error_message(i,angs,flag,instance="ang")

     If( KeyLeft == KeyRight ) call error_message(i,angs,flag,instance="ang")

end do

! checking diheds topology ...
allocate( AngKeys(Nangs) )
do i = 1 , Nangs

     x = angs(i,1)  ;  y = angs(i,2)   ;  z = angs(i,3) 
     AngKeys(i) = CantorPairing( x , y , z ) 

end do

do i = 1 , Ndiheds

     flag = .false.

     x = diheds(i,1)  ;  y = diheds(i,2)   ;  z = diheds(i,3) 
     KeyLeft = CantorPairing( x , y , z ) 
     If( .not. any(KeyLeft == AngKeys) ) call error_message(i,diheds,flag,instance="dihed")

     x = diheds(i,2)  ;  y = diheds(i,3)   ;  z = diheds(i,4) 
     KeyRight = CantorPairing( x , y , z ) 
     If( .not. any(KeyRight == AngKeys) ) call error_message(i,diheds,flag,instance="dihed")

end do

! prepare to leave ...
if( done ) then  
    TorF = .true.     ! <==  error detected
    close(10)
else
    TorF = .false.    ! <==  NO error detected
end If

end function Checking_Topology
!
!
!
!
!========================================
 function CantorPairing(i,j,k) result(R)
! 3-tupling Cantor Function ...
! f(i,j,k) = f(k,j,i)
!========================================
implicit none
integer            , intent(in) :: i,j,k

! local variables ... 
integer :: R , L , a , b

! Symmetric Pairing for (i,k)-tuple ...
a = max(i,k)  ;  b = min(i,k)
L = a*(a+1)/2 + b 

! Cantor pairing with the center pairing ...
R = (L+j)*(L+j+1)/2 + L 

end function CantorPairing
!
!
!
!===============================================
 function PairingFunction(i,j,verbose) result(k)
!===============================================
implicit none
integer            , intent(in) :: i,j
logical , optional , intent(in) :: verbose

! local variables ... 
integer :: k , a , b

If( (i == j) .and. present(verbose) ) then
    Print 232, i , j
    stop
end If

! the symmetric pairing satisfies f(i,j)=f(j,i) ...

! Symmetric Cantor Pairing ...
!k = (i+j)*(i+j+1)/2 + (i*j) 

! Symmetric Pairing ...
a = max(i,j)  ;  b = min(i,j)
k = a*(a+1)/2 + b 

232 format(/,1x,'>>>  Degenerate Pairing Function in Topology file.....: ',I4,I4)

end function PairingFunction
!
!
!
!==================================================
 subroutine error_message(i , a , flag , instance ) 
!==================================================
implicit none
integer          , intent(in) :: i
integer          , intent(in) :: a(:,:)
logical          , intent(inout) :: flag
character(len=*) , intent(in) :: instance

If( .not. done ) open (10, file='Topology.log', status='unknown')
done = .true.

If( flag == .true. ) return
select case (instance)

       case("ang")
       write(10,231) a(i,1) , a(i,2) , a(i,3) 

       case("dihed")
       write(10,233) a(i,1) , a(i,2) , a(i,3)  , a(i,4) 

end select
flag = .true.

231 format(/,1x,'>>> Error detected in Toplogy file .....: Angle (',I4,',',I4,',',I4,')' )
233 format(/,1x,'>>> Error detected in Toplogy file .....: Dihedral (',I4,',',I4,',',I4,',',I4,')' )

end subroutine error_message
!
!
!

end module Topology_routines

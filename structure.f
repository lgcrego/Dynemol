 module Structure_Builder

    use type_m
    use Babel_m                   , only : Read_from_XYZ , Read_from_Poscar , Read_from_PDB , Read_PDB , Read_VASP  &
                                         , Identify_Fragments , trj , System_Characteristics
    use Allocation_m              , only : Allocate_Structures
    use Semi_Empirical_Parms      , only : atom

    type(structure)               , public  :: Unit_Cell , Extended_Cell 
    type(STO_basis) , allocatable , public  :: ExCell_basis(:)

    public :: Read_Structure , Generate_Structure , Basis_Builder 

    private

 contains
!
!
!
!=========================
subroutine Read_Structure
!=========================

select case( file_type )

    case( "structure" )

        select case( file_format )
            case( "xyz" )
                CALL Read_from_XYZ      ( Unit_Cell ) 
            case( "pdb" )
                CALL Read_from_PDB      ( Unit_Cell )
            case( "vasp" )
                CALL Read_from_POSCAR   ( Unit_Cell )
            end select

    case( "trajectory" )
    
        select case( file_format ) 
            case( "pdb" ) 
                CALL Read_PDB   ( trj ) 
            case( "vasp" )
                CALL Read_VASP  ( trj )
            end select

end select

Print 70, System_Characteristics

include 'formats.h'

end subroutine Read_Structure
!
!
!
!======================================
 subroutine Generate_Structure( frame )
!======================================
integer , intent(in) :: frame

! local variables ...
integer :: copy

!----------------------------------------------------------
! GENERATES   THE   EXTENDED-STRUCTURE (REAL,not periodic)
!----------------------------------------------------------

If( .NOT. allocated(Extended_Cell%coord) ) CALL Allocate_Structures( (2*nnx+1)*(2*nny+1)*unit_cell%atoms , extended_cell )

 k = 0
 copy = 0

 DO ix = -nnx , nnx 
 DO iy = -nny , nny


     If( (ix /= 0) .OR. (iy /= 0) ) THEN 

        copy = copy + 1
        FORALL( n=1:unit_cell%atoms )

            extended_cell % coord    (k+n,1) =  unit_cell % coord    (n,1) + ix * unit_cell%T_xyz(1)
            extended_cell % coord    (k+n,2) =  unit_cell % coord    (n,2) + iy * unit_cell%T_xyz(2)
            extended_cell % coord    (k+n,3) =  unit_cell % coord    (n,3) 
            extended_cell % AtNo     (k+n)   =  unit_cell % AtNo     (n)
            extended_cell % k_WH     (k+n)   =  unit_cell % k_WH     (n)
            extended_cell % fragment (k+n)   =  unit_cell % fragment (n)
            extended_cell % Symbol   (k+n)   =  unit_cell % Symbol   (n)
            extended_cell % MMSymbol (k+n)   =  unit_cell % MMSymbol (n)
            extended_cell % residue  (k+n)   =  unit_cell % residue  (n)
            extended_cell % copy_No  (k+n)   =  copy
        
        END FORALL

        k = k + unit_cell%atoms

     END IF

 END DO
 END DO

 FORALL( n=1:unit_cell%atoms )     ! <== the DONOR CELL is at the end (extended_cell%copy_No = 0)

    extended_cell % coord    (k+n,1:3)  =  unit_cell % coord    (n,1:3)
    extended_cell % AtNo     (k+n)      =  unit_cell % AtNo     (n)
    extended_cell % k_WH     (k+n)      =  unit_cell % k_WH     (n)
    extended_cell % fragment (k+n)      =  unit_cell % fragment (n)
    extended_cell % symbol   (k+n)      =  unit_cell % Symbol   (n)
    extended_cell % MMSymbol (k+n)      =  unit_cell % MMSymbol (n)
    extended_cell % residue  (k+n)      =  unit_cell % residue  (n)
    extended_cell % copy_No  (k+n)      =  0

 END FORALL

 ! . define the DONOR fragment 
 where( (extended_cell%fragment == 'M') .AND. (extended_cell%copy_No == 0) ) extended_cell%fragment = 'D' 
 CALL Identify_Fragments( Extended_Cell )    

 extended_cell%T_xyz(1) = (2*nnx+1)*unit_cell%T_xyz(1)
 extended_cell%T_xyz(2) = (2*nny+1)*unit_cell%T_xyz(2)
 extended_cell%T_xyz(3) = unit_cell%T_xyz(3)

 if( frame == 1 ) CALL diagnosis( Extended_Cell )
    
 CALL BoundingBox(unit_cell)

 include 'formats.h'

 end subroutine generate_structure
! 
!
! 
!----------------------------------------------------------
!  the order orbitals are stored
! 
!       S      -->  1   --> l = 0  ,  m =  0           
!       Py     -->  2   --> l = 1  ,  m = -1    
!       Pz     -->  3   --> l = 1  ,  m =  0         
!       Px     -->  4   --> l = 1  ,  m = +1
!       Dxy    -->  5   --> l = 2  ,  m = -2      
!       Dyz    -->  6   --> l = 2  ,  m = -1
!       Dz2    -->  7   --> l = 2  ,  m =  0     
!       Dxz    -->  8   --> l = 2  ,  m = +1        
!       Dx2y2  -->  9   --> l = 2  ,  m = +2        
!----------------------------------------------------------
!
!
!
!==========================================
 subroutine Basis_Builder( system , basis )
!==========================================
 type(structure)               , intent(inout) :: system
 type(STO_basis) , allocatable , intent(out)   :: basis(:)

! local variables 
 integer :: k , i , l , m

! total number of orbitals
 N_of_orbitals = sum(atom(system%AtNo)%DOS)

! building AO basis  
 allocate(basis(N_of_orbitals))
 
 k = 1
 do i = 1 , system%atoms

    AtNo = system%AtNo(i)

    system%BasisPointer(i) = k-1  ! <== BasisPointer + {DOS} = {atom subspace}

    do l = 0 , atom(AtNo)%AngMax

        do m = -l , +l

            basis(k)%atom      =  i
            basis(k)%AtNo      =  AtNo
            basis(k)%copy_No   =  system%copy_No  (i)
            basis(k)%symbol    =  system%symbol   (i)
            basis(k)%fragment  =  system%fragment (i)
            basis(k)%EHSymbol  =  system%MMSymbol (i)
            basis(k)%residue   =  system%residue  (i)

            basis(k)%n         =  atom(AtNo)%Nquant(l)
            basis(k)%l         =  l
            basis(k)%m         =  m

            basis(k)%IP        =  atom(AtNo)%IP    (l)
            basis(k)%Nzeta     =  atom(AtNo)%Nzeta (l)
            basis(k)%coef(1)   =  atom(AtNo)%coef  (l,1)
            basis(k)%coef(2)   =  atom(AtNo)%coef  (l,2)
            basis(k)%zeta(1)   =  atom(AtNo)%zeta  (l,1)
            basis(k)%zeta(2)   =  atom(AtNo)%zeta  (l,2)
            basis(k)%k_WH      =  system%k_WH(i)

            basis(k)%x         =  system%coord (i,1)
            basis(k)%y         =  system%coord (i,2)
            basis(k)%z         =  system%coord (i,3)

            k = k + 1

        end do
    end do
 end do

 end subroutine Basis_Builder
!
!
!
!
 subroutine BoundingBox(a)

 type(structure) :: a

!  size of the box
 forall(i=1:3) &
 a%BoxSides(i) = maxval(a%coord(:,i)) - minval(a%coord(:,i)) + 4.0

!  find the center of mass
 forall(i=1:3) &
 a%Center_of_Mass(i) = sum(a%coord(:,i)) / a%atoms

 end subroutine BoundingBox
!
!
!
!=========================
 subroutine Diagnosis( a )
!=========================
 implicit none
 type(structure)    , intent(inout)    :: a

! local variables ...
integer :: N_of_orbitals, N_of_electron, N_of_atom_type, AtNo , residue , N_of_residue_type , fragment , N_of_fragment_type

! total number of orbitals ...
N_of_orbitals = sum(atom(a%AtNo)%DOS)
Print 120 , N_of_orbitals                       

! total number of electrons ...
a%N_of_electrons = sum(atom(a%AtNo)%Nvalen)
Print 140 , a%N_of_electrons

! total number of atoms ...
Print 141 , a%atoms

! total number of atoms of given type ...
do AtNo = 1 , size(atom)

    N_of_atom_type = count( a%AtNo == AtNo )
   
    If( N_of_atom_type /= 0 )       Print 121 , atom(AtNo)%symbol , N_of_atom_type

end do

! total number of residues ...
do residue = 1 , size(Unit_Cell%list_of_residues)

    N_of_residue_type = count( a%residue == Unit_Cell%list_of_residues(residue) )
   
    If( N_of_residue_type /= 0 )    Print 122 , Unit_Cell%list_of_residues(residue) , N_of_residue_type

end do

! total number of fragments ...
do fragment = 1 , size(a%list_of_fragments)

    N_of_fragment_type = count( a%fragment == a%list_of_fragments(fragment) )
    
    If( N_of_fragment_type /= 0 )   Print 123 , a%list_of_fragments(fragment) , N_of_fragment_type

end do

print * , " "

include 'formats.h'

end subroutine diagnosis
!
!
end module Structure_Builder

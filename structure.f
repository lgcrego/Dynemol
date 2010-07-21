 module Structure_Builder

    use type_m
    use Babel_m                     , only : Read_from_XYZ ,            &
                                             Read_from_Poscar ,         &
                                             Read_from_PDB ,            &
                                             Read_PDB ,                 &
                                             Read_VASP  ,               &
                                             Identify_Fragments ,       &
                                             System_Characteristics ,   & 
                                             trj 
    use Allocation_m                , only : Allocate_Structures
    use Semi_Empirical_Parms        , only : atom ,                     &
                                             Include_OPT_parameters

    type(structure)                 , public  :: Unit_Cell , Extended_Cell 
    type(STO_basis) , allocatable   , public  :: ExCell_basis(:)

    public :: Read_Structure , Generate_Structure , Basis_Builder 

    private

 contains
!
!
!=========================
subroutine Read_Structure
!=========================
implicit none

select case( file_type )

    case( "structure" )

        select case( file_format )
            case( "xyz" )
                CALL Read_from_XYZ      ( Unit_Cell ) 
            case( "pdb" )
                CALL Read_from_PDB      ( Unit_Cell )
            case( "vasp" )
                CALL Read_from_POSCAR   ( Unit_Cell )
            case default
                print*, ">>> check file type selection <<< : " , file_format
            end select

    case( "trajectory" )
    
        select case( file_format ) 
            case( "pdb" ) 
                CALL Read_PDB   ( trj ) 
            case( "vasp" )
                CALL Read_VASP  ( trj )
            case default
                print*, ">>> check file type selection <<< : " , file_format
            end select

    case default
        print*, ">>> check file type selection <<< : " , file_type

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
implicit none
integer , intent(in) :: frame

! local variables ...
integer :: copy , nr_sum , ix , iy , k , n
! local parameter: offset for ASCII numerical characters ...
integer :: ASC_offset = 48      

!----------------------------------------------------------
! GENERATES   THE   EXTENDED-STRUCTURE (REAL,not periodic)
!----------------------------------------------------------

 If( .NOT. allocated(Extended_Cell%coord) ) CALL Allocate_Structures( (2*nnx+1)*(2*nny+1)*unit_cell%atoms , extended_cell )

 k      =  0
 copy   =  0
 nr_sum =  0

 DO ix = -nnx , nnx 
 DO iy = -nny , nny

     If( (ix /= 0) .OR. (iy /= 0) ) THEN 

        copy    =  copy   + 1
        nr_sum  =  nr_sum + maxval( unit_cell%nr )

        FORALL( n=1:unit_cell%atoms )

            extended_cell % coord    (k+n,1) =  unit_cell % coord    (n,1) + ix * unit_cell%T_xyz(1)
            extended_cell % coord    (k+n,2) =  unit_cell % coord    (n,2) + iy * unit_cell%T_xyz(2)
            extended_cell % coord    (k+n,3) =  unit_cell % coord    (n,3) 
            extended_cell % AtNo     (k+n)   =  unit_cell % AtNo     (n)
            extended_cell % k_WH     (k+n)   =  unit_cell % k_WH     (n)
            extended_cell % fragment (k+n)   =  unit_cell % fragment (n)
            extended_cell % Symbol   (k+n)   =  unit_cell % Symbol   (n)
            extended_cell % MMSymbol (k+n)   =  unit_cell % MMSymbol (n)
            extended_cell % nr       (k+n)   =  unit_cell % nr       (n)   + nr_sum
            extended_cell % residue  (k+n)   =  unit_cell % residue  (n)
            extended_cell % copy_No  (k+n)   =  copy
        
        END FORALL

        k      =  k      + unit_cell%atoms
        nr_sum =  nr_sum + maxval( unit_cell%nr )

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
    extended_cell % nr       (k+n)      =  unit_cell % nr       (n)
    extended_cell % copy_No  (k+n)      =  0

 END FORALL

! define the DONOR fragment ...
 where( (extended_cell%fragment == 'M') .AND. (extended_cell%copy_No == 0) ) extended_cell%fragment = 'D' 
 where( (extended_cell%fragment == 'D') .AND. (extended_cell%copy_No /= 0) ) extended_cell%fragment = 'M' 
 where(  extended_cell%fragment == 'M'                                     ) extended_cell%fragment = achar( ASC_offset + extended_cell%copy_No) 

! create_&_allocate Extended_Cell%list_of_fragments ...     
 CALL Identify_Fragments( Extended_Cell )    

 extended_cell % N_of_Solvent_Molecules = (2*nnx+1) * (2*nny+1) * unit_cell % N_of_Solvent_Molecules

 extended_cell%T_xyz(1) = (2*nnx+1)*unit_cell%T_xyz(1)
 extended_cell%T_xyz(2) = (2*nny+1)*unit_cell%T_xyz(2)
 extended_cell%T_xyz(3) = unit_cell%T_xyz(3)

 if( frame == 1 ) CALL diagnosis( Extended_Cell )

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
 implicit none
 type(structure)               , intent(inout) :: system
 type(STO_basis) , allocatable , intent(out)   :: basis(:)

! local variables ...
 integer :: k , i , l , m , AtNo , N_of_orbitals

! total number of orbitals ...
 N_of_orbitals = sum(atom(system%AtNo)%DOS)

! building AO basis ...  
 allocate( basis(N_of_orbitals) )
 
 k = 1
 do i = 1 , system%atoms

    AtNo = system%AtNo(i)

    system%BasisPointer(i) = k-1  ! <== BasisPointer + {DOS} = {atom subspace}

    do l = 0 , atom(AtNo)%AngMax

        do m = -l , +l

            basis(k) % atom      =  i
            basis(k) % AtNo      =  AtNo
            basis(k) % nr        =  system % nr       (i)
            basis(k) % copy_No   =  system % copy_No  (i)
            basis(k) % symbol    =  system % symbol   (i)
            basis(k) % fragment  =  system % fragment (i)
            basis(k) % EHSymbol  =  system % MMSymbol (i)
            basis(k) % residue   =  system % residue  (i)

            basis(k) % n         =  atom(AtNo) % Nquant(l)
            basis(k) % l         =  l
            basis(k) % m         =  m

            basis(k) % IP        =  atom(AtNo) % IP    (l)
            basis(k) % Nzeta     =  atom(AtNo) % Nzeta (l)
            basis(k) % coef(1)   =  atom(AtNo) % coef  (l,1)
            basis(k) % coef(2)   =  atom(AtNo) % coef  (l,2)
            basis(k) % zeta(1)   =  atom(AtNo) % zeta  (l,1)
            basis(k) % zeta(2)   =  atom(AtNo) % zeta  (l,2)
            basis(k) % k_WH      =  system % k_WH(i)

            basis(k) % x         =  system % coord (i,1)
            basis(k) % y         =  system % coord (i,2)
            basis(k) % z         =  system % coord (i,3)

            basis(k) % indx      = k

            k = k + 1

        end do
    end do
 end do

 If( OPT_basis ) CALL Include_OPT_parameters( basis )

 end subroutine Basis_Builder
!
!
!
!
!=========================
 subroutine Diagnosis( a )
!=========================
 implicit none
 type(structure) , intent(inout) :: a

! local variables ...
integer :: N_of_orbitals, N_of_electron, N_of_atom_type, AtNo , residue , N_of_residue_type , fragment , N_of_fragment_type
integer :: first_nr , last_nr , N_of_residue_members 
integer :: i


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
   
!   finding positions of residues in structure and number of residue members ; for continuous blocks ...
    first_nr = minval( a%nr , (a%residue == Unit_Cell%list_of_residues(residue)) .AND. (a%copy_No == 0) )
    last_nr  = maxval( a%nr , (a%residue == Unit_Cell%list_of_residues(residue)) .AND. (a%copy_No == 0) )

    N_of_residue_members = last_nr - first_nr + 1

    If( N_of_residue_type /= 0 )    Print 122 , Unit_Cell%list_of_residues(residue) , N_of_residue_members , N_of_residue_type

end do

! total number of fragments ...
do fragment = 1 , size(a%list_of_fragments)

    N_of_fragment_type = count( a%fragment == a%list_of_fragments(fragment) )
    
    If( N_of_fragment_type /= 0 )   Print 123 , a%list_of_fragments(fragment) , N_of_fragment_type

end do

! dumping information about structure ...
CALL dump_diagnosis( a )

Print 47

include 'formats.h'

end subroutine diagnosis
!
!
!
!==============================
 subroutine dump_diagnosis( a )
!==============================
implicit none
type(structure) , intent(inout) :: a

! local variables ...
integer :: i , j

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
OPEN( unit=3 , file='structure.log' , status='unknown' )

write(3,*) System_Characteristics

do i = 1 , a%atoms

    write(3,100) i , a%Symbol(i) , a%MMSymbol(i) , a%nr(i) , a%residue(i) , a%fragment(i) , (a%coord(i,j) , j=1,3) 

end do

close(3)

100 format(I7,A6,A7,I5,A7,A4,3F10.4)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
OPEN(unit=4,file='structure.pdb',status='unknown')
write(4,6) 'COMPND' , System_Characteristics

write(4,1) 'CRYST1' , a%T_xyz(1) , a%T_xyz(2) , a%T_xyz(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

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
write(4,*) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a6,a72)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

end subroutine dump_diagnosis
!
!
end module Structure_Builder

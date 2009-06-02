 module Structure_Builder

    use type_m
    use Babel_m
    use Allocation_m
    use EHT_parameters

    type(structure)               , public             :: Unit_Cell , Extended_Cell 
    type(STO_basis) , allocatable , public , protected :: ExCell_basis(:)

    public :: Read_Structure , Generate_Structure , Basis_Builder 

    private

 contains
!
!
!
 subroutine Read_Structure

 If (POSCAR) then

    CALL Read_from_Poscar(unit_cell)

 else

    CALL Read_from_XYZ(unit_cell)

 end If
 
 end subroutine Read_Structure
!
!
!
 subroutine Generate_Structure(t)

 real*8 , intent(in) :: t

 integer :: copy, N_of_orbitals, N_of_electron, N_of_atom_type, AtNo

!----------------------------------------------------------
!           GENERATES   THE   STRUCTURE
!----------------------------------------------------------
!                  ORIGINAL  CELL 
!     places the molecule at the end of the list
!----------------------------------------------------------
! DEFINES  THE  EXTENDED-STRUCTURE (REAL,not periodic)

 CALL Allocate_Structures( (2*nnx+1)*(2*nny+1)*unit_cell%atoms , extended_cell )

 k = 0
 copy = 0

 DO ix = -nnx , nnx 
 DO iy = -nny , nny

     If( (ix /= 0) .OR. (iy /= 0) ) THEN 

        copy = copy + 1
        DO n = 1 , unit_cell%atoms

            k = k + 1

            extended_cell%coord    (k,1) =  unit_cell%coord(n,1) + ix * unit_cell%T_xyz(1)
            extended_cell%coord    (k,2) =  unit_cell%coord(n,2) + iy * unit_cell%T_xyz(2)
            extended_cell%coord    (k,3) =  unit_cell%coord(n,3)
            extended_cell%AtNo     (k)   =  unit_cell%AtNo(n)
            extended_cell%k_WH     (k)   =  unit_cell%k_WH(n)
            extended_cell%fragment (k)   =  unit_cell%fragment(n)
            extended_cell%symbol   (k)   =  atom( unit_cell%AtNo(n) )%symbol
            extended_cell%copy_No  (k)   =  copy
        
        END DO
     END IF

 END DO
 END DO

 extended_cell%molecule = k + unit_cell%molecule

 DO n = 1 , unit_cell%atoms     ! <== the DONOR CELL is at the end

    k = k + 1

    extended_cell%coord    (k,1:3)  =  unit_cell%coord    (n,1:3)
    extended_cell%AtNo     (k)      =  unit_cell%AtNo     (n)
    extended_cell%k_WH     (k)      =  unit_cell%k_WH     (n)
    extended_cell%fragment (k)      =  unit_cell%fragment (n)
    extended_cell%symbol   (k)      =  atom( unit_cell%AtNo(n) )%symbol
    extended_cell%copy_No  (k)      =  0

 END DO    

 extended_cell%T_xyz(1) = (2*nnx+1)*unit_cell%T_xyz(1)
 extended_cell%T_xyz(2) = (2*nny+1)*unit_cell%T_xyz(2)
 extended_cell%T_xyz(3) = unit_cell%T_xyz(3)

!------------------------------------------------------------

! total number of orbitals
 N_of_orbitals = sum(atom(extended_cell%AtNo)%DOS)
 Print 120 , N_of_orbitals                       

! total number of electrons
 extended_cell%N_of_electrons = sum(atom(extended_cell%AtNo)%Nvalen)
 Print 140 , extended_cell%N_of_electrons

! total number of atoms
 Print 141 , extended_cell%atoms

! total number of atoms of given type 
 do AtNo = 1 , size(atom)

    N_of_atom_type = count(extended_cell%AtNo == AtNo)
    
    If( N_of_atom_type /= 0 ) Print 121 , atom(AtNo)%symbol , N_of_atom_type

 end do
    
 print * , ' '
!------------------------------------------------------------

 CALL seed_Yaehmop(extended_cell)

 CALL seed_VASP(extended_cell)

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
 subroutine Basis_Builder(system,basis)

 type(structure)               , intent(inout) :: system
 type(STO_basis) , allocatable , intent(out)   :: basis(:)

 integer :: k , i , l , m

! total number of orbitals
 N_of_orbitals = sum(atom(system%AtNo)%DOS)

! => building AO basis <= 
 allocate(basis(N_of_orbitals))

 k = 1
 do i = 1 , system%atoms

    AtNo = system%AtNo(i)

    system%BasisPointer(i) = k-1  ! <== BasisPointer + {DOS} = {atom subspace}

    do l = 0 , atom(AtNo)%AngMax

        do m = -l , +l

            basis(k)%atom     =  i
            basis(k)%copy_No  = system%copy_No(i)
            basis(k)%AtNo     =  AtNo
            basis(k)%symbol   =  system%symbol(i)

            basis(k)%n        =  atom(AtNo)%Nquant(l)
            basis(k)%l        =  l
            basis(k)%m        =  m

            basis(k)%IP       =  atom(AtNo)%IP(l)
            basis(k)%Nzeta    =  atom(AtNo)%Nzeta(l)
            basis(k)%coef(1)  =  atom(AtNo)%coef(l,1)
            basis(k)%coef(2)  =  atom(AtNo)%coef(l,2)
            basis(k)%zeta(1)  =  atom(AtNo)%zeta(l,1)
            basis(k)%zeta(2)  =  atom(AtNo)%zeta(l,2)
            basis(k)%k_WH     =  system%k_WH(i)

            basis(k)%x        =  system%coord(i,1)
            basis(k)%y        =  system%coord(i,2)
            basis(k)%z        =  system%coord(i,3)

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
 end module Structure_Builder

















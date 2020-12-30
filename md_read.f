module MD_read_m

    use constants_m
    use atomicmass
    use MM_input       
    use parameters_m            , only : restart , ad_hoc , driver , preview , resume
    use MM_types                , only : MM_molecular, MM_atomic, debug_MM, DefinePairs
    use syst                    , only : bath_T, press, talt, talp, initial_density 
    use for_force               , only : KAPPA, Dihedral_potential_type, rcut, forcefield
    use MM_tuning_routines      , only : ad_hoc_MM_tuning 
    use gmx2mdflex              , only : itp2mdflex, top2mdflex, SpecialPairs
    use namd2mdflex             , only : psf2mdflex, prm2mdflex, convert_NAMD_velocities
    use Babel_m                 , only : QMMM_key
    use Structure_Builder       , only : Unit_Cell

    type(MM_molecular) , allocatable   :: molecule(:)
    type(MM_atomic)    , allocatable   :: atom(:) , FF(:)

    ! module variables ...
    public :: Build_MM_Environment , MM , atom , molecule , species , FF 

    private

contains
!
!
!
!================================
subroutine Build_MM_Environment
!================================
implicit none

! local variables ...
integer         :: i , j , k , l , atmax , Total_N_of_atoms_of_species_i , nresid , i1, i2, i3, sp, nr 


!=======================  setting up system  ============================= 

CALL Define_MM_Environment

bath_T        = temperature                     !  Temperature (K)
press         = pressure                        !  Pressure (atm)
rcut          = cutoff_radius                   !  Cutoff radius (Angs.)
talt          = thermal_relaxation_time         !  Temperature coupling
talp          = pressure_relaxation_time        !  Pressure coupling 
KAPPA         = damping_Wolf                    !  Wolf's method damping paramenter (length^{-1}) ; &
                                                !  Ref: J. Chem. Phys. 1999; 110(17):8254
atmax = sum( species(:) % N_of_atoms )

If( any( species % N_of_atoms == 0 ) ) stop ' >> you forgot to define a MM species ; check parameters_MM.f << '

! =====================================================================================
! types of molecules ...
CALL allocate_molecule( MM % N_of_molecules )

MM % N_of_atoms = 0
l = 1
do i = 1 , MM % N_of_species
    Total_N_of_atoms_of_species_i = species(i) % N_of_molecules * species(i) % N_of_atoms
    do j = 1 , species(i) % N_of_molecules
        molecule(l) % my_species = i
        l = l + 1
    end do
    MM % N_of_atoms = MM % N_of_atoms + Total_N_of_atoms_of_species_i 
end do

allocate ( atom ( MM % N_of_atoms ) )
k = 1
do i = 1 , MM % N_of_species
    Total_N_of_atoms_of_species_i = species(i) % N_of_molecules * species(i) % N_of_atoms
    do j = 1 , Total_N_of_atoms_of_species_i
        atom(k) % my_species = i
        k = k + 1
    end do
end do

! pass information from structure to molecular dynamics ...
CALL Structure_2_MD

!----------------------
! finished reading ...
!----------------------

! get Atomic Number (AtNo) ...
CALL Symbol_2_AtNo( atom )

! Define atomic mass ...
atom % mass = Atomic_mass( atom % AtNo )

If( Selective_Dynamics ) CALL ad_hoc_MM_tuning( atom , instance="MegaMass")

do i = 1 , MM % N_of_species
    species(i) % mass = sum( atom % mass , atom % my_species == i ) / species(i) % N_of_molecules
    where( molecule % my_species == i ) molecule % mass = species(i) % mass
end do

! mass of molecule species in kg ...
species % mass = species % mass * imol
molecule % mass = molecule % mass * imol

! initial density of the box ... 
initial_density = sum( molecule % mass ) * MM % ibox(1) * MM % ibox(2) * MM % ibox(3) * (mts_2_nano**3)

99  format(a72)
100 format(t10, f6.3, t19, f6.3, t28, f6.3)
105 format(a6)
115 FORMAT(t8,i4,t14,a3,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3)

!=======================  reading  potential.inpt  ============================= 
CALL allocate_FF( atmax )

select case ( MM_input_format )

    case( "GMX" )  ! <== reads FF data in GMX format ... 
     
         If( ad_hoc ) CALL ad_hoc_MM_tuning(instance="SpecialBonds")

         CALL itp2mdflex( MM , atom , species , FF)

         CALL top2mdflex( MM , species , FF )

         do i = 1 , size(species)
             CALL MMSymbol_2_Symbol( species(i) % atom )
         end do

         CALL MMSymbol_2_Symbol ( FF )
         CALL Symbol_2_AtNo     ( FF )

         do i = 1 , MM % N_of_species
             do j = 1 , species(i) % N_of_atoms
                 where( FF % residue == species(i) % atom(j) % residue )
                     FF % nr         = species(i) % atom(j) % nr
                     FF % flex       = species(i) % atom(j) % flex
                     FF % my_species = species(i) % my_species
                 end where
             end do
         end do

    case( "NAMD" , "GAFF" )  ! <== reads FF data in NAMD format ... 

         If( ad_hoc ) CALL ad_hoc_MM_tuning(instance="SpecialBonds")
  
         CALL psf2mdflex( MM , atom , species , FF)
 
         CALL prm2mdflex( MM , species , FF)

         do i = 1 , size(species)
             CALL MMSymbol_2_Symbol( species(i) % atom )
         end do

        CALL MMSymbol_2_Symbol ( FF )
        CALL Symbol_2_AtNo     ( FF )

         do i = 1 , MM % N_of_species
             do j = 1 , species(i) % N_of_atoms
                 where( FF % residue == species(i) % atom(j) % residue )
                     FF % nr         = species(i) % atom(j) % nr
                     FF % flex       = species(i) % atom(j) % flex
                     FF % my_species = species(i) % my_species
                 end where
             end do
         end do

    case default
         Print*, " >>> Check your MM_input_format option in parameters_MM.f <<< :" , MM_input_format
         stop

end select 

If( ad_hoc ) CALL ad_hoc_MM_tuning( atom , instance = "General" )

Unit_Cell% flex(:) = atom(:)% flex
!=======================  finished  reading  potential.inpt  ============================= 

!=======================         set-up molecule(:)          ============================= 
do i = 1 , MM % N_of_species
    where( molecule % my_species == i ) molecule % N_of_atoms = species(i) % N_of_atoms
    where( molecule % my_species == i ) molecule % Nbonds     = species(i) % Nbonds
    where( molecule % my_species == i ) molecule % Nangs      = species(i) % Nangs
    where( molecule % my_species == i ) molecule % Ndiheds    = species(i) % Ndiheds
    where( molecule % my_species == i ) molecule % Nbonds14   = species(i) % Nbonds14
    where( molecule % my_species == i ) molecule % NintraLJ   = species(i) % NintraLJ
end do

!=======================         set-up atom(:) <=> FF(:)    ============================= 

! defining atom % my_intra_id ...
do i = 1 , size(atom)
    atom(i) % my_intra_id = i + molecule( atom(i) % nr ) % N_of_atoms - sum( molecule(1:atom(i) % nr) % N_of_atoms )
end do

! defining molecule % nr from atom % nr ...
do i = 1 , MM % N_of_atoms
    nresid = atom(i) % nr
    molecule(nresid) % nr = nresid
end do

! passing MMSymbol from FF to atom ...
i1 = 1
do nr = 1 , atom(MM%N_of_atoms) % nr
    sp = atom(i1) % my_species
    i3 = count(atom%nr==nr)
    i2 = i1 + (i3-1)
    atom(i1:i2)%MMSymbol = pack(FF % MMSymbol , FF % my_species == sp)
    i1 = i2+1
end do

atom % MMSymbol = adjustl(atom % MMSymbol)

do i = 1 , size(FF)
    where( atom % MMSymbol == FF(i) % MMSymbol )
        atom % eps   = FF(i) % eps
        atom % eps14 = FF(i) % eps14
        atom % sig   = FF(i) % sig
        atom % sig14 = FF(i) % sig14
    end where
end do

!=======================         set-up molecule(:)          ============================= 
do i = 1 , MM % N_of_molecules

    k = size( species(molecule(i) % my_species) % kdihed0(1,:) )

    If( molecule(i)%Nbonds14 > 0 ) & 
    allocate( molecule(i) % bonds14       ( molecule(i) % Nbonds14 , 2 ) )

    If( molecule(i)%Nbonds > 0 ) then
    allocate( molecule(i) % bonds         ( molecule(i) % Nbonds   , 2 ) )
    allocate( molecule(i) % kbond0        ( molecule(i) % Nbonds   , 3 ) )
    allocate( molecule(i) % bond_type     ( molecule(i) % Nbonds       ) )
    allocate( molecule(i) % funct_bond    ( molecule(i) % Nbonds       ) )
    End If

    If( molecule(i)%Nangs > 0 ) then
    allocate( molecule(i) % angs          ( molecule(i) % Nangs    , 3 ) )
    allocate( molecule(i) % kang0         ( molecule(i) % Nangs    , 4 ) )
    allocate( molecule(i) % angle_type    ( molecule(i) % Nangs        ) )
    allocate( molecule(i) % funct_angle   ( molecule(i) % Nangs        ) )
    End If

    If( molecule(i)%Ndiheds > 0 ) then
    allocate( molecule(i) % diheds        ( molecule(i) % Ndiheds  , 4 ) )
    allocate( molecule(i) % kdihed0       ( molecule(i) % Ndiheds  , k ) )
    allocate( molecule(i) % Dihedral_Type ( molecule(i) % Ndiheds      ) )
    allocate( molecule(i) % funct_dihed   ( molecule(i) % Ndiheds      ) )
    End If

    If( molecule(i)%NintraLJ > 0 ) & 
    allocate( molecule(i) % IntraLJ       ( molecule(i) % NintraLJ , 2 ) )

end do

k = 0
do i = 1 , MM % N_of_molecules

    If( molecule(i)%Nbonds14 > 0 ) & 
    molecule(i) % bonds14       = species(molecule(i) % my_species) % bonds14 + k

    If( molecule(i)%Nbonds > 0 ) then
    molecule(i) % kbond0        = species(molecule(i) % my_species) % kbond0
    molecule(i) % bonds         = species(molecule(i) % my_species) % bonds + k 
    molecule(i) % bond_type     = species(molecule(i) % my_species) % bond_type
    molecule(i) % funct_bond    = species(molecule(i) % my_species) % funct_bond
    End If

    If( molecule(i)%Nangs > 0 ) then
    molecule(i) % kang0         = species(molecule(i) % my_species) % kang0
    molecule(i) % angs          = species(molecule(i) % my_species) % angs + k
    molecule(i) % angle_type    = species(molecule(i) % my_species) % angle_type
    molecule(i) % funct_angle   = species(molecule(i) % my_species) % funct_angle
    End If

    If( molecule(i)%Ndiheds > 0 ) then
    molecule(i) % kdihed0       = species(molecule(i) % my_species) % kdihed0
    molecule(i) % diheds        = species(molecule(i) % my_species) % diheds + k
    molecule(i) % Dihedral_Type = species(molecule(i) % my_species) % Dihedral_Type
    molecule(i) % funct_dihed   = species(molecule(i) % my_species) % funct_dihed
    End If

    If( molecule(i)%NintraLJ > 0 ) & 
    molecule(i) % IntraLJ       = species(molecule(i) % my_species) % IntraLJ + k

    k = k + molecule(i) % N_of_atoms

end do

do i = 1 , MM % N_of_species
    k = count( molecule % my_species == species(i) % my_species )
    where( molecule % my_species == i ) 
        molecule % N_of_molecules = k
        molecule % residue        = species(i) % residue
        molecule % flex           = species(i) % flex
    end where
end do

!call debug_MM( atom )
!========================================================================================= 

CALL MM_diagnosis( )

end subroutine Build_MM_Environment
!
!
!
!================================
subroutine allocate_molecule( N )
!================================
implicit none
integer , intent(in)    :: N

! local variables ...
integer :: i

allocate( molecule ( N ) )

do i = 1 , N
    molecule(i) % my_species     = 0
    molecule(i) % N_of_atoms     = 0
    molecule(i) % N_of_molecules = 0
    molecule(i) % cm(3)          = 0.0d0
    molecule(i) % mass           = 0.0d0
    molecule(i) % flex           = .false.
    molecule(i) % residue        = "XXX"
    molecule(i) % nr             = 0
    molecule(i) % Nbonds         = 0
    molecule(i) % Nangs          = 0
    molecule(i) % Ndiheds        = 0
    molecule(i) % Nbonds14       = 0
end do

end subroutine allocate_molecule
!
!
!
!==========================
subroutine allocate_FF( N )
!==========================
implicit none
integer , intent(in)    :: N

! local variables ...
integer :: i

allocate( FF ( N ) )

do i = 1 , N
    FF(i) % AtNo         = 0
    FF(i) % my_id        = 0
    FF(i) % my_species   = 0
    FF(i) % nr           = 0
    FF(i) % residue      = "XXX"
    FF(i) % Symbol       = "XXX"
    FF(i) % MMSymbol     = "XXX"
    FF(i) % EHSymbol     = "XXX"
    FF(i) % xyz(:)       = 0.0d0
    FF(i) % vel(:)       = 0.0d0
    FF(i) % fbond(:)     = 0.0d0
    FF(i) % fang(:)      = 0.0d0
    FF(i) % fdihed(:)    = 0.0d0
    FF(i) % fnonbd14(:)  = 0.0d0
    FF(i) % fnonch14(:)  = 0.0d0
    FF(i) % f_MM(:)      = 0.0d0
    FF(i) % Ehrenfest(:) = 0.0d0
    FF(i) % ftotal(:)    = 0.0d0
    FF(i) % fch(:)       = 0.0d0
    FF(i) % fsr(:)       = 0.0d0
    FF(i) % mass         = 0.0d0
    FF(i) % charge       = 0.0d0
    FF(i) % MM_charge    = 0.0d0
    FF(i) % eps          = 0.0d0
    FF(i) % eps14        = 0.0d0
    FF(i) % sig          = 0.0d0
    FF(i) % sig14        = 0.0d0
    FF(i) % flex         = .false.
end do

end subroutine allocate_FF
!
!
!
!===========================
 subroutine Symbol_2_AtNo(a)
!===========================
type( MM_atomic ) , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(adjustl(a(i)%Symbol))
        case( 'H')
            a(i)%AtNo = 1
        case( 'LI','Li')
            a(i)%AtNo = 3
        case( 'BE','Be')
            a(i)%AtNo = 4
        case( 'B' )
            a(i)%AtNo = 5
        case( 'C')
            a(i)%AtNo = 6
        case( 'N')
            a(i)%AtNo = 7
        case( 'O')
            a(i)%AtNo = 8
        case( 'F')
            a(i)%AtNo = 9
        case( 'AL','Al')
            a(i)%AtNo = 13
        case( 'SI','Si')
            a(i)%AtNo = 14
        case( 'P' )
            a(i)%AtNo = 15
        case( 'S','s')
            a(i)%AtNo = 16
        case( 'CL','Cl')
            a(i)%AtNo = 17
        case( 'TI','Ti')
            a(i)%AtNo = 22
        case( 'MN','Mn')
            a(i)%AtNo = 25
        case( 'FE','Fe')
            a(i)%AtNo = 26
        case( 'RU','Ru')
            a(i)%AtNo = 44
        case( 'I')
            a(i)%AtNo = 53
        case( 'Pb')
            a(i)%AtNo = 82
        case default
            print*, ' >> Symbol_2_AtNo ; unknown atom found <<' , '[',a(i)%Symbol,']' , i
            stop
    end select

 END DO

 end subroutine Symbol_2_AtNo
!
!
!
!=========================
 subroutine Structure_2_MD
!=========================
implicit none

! local variables ...
logical :: exist
integer :: i , j , nr , indx

MM % box  = Unit_Cell % T_xyz
MM % ibox =  D_one / MM % box

do j = 1 , MM % N_of_atoms

    i = QMMM_key(j)

    atom(i) % AtNo         = Unit_Cell % AtNo(j)
    atom(i) % my_id        = i
    atom(i) % nr           = Unit_Cell % nr(j)
    atom(i) % residue      = adjustl(Unit_Cell % residue(j))
    atom(i) % Symbol       = adjustl(Unit_Cell % Symbol(j))
    atom(i) % EHSymbol     = adjustl(Unit_Cell % MMSymbol(j))
    atom(i) % xyz(:)       = Unit_Cell % coord(j,:)
    atom(i) % vel(:)       = 0.d0
    atom(i) % fbond(:)     = 0.d0
    atom(i) % fang(:)      = 0.d0
    atom(i) % fdihed(:)    = 0.d0
    atom(i) % fnonbd14(:)  = 0.d0
    atom(i) % fnonch14(:)  = 0.d0
    atom(i) % f_MM(:)      = 0.d0
    atom(i) % Ehrenfest(:) = 0.d0
    atom(i) % ftotal(:)    = 0.d0
    atom(i) % fch(:)       = 0.d0
    atom(i) % fsr(:)       = 0.d0
    atom(i) % mass         = 0.d0
    atom(i) % charge       = 0.d0
    atom(i) % MM_charge    = 0.d0
    atom(i) % eps          = 0.d0
    atom(i) % eps14        = 0.d0
    atom(i) % sig          = 0.d0
    atom(i) % sig14        = 0.d0
    atom(i) % flex         = Unit_Cell % flex(i)
end do

! this is assumed a priori , but can be changed otherwise if required by the Force Field ...
atom % MMSymbol = atom % EHSymbol  

If( ad_hoc ) CALL ad_hoc_MM_tuning( atom , instance = "General" )

indx = 1
do nr = 1 , atom( MM % N_of_atoms ) % nr
    do i = 1 , count( atom%nr == nr ) 
        atom(indx) % my_intra_id = i
        indx = indx + 1
    end do
end do

! setting up initial velocities ...
do i = 1 , MM % N_of_atoms
    atom(i) % vel(1:3) = 0.d0
end do

if( read_velocities ) then

    if( restart ) STOP ' Should NOT read velocity_MM.inpt if restart = .true. !'

    if( preview ) CALL convert_NAMD_velocities( MM% N_of_atoms )

    inquire(file="velocity_MM.out", EXIST=exist)
    if (exist .AND. resume) then
        CALL system("sed '11i >> must update inpt file:   mv[velocity_MM.out ==> velocity.inpt]   or   rm velocity_MM.out << ' warning.signal |cat")
        STOP 
    end If

    inquire(file="velocity_MM.inpt", EXIST=exist)
    if (exist) then
        open(unit=33 , file='velocity_MM.inpt' , status='old' , action='read')
        do i = 1 , size(atom)
            read(33,*) atom(i) % vel(1) , atom(i) % vel(2) , atom(i) % vel(3)
        end do
        close(33)
    else
        STOP ' >> read_velocity = .true. but file velocity_MM.inpt was not found ! << '
    endif

elseif( resume ) then 

       ! read_velocity flag = F_
        CALL system("sed '11i >> read_velocity = .false. in parametes_MM.f, must be true in resume simulations ! << ' warning.signal |cat")
        STOP 

end if

! check list of input data ...
If( sum(species%N_of_Molecules * species%N_of_atoms) /= Unit_Cell%atoms ) then
    CALL system("sed '11i >>> error: sum(species%N_of_Molecules * species%N_of_atoms) /= Unit_Cell%atoms ; check parameters_MM.f <<<' warning.signal |cat")
    STOP 
end If

If( Unit_Cell%atoms /= MM% N_of_atoms ) then
    CALL system("sed '11i >>> error: Unit_Cell%atoms /= MM% N_of_atoms <<<' warning.signal |cat")
    STOP 
end If

If( maxval(atom%nr) < MM%N_of_species ) then
    CALL system("sed '11i >>> # of residues must be (>=) # of species; check input.pdb and parameters_MM.f <<<' warning.signal |cat")
    STOP 
end If

end subroutine Structure_2_MD
!
!
!
!===============================
 subroutine MMSymbol_2_Symbol(a)
!===============================
type(MM_atomic) , intent(inout) :: a(:)

! local variables ...
integer             :: i
character(len=1)    :: element

DO i = 1 , size(a)

    write( element,'(A1)' ) adjustl( a(i)%MMSymbol )

    select case( element )
        case( 'I' )
            a(i)%Symbol = 'I'
        case( 'P' )
            a(i)%Symbol = 'P'
        case( 'S' )
            a(i)%Symbol = 'S'
        case( 'C' )
            a(i)%Symbol = 'C'
        case( 'N' )
            a(i)%Symbol = 'N'
        case( 'O' )
            a(i)%Symbol = 'O'
        case( 'H' )
            a(i)%Symbol = 'H'
    end select

    select case( adjustl( a(i)%MMSymbol) )
        case( 'Ox' , 'OxH' , 'OS' )
            a(i)%Symbol = 'O'
        case( 'YN' , 'NTr' , 'Nx' , 'N2' , 'N32' , 'N42' , 'N52' )
            a(i)%Symbol = 'N'
        case( 'Al' )
            a(i)%Symbol = 'Al'
        case( 'Ti' , 'TI' )
            a(i)%Symbol = 'Ti'
        case( 'Li' )
            a(i)%Symbol = 'Li'
        case( 'Ru' )
            a(i)%Symbol = 'Ru'
        case( 'HC' , 'HA' )
            a(i)%Symbol = 'H'
        case( 'C=' , 'CTr' , 'CS' , 'CC' , 'CM' , 'YC' , 'CT' , 'CA' , 'CB' )
            a(i)%Symbol = 'C'
        case( 'SS' )
            a(i)%Symbol = 'S'
        case( 'Pb' , 'PB' )
            a(i)%Symbol = 'Pb'
        case( 'P' )
            a(i)%Symbol = 'P'
    end select

END DO

a % Symbol = adjustl(a % Symbol)

end subroutine MMSymbol_2_Symbol
!
!
!
!
!==========================
 subroutine MM_diagnosis( ) 
!==========================
implicit none

! local variabbles ...
integer                         :: i , j , m , at1 , at2 , at3 , at4 , funct_dih , multiples , prototype
real*8                          :: factor , factor_1 , factor_2 
character(3)                    :: funct_type , flag
character(len=:) , allocatable  :: string(:)

 open( unit = 51 , file = "log.trunk/MM_parms_log.out" , status = "replace", action = "write" , position = "append" )

 !========================================================================================================
 write(51, *) " "
 ! Print MM_input_format
 write (51, 215) MM_input_format
                
 write(51, *) " "
 ! Print # of atoms-types
 write (51, 205) MM % N_of_AtomTypes , &
                 count( [ ( .NOT. any(FF(1:i-1)% MMSymbol == FF(i)% MMSymbol) , i=1,size(FF)) ] ) 

 write(51, *) " "
 do i = 1 , MM % N_of_species
 ! Print # of atoms
     write (51, 201) species(i) % residue, species(i) % N_of_atoms
 ! Print # of bonds
     write(51,202 )  species(i) % residue, species(i) % Nbonds
 ! Print # of angles
     write(51, 203) species(i) % residue, species(i) % Nangs
 ! Print # of dihedrals
     write(51, 204) species(i) % residue, species(i) % Ndiheds
 ! Print # of Torsion DHDs
     write(51, 214) species(i) % residue, species(i) % NTorsions
 ! Print # of Improper DHSs
     write(51, 224) species(i) % residue, species(i) % NImpropers
 ! Print total MM_charge
     write(51, 225) species(i) % residue, sum(species(i)%atom%MM_charge)
     write(51, *) " "
 end do
 !========================================================================================================
 ! Force Field Parameters ...
 write(51,"(A)") "Force Field Parameters:"               

 write(51, 206) forcefield
 write(51, 207) MM % CombinationRule
 write(51, 208) MM % fudgeQQ
 write(51, 209) MM % fudgeLJ

 !========================================================================================================
 ! atom types saving ...
 write(51, *) " "
 write(51,"(A)") "[ atomtypes ]"               

 ! General NonBonded parms ...
 do i = 1 , size(FF)

    if( .NOT. any(FF(1:i-1)% MMSymbol == FF(i)% MMSymbol) ) then 

        ! warns if NB paramater was not assigned to this atom  ...
        flag = merge( "<==" , "   " , FF(i)% sig * FF(i)%eps == 0 )

        write(51,'(I5,A5,2F12.5,A4)') count(FF% MMSymbol == FF(i)% MMSymbol) , &
                                      FF(i)% MMSymbol                        , & 
                                      FF(i)% sig                             , &
                                      FF(i)% eps                             , &
                                      flag

    end if

 end do

 write(51, *) " "
 write(51,"(A)") "[ SpecialPairs ]"               

 ! NonBonded SpecialPairs parms ...
 do i = 1 , size(SpecialPairs)

    ! warns if NB paramater was not assigned to this atom  ...
    flag = merge( "<==" , "   " , SpecialPairs(i)% Parms(1) * SpecialPairs(i)% Parms(2) == 0 )

    write(51,'(2A5,2F12.5,A4)')   SpecialPairs(i)% MMSymbols(1)          , & 
                                  SpecialPairs(i)% MMSymbols(2)          , & 
                                  SpecialPairs(i)% Parms(1)              , &
                                  SpecialPairs(i)% Parms(2)              , &
                                  flag

 end do

 !========================================================================================================
 ! bond parms saving ...
 write(51, *) " "
 write(51,"(A)") "[ bondtypes ]"               

 prototype = 1
 do m = 1 , MM % N_of_Molecules
 
    if( molecule(m)%my_species /= prototype ) cycle
 
    allocate( character(len=2*len(atom(at1)%MMSymbol)) :: string(molecule(m)%Nbonds) )
 
    do i = 1 , molecule(m)%Nbonds
 
       at1 = molecule(m)%bonds(i,1)
       at2 = molecule(m)%bonds(i,2)
 
       string(i) = atom(at1)%MMSymbol//atom(at2)%MMSymbol
 
       if( .NOT. any(string(1:i-1) == string(i)) ) then 
 
           ! warns if paramater was not assigned to this bond ...
           flag = merge( "<==" , "   " , sum(molecule(m)%kbond0(i,:)) == 0 )
 
           funct_type = molecule(m) % funct_bond(i) 
    
           factor = factor2 * imol  
           if( funct_type == "3" ) factor = factor1 * imol
 
           write(51,'(3A4,F15.5,2F15.3,A3)')  atom(at1)%MMSymbol                      , &
                                              atom(at2)%MMSymbol                      , &
                                              funct_type                              , &
                                              molecule(m)%kbond0(i,2) / nano_2_angs   , &
                                              molecule(m)%kbond0(i,1) / factor        , &
                                              molecule(m)%kbond0(i,3) * nano_2_angs   , &
                                              flag
       end if
    end do
    deallocate(string)
 
    prototype = prototype + 1 
 
 end do

 !========================================================================================================
 ! angle parms saving ...
 write(51,*) " "
 write(51,"(A)") "[ angletypes ]"

 prototype = 1
 do m = 1 , MM % N_of_Molecules
 
    if( molecule(m)%my_species /= prototype ) cycle
 
    allocate( character(len=3*len(atom(at1)%MMSymbol)) :: string(molecule(m)%Nangs) )

    do i = 1 , molecule(m)%Nangs

       at1 = molecule(m)%angs(i,1)
       at2 = molecule(m)%angs(i,2)
       at3 = molecule(m)%angs(i,3)

       string(i) = atom(at1)%MMSymbol//atom(at2)%MMSymbol//atom(at3)%MMSymbol

       if( .NOT. any(string(1:i-1) == string(i)) ) then 

           ! warns if paramater was not assigned to this angle ...
           flag = merge( "<==" , "   " , sum(molecule(m)%kang0(i,:)) == 0 )

           funct_type = molecule(m) % funct_angle(i)

           factor_1 = factor1 * imol
           factor_2 = factor2 * imol

           write(51,'(4A4,2F15.3)',advance="no") atom(at1)%MMSymbol  , &
                                                 atom(at2)%MMSymbol  , &
                                                 atom(at3)%MMSymbol  , &
                                                 funct_type       

           select case( adjustl(molecule(m) % Angle_Type(i)) )

               case ('harm') 

                   write(51,'(2F15.3,A3)') molecule(m)%kang0(i,2) / deg_2_rad   , &
                                           molecule(m)%kang0(i,1) / factor_1    , &
                                           flag 

               case('urba')

                   write(51,'(4F15.3,A3)') molecule(m)%kang0(i,2) / deg_2_rad   , &
                                           molecule(m)%kang0(i,1) / factor_1    , &
                                           molecule(m)%kang0(i,4) / nano_2_angs , &
                                           molecule(m)%kang0(i,3) / factor_2    , &
                                           flag

               case default

                   write(*,'(A5)',advance="no") adjustl(molecule(m) % Angle_Type(i))
                   stop " <== angle FF not supported in FF_OPT_class%output"

           end select
       end if
    end do
    deallocate(string)
 
    prototype = prototype + 1 
 
 end do

 !========================================================================================================
 ! dihedral parms saving ...
 write(51,*) " "
 write(51,"(A)") "[ dihedraltypes ]"

 prototype = 1
 do m = 1 , MM % N_of_Molecules
 
    if( molecule(m)%my_species /= prototype ) cycle

    allocate( character(len=4*len(atom(at1)%MMSymbol)+len(molecule(m)%Dihedral_Type)) :: string(molecule(m)%Ndiheds) )
    do i = 1 , molecule(m)%Ndiheds

       at1 = molecule(m)%diheds(i,1)
       at2 = molecule(m)%diheds(i,2)
       at3 = molecule(m)%diheds(i,3)
       at4 = molecule(m)%diheds(i,4)

       string(i) = atom(at1)%MMSymbol//atom(at2)%MMSymbol//atom(at3)%MMSymbol//atom(at4)%MMSymbol//molecule(m)%Dihedral_Type(i)

       if( (.NOT. any(string(1:i-1) == string(i))) .OR. (.NOT. any(molecule(m)%kdihed0(1:i-1,1) == molecule(m)%kdihed0(i,1))) ) then 

           ! warns if paramater was not assigned to this dihedral ...
           flag = merge( "<==" , "   " , sum(abs(molecule(m)%kdihed0(i,:))) == 0 )

           funct_dih = molecule(m) % funct_dihed(i)

           factor = factor1 * imol

           write(51,'(4A4,I5)',advance="no") atom(at1)%MMSymbol                  , &
                                             atom(at2)%MMSymbol                  , &
                                             atom(at3)%MMSymbol                  , &
                                             atom(at4)%MMSymbol                  , &
                                             funct_dih      

           select case( adjustl(molecule(m) % Dihedral_Type(i)) )

               case ('cos' , 'imp')  ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] ; Eq. 4.60 (GMX 5.0.5 manual)

                   write(51,'(3F12.5,A3)') molecule(m)%kdihed0(i,1) / deg_2_rad  , &
                                           molecule(m)%kdihed0(i,2) / factor     , &  
                                           molecule(m)%kdihed0(i,3)              , &
                                           flag

               case ('cos3') ! V = C0 + C1*cos(phi-180) + C2*cos^2(phi-180) + C3*cos^3(phi-180) + C4*cos^4(phi-180) + C5*cos(phi-180)  
                             ! Eq. 4.61 (GMX 5.0.5 manual)

                   write(51,'(6F12.5,A3)') molecule(m)%kdihed0(i,1) / factor     , &
                                           molecule(m)%kdihed0(i,2) / factor     , &  
                                           molecule(m)%kdihed0(i,3) / factor     , &
                                           molecule(m)%kdihed0(i,4) / factor     , & 
                                           molecule(m)%kdihed0(i,5) / factor     , &
                                           molecule(m)%kdihed0(i,6) / factor     , &
                                           flag

               case ('harm') ! V = 1/2.k[cos(phi) - cos(phi0)]²
                        ! factor1 = 1.0d26      <== Factor used to correct the units 
                        ! kdihed0(:,1) = xi_0   ==> angle (deg) * deg_2_rad
                        ! kdihed0(:,2) = K_(xi) ==> force constant (kcal.mol⁻¹.rad⁻²) * factor1 * imol * cal_2_J

                   write(51,'(6F12.5,A3)') molecule(m)%kdihed0(i,1) / factor     , &
                                           molecule(m)%kdihed0(i,2) / factor     , &  
                                           molecule(m)%kdihed0(i,3) / factor     , &
                                           flag

               case ('chrm')  ! V = k_phi * [ 1 + cos( n * phi - phi_s ) ] (multiple) ; Eq. 4.60 (GMX 5.0.5 manual)


                   multiples = count(molecule(m)%kdihed0(i,:) /= 0) 

                       If( multiples <= 9 ) then   

                            write(51,'(3F12.5)',advance='no') molecule(m)%kdihed0(i,1) / deg_2_rad  , &
                                                              molecule(m)%kdihed0(i,2) / factor     , &  
                                                              molecule(m)%kdihed0(i,3)              
                       if( multiples >= 4 ) then    

                            write(51,'(3F12.5)',advance='no') molecule(m)%kdihed0(i,4) / deg_2_rad  , &
                                                              molecule(m)%kdihed0(i,5) / factor     , &  
                                                              molecule(m)%kdihed0(i,6)              
                       if( multiples >= 7 ) then    

                            write(51,'(3F12.5)',advance='no') molecule(m)%kdihed0(i,7) / deg_2_rad  , &
                                                              molecule(m)%kdihed0(i,8) / factor     , &  
                                                              molecule(m)%kdihed0(i,9)              
                       endif; endif; EndIf
                   write(51,'(A3)') flag

               case default

                   write(*,'(A5)',advance="no") adjustl(molecule(m) % Dihedral_Type(i))
                   stop " <== dihedral FF not supported in FF_OPT_class%output"

           end select
       end if
    end  do
    deallocate(string)
    
    prototype = prototype + 1 
 
 end do

 !========================================================================================================
 ! charge parms saving ...
 write(51, *) " "
 write(51,"(A)") "[ charges ]"               

 do j = 1 , size(species)
    
    write(51,"(3A)") "( ",species(j)%residue," )"               
    write(51,"(A5,F8.4,A5)") (species(j)% atom(i)% MMSymbol , species(j)% atom(i) % MM_Charge , &
                              merge("<==" , "   " , species(j)% atom(i) % MM_Charge == 0), i = 1,species(j)% N_of_atoms)
 end do
!========================================================================================================

 close(51)

include 'formats.h'

end subroutine MM_diagnosis
!
!
!
end  module MD_read_m

module MD_read_m
 
    use constants_m
    use atomicmass
    use MM_input                
    use parameters_m            , only : restart , ad_hoc , driver , preview
    use MM_types                , only : MM_molecular, MM_atomic, debug_MM, DefinePairs
    use syst                    , only : bath_T, press, talt, talp, initial_density 
    use for_force               , only : KAPPA, Dihedral_potential_type, rcut
    use MM_tuning_routines      , only : ad_hoc_MM_tuning 
    use gmx2mdflex              , only : itp2mdflex, top2mdflex
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

if ( driver == "diagnostic" ) return

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

    case( "NAMD" )  ! <== reads FF data in NAMD format ... 

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
    ! redefining atom % nr for NAMD input format ... 
    if ( MM_input_format == "NAMD" ) then
            if( atom(i) % my_species == 1 ) then 
                    atom(i) % nr = ceiling( real(i) / real( species(1) % N_of_atoms ) )
            else 
                    atom(i) % nr = atom(i) % nr + sum( species( 1:atom(i)%my_species-1 ) % N_of_molecules ) 
            end if 
    end if
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

end if

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
end  module MD_read_m

module MD_read_m
 
    use constants_m
    use atomicmass
    use MM_input                
    use parameters_m            , only : restart
    use MM_types                , only : MM_system , MM_molecular , MM_atomic , debug_MM
    use syst                    , only : bath_T, press, talt, talp, initial_density 
    use for_force               , only : KAPPA, Dihedral_potential_type, forcefield, rcut
    use MM_tuning_routines      , only : ad_hoc_MM_tuning 
    use gmx2mdflex              , only : itp2mdflex, top2mdflex
    use Babel_m                 , only : QMMM_key
    use Structure_Builder       , only : Unit_Cell

    type(MM_molecular) , allocatable   :: molecule(:)
    type(MM_atomic)    , allocatable   :: atom(:) , FF(:)

    ! module variables ...
    logical :: read_from_gmx

    public :: Build_MM_Environment , MM , atom , molecule , species , FF , read_from_gmx

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
real*8          :: massa 
integer         :: i , j , k , l , a , b , atmax , Total_N_of_atoms_of_species_i , ioerr , nresid 
logical         :: exist , read_vel
character(10)   :: string


!=======================  setting up system  ============================= 

CALL Define_MM_Environment

bath_T        = temperature                     !  Temperature (K)
press         = pressure                        !  Pressure
rcut          = cutoff_radius                   !  Cutoff radius (Angs.)
talt          = thermal_relaxation_time         !  Temperature coupling
talp          = pressure_relaxation_time        !  Pressure coupling 
KAPPA         = damping_Wolf                    !  Wolf's method damping paramenter (length^{-1}) ; (J. Chem. Phys. 1999; 110(17):8254)
read_vel      = read_velocities                 ! .T. , .F.
read_from_gmx = gmx_input_format                ! .T. , .F.

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

! ==============  pass information from structure to molecular dynamics  ==============
CALL Structure_2_MD( read_vel )

do i = 1 , MM % N_of_atoms
    nresid = atom(i) % nr
    molecule(nresid) % nr = nresid
end do

!----------------------
! finished reading ...
!----------------------

! get Atomic Number (AtNo) ...
CALL Symbol_2_AtNo( atom )

! Define atomic mass ...
atom % mass = Atomic_mass( atom % AtNo )
do i = 1 , MM % N_of_species
    species(i) % mass = sum( atom % mass , atom % my_species == i ) / species(i) % N_of_molecules
    where( molecule % my_species == i ) molecule % mass = species(i) % mass
end do

! mass of molecule species in Kg ...
species % mass = species % mass * imol
molecule % mass = molecule % mass * imol

! initial density of the box ... 
initial_density = sum( molecule % mass ) * MM % ibox(1) * MM % ibox(2) * MM % ibox(3) * (mts_2_nano**3)

10 if( ioerr > 0 ) stop "file_name file not found; terminating execution"
11 if( ioerr > 0 ) stop "input.pdb file not found; terminating execution"
12 if( ioerr > 0 ) stop "input.gro file not found; terminating execution"

99  format(a72)
100 format(t10, f6.3, t19, f6.3, t28, f6.3)
105 format(a6)
115 FORMAT(t8,i4,t14,a3,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3)

!=======================  reading  potential.inpt  ============================= 
CALL allocate_FF( atmax )

If( read_from_gmx ) then

    CALL ad_hoc_MM_tuning(instance="SpecialBonds")

    CALL itp2mdflex( MM , atom , species , FF)

    CALL top2mdflex( MM , atom , species , FF )

    do i = 1 , size(species)
        CALL MMSymbol_2_Symbol( species(i) % atom )
    end do

    CALL MMSymbol_2_Symbol ( FF )
    CALL Symbol_2_AtNo     ( FF )

    do i = 1 , size( Atomic_Mass )
        where( atom % AtNo == i ) atom % mass = Atomic_Mass(i)
        where( FF % AtNo == i ) FF % mass = Atomic_Mass(i)
    end do

    do i = 1 , MM % N_of_species
        do j = 1 , species(i) % N_of_atoms
            where( FF % residue == species(i) % atom(j) % residue )
                FF % nr         = species(i) % atom(j) % nr
                FF % flex       = species(i) % atom(j) % flex
                FF % my_species = species(i) % my_species
            end where
        end do
    end do

else

    open (30, file='potential.inpt', status='old')

    read(30,*) forcefield
    read(30,*) MM % CombinationRule

    select case( forcefield )

        case ( 1 )    
        !Born-Mayer potential ...

        case( 2 ) 
        ! L-J potential ...   
            do i = 1, atmax
                read (30,*) FF(i) % MMSymbol, FF(i) % MM_charge, FF(i) % sig, FF(i) % eps
                FF(i) % MMSymbol = adjustl( FF(i) % MMSymbol )
                ! factor1 = 1.0d26       <== factor used to not work with small numbers
                FF(i) % eps = FF(i) % eps * factor1 * imol
                FF(i) % eps = SQRT( FF(i) % eps )
                FF(i) % sig = ( FF(i) % sig * nano_2_angs )
                select case ( MM % CombinationRule )
                     case (2) 
                     FF(i) % sig = FF(i) % sig / TWO
                     case (3)
                     FF(i) % sig = sqrt ( FF(i) % sig )
                end select
              
                where( atom % MMSymbol == FF(i) % MMSymbol ) atom % eps       = FF(i) % eps
                where( atom % MMSymbol == FF(i) % MMSymbol ) atom % sig       = FF(i) % sig
                where( atom % MMSymbol == FF(i) % MMSymbol ) atom % MM_charge = FF(i) % MM_charge
                 
            end do
    end select

    CALL MMSymbol_2_Symbol( FF )

    CALL Symbol_2_AtNo( FF )

    do i = 1 , size( Atomic_Mass )
        where( atom % AtNo == i ) atom % mass = Atomic_Mass(i)
        where( FF % AtNo == i ) FF % mass = Atomic_Mass(i)
    end do

    do i = 1 , MM % N_of_species

        allocate( species(i) % atom( species(i) % N_of_atoms ) )

        do j = 1 , size(atom)

            if( species(i) % residue == atom(j) % residue ) then
                species(i) % atom(1:species(i) % N_of_atoms) % residue    = species(i) % residue
                species(i) % atom(1:species(i) % N_of_atoms) % my_species = species(i) % my_species
                species(i) % atom(1:species(i) % N_of_atoms) % my_id      = atom(j:j + species(i) % N_of_atoms) % my_id
                species(i) % atom(1:species(i) % N_of_atoms) % nr         = atom(j:j + species(i) % N_of_atoms) % nr
                species(i) % atom(1:species(i) % N_of_atoms) % flex       = atom(j:j + species(i) % N_of_atoms) % flex
                exit
            end if
            
        end do

    end do

    do k = 1 , MM % N_of_species
        do i = 1 , size(FF) - species(k) % N_of_atoms + 1
            do j = 1 , size(atom) - species(k) % N_of_atoms + 1

                l = 0
                do
                    if( i + l > size(FF) .OR. j + l > size(atom) ) exit
                    if( FF(i+l) % Symbol == atom(j+l) % Symbol ) then
                        l = l + 1
                    else
                        exit
                    end if
                end do

                if( l == species(k) % N_of_atoms ) then
                    FF(i:i+l-1) % my_species = atom(j:j+l-1) % my_species
                    FF(i:i+l-1) % residue    = atom(j:j+l-1) % residue
                    FF(i:i+l-1) % my_id      = atom(j:j+l-1) % my_id
                    FF(i:i+l-1) % nr         = atom(j:j+l-1) % nr
                    FF(i:i+l-1) % flex       = atom(j:j+l-1) % flex
                    exit
                end if

            end do

            if( l == species(k) % N_of_atoms ) exit

        end do
    end do

    ! Internal bondig parameters ...
    do a = 1 , MM % N_of_species

        if( species(a) % flex ) then
        ! species(a) is flexible ...

            ! Intramolecular nonbonding pairs 1-4 for species(a) ...
            read(30,*) species(a) % residue
            read(30,*,iostat = ioerr) species(a) % Nbonds14
            If( ioerr > 0 ) then
               string = species(a) % residue // '.inpt14'
               open (70, file=string, status='old',iostat=ioerr,err=101)
               read (70,*) species(a) % Nbonds14 
               allocate( species(a) % bonds14 ( species(a) % Nbonds14,2 ) )
               do b = 1 , species(a) % Nbonds14
                   read(70,*) species(a) % bonds14(b,1:2) , MM % fudgeLJ , MM % fudgeQQ
               end do
               close (70)
            else
               allocate( species(a) % bonds14 ( species(a) % Nbonds14,2 ) )
               do b = 1 , species(a) % Nbonds14
                   read(30,*) species(a) % bonds14(b,1:2) , MM % fudgeLJ , MM % fudgeQQ
               end do        
            end if

            101 if( ioerr > 0 )then
                print*, string,' file not found; terminating execution'; stop
            end if

            ! bond stretching pairs for species(a) ...
            read(30,*) species(a) % Nbonds
            allocate( species(a) % bonds  ( species(a) % Nbonds,2 ) )
            allocate( species(a) % kbond0 ( species(a) % Nbonds,2 ) )
            do b = 1 , species(a) % Nbonds
                read(30,*) species(a) % bonds(b,1:2) , (species(a) % kbond0(b,k) ,k=2,1,-1)
            end do
            ! factor used to not work with small numbers ...
            ! factor2 = 1.0d24
            species(a) % kbond0(:,1) = species(a) % kbond0(:,1) * factor2 * imol
            species(a) % kbond0(:,2) = species(a) % kbond0(:,2) * nano_2_angs

            ! bond angle pairs for species(a) ...
            read(30,*) species(a) % Nangs
            allocate( species(a) % angs  ( species(a) % Nangs,3 ) )
            allocate( species(a) % kang0 ( species(a) % Nangs,2 ) )
            do b = 1 , species(a) % Nangs
                read(30,*) species(a) % angs(b,1:3), (species(a) % kang0(b,k) ,k=2,1,-1)
            end do
            species(a) % kang0(:,2) = species(a) % kang0(:,2) * deg_2_rad
            ! factor used to not work with small numbers ...
            ! factor1 = 1.0d26
            species(a) % kang0(:,1) = species(a) % kang0(:,1) * factor1 * imol

            ! Dihedral Angle Potential for species(a) : general form ...
            read(30,*) Dihedral_Potential_Type
            select case( adjustl(Dihedral_Potential_Type) ) 
          
               ! factor1 = 1.0d26       <== factor used to not work with small numbers

               case ('cos') 
                    read(30,*) species(a) % Ndiheds
                    allocate( species(a) % diheds ( species(a) % Ndiheds,4 ) )
                    allocate( species(a) % kdihed0( species(a) % Ndiheds,2 ) )                  
                    do b = 1 , species(a) % Ndiheds
                       read(30,*) species(a) % diheds(b,1:4), species(a) % kdihed0(b,1:2), species(a) % Nharm
                    end do
                    species(a) % kdihed0(:,2) = species(a) % kdihed0(:,2) * deg_2_rad
                    species(a) % kdihed0(:,1) = species(a) % kdihed0(:,1) * factor1 * imol

               case ('harm')
                    read(30,*) species(a) % Ndiheds
                    allocate( species(a) % diheds( species(a) % Ndiheds,4  ) )
                    allocate( species(a) % kdihed0( species(a) % Ndiheds,2 ) )
                    do b = 1 , species(a) % Ndiheds
                       read(30,*) species(a) % diheds(b,1:4), species(a) % kdihed0(b,1:2)
                    end do
                    species(a) % kdihed0(:,2) = species(a) % kdihed0(:,2) * deg_2_rad
                    species(a) % kdihed0(:,1) = species(a) % kdihed0(:,1) * factor1 * imol
           
               case ('hcos')
                    read(30,*) species(a) % Ndiheds
                    allocate( species(a) % diheds( species(a) % Ndiheds,4  ) )
                    allocate( species(a) % kdihed0( species(a) % Ndiheds,2 ) )
                    do b = 1 , species(a) % Ndiheds
                       read(30,*) species(a) % diheds(b,1:4), species(a) % kdihed0(b,1:2)
                    end do
                    species(a) % kdihed0(:,2) = species(a) % kdihed0(:,2) * deg_2_rad
                    species(a) % kdihed0(:,1) = species(a) % kdihed0(:,1) * factor1 * imol
           
               case ('cos3')
                    read(30,*) species(a) % Ndiheds
                    allocate( species(a) % diheds( species(a) % Ndiheds,4  ) )
                    allocate( species(a) % kdihed0( species(a) % Ndiheds,3 ) )
                    do b = 1 , species(a) % Ndiheds
                       read(30,*) species(a) % diheds(b,1:4), species(a) % kdihed0(b,1:3)
                    end do
                    species(a) % kdihed0(:,1:3) = species(a) % kdihed0(:,1:3) * factor1 * imol
           
               case ('ryck')
                    read(30,*) species(a) % Ndiheds
                    allocate( species(a) % diheds( species(a) % Ndiheds,4  ) )
                    allocate( species(a) % kdihed0( species(a) % Ndiheds,6 ) )
                    do b = 1 , species(a) % Ndiheds
                       read(30,*) species(a) % diheds(b,1:4), species(a) % kdihed0(b,1:6)
                    end do
                    species(a) % kdihed0(:,1:6) = species(a) % kdihed0(:,1:6) * factor1 * imol
           
               case ('opls')
                    read(30,*) species(a) % Ndiheds
                    allocate( species(a) % diheds( species(a) % Ndiheds,4  ) )
                    allocate( species(a) % kdihed0( species(a) % Ndiheds,4 ) )
                    do b = 1 , species(a) % Ndiheds
                       read(30,*) species(a) % diheds(b,1:4), species(a) % kdihed0(b,1:4)
                    end do
                    species(a) % kdihed0(:,1:4) = species(a) % kdihed0(:,1:4) * factor1 * imol
            end select 

            allocate( species(a) % Dihedral_Type(species(a)%Ndiheds) , source = Dihedral_Potential_Type )

        else
            ! species(a) is NOT flexible ...    
            species(a) % Nbonds14  =  0
            species(a) % Nbonds    =  0
            species(a) % Nangs     =  0
            species(a) % Ndiheds   =  0
        endif
    end do

    close (30)

end If

!=======================  finished  reading  potential.inpt  ============================= 
do i = 1 , MM % N_of_species
    where( molecule % my_species == i ) molecule % N_of_atoms = species(i) % N_of_atoms
    where( molecule % my_species == i ) molecule % Nbonds     = species(i) % Nbonds
    where( molecule % my_species == i ) molecule % Nangs      = species(i) % Nangs
    where( molecule % my_species == i ) molecule % Ndiheds    = species(i) % Ndiheds
    where( molecule % my_species == i ) molecule % Nharm      = species(i) % Nharm
    where( molecule % my_species == i ) molecule % Nbonds14   = species(i) % Nbonds14
    where( molecule % my_species == i ) molecule % NintraLJ   = species(i) % NintraLJ
end do

do i = 1 , MM % N_of_molecules
    allocate( molecule(i) % bonds14       ( molecule(i) % Nbonds14 , 2 ) )
    allocate( molecule(i) % bonds         ( molecule(i) % Nbonds   , 2 ) )
    allocate( molecule(i) % kbond0        ( molecule(i) % Nbonds   , 2 ) )
    allocate( molecule(i) % angs          ( molecule(i) % Nangs    , 3 ) )
    allocate( molecule(i) % kang0         ( molecule(i) % Nangs    , 2 ) )
    allocate( molecule(i) % diheds        ( molecule(i) % Ndiheds  , 4 ) )
    allocate( molecule(i) % kdihed0       ( molecule(i) % Ndiheds  , 7 ) )
    allocate( molecule(i) % harm          ( molecule(i) % Ndiheds      ) )
    allocate( molecule(i) % Dihedral_Type ( molecule(i) % Ndiheds      ) )
    allocate( molecule(i) % IntraLJ       ( molecule(i) % NintraLJ , 2 ) )
end do

k = 0
do i = 1 , MM % N_of_molecules

    molecule(i) % kbond0        = species(molecule(i) % my_species) % kbond0
    molecule(i) % bonds         = species(molecule(i) % my_species) % bonds + k 
    molecule(i) % kang0         = species(molecule(i) % my_species) % kang0
    molecule(i) % angs          = species(molecule(i) % my_species) % angs + k
    molecule(i) % kdihed0       = species(molecule(i) % my_species) % kdihed0
    molecule(i) % harm          = species(molecule(i) % my_species) % harm
    molecule(i) % diheds        = species(molecule(i) % my_species) % diheds + k
    molecule(i) % bonds14       = species(molecule(i) % my_species) % bonds14 + k
    molecule(i) % IntraLJ       = species(molecule(i) % my_species) % IntraLJ + k
    molecule(i) % Dihedral_Type = species(molecule(i) % my_species) % Dihedral_Type

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

do i = 1 , size(atom)
    atom(i) % my_intra_id = i + molecule( atom(i) % nr ) % N_of_atoms - sum( molecule(1:atom(i) % nr) % N_of_atoms )
end do

CALL ad_hoc_MM_tuning( atom , instance = "General" )

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
    molecule(i) % cm(:)          = 0.0d0
    molecule(i) % mass           = 0.0d0
    molecule(i) % flex           = .false.
    molecule(i) % residue        = "XXX"
    molecule(i) % nr             = 0
    molecule(i) % Nbonds         = 0
    molecule(i) % bonds          = 0.0d0
    molecule(i) % Nangs          = 0
    molecule(i) % Ndiheds        = 0
    molecule(i) % Nharm          = 0
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
    FF(i) % ftotal(:)    = 0.0d0
    FF(i) % fch(:)       = 0.0d0
    FF(i) % fsr(:)       = 0.0d0
    FF(i) % mass         = 0.0d0
    FF(i) % charge       = 0.0d0
    FF(i) % MM_charge    = 0.0d0
    FF(i) % eps          = 0.0d0
    FF(i) % sig          = 0.0d0
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
        case default
            print*, ' >> Symbol_2_AtNo ; unknown atom found <<' , '[',a(i)%Symbol,']' , i
            stop
    end select

 END DO

 end subroutine Symbol_2_AtNo
!
!
!
!====================================
subroutine Structure_2_MD( read_vel )
!====================================
implicit none
logical , intent(in) :: read_vel

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
    atom(i) % ftotal(:)    = 0.d0
    atom(i) % fch(:)       = 0.d0
    atom(i) % fsr(:)       = 0.d0
    atom(i) % mass         = 0.d0
    atom(i) % charge       = 0.d0
    atom(i) % MM_charge    = 0.d0
    atom(i) % eps          = 0.d0
    atom(i) % sig          = 0.d0
    atom(i) % flex         = .true.
end do

! this is assumed a priori , but can be changed otherwise if required by the Force Field ...
atom % MMSymbol = atom % EHSymbol  

CALL ad_hoc_MM_tuning( atom , instance = "General" )

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

if( read_vel ) then

    if( restart ) STOP ' Should NOT read velocity_MM.inpt if restart = .true. !'

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
    end select

END DO

a % Symbol = adjustl(a % Symbol)

end subroutine MMSymbol_2_Symbol
!
!
!
end  module MD_read_m

module tuning_m

    use type_m
    use constants_m
    use card_reading       , only : ReadInputCard_ADHOC , electron_fragment , hole_fragment , solvent_QM_droplet_radius
    use parameters_m       , only : static , electron_state , hole_state , n_part , Survival , ad_hoc_droplet

    public :: ad_hoc_tuning , eh_tag , orbital 

    private

    ! module variables ...
    integer      , allocatable :: orbital(:)
    character(2) , allocatable :: eh_tag(:)
    logical                    :: done = .false.

    ! module parameters ...
    logical, parameter :: T_ = .true. , F_ = .false.

    contains
!
!
!
!================================
 subroutine ad_hoc_tuning( univ )
!================================
implicit none
type(universe) , intent(inout) :: univ

!local variables ...
integer :: i
logical :: exist

! edit structure  ...

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      univ % atom(:) % etc
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inquire( file=dynemolworkdir//"makefile" , EXIST = exist )

if( (.NOT. exist) .AND. (.NOT. done) ) then

    call ReadInputCard_ADHOC( structure=univ )
    done = .true.

else

    !-----------------------------------
    !      define %atom
    !-----------------------------------
    
    !-----------------------------------
    !      define %residue
    !-----------------------------------
    
    !-----------------------------------
    !      define %nr
    !-----------------------------------
    
    !---------------------------------------------------
    !      define %DPF     (Dipole Fragment) 
    !      define %V_shift (FMO offset shift)
    !---------------------------------------------------
    
    !---------------------------------------------------
    !      define %QMMM  
    !      default is QMMM = QM; set QMMM = MM for classical atoms ... 
    !---------------------------------------------------

end if

!---------------------------------------------------
!      define %El   : mandatory !!
!---------------------------------------------------
where(univ % atom % residue == electron_fragment ) univ % atom % El = .true.
!---------------------------------------------------
!      define %Hl   : must be T_ for El/Hl calcs ...
!---------------------------------------------------
where(univ % atom % residue == hole_fragment) univ % atom % Hl = .true.
!----------------------------------------------------
!      define %fragment 
!----------------------------------------------------

! ---------- Default fragments -------------
!   Acceptor    =   A       
!   Donor       =   D 
!   Molecule    =   M
!   Bridge      =   B
!   Solvent     =   S  
!-------------------------------------------

DO i = 1 , size(univ%atom)

    select case(univ%atom(i)%residue)

        ! this is water, this is water, this is water ...
        ! consider RDF(H-O) ~ 1.85A , RDF(H-H) ~ 2.35A, RDF(O-O) ~ 2.75A
        ! radius of solvation hardcore = radius of first solvation shell; RDF_2(H-O) ~ 3.25A
        ! for r < R_hardcore, eHuckel applies; for r > R_hardcore nonbonding hamiltonian applies
        case( 'H2O' , 'WAT' , 'TIP' )
            univ%atom(i)%fragment = 'S'
            univ%atom(i)%solvation_hardcore = 3.d0

        ! this is acetonitrile ...
        case( 'ACN' )
            univ%atom(i)%fragment = 'S'
            univ%atom(i)%solvation_hardcore = 3.d0

        ! case default ...
        case( 'SOL' )
            univ%atom(i)%fragment = 'S'
            univ%atom(i)%solvation_hardcore = 3.d0

    end select

END DO

if(ad_hoc_droplet) call QM_droplet( univ )

call warnings( univ%atom ) 

Print 46

include 'formats.h'

end subroutine ad_hoc_tuning
!
!
!
!=============================
 subroutine QM_droplet( univ )
!=============================
implicit none
type(universe) , intent(inout) :: univ

!local variables ...
integer :: i, nr , nr_max
real*8  :: distance
real*8  :: solvent_CG(3) , solute_CG(3)

call ReGroup_Molecules(univ)

nr_max =  maxval(univ%atom(:)%nr)

if( any(univ % atom % El) ) &
then
       univ % atom % solute = univ % atom % El 
else
       where( univ % atom % nr == 1 ) univ % atom % solute = .true. 
end if

! identify the CG of the solute ...
forall( i=1:3 ) solute_CG(i) = sum( univ%atom%xyz(i) , univ%atom%solute == .true. ) / count(univ%atom%solute)

do nr = 1 , nr_max

      forall( i=1:3 ) solvent_CG(i) = sum( univ%atom%xyz(i) , univ%atom%nr == nr ) / count(univ%atom%nr==nr)
      
      distance = sqrt( (solute_CG(1) - solvent_CG(1))**2 + &
                       (solute_CG(2) - solvent_CG(2))**2 + &
                       (solute_CG(3) - solvent_CG(3))**2   )

      if( distance <= solvent_QM_droplet_radius ) &
      then
          where( univ%atom%nr == nr ) 
               univ%atom%QMMM    = "QM"
               univ%atom%fragment = "Q"
          end where
      end if
end do

end subroutine QM_droplet
!
!=================================
subroutine ReGroup_Molecules(univ)
!=================================
implicit none
type(universe) , intent(inout) :: univ

!local variables
integer :: nr , i , j , indx1 , indx2 
real*8  :: delta(3) 

do nr = minval(univ%atom%nr) , maxval(univ%atom%nr)

    ! atomic pointers for molecule with nresidue = nr
    indx1 = minval( [(i , i=1,size(univ%atom))] , (univ%atom%nr == nr) )
    indx2 = maxval( [(i , i=1,size(univ%atom))] , (univ%atom%nr == nr) )

    do i = indx1+1 , indx2

        delta = univ%atom(i)%xyz - univ%atom(indx1)%xyz

        do j = 1 , 3
             if( abs(delta(j)) > univ%box(j)*HALF ) univ%atom(i)%xyz(j) = univ%atom(i)%xyz(j) - sign( univ%box(j) , delta(j) )
        end do 
    end do
end do

end subroutine ReGroup_Molecules
!
!
!
!
!===========================
 subroutine warnings( a )
!===========================
implicit none
type(atomic) , intent(inout) :: a(:)

!local variables ...
logical :: propagate_el , propagate_hl

! orbitals to be propagated ...
If( Survival ) then

   allocate( orbital(n_part) , source = 0    )
   allocate( eh_tag(n_part)  , source = "XX" )

   If( any(a% El) ) then
       orbital(1) = electron_state ; eh_tag(1) = "el"
   End If
   If( any(a% Hl) ) then
       orbital(2) = hole_state     ; eh_tag(2) = "hl"
   End If

end If

!default: %El => DONOR
If( any(a% El) ) then      ! <== first priority ...
    where( a% El ) a% fragment = "D"
else If( any(a% Hl) ) then ! <== second priority, only for cationic systems ...
    where( a% Hl ) a% fragment = "A"
else If(.NOT. static .AND. electron_state /= 0 ) then
        CALL warning("execution stopped, must define eletron ...%El in ad_hoc_tuning; is ad_hoc = T_?")
        stop 
end If

propagate_el = any(a% El) .EQV. (electron_state /= 0)
If( .not. propagate_el ) then
        CALL warning("execution stopped, ELECTRON wavepacket setup is not consistent: check electron_state (parameters.f) and %El (tuning.f)")
        stop 
end If

propagate_hl = any(a% Hl) .EQV. (hole_state /= 0)
If( .not. propagate_hl ) then
        CALL warning("execution stopped, HOLE wavepacket setup is not consistent: check hole_state (parameters.f) and %Hl (tuning.f)")
        stop 
end If

end subroutine warnings
!
!
!
end module tuning_m
!
!
!
!
module MM_tuning_routines

    use atomicmass
    use constants_m     , only: large
    use type_m          , only: dynemolworkdir
    use parameters_m    , only: static
    use card_reading    , only: ReadInputCard_ADHOC 
    use MM_types        , only: MM_atomic, DefineBonds, DefineAngles

    private

    public :: ad_hoc_MM_tuning , SpecialBonds, SpecialAngs

    ! module variables ...
    type(DefineBonds) , allocatable :: SpecialBonds(:)
    type(DefineAngles), allocatable :: SpecialAngs(:)
    logical :: done = .false.

    ! module parameters ...
    logical, parameter :: T_ = .true. , F_ = .false.

    contains

!================================================
 subroutine ad_hoc_MM_tuning( atom , instance )
!================================================
implicit none
type(MM_atomic) , optional , intent(inout) :: atom(:)
character(*)               , intent(in)    :: instance

! local variables ...
logical :: exist 

inquire( file=dynemolworkdir//"makefile" , EXIST = exist )
if( (.NOT. exist) .AND. (.NOT. done) ) then
    call ReadInputCard_ADHOC( atom=atom )
    done = .true.
    return
endif

! edit structure  ...

select case ( instance ) 

!==================================
    case ("General")

!----------------------------------
!      define SPECIAL atoms 
!----------------------------------

!----------------------------------
!      define flex
!----------------------------------

!----------------------------------
!      define MM atom types 
!----------------------------------

!----------------------------------
!      define my_species
!----------------------------------

!----------------------------------
!      define residues
!----------------------------------

!----------------------------------
!      define nr
!----------------------------------

!----------------------------------
!       charge of the atoms 
!----------------------------------


!----------------------------------
!     Selective_Dynamics 
!----------------------------------
    case ("MegaMass")
         where( atom% mass < 0.0 ) 
             atom% mass = Atomic_mass( atom % AtNo ) * 1.d3
         elsewhere
             atom% mass = Atomic_mass( atom % AtNo )
         end where

!----------------------------------
!      define SPECIAL bonds
!----------------------------------
    case( 'SpecialBonds' )

!----------------------------------
!      define SPECIAL angles 
!----------------------------------

!=====================================

end select

end subroutine ad_hoc_MM_tuning

end module MM_tuning_routines
!
!
!
module syst
implicit none
real*8  :: bath_T, press, tau_temp, tau_p, Initial_density

type kind_of_ensemble
     logical :: inter = .false.
     logical :: intra = .false.
     logical :: anyone = .false.
end  type kind_of_ensemble

type(kind_of_ensemble) :: using_barostat

end module syst
!
!
!
module for_force
    implicit none

    integer                               :: forcefield
    real*8, dimension(:,:)  , allocatable :: vscut, fscut
    real*8                                :: rcut, vrecut, frecut, rcut2, KAPPA
    character(4)                          :: Dihedral_Potential_Type
 
    real*8 :: Coul_inter = 0.d0 
    real*8 :: Vself      = 0.d0
    real*8 :: evdw       = 0.d0
    real*8 :: bdpot      = 0.d0
    real*8 :: harm_bond  = 0.d0
    real*8 :: morse_bond = 0.d0
    real*8 :: Morspot    = 0.d0
    real*8 :: angpot     = 0.d0
    real*8 :: dihpot     = 0.d0
    real*8 :: proper_dih = 0.d0
    real*8 :: ryck_dih   = 0.d0
    real*8 :: harm_dih   = 0.d0
    real*8 :: imp_dih    = 0.d0
    real*8 :: LJ_14      = 0.d0
    real*8 :: LJ_intra   = 0.d0
    real*8 :: Coul_14    = 0.d0
    real*8 :: Coul_intra = 0.d0
    real*8 :: DWFF_intra = 0.d0
    real*8 :: DWFF_inter = 0.d0
    real*8 :: pot_INTER  = 0.d0
    real*8 :: pot_total  = 0.d0

end module for_force

module FF_OPT_class_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use parameters_m            , only : PBC
    use MD_read_m               , only : atom , molecule , MM 
    use MM_types                , only : MM_system , MM_atomic , MM_molecular , debug_MM
    use MM_input                , only : driver_MM
    use F_intra_m               , only : ForceIntra, Pot_Intra                                     
    use OPT_Parent_class_m      , only : OPT_Parent
    use NonlinearSidekick_m     , only : Fletcher_Reeves_Polak_Ribiere_minimization

    public :: FF_OPT , LogicalKey

    private

    type, extends(OPT_Parent)   :: FF_OPT
        integer                 :: ITMAX_FF = 100           ! <== 100-300 is a good compromise of accuracy and safety
        real*8                  :: BracketSize_FF = 1.d-4   ! <== this value may vary between 1.0d-3 and 1.0d-4
        logical                 :: profiling_FF = .FALSE.
    contains
        procedure :: cost => cost_evaluator
        procedure :: cost_variation => Generalized_Forces
        procedure :: output => dump_parameters
    end type FF_OPT

    interface FF_OPT
        module procedure  constructor
    end interface

    ! module types ...
    type LogicalKey
        logical :: bonds(3)
        logical :: angs(2)
        logical :: diheds(7)
    end type

    ! module variables ...
    integer                                     :: bonds , angs , diheds
    logical             , allocatable           :: bonds_mask(:,:) , angs_mask(:,:) , diheds_mask(:,:)
    real*8              , allocatable , target  :: bond_target(:,:) , ang_target(:,:) , dihed_target(:,:)
    type(real_pointer)  , allocatable           :: bond(:,:) , ang(:,:) , dihed(:,:)
    type(MM_atomic)     , allocatable           :: atom0(:)
    character(len=11)                           :: method
    logical                                     :: done = .false.

contains
!
!
!
!=========================================================
 function constructor( key , kernel , ref ) result( me )
!=========================================================
implicit none
type(LogicalKey)            , intent(in) :: key
character(*)                , intent(in) :: kernel
type(MM_atomic)  , optional , intent(in) :: ref(:)
type(FF_OPT)                             :: me 

!local variable ...
integer :: i , j 

method = kernel
if( present(ref) ) allocate( atom0 , source = ref )

! Cause you are my kind / You're all that I want ...
me % ITMAX       = me % ITMAX_FF
me % BracketSize = me % BracketSize_FF
me % profiling   = me % profiling_FF
if( kernel == "NormalModes" ) me % BracketSize = 100000*me % BracketSize
if( kernel == "NormalModes" ) me % itmax = 20
me % driver = driver_MM

If( driver_MM == "Parametrize" ) me % profiling = .true.

if( .NOT. done ) then
    ! catch the distinct parameters and make (bond, ang, dihed) point to them ...
    allocate(bond_target , source=Select_Different( molecule(1)%kbond0                 ) )
    allocate(bond        , source=Define_Pointer  ( molecule(1)%kbond0  , bond_target  ) )

    allocate(ang_target  , source=Select_Different( molecule(1)%kang0                  ) )
    allocate(ang         , source=Define_Pointer  ( molecule(1)%kang0   , ang_target   ) )

    allocate(dihed_target, source=Select_Different( molecule(1)%kdihed0                ) )
    allocate(dihed       , source=Define_Pointer  ( molecule(1)%kdihed0 , dihed_target ) )

    ! the parameters = zero are not optimized ...
    allocate( bonds_mask  ( size(bond_target (:,1)) , size(bond_target (1,:)) ) , source = (abs(bond_target)  > low_prec) )
    allocate( angs_mask   ( size(ang_target  (:,1)) , size(ang_target  (1,:)) ) , source = (abs(ang_target)   > low_prec) )
    allocate( diheds_mask ( size(dihed_target(:,1)) , size(dihed_target(1,:)) ) , source = (abs(dihed_target) > low_prec) )
    done = .true.
end if

! logical key to set the optimization routine ...
forall(j=1:3) bonds_mask(:,j)  =  key%bonds(j)
forall(j=1:2) angs_mask(:,j)   =  key%angs(j)
forall(j=1:7) diheds_mask(:,j) =  key%diheds(j)

! number of degrees of freedom in optimization space ...
bonds  = count( bonds_mask  ) 
angs   = count( angs_mask   ) 
diheds = count( diheds_mask )

me % N_of_Freedom = bonds + angs + diheds

allocate( me % p( me % N_of_Freedom ) , source = D_zero )

me % p( 1 :            ) = pack( bond_target  , bonds_mask  , me%p )
me % p( bonds+1 :      ) = pack( ang_target   , angs_mask   , me%p ) 
me % p( bonds+angs+1 : ) = pack( dihed_target , diheds_mask , me%p ) 

end function constructor
!
!
!
!==========================================
 function cost_evaluator( me ) result(cost)
!==========================================
implicit none
class(FF_OPT) , intent(inout)  :: me
real*8                         :: cost

! local variables ...
integer         :: i , j  
real*8          :: energy 
type(R_eigen)   :: Hesse
    
! reset forces ...
forall( i=1:3 ) atom(:) % ftotal(i) = D_zero

bond_target  = unpack( me%p(:bonds)             , bonds_mask  , bond_target  )
ang_target   = unpack( me%p(bonds+1:bonds+angs) , angs_mask   , ang_target   )
dihed_target = unpack( me%p(bonds+angs+1:)      , diheds_mask , dihed_target )

forall( i=1:molecule(1)%nbonds  , j=1:size(molecule(1)%kbond0 (1,:)) ) molecule(1)%kbond0(i,j)  =  bond(i,j)%PTR

forall( i=1:molecule(1)%nangs   , j=1:size(molecule(1)%kang0  (1,:)) ) molecule(1)%kang0(i,j)   =  ang (i,j)%PTR

forall( i=1:molecule(1)%ndiheds , j=1:size(molecule(1)%kdihed0(1,:)) ) molecule(1)%kdihed0(i,j) =  dihed(i,j)%PTR

select case ( method )

    case( "energy" )

        CALL ForceIntra
        Energy = Pot_Intra * mol * micro / MM % N_of_molecules

        cost   = Energy

    case( "NormalModes" )

        CALL normal_modes( Hesse )

        cost = sqrt( (hesse%erg(10)-116.d0)**2 + (hesse%erg(14)-258.d0)**2 ) 

    case default 

        stop

end select

end function cost_evaluator
!
!
!
!
!============================================
 subroutine Generalized_Forces( me , vector )
!============================================
implicit none
class(FF_OPT)   , intent(in)    :: me
real*8          , intent(inout) :: vector(:)

! local parameters ...
real*8  , parameter :: small = 1.d-4

! local variables ...
integer         :: i 
type(FF_OPT)    :: before , after

do i = 1 , me % N_of_Freedom

    after  = me
    before = me

    after  % p(i) = me % p(i) * (D_one + small)
    before % p(i) = me % p(i) * (D_one - small)

    vector(i) = ( after%cost() - before%cost() ) / (two*small*me%p(i))

end do
 
end subroutine Generalized_Forces
!
!
!
!======================================
 subroutine dump_parameters( me , iter)
!======================================
implicit none
class(FF_OPT)   , intent(in) :: me
integer         , optional , intent(in) :: iter

! local variables ...
integer                         :: i , at1 , at2 , at3 , at4 , funct_dih
character(3)                    :: funct_type
real*8                          :: factor
character(len=:) , allocatable  :: string(:)

 open( unit = 51 , file = "topol-OPT.top" , status = "replace", action = "write" , position = "append" )

 !========================================================================================================
 ! bond parms saving ...
 write(51,"(A)") "[ bondtypes ]"               
 write(51,"(A)") "; Optimized by OOP: flexible inheritance of objects for nonlinear optimization"

 allocate( character(len=2*len(atom(at1)%MMSymbol)) :: string(molecule(1)%Nbonds) )
 do i = 1 , molecule(1)%Nbonds

    at1 = molecule(1)%bonds(i,1)
    at2 = molecule(1)%bonds(i,2)

    string(i) = atom(at1)%MMSymbol//atom(at2)%MMSymbol

    if( .NOT. any(string(1:i-1) == string(i)) ) then 

        funct_type = molecule(1) % funct_bond(i) 
 
        factor = factor2 * imol  
        if( funct_type == "3" ) factor = factor1 * imol

        write(51,'(3A4,F15.5,2F15.3)')  atom(at1)%MMSymbol                      , &
                                        atom(at2)%MMSymbol                      , &
                                        funct_type                              , &
                                        molecule(1)%kbond0(i,2) / nano_2_angs   , &
                                        molecule(1)%kbond0(i,1) / factor        , &
                                        molecule(1)%kbond0(i,3) * nano_2_angs
    end if
 end do
 deallocate(string)
 !========================================================================================================
 ! angle parms saving ...
 write(51,*) " "
 write(51,"(A)") "[ angletypes ]"
 write(51,"(A)") "; Optimized by OOP: flexible inheritance of objects for nonlinear optimization"

 allocate( character(len=3*len(atom(at1)%MMSymbol)) :: string(molecule(1)%Nangs) )
 do i = 1 , molecule(1)%Nangs

    at1 = molecule(1)%angs(i,1)
    at2 = molecule(1)%angs(i,2)
    at3 = molecule(1)%angs(i,3)

    string(i) = atom(at1)%MMSymbol//atom(at2)%MMSymbol//atom(at3)%MMSymbol

    if( .NOT. any(string(1:i-1) == string(i)) ) then 

        funct_type = molecule(1) % funct_angle(i)

        factor = factor1 * imol

        write(51,'(4A4,2F15.3)') atom(at1)%MMSymbol                 , &
                                atom(at2)%MMSymbol                  , &
                                atom(at3)%MMSymbol                  , &
                                funct_type                          , &
                                molecule(1)%kang0(i,2) / deg_2_rad  , &
                                molecule(1)%kang0(i,1) / factor      
    end if
 end do
 deallocate(string)
 !========================================================================================================
 ! dihedral parms saving ...
 write(51,*) " "
 write(51,"(A)") "[ dihedraltypes ]"
 write(51,"(A)") "; Optimized by OOP: flexible inheritance of objects for nonlinear optimization"

 allocate( character(len=4*len(atom(at1)%MMSymbol)) :: string(molecule(1)%Ndiheds) )
 do i = 1 , molecule(1)%Ndiheds

    at1 = molecule(1)%diheds(i,1)
    at2 = molecule(1)%diheds(i,2)
    at3 = molecule(1)%diheds(i,3)
    at4 = molecule(1)%diheds(i,4)

    string(i) = atom(at1)%MMSymbol//atom(at2)%MMSymbol//atom(at3)%MMSymbol//atom(at4)%MMSymbol

    if( .NOT. any(string(1:i-1) == string(i)) ) then 

        funct_dih = molecule(1) % funct_dihed(i)

        factor = factor1 * imol

        write(51,'(4A4,I5,2F15.5)') atom(at1)%MMSymbol                  , &
                                    atom(at2)%MMSymbol                  , &
                                    atom(at3)%MMSymbol                  , &
                                    atom(at4)%MMSymbol                  , &
                                    funct_dih                           , &
                                    molecule(1)%kdihed0(i,1) / factor   , &
                                    molecule(1)%kdihed0(i,2) / factor   
    end if
 end  do
 deallocate(string)
 !========================================================================================================

 close(51)

end subroutine dump_parameters
!
!
!
!==========================================
 function Select_Different( a ) result( b )
!==========================================
implicit none
real*8        , intent(in)  :: a(:,:)
real*8        , allocatable :: b(:,:)

! local variables ...
integer                 :: j , k , diff , diff_size
integer , allocatable   ::  indx(:)
real*8  , allocatable   ::  tmp(:,:)

allocate( tmp  (size(a(:,1)),size(a(1,:))), source=D_zero )
allocate( indx (size(a(:,1))             ), source=I_zero )

! select different elements of a ...
diff_size = 1
indx(1)   = 1
tmp(1,:)  = a(1,:)

do j = 2 , size(a(:,1))
    
    do k = 1 , diff_size
        if( all(tmp(k,:) == a(j,:)) )  exit
    end do

    if( k > diff_size ) then
        tmp(k,:)  = a(j,:)
        indx(k)   = j
        diff_size = k
    end if

end do

allocate( b , source=tmp(1:diff_size,:) )

deallocate( tmp , indx )

end function Select_Different
!
!
!
!============================================
 function Define_Pointer( a , b ) result( c )
!============================================
implicit none
real*8 , intent(in)                 :: a(:,:)
real*8 , intent(in) , target        :: b(:,:)
type(real_pointer)  , allocatable   :: c(:,:) 

! local variables ...
integer :: i , j , k 

allocate( c(size(a(:,1)),size(a(1,:))) )

! creates a copy of a that points to b, i.e., c=a and c => b ...
do i = 1 , size(b(:,1))
    do j = 1  , size( a(:,1) )
        if( all(a(j,:) == b(i,:)) ) forall( k=1:size(a(1,:)) ) c(j,k)%PTR => b(i,k)
    end do
end do

end function Define_Pointer
!
!
!
!================================
 subroutine normal_modes( Hesse )
!================================
use MM_ERG_class_m  , only    : MM_OPT

implicit none
type(R_eigen)   , intent(out) :: Hesse

! local variables ...
integer                       :: i , j , k , l , column , info
real*8                        :: local_energy_minimum
real*8          , allocatable :: Hessian(:,:)
type(MM_atomic) , allocatable :: equilibrium(:) , atom_fwd(:) , atom_bwd(:)
type(MM_OPT)                  :: MM_erg

!local parameters ...
real*8 , parameter :: delta         = 1.d-5             ! <== displacement in Angs.
real*8 , parameter :: eV_2_cm_inv   = 1.d-12*8065.73    ! <== displacement in Angs.

! start the normal mode calculations from an energy minimum ...

! instantiating MM ...
MM_erg = MM_OPT( )
CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_erg , MM_erg%N_of_Freedom , local_energy_minimum )

! allocate all the stuff once ...
if( .not. allocated(Hessian) ) then
    allocate( equilibrium (   MM % N_of_atoms                     ) )
    allocate( atom_fwd    (   MM % N_of_atoms                     ) )
    allocate( atom_bwd    (   MM % N_of_atoms                     ) )
    allocate( Hessian     ( 3*MM % N_of_atoms , 3*MM % N_of_atoms ) )
    allocate( Hesse%erg   ( 3*MM % N_of_atoms                     ) )
end if

do i = 1 , 3
    equilibrium % xyz(i) = MM_erg % p( (i-1) * MM%N_of_atoms + 1 : i * MM%N_of_atoms ) 
end do

! reset the atomic forces ...
forall(i=1:3) atom % ftotal(i) = D_zero

! start build up of Hessian matrix ...
column = 1
do k = 1 , MM%N_of_atoms

    do j = 1 , 3

        atom(k) % xyz(j) = equilibrium(k) % xyz(j) + delta
        forall(i=1:3) atom % ftotal(i) = D_zero
        CALL ForceIntra
        forall(i=1:3) atom_fwd % ftotal(i) = atom % ftotal(i)

        atom(k) % xyz(j) = equilibrium(k) % xyz(j) - delta
        forall(i=1:3) atom % ftotal(i) = D_zero
        CALL ForceIntra
        forall(i=1:3) atom_bwd % ftotal(i) = atom % ftotal(i)

        do i = 1 , MM%N_of_atoms
            do l = 1 , 3
                Hessian( (i-1)*3+l , column ) = - (atom_fwd(i) % ftotal(l) - atom_bwd(i) % ftotal(l)) / (two*delta)
            end do
        end do

        column = column + 1

    end do

end do

! atom%mass * imol = atomic mass in kg ; this is the unit of mass used in the verlet equations ...
forall( i=1:MM%N_of_atoms , l=1:3 ) Hessian( (i-1)*3 + l , : ) = Hessian( (i-1)*3 + l , : ) / sqrt( atom(i)%mass*imol )
forall( j=1:MM%N_of_atoms , k=1:3 ) Hessian( : , (j-1)*3 + k ) = Hessian( : , (j-1)*3 + k ) / sqrt( atom(j)%mass*imol )

! fixing correct units ...
Hessian = Hessian / Angs_2_mts

CALL SYEV( Hessian , Hesse % erg , 'V' , 'U' , info )
If ( info /= 0 ) write(*,*) 'info = ',info,' in SYEV in vibes normal modes '

! transforming back the normal modes: A --> M^{-1/2}*A ...
forall( i=1:MM%N_of_atoms , l=1:3 ) Hessian( (i-1)*3 + l , : ) = Hessian( (i-1)*3 + l , : ) / sqrt(atom(i)%mass)

! convert units to cm^{-1} ...
hesse%erg = sqrt( abs(hesse%erg) ) * h_bar*ev_2_cm_inv

end subroutine normal_modes
!
!
!
!=============
 function RMSD 
!=============
implicit none
real*8  :: RMSD

!local variable ...
integer  :: i
real*8   :: rij(3) , distance

distance = D_zero
do i = 1 , size(atom)
    rij(:)   = atom(i) % xyz(:) - atom0(i) % xyz(:)
    rij(:)   = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
    RMSD = RMSD + SQRT( rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3) )
end do

end function RMSD
!
!
!
end module FF_OPT_class_m

module diagnostic_m

 use type_m
 use omp_lib
 use constants_m
 use DOS_m
 use parameters_m            , only : spectrum , DP_Moment , &
                                      survival , EnvField_ , &
                                      Alpha_Tensor ,         &
                                      GaussianCube ,         &
                                      HFP_Forces ,           &
                                      DOS_range , sigma
 use Solvated_M              , only : DeAllocate_TDOS ,      &
                                      DeAllocate_PDOS ,      &
                                      DeAllocate_SPEC 
 use QCModel_Huckel          , only : EigenSystem
 use Babel_routines_m        , only : TO_UPPER_CASE
 use Structure_Builder       , only : Unit_Cell ,            &
                                      Extended_Cell ,        &
                                      Generate_Structure ,   &
                                      Basis_Builder 
 use GA_QCModel_m            , only : Mulliken
 use DP_main_m               , only : Dipole_Matrix
 use Dielectric_Potential    , only : Environment_SetUp
 use Oscillator_m            , only : Optical_Transitions
 use Psi_squared_cube_format , only : Gaussian_Cube_Format
 use Data_Output             , only : Dump_stuff
 use Embedded_FF_Alpha       , only : AlphaPolar
 use HuckelForces_m          , only : HuckelForces

 public :: diagnostic

 private

 !module variables ...
 character(8) :: token
 logical      :: atomicPDOS_called  = .false.
 logical      :: orbitalPDOS_called = .false.

 contains
!
!
!
!====================
subroutine diagnostic
!====================
implicit none 

! local variables ...
 integer                        :: i , nr , N_of_residues
 integer         , allocatable  :: MOnum(:)
 real*8                         :: DP(3) 
 type(R_eigen)                  :: UNI
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(STO_basis) , allocatable  :: ExCell_basis(:)

 
! preprocessing stuff ...................................

IF ( survival ) pause " >>> quit: diagnostic driver does not carry q_dynamics calculations <<< "

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

! reading command line arguments for plotting MO cube files ...
CALL Read_Command_Lines_Arguments( MOnum )

!.........................................................

 CALL Generate_Structure(1)

 CALL Basis_Builder( Extended_Cell, ExCell_basis )

 If( any([DP_Moment,Spectrum,EnvField_]) ) CALL Dipole_Matrix( Extended_Cell, ExCell_basis, UNI%L, UNI%R , DP )

 If( EnvField_ ) CALL Environment_SetUp( Extended_Cell )

 If( Alpha_Tensor .AND. DP_Moment ) CALL AlphaPolar( Extended_Cell, ExCell_basis ) 

 CALL EigenSystem( Extended_Cell, ExCell_basis, UNI )

 CALL Total_DOS( UNI%erg , TDOS )

 if( atomicPDOS_called  ) CALL atomicPDOS ( Extended_Cell, ExCell_basis, UNI )   
 if( orbitalPDOS_called ) CALL orbitalPDOS( Extended_Cell, ExCell_basis, UNI , MOnum)   
 
 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( Spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( HFP_Forces ) CALL HuckelForces( Extended_Cell, ExCell_basis, UNI )

 Print*, " " 
 Print*, "dE1 = ", UNI%erg(32) - UNI%erg(31) ,   "L-H"
 Print*, "dE2 = ", UNI%erg(32) - UNI%erg(30) ,   "L-H-1"
 Print*, "dE3 = ", UNI%erg(32) - UNI%erg(29) ,   "L-H-2"
 Print*, "dE4 = ", UNI%erg(32) - UNI%erg(28) ,   "L-H-3"
 Print*, "dE5 = ", UNI%erg(33) - UNI%erg(31) ,   "L+1-H"
 Print*, "dE6 = ", UNI%erg(34) - UNI%erg(31) ,   "L+2-H"
 
 If( GaussianCube .AND. (size(MOnum) > 0) ) then

     ! use this to check the orbitals separately ... 
     ! Extended_Cell% coord = Extended_Cell% coord * TWO

     do i = 1 , size(MOnum)
         CALL Gaussian_Cube_Format( UNI%L(MOnum(i),:) , UNI%R(:,MOnum(i)) , MOnum(i) , 0.d0 )
     end do
     Print 220 , MOnum(:)

 end if

 CALL Dump_stuff( TDOS , PDOS , SPEC )

 CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
 CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
 CALL DeAllocate_SPEC( SPEC , flag="dealloc" )

 include 'formats.h'

end subroutine diagnostic
!
!
!
!=========================================
subroutine atomicPDOS( sys , basis , UNI )
!=========================================
implicit none
type(structure) , intent(in) :: sys
type(STO_basis) , intent(in) :: basis(:)
type(R_eigen)   , intent(in) :: UNI

! local variables ...
type(f_grid)         :: PDOS_atomic
integer              :: i , j , k , i1 , i2 , n_of_DOS_states
real*8               :: gauss_norm , two_sigma2 , step , delta_E
real*8 , allocatable :: PDOS(:,:) , erg_MO(:) , projection(:,:)
character(len=21)    :: string
character(len=2) , allocatable :: list(:)
! local parameters ...
integer , parameter  :: grid_size = 1500


allocate( list , source = fetch_names(sys , basis , instance=token) )

PDOS = evaluate_PDOS(UNI,basis,list)

!----------------------------------------------------------------------------------------------
!                             Guassian broadden of PDOS peaks

gauss_norm = 1.d0  ! 1.d0 / (sigma*sqrt2PI) ! <== for gauss_norm = 1 the gaussians are not normalized ...
two_sigma2 = 2.d0 * sigma*sigma
step = (DOS_range%fim-DOS_range%inicio) / float(grid_size-1)

! states in the range [DOS_range% inicio , DOS_range% fim]
i1 = maxloc(UNI%erg , 1 , UNI%erg <  DOS_range%inicio)
i2 = maxloc(UNI%erg , 1 , UNI%erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1 + 1

call allocate_stuff( n_of_DOS_states , size(list) , grid_size , PDOS_atomic , erg_MO , projection )

erg_MO = UNI%erg  (i1:i2)
projection = PDOS (i1:i2 , :)

forall(k=1:grid_size) PDOS_atomic%grid(k) = (k-1)*step + DOS_range%inicio

do i = 1 , size(list)

     PDOS_atomic% peaks = 0.d0
     PDOS_atomic% func  = 0.d0

     do j = 1 , n_of_DOS_states
          do k = 1 , grid_size
     
              delta_E = PDOS_atomic% grid(k)-erg_MO(j)

              if(dabs(delta_E) < (step/two) )&
              then
                    PDOS_atomic% peaks(k) = PDOS_atomic% peaks(k) + projection(j,i)                                                              
              end if
     
              if( delta_E**2/two_sigma2 < 25.d0 )&
              then
                    PDOS_atomic% func(k) = PDOS_atomic% func(k) + projection(j,i)*gauss_norm*exp( -delta_E**2/two_sigma2 )
              end if

          end do
     end do

     string = "dos.trunk/PDOS-"//trim(list(i))//".dat"

     OPEN( unit=31 , file=string , status='unknown' )
     do k = 1 , grid_size
          write(31,"(3F12.5)") PDOS_atomic%grid(k) , PDOS_atomic%func(k) , PDOS_atomic%peaks(k) 
     end do
     CLOSE(31)

end do

end subroutine atomicPDOS
!
!
!
!==========================================
subroutine orbitalPDOS( sys , basis , UNI , MOnum )
!==========================================
implicit none
type(structure)       , intent(in) :: sys
type(STO_basis)       , intent(in) :: basis(:)
type(R_eigen)         , intent(in) :: UNI
integer , allocatable , intent(in) :: MOnum(:)

! local variables  ...
integer              :: i , j , k
real*8               :: pop(9)
character(len=2) , allocatable :: list(:)
character(len=5)     :: AO_orbs(9)=["s" , "py" , "pz" , "px" , "dxy" , "dyz" , "dz2" , "dxz" , "dx2y2"]

!----------------------------------------------------------------------------------------------
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {Symbol} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!----------------------------------------------------------------------------------------------

allocate( list , source = fetch_names(sys , basis , instance=token) )

open( unit = 18 , file = "dos.trunk/orbitalPDOS.dat" , status = "unknown", action = "write" )

do i = 1 , size(MOnum)
     write(18,"(/a5,i4)") "MO = ",MOnum(i)
     write(18,"(t13,a83)") "S        Py        Pz        Px        Dxy       Dyz       Dz2       Dxz      Dx2y2"
     do k = 1 , size(list)
          write(18,"(a5)",advance='no') list(k) 
          do j = 1 , 9     
               pop(j) = Mulliken( UNI , basis , MO=MOnum(i) , AO = AO_orbs(j) , EHSymbol=list(k) )
          end do
     write(18,"(9F10.4)") (pop(j) , j=1,9)
     end do
end do

close(18)

end subroutine orbitalPDOS
!
!
!
!===============================================
subroutine Read_Command_Lines_Arguments( MOnum )
!===============================================
implicit none
integer , allocatable  , intent(out) :: MOnum(:)

! local variables ...
integer      :: i , total , MO_first , MO_last , arg_i , arg_f
character(4) :: str
character(6) :: MOstr

total = COMMAND_ARGUMENT_COUNT()

select case (total)

       case (0) 
            !fine, carry on

       case (2)
            call get_command_argument(1,str)
            str = TO_UPPER_CASE(str)           
            if( str == "PDOS")& 
            then 
                 call get_command_argument(2,token)
                 token = TO_UPPER_CASE(token)           
                 select case (token)
                     case( "EHSYMBOL", "SYMBOL" )
                         if(token=="EHSYMBOL") token = "EHSymbol"
                         if(token=="SYMBOL")   token = "Symbol"
                         !fine, carry on
                         atomicPDOS_called = .true.
                     case default
                         print*, "usage: dynemol PDOS {EHSymbol,Symbol}"
                         stop
                 end select
           end if  

       case (4:5) 

            call get_command_argument(1,str)
            str = TO_UPPER_CASE(str)           
            if( str == "PDOS")& 
            then 
                 token = "EHSymbol"
                 orbitalPDOS_called = .true.
            end if

            arg_i = total - 2
            arg_f = total

            call get_command_argument(arg_i+1,MOstr)

            select case (MOstr)
              case( "-" )  ! <== orbitals within a range ...

                  CALL GET_COMMAND_ARGUMENT(arg_i, MOstr)
                  read( MOstr,*) MO_first
                  CALL GET_COMMAND_ARGUMENT(arg_f, MOstr)
                  read( MOstr,*) MO_last

                  total  = MO_last - MO_first + 1
                  allocate( MOnum(total) )

                  do i = 1 , total
                     MOnum(i) = MO_first + (i-1)
                  end do

              case default  ! <== list of  3 orbitals ...     

                  allocate( MOnum(total) )
                  do i = 1 , total
                      CALL GET_COMMAND_ARGUMENT(i, MOstr)
                      read( MOstr,*) MOnum(i)
                  end do
            end select

       case default
            print*, "please check dynemol -h for usage"
            stop

end select

end subroutine Read_Command_Lines_Arguments
!
!
!
!=====================================================
 function evaluate_PDOS (UNI,basis,list)  result(PDOS)
!=====================================================
implicit none
type(R_eigen)   , intent(in) :: UNI
type(STO_basis) , intent(in) :: basis(:)
character(*)    , intent(in) :: list(:)

! local variables ...
integer              :: i , j , k
real*8 , allocatable :: PDOS(:,:)
character(len=5)     :: AO_orbs(9)=["s" , "py" , "pz" , "px" , "dxy" , "dyz" , "dz2" , "dxz" , "dx2y2"]

!----------------------------------------------------------------------------------------------
! ==> Mulliken( OPT_UNI , basis , MO , {atom}=[.,.,.] , {AO} , {EHSymbol} , {residue} , {Symbol} )
! Population analysis ...
! {...} terms are optional  
! AO = s , py , pz , px , dxy , dyz , dz2 , dxz , dx2y2
!----------------------------------------------------------------------------------------------

allocate( PDOS (size(UNI%erg),size(list)) , source = 0.d0 )

! uncomment for debugging
!write(*,"(t1,i3,t12,a6,t22,a5,t35,a2,t45,a2,t55,a2,t65,a2)") i , "MO erg" , "total" , list(:)

select case (token)

       case("EHSymbol")
           do i = 1 , size(UNI%erg)
                do j = 1 , 9     
                     do k = 1 , size(list)
                          PDOS(i,k) = PDOS(i,k) + Mulliken( UNI , basis , MO=i , AO = AO_orbs(j) , EHSymbol=list(k) )
                     end do
                end do
                !write(*,"(t1,i3,t6,F12.4,F10.5,5F10.3)") i , UNI%erg(i) , sum(PDOS(i,:)) , PDOS(i,:)
           end do

       case("Symbol")
           do i = 1 , size(UNI%erg)
                do j = 1 , 9     
                     do k = 1 , size(list)
                          PDOS(i,k) = PDOS(i,k) + Mulliken( UNI , basis , MO=i , AO = AO_orbs(j) , Symbol=list(k) )
                     end do
                end do
                !write(*,"(t1,i3,t6,F12.4,F10.5,5F10.3)") i , UNI%erg(i) , sum(PDOS(i,:)) , PDOS(i,:)
           end do

end select

end function evaluate_PDOS
!
!
!
!=======================================================
 function fetch_names (sys,basis,instance)  result(list)
!=======================================================
type(structure) , intent(in) :: sys
type(STO_basis) , intent(in) :: basis(:)
character(*)    , intent(in) :: instance

!local variables ...
integer :: i
integer , allocatable :: temp(:)
character(len=2) , allocatable :: list(:), string(:)

allocate( string( size(basis) ) )
allocate( temp  ( size(basis) ) , source = 0 )

select case (instance)

       case( "Symbol" ) 
       string(:) = basis(:)% Symbol

       case( "EHSymbol" )
       string(:) = basis(:)% EHSymbol

end select

do i = 1 , size(basis)

    ! find different Symbols ... 
    if( .NOT. any( string(1:i-1) == string(i) ) ) then
        temp(i) = i
    end if

end do

allocate( list( count(temp/=0) ) )
list = pack(string , temp/=0) 

deallocate( string , temp )

end function fetch_names
!
!
!
!=========================================================================================================
 subroutine allocate_stuff( n_of_DOS_states , list_size , grid_size , PDOS_atomic , erg_MO , projection )
!=========================================================================================================
implicit none
integer               , intent(in)  :: n_of_DOS_states
integer               , intent(in)  :: list_size
integer               , intent(in)  :: grid_size
type(f_grid)          , intent(out) :: PDOS_atomic
real*8  , allocatable , intent(out) :: erg_MO(:)
real*8  , allocatable , intent(out) :: projection(:,:)

allocate( erg_MO           (n_of_DOS_states)          )
allocate( projection       (n_of_DOS_states,list_size))
allocate(PDOS_atomic% grid (grid_size))
allocate(PDOS_atomic% peaks(grid_size) , source = 0.d0)
allocate(PDOS_atomic% func (grid_size) , source = 0.d0)

end  subroutine allocate_stuff
!
!
!
end module diagnostic_m

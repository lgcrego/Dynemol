module diagnostic_m

 use type_m
 use omp_lib
 use constants_m
 use DOS_m
 use parameters_m               , only : spectrum , DP_Moment , &
                                         survival , EnvField_ , &
                                         Alpha_Tensor ,         &
                                         GaussianCube ,         &
                                         HFP_Forces
 use Solvated_M                 , only : DeAllocate_TDOS ,      &
                                         DeAllocate_PDOS ,      &
                                         DeAllocate_SPEC 
 use QCModel_Huckel             , only : EigenSystem
 use Structure_Builder          , only : Unit_Cell ,            &
                                         Extended_Cell ,        &
                                         Generate_Structure ,   &
                                         Basis_Builder ,        &
                                         ExCell_basis
 use GA_QCModel_m               , only : Mulliken
 use DP_main_m                  , only : Dipole_Matrix
 use Dielectric_Potential       , only : Environment_SetUp
 use Oscillator_m               , only : Optical_Transitions
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format
 use Data_Output                , only : Dump_stuff
 use Embedded_FF_Alpha          , only : AlphaPolar
 use HuckelForces_m             , only : HuckelForces

 public :: diagnostic

 private

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

 do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
 end do

 If( Spectrum ) CALL Optical_Transitions( Extended_Cell, ExCell_basis, UNI , SPEC )

 If( HFP_Forces ) CALL HuckelForces( Extended_Cell, ExCell_basis, UNI )

 Print*, " " 
 Print*, "dE1 = ", UNI%erg(115) - UNI%erg(114) ,   2.6470
 Print*, "dE2 = ", UNI%erg(114) - UNI%erg(113) ,   0.3040
 Print*, "dE3 = ", UNI%erg(115) - UNI%erg(113) ,   2.9510
 Print*, "dE4 = ", UNI%erg(113) - UNI%erg(112) ,   0.8950
 Print*, "dE5 = ", UNI%erg(112) - UNI%erg(111) ,   0.4360
 Print*, "dE6 = ", UNI%erg(117) - UNI%erg(116) ,   1.6000
 
 If( GaussianCube .AND. (size(MOnum) > 0) ) then
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
!===============================================
subroutine Read_Command_Lines_Arguments( MOnum )
!===============================================
implicit none
integer , allocatable  , intent(out) :: MOnum(:)

! local variables ...
integer      :: i , MO_total , MO_first , MO_last
character(6) :: MOstr

MO_total  = COMMAND_ARGUMENT_COUNT()

IF( MO_total == 3) then
    call get_command_argument(2,MOstr)

    select case (MOstr)

      case( ":" )  ! <== orbitals within a range ...

          CALL GET_COMMAND_ARGUMENT(1, MOstr)
          read( MOstr,*) MO_first
          CALL GET_COMMAND_ARGUMENT(3, MOstr)
          read( MOstr,*) MO_last

          MO_total  = MO_last - MO_first + 1
          allocate( MOnum(MO_total) )

          do i = 1 , MO_total
             MOnum(i) = MO_first + (i-1)
          end do

      case default  ! <== list of  3 orbitals ...     

          allocate( MOnum(MO_total) )
          do i = 1 , MO_total
              CALL GET_COMMAND_ARGUMENT(i, MOstr)
              read( MOstr,*) MOnum(i)
          end do

   end select

ELSE
   ! arbitrary (/= 3) list of orbitals ...     
   allocate( MOnum(MO_total) )
   do i = 1 , MO_total
       CALL GET_COMMAND_ARGUMENT(i, MOstr)
       read( MOstr,*) MOnum(i)
   end do

end IF

end subroutine Read_Command_Lines_Arguments
!
!
!
end module diagnostic_m

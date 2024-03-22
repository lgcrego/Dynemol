module GA_driver_m

 use type_m
 use constants_m
 use MPI_definitions_m          , only : slave 
 use parameters_m               , only : spectrum , DP_Moment , GaussianCube , Alpha_Tensor , OPT_parms
 use Solvated_m                 , only : DeAllocate_TDOS , DeAllocate_PDOS , DeAllocate_SPEC 
 use GA_m                       , only : Genetic_Algorithm , Dump_OPT_parameters
 use GA_QCModel_m               , only : GA_eigen , Mulliken , GA_DP_Analysis , AlphaPolar 
 use cost_EH                    , only : evaluate_cost , REF_DP , REF_Alpha 
 use Semi_Empirical_Parms       , only : EH_atom
 use Multipole_Routines_m       , only : Util_multipoles
 use Structure_Builder          , only : Generate_Structure , Extended_Cell , Unit_Cell , Basis_Builder , ExCell_basis , Cube_Coef , Cube_Zeta
 use Oscillator_m               , only : Optical_Transitions
 use Data_Output                , only : Dump_stuff
 use Psi_squared_cube_format    , only : Gaussian_Cube_Format
 use DOS_m 

 public :: GA_driver

 private

 contains
!
!
!
!====================
 subroutine GA_driver
!====================
implicit none 

! local variables ...
 integer                        :: i , nr , N_of_residues 
 integer         , allocatable  :: MOnum(:)
 real*8                         :: first_cost , Alpha_ii(3)
 logical                        :: DIPOLE_
 type(R_eigen)                  :: UNI
 type(f_grid)                   :: TDOS , SPEC
 type(f_grid)    , allocatable  :: PDOS(:) 
 type(STO_basis) , allocatable  :: OPT_basis(:)
 real*8                         :: DP(3) = d_zero 

! preprocessing stuff ...................................

DIPOLE_ = ( spectrum .OR. DP_Moment )

CALL DeAllocate_TDOS( TDOS , flag="alloc" )
CALL DeAllocate_PDOS( PDOS , flag="alloc" )
CALL DeAllocate_SPEC( SPEC , flag="alloc" )

N_of_residues = size( Unit_Cell%list_of_residues )

! reading command line arguments for plotting MO cube files ...
CALL Read_Command_Lines_Arguments( MOnum )

!.........................................................

! setting up the system ...

CALL Generate_Structure(1)

CALL Basis_Builder( Extended_Cell, ExCell_basis )

! setting up constraints ...
CALL GA_eigen( Extended_Cell, ExCell_basis, UNI )

If( DIPOLE_ ) &
then
    CALL Util_multipoles
    CALL GA_DP_Analysis( Extended_Cell, ExCell_basis, UNI%L, UNI%R, DP )  
end if

! calculates the cost on input, for future comparison ...
first_cost = evaluate_cost(Extended_Cell, UNI, ExCell_basis, DP ) 

! Optimization of Huckel parameters ... 
CALL Genetic_Algorithm( ExCell_basis, OPT_basis )

! calculations with new parameters ...
CALL GA_eigen( Extended_Cell, OPT_basis, UNI , flag=2 )

CALL Total_DOS( UNI%erg, TDOS )

do nr = 1 , N_of_residues
    CALL Partial_DOS( Extended_Cell , UNI , PDOS , nr )            
end do

If( DIPOLE_ ) CALL GA_DP_Analysis( Extended_Cell, OPT_basis, UNI%L, UNI%R, DP )  

If( Alpha_Tensor .AND. DP_Moment ) CALL AlphaPolar( Extended_Cell, OPT_basis , Alpha_ii )

If( spectrum ) CALL Optical_Transitions( Extended_Cell, OPT_basis, UNI , SPEC )

!----------------------------------------------
! printing zone ...

! compare costs to evalualte otimization ...
Print*, " " 
If( DIPOLE_ ) &
then
    Print 210 , evaluate_cost( Extended_Cell, UNI, OPT_basis, DP ,ShowCost=.true. ) , first_cost
    Print 154, DP, sqrt( dot_product(DP,DP) )
   !Print 189 , Alpha_ii , sum( Alpha_ii ) / three 
else
    Print 210 , evaluate_cost( Extended_Cell, UNI, OPT_basis, ShowCost=.true. ) , first_cost
end if

Print*, " " 
Print 10, "dE1 = ",UNI%erg(5) - UNI%erg(4) , "  vs " , 8.00d0 , "  => error = ", ( UNI%erg(5) - UNI%erg(4)) - 8.0d0
Print 10, "dE2 = ",UNI%erg(4) - UNI%erg(3) , "  vs " , 2.10d0 , "  => error = ", ( UNI%erg(4) - UNI%erg(3)) - 2.10d0
Print 10, "dE3 = ",UNI%erg(3) - UNI%erg(2) , "  vs " , 4.d0 , "  => error = ", ( UNI%erg(3) - UNI%erg(2)) - 4.d0

10 format(A6,F9.5,A5,F9.5,A13,F9.5)

CALL Dump_OPT_parameters( OPT_basis , output='STDOUT' )

! Population analysis ...
If( GaussianCube ) &
then
    ! do it again because the STO have been optimized
    do i = 1 , size(OPT_basis)  
        Cube_Coef(i,:) = OPT_basis(i) % coef(:)
        Cube_Zeta(i,:) = OPT_basis(i) % zeta(:) * a_Bohr
    end do

    do i = 1 , size(MOnum)
        CALL Gaussian_Cube_Format( UNI%L(MOnum(i),:) , UNI%R(:,MOnum(i)) , MOnum(i) , 0.d0 )
    end do

    Print 220 , MOnum(:)
end if

If( OPT_parms ) then    
     Print 445
     Print 45 , ( EH_atom(i)% EHSymbol , EH_atom(i)% residue , i = 1,size(EH_atom) )
else 
     Print*, ">> OPT_parms were not used <<"
end if

!----------------------------------------------

CALL Dump_stuff( TDOS , PDOS , SPEC )

CALL DeAllocate_TDOS( TDOS , flag="dealloc" )
CALL DeAllocate_PDOS( PDOS , flag="dealloc" )
CALL DeAllocate_SPEC( SPEC , flag="dealloc" )

include 'formats.h'

end subroutine GA_driver
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

IF( MO_total == 4) then

    call get_command_argument(3,MOstr)

    select case (MOstr)
      
      case( "-" )  ! <== orbitals within a range ...
         
          CALL GET_COMMAND_ARGUMENT(2, MOstr)
          read( MOstr,*) MO_first
          CALL GET_COMMAND_ARGUMENT(4, MOstr)
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
end module GA_driver_m

 module Schroedinger_m

 use type_m
 use constants_m
 use f95_precision
 use blas95
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube , preview, &
                                         GaussianCube_step ,  DP_Moment , electron_state ,  &
                                         Coulomb_ , DensityMatrix , driver
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use FMO_m                      , only : FMO_analysis , orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations , Net_Charge
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format
 use Backup_m                   , only : Security_Copy , Restart_state
 use Auto_Correlation_m         , only : MO_Occupation


    public :: Simple_dynamics , DeAllocate_QDyn , RunningStat

    private

    ! module variables ...
    Complex*16   , ALLOCATABLE , dimension(:,:)   :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Real*8       , ALLOCATABLE , dimension(:,:,:) :: Pops
    type(f_time)                                  :: Mean_QDyn , aux
    integer                                       :: iter = 0 

 contains
!
!
!=====================================================
 subroutine Simple_dynamics(system, basis, UNI, QDyn )
!=====================================================
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(inout) :: basis(:)
 type(R_eigen)   , intent(in)    :: UNI
 type(f_time)    , intent(inout) :: QDyn

! local variables ...
integer                          :: j , nn , mm
integer                          :: it , n , it_init
real*8                           :: t , t_rate
real*8                           :: Total_DP(3)
complex*16      , ALLOCATABLE    :: phase(:)
type(R_eigen)                    :: el_FMO , hl_FMO

! ------------------ preprocess stuff --------------------

allocate( Pops( n_t , 0:size(system%list_of_fragments)+1 , n_part ) )

mm = size(basis) ; nn = n_part

CALL Allocate_Brackets( mm , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

t  = t_i
it = 1

!   building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
do n = 1 , n_part
    select case( eh_tag(n) )

        case( "el" )

            CALL FMO_analysis ( system , basis , UNI , el_FMO , instance="E" )

            MO_bra( : , n ) = el_FMO%L( orbital(n) , : )
            MO_ket( : , n ) = el_FMO%R( : , orbital(n) )

            Print 591, orbital(n) , el_FMO%erg(orbital(n))

        case( "hl" )

            CALL FMO_analysis ( system , basis , UNI , hl_FMO , instance="H" )

            MO_bra( : , n ) = hl_FMO%L( orbital(n) , : )
            MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )

            Print 592, orbital(n) , hl_FMO%erg(orbital(n))
            If( (orbital(n) > hl_FMO%Fermi_State) ) write(*,"(/a)") '>>> warning: hole state above the Fermi level <<<'

        end select
end do

! deallocate after use ...
deallocate( el_FMO%L , el_FMO%R , el_FMO%erg , hl_FMO%L , hl_FMO%R , hl_FMO%erg )

! stop here to preview and check input and system info ...
If( preview ) stop

! DUAL representation for efficient calculation of survival probabilities ...
CALL DZgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%R , mm , MO_ket , mm , C_zero , DUAL_ket , mm )
CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_bra , mm , C_zero , DUAL_bra , mm )

! save populations ...
Pops(1,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )

QDyn%dyn(1,:,:) = Pops(1,:,:)

If( DensityMatrix ) then
    If( n_part == 1 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI )
    If( n_part == 2 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI, UNI )
End If

!   save the initial GaussianCube file ...
If( GaussianCube ) then

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra
    CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

    do n = 1 , n_part
        if( eh_tag(n) == "XX" ) cycle
        CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
    end do

end If
!-------------------------------------------------------------
!                       Q-DYNAMICS

it_init = it + 1

t_rate = (t_f - t_i) / float(n_t)

DO it = it_init , n_t

    t = t + t_rate

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    forall( j=1:n_part)
        MO_bra(:,j) = merge( conjg(phase(:)) * MO_bra(:,j) , C_zero , eh_tag(j) /= "XX" )
        MO_ket(:,j) = merge(       phase(:)  * MO_ket(:,j) , C_zero , eh_tag(j) /= "XX" )
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    CALL dzgemm( 'N' , 'N' , mm , nn , mm , C_one , UNI%R , mm , MO_ket , mm , C_zero , DUAL_ket , mm )
    CALL dzgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_bra , mm , C_zero , DUAL_bra , mm )

    ! get populations ...
    Pops(it,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )
    QDyn%dyn(it,:,:) = Pops(it,:,:)

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra
    CALL DZgemm( 'T' , 'N' , mm , nn , mm , C_one , UNI%L , mm , MO_ket , mm , C_zero , AO_ket , mm )

    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then

        do n = 1 , n_part
            if( eh_tag(n) == "XX" ) cycle
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
        end do

    end If

    if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

END DO

! saving populations & all that jazz ...
select case (driver) 

       case( "q_dynamics" ) 
           ! dump final populations ...
           CALL dump_QDyn( QDyn , UNI )

       case( "avrg_confgs" ) 

           CALL RunningStat( Qdyn )

           ! deliver on-the-fly averages back to calling routine ...
           Qdyn% dyn = Mean_Qdyn% dyn

           ! dump on-the-fly averages in here ...
           CALL dump_QDyn( QDyn , UNI )

       case default
           Print*, " >>> Check your driver options <<< :" , driver
           stop

end select 

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

include 'formats.h'

end subroutine Simple_dynamics
!
!
!
!
!=========================================
 subroutine RunningStat( Qdyn , instance )
!=========================================
implicit none
type(f_time)            , intent(inout) :: QDyn
character    , optional , intent(in)    :: instance

! Donald Knuth’s Art of Computer Programming, Vol 2, page 232, 3rd edition 
! M[1] = x[1] ; S[1] = 0
! M[k] = ( (k-1)*M[k-1] + x[k] ) / k
! S[k] = S[k-1] + (x[k] – M[k-1]) * (x[k] – M[k])
! M[] = mean value
! S[n]/(n-1)  = S^2 = variance

!===============================================
If( present(instance) ) then
   Qdyn% std = sqrt( Qdyn% std/float(iter-1) )
   ! preserving time; does not undergo std ...
   Qdyn% std(:,0,:) = Qdyn% dyn(:,0,:)  
  
   return
end If
!===============================================

iter = iter + 1

If( .not. allocated(Mean_Qdyn% dyn) ) then
    allocate( Mean_Qdyn% dyn , source = QDyn% dyn )
    allocate( aux% dyn       , source = QDyn% dyn )
End If

aux% dyn = Mean_Qdyn% dyn

Mean_Qdyn% dyn = ( (iter-1)*aux% dyn + Qdyn% dyn ) / iter

If(iter == 1 ) then

   Qdyn% std = D_zero

else

   Qdyn% std = Qdyn% std + (Qdyn% dyn - aux% dyn)*(Qdyn% dyn - Mean_Qdyn% dyn)

End If

end subroutine RunningStat
!
!
!
!
!==================================
 subroutine dump_Qdyn( Qdyn , UNI )
!==================================
implicit none
type(f_time)    , intent(in) :: QDyn
type(R_eigen)   , intent(in) :: UNI

! local variables ...
integer     :: nf , n , it
complex*16  :: wp_energy

do n = 1 , n_part

      if( eh_tag(n) == "XX" ) cycle
      
      open( unit = 52 , file = "dyn.trunk/"//eh_tag(n)//"_survival.dat" , status = "replace" , action = "write" , position = "append" )
      write(52,15) "#" ,( nf+1 , nf=0,size(QDyn%fragments)+1 )  ! <== numbered columns for your eyes only ...
      write(52,12) "#" , QDyn%fragments , "total"
      
      open( unit = 53 , file = "dyn.trunk/"//eh_tag(n)//"_wp_energy.dat" , status = "replace" , action = "write" , position = "append" )

      wp_energy = sum(MO_bra(:,n) * UNI%erg(:) * MO_ket(:,n))
      
      DO it = 1 , n_t

      
          ! dumps el-&-hl populations ...
          write(52,13) ( QDyn%dyn(it,nf,n) , nf=0,size(QDyn%fragments)+1 )
      
          ! dumps el-&-hl wavepachet energies ...
          write(53,14) QDyn%dyn(it,0,n) , real(wp_energy) , dimag(wp_energy)
      
      end do

      close(52)
      close(53)
      
end do

12 FORMAT(/15A10)
13 FORMAT(F11.6,14F10.5)
14 FORMAT(3F12.6)
15 FORMAT(A,I9,14I10)

end subroutine dump_Qdyn
!
!
!
!
!=========================================
 subroutine DeAllocate_QDyn( QDyn , flag )
!=========================================
implicit none
type(f_time)  , intent(inout) :: QDyn
character(*)  , intent(in)    :: flag

! local variable ...
integer      :: N_of_fragments
character(1) :: first_in_line
logical      :: A_flag

select case( flag )

    case( "alloc" )

        if( allocated(trj) ) then

            CALL Coords_from_Universe( Unit_Cell, trj(2) )          ! <== use number 2 to avoid verbose
            CALL Generate_Structure( 2 )
            N_of_fragments = size( Extended_Cell%list_of_fragments )

        else

            CALL Generate_Structure( 2 )                            ! <== use number 2 to avoid verbose
            N_of_fragments = size( Extended_Cell%list_of_fragments )

        end if

        ! for the sake of having the DONOR or ACCEPTOR survival probability in the first column at output ...
        A_flag = any(Extended_Cell%list_of_fragments == "A")
        first_in_line = Extended_Cell%list_of_fragments(1)
        If( A_flag ) then
            where( Extended_Cell%list_of_fragments == "A" ) Extended_Cell%list_of_fragments = first_in_line
        else
            where( Extended_Cell%list_of_fragments == "D" ) Extended_Cell%list_of_fragments = first_in_line
        end If
        Extended_Cell%list_of_fragments(1) = merge( "A" , "D" , A_flag )

        ! QDyn%dyn = ( time ; fragments ; all fragments ) ...
        allocate( QDyn%fragments( size(Extended_Cell % list_of_fragments) ) , source = Extended_Cell % list_of_fragments )
        allocate( QDyn%dyn      ( n_t , 0:N_of_fragments+1 , n_part       ) , source = 0.d0                              )
        If( driver == "avrg_confgs" ) then
              allocate( QDyn%std( n_t , 0:N_of_fragments+1 , n_part       ) , source = 0.d0                              )
        End If
 
        ! allocatating Net_Charte for future use ...
        allocate( Net_Charge(Extended_Cell%atoms) , source = D_zero )

        ! cleaning the mess ...
        CALL DeAllocate_Structures( Extended_Cell )

    case( "dealloc" )

        deallocate( QDyn%dyn , QDyn%fragments )
        
        If( allocated(QDyn% std) ) deallocate(QDyn% std)

end select

end subroutine DeAllocate_QDyn
!
!
!
end module Schroedinger_m

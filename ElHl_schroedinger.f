 module Schroedinger_m

 use type_m
 use constants_m
 use f95_precision
 use blas95
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube , preview, &
                                         GaussianCube_step ,  DP_Moment , electron_state ,  &
                                         Coulomb_ , DensityMatrix , driver , SOC , comb
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use FMO_m                      , only : FMO_analysis , orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations , Net_Charge , FileName
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format
 use Backup_m                   , only : Security_Copy , Restart_state
 use Auto_Correlation_m         , only : MO_Occupation

 use Overlap_Builder  , only : Overlap_Matrix


    public :: Simple_dynamics , DeAllocate_QDyn , RunningStat

    private

    ! module variables ...
    Complex*16   , ALLOCATABLE , dimension(:,:)     :: MO_bra , MO_ket , AO_bra , AO_ket , DUAL_ket , DUAL_bra
    Real*8       , ALLOCATABLE , dimension(:,:,:,:) :: Pops
    type(f_time)                                    :: Mean_QDyn , aux
    integer                                         :: n_spin
    integer                                         :: iter = 0 

 contains
!
!
!=====================================================
 subroutine Simple_dynamics(system, basis, UNI, QDyn )
!=====================================================
 implicit none
 type(structure) , intent(inout) :: system
 type(STO_basis) , intent(inout) :: basis(:)
 type(C_eigen)   , intent(in)    :: UNI
 type(f_time)    , intent(inout) :: QDyn

! local variables ...
integer                          :: i , j , nn , mm
integer                          :: it , n , it_init
real*8                           :: t , t_rate
real*8                           :: Total_DP(3)
complex*16      , ALLOCATABLE    :: phase(:)
type(C_eigen)                    :: el_FMO , hl_FMO

! ------------------ preprocess stuff --------------------

mm = merge(size(system%list_of_fragments)+2,size(system%list_of_fragments)+1,SOC)

allocate( Pops(n_t , 0:mm , n_part , n_spin) )

mm = size(basis) ; nn = n_part

CALL Allocate_Brackets( mm , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

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

if( comb ) then
    MO_bra(:,1) = ( MO_bra(:,1) + MO_bra(:,2) ) / dsqrt(TWO)
    MO_ket(:,1) = ( MO_ket(:,1) + MO_ket(:,2) ) / dsqrt(TWO)
    MO_bra(:,2) = MO_bra(:,1)
    MO_ket(:,2) = MO_ket(:,1)
end if

! deallocate after use ...
if( eh_tag(1) == "el" ) deallocate( el_FMO%L , el_FMO%R , el_FMO%erg )
if( eh_tag(2) == "hl" ) deallocate( hl_FMO%L , hl_FMO%R , hl_FMO%erg )

! stop here to preview and check input and system info ...
If( preview ) stop

! DUAL representation for efficient calculation of survival probabilities ...
CALL gemm( UNI%R , MO_ket , DUAL_ket , 'N' , 'N' )
CALL gemm( UNI%L , MO_bra , DUAL_bra , 'T' , 'N' )

! save populations ...
t  = t_i
it = 1
Pops(it,:,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )
QDyn%dyn(it,:,:,:) = Pops(it,:,:,:)

!If( DensityMatrix ) then
!    If( n_part == 1 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI )
!    If( n_part == 2 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI, UNI )
!End If

!   save the initial GaussianCube file ...
If( GaussianCube ) then

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra
    CALL gemm( UNI%L , MO_ket , AO_ket , 'T' , 'N' )

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
    CALL gemm( UNI%R , MO_ket , DUAL_ket , 'N' , 'N' )
    CALL gemm( UNI%L , MO_bra , DUAL_bra , 'T' , 'N' )

    ! get populations ...
    Pops(it,:,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )
    QDyn%dyn(it,:,:,:) = Pops(it,:,:,:)

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra
    CALL gemm( UNI%L , MO_ket , AO_ket , 'T' , 'N' )

    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then

        do n = 1 , n_part
            if( eh_tag(n) == "XX" ) cycle
            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it , t , eh_tag(n) )
        end do

    end If

!    if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

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
   Qdyn% std(:,0,:,:) = Qdyn% dyn(:,0,:,:)  
  
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
!==================================
 subroutine dump_Qdyn( Qdyn , UNI )
!==================================
implicit none
type(f_time)    , intent(in) :: QDyn
type(C_eigen)   , intent(in) :: UNI

! local variables ...
logical                       :: erg_done
integer                       :: nf , n , nc , it , spin
complex*16                    :: wp_energy
character(len=:) ,allocatable :: f_name

real*8 , allocatable :: vec_up(:) , vec_down(:) , p(:) , vec(:)
integer :: i

nc = merge(size(QDyn%fragments)+2,size(QDyn%fragments)+1,SOC)

erg_done = no
do spin = 1 , n_spin
    
    do n = 1 , n_part

        if( eh_tag(n) == "XX" ) cycle

        call FileName( f_name , n , spin , instance="dens" )
        open( unit = 52 , file = f_name , status = "replace" , action = "write" , position = "append" )
        write(52,15) "#" , ( nf+1 , nf=0,nc )  ! <== numbered columns for your eyes only ...

        if( SOC ) then
            write(52,12) "#" , QDyn%fragments , "total" , "both spins"
        else
            write(52,12) "#" , QDyn%fragments , "total"
        end if

        do it = 1 , n_t

            ! dumps el-&-hl populations ...
            write(52,13) ( QDyn%dyn(it,nf,n,spin) , nf=0,nc )
        
        end do

        close(52)
        
        if( not(erg_done) ) then

            call FileName( f_name , n , spin , instance="erg" )
            open( unit = 53 , file = f_name , status = "replace" , action = "write" , position = "append" )
            wp_energy = sum(MO_bra(:,n) * UNI%erg(:) * MO_ket(:,n))
        
            do it = 1 , n_t

                ! dumps el-&-hl wavepachet energies ...
                write(53,14) QDyn%dyn(it,0,n,spin) , dreal(wp_energy) , dimag(wp_energy)
        
            end do
        
            close(53)
        
            erg_done = yes

        end if
         
    end do
end do

if( SOC ) then

    open( unit = 54 , file = "dyn.trunk/"//"el"//"_dens.dat" , status = "replace" , action = "write" , position = "append" )

    write(54,15) "#" , ( nf+1 , nf=0,nc-1 )  ! <== numbered columns for your eyes only ...
    write(54,12) "#" , QDyn%fragments , "total"

    do it = 1 , n_t

        write(54,13) QDyn%dyn(it,0,1,1) , ( QDyn%dyn(it,nf,1,1) + QDyn%dyn(it,nf,1,2) , nf=1,nc-1 )

    end do

    close(54)

    if( n_part == 2 ) then

        open( unit = 54 , file = "dyn.trunk/"//"hl"//"_dens.dat" , status = "replace" , action = "write" , position = "append" )

        write(54,15) "#" , ( nf+1 , nf=0,nc-1 )  ! <== numbered columns for your eyes only ...
        write(54,12) "#" , QDyn%fragments , "total"

        do it = 1 , n_t

            write(54,13) QDyn%dyn(it,0,2,1) , ( QDyn%dyn(it,nf,2,1) + QDyn%dyn(it,nf,2,2) , nf=1,nc-1 )

        end do

        close(54)

    end if
end if

if( SOC ) then

    allocate( vec ( nc - 1 ) , source = D_zero )

    if( not(comb) ) then

        allocate( vec_up   ( nc - 1 ) , source = D_zero )
        allocate( vec_down ( nc - 1 ) , source = D_zero )

        open( unit = 54 , file = "dyn.trunk/"//"Polarization.dat" , status = "replace" , action = "write" , position = "append" )

        write(54,15) "#" , ( nf+1 , nf=0,nc-1 )  ! <== numbered columns for your eyes only ...
        write(54,12) "#" , QDyn%fragments , "total"

        ! el == py up
        ! hl == py down
        vec = D_zero
        do it = 2 , n_t

            do nf = 1 , nc-1
                vec_up(nf)   = QDyn%dyn(it,nf,1,1) + QDyn%dyn(it,nf,1,2)
                vec_up(nf)   = QDyn%dyn(it,nf,1,1) / vec_up(nf)
                vec_down(nf) = QDyn%dyn(it,nf,2,1) + QDyn%dyn(it,nf,2,2)
                vec_down(nf) = QDyn%dyn(it,nf,2,2) / vec_down(nf)
                vec(nf) = vec(nf) + vec_up(nf) - vec_down(nf)
            end do

            write(54,13) QDyn%dyn(it,0,1,1) , ( vec_up(nf) - vec_down(nf) , nf=1,nc-1 )

        end do

        vec = vec / dfloat( n_t - 1 )
        write(54,fmt='(a4,10ES17.7)') "#" , ( vec(nf) , nf=1,nc-1 )

        close(54)

        deallocate( vec_up , vec_down )

    else

        allocate( p ( nc - 1 ) , source = D_zero )

        open(unit=54 , file = "dyn.trunk/"//"Polarization_comb_+.dat" , status = "replace" , action = "write" , position = "append" )

        write(54,15) "#" , ( nf+1 , nf=0,nc-1 )  ! <== numbered columns for your eyes only ...
        write(54,12) "#" , QDyn%fragments , "total"

        ! el == py up
        ! hl == py down
        vec = D_zero
        do it = 2 , n_t

            do nf = 1 , nc-1
                p(nf) = QDyn%dyn(it,nf,1,1) + QDyn%dyn(it,nf,1,2)
                p(nf) = ( QDyn%dyn(it,nf,1,1) - QDyn%dyn(it,nf,1,2) ) / p(nf)
                vec(nf) = vec(nf) + p(nf)
            end do

            write(54,13) QDyn%dyn(it,0,1,1) , ( p(nf) , nf=1,nc-1 )

        end do

        vec = vec / dfloat( n_t - 1 )
        write(54,fmt='(a4,10ES17.7)') "#" , ( vec(nf) , nf=1,nc-1 )

        close(54)

        deallocate( p )

        allocate( p ( nc - 1 ) , source = D_zero )

        open(unit=54 , file = "dyn.trunk/"//"Polarization_comb_-.dat" , status = "replace" , action = "write" , position = "append" )

        write(54,15) "#" , ( nf+1 , nf=0,nc-1 )  ! <== numbered columns for your eyes only ...
        write(54,12) "#" , QDyn%fragments , "total"

        ! el == py up
        ! hl == py down
        vec = D_zero
        do it = 2 , n_t

            do nf = 1 , nc-1
                p(nf) = QDyn%dyn(it,nf,2,1) + QDyn%dyn(it,nf,2,2)
                p(nf) = ( QDyn%dyn(it,nf,2,1) - QDyn%dyn(it,nf,2,2) ) / p(nf)
                vec(nf) = vec(nf) + p(nf)
            end do

            write(54,13) QDyn%dyn(it,0,1,1) , ( p(nf) , nf=1,nc-1 )

        end do

        vec = vec / dfloat( n_t - 1 )
        write(54,fmt='(a4,10ES17.7)') "#" , ( vec(nf) , nf=1,nc-1 )

        close(54)

        deallocate( p )
      
    end if

    deallocate( vec )

end if

12 FORMAT(/15A11)
!13 FORMAT(F11.6,14F11.5)
13 FORMAT(F11.6,7F16.10)
14 FORMAT(3F12.6)
15 FORMAT(A,I10,14I11)

end subroutine dump_Qdyn
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
integer      :: n , N_of_fragments 
character(1) :: first_in_line
logical      :: A_flag

n_spin = merge(2,1,SOC)

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
        allocate( QDyn%fragments( size(Extended_Cell % list_of_fragments)    ) , source = Extended_Cell % list_of_fragments )
        n = merge(N_of_fragments+2,N_of_fragments+1,SOC)
        allocate( QDyn%dyn      ( n_t , 0:n , n_part , n_spin ) , source = 0.d0 )
        If( driver == "avrg_confgs" ) then
              allocate( QDyn%std( n_t , 0:N_of_fragments+1 , n_part , n_spin ) , source = 0.d0                              )
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

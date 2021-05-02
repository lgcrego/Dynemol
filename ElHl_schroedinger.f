 module Schroedinger_m

 use type_m
 use constants_m
 use f95_precision
 use blas95
 use parameters_m               , only : t_i , t_f , n_t , n_part , GaussianCube , preview, &
                                         GaussianCube_step ,  DP_Moment , electron_state ,  &
                                         Coulomb_ , DensityMatrix , driver , SOC , &
                                         comb , tp_comb , h1_state , h2_state , e1_state , e2_state , hole_state , electron_state
 use Allocation_m               , only : Allocate_Brackets , DeAllocate_Structures
 use Babel_m                    , only : trj , Coords_from_Universe
 use Structure_Builder          , only : Unit_Cell , Extended_Cell , Generate_Structure
 use FMO_m                      , only : FMO_analysis , orbital , eh_tag
 use DP_main_m                  , only : Dipole_Moment
 use Data_Output                , only : Populations , Net_Charge , FileName
 use Psi_Squared_Cube_Format    , only : Gaussian_Cube_Format , probability_flux
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

real*8     , allocatable :: i_up(:,:) , i_down(:,:) , i_total(:,:) , i_SO(:,:) , Sup(:) , Sdown(:) , S_real(:,:)
real*8                :: z_LFT , z_RGT , a , b , c , d
integer               :: iLFT , iRGT

! ------------------ preprocess stuff --------------------

mm = merge(size(system%list_of_fragments)+2,size(system%list_of_fragments)+1,SOC)

allocate( Pops(n_t , 0:mm , n_part , n_spin) )

mm = size(basis) ; nn = n_part

allocate( i_up    ( n_t , 4 ) , source = D_zero )
allocate( i_down  ( n_t , 4 ) , source = D_zero )
allocate( i_SO    ( n_t , 4 ) , source = D_zero )
allocate( i_total ( n_t , 4 ) , source = D_zero )
allocate( Sup     ( n_t )     , source = D_zero )
allocate( Sdown   ( n_t )     , source = D_zero )

do i = 1 , mm
    if( basis(i) % residue == "LFT" ) then
!    if( basis(i) % residue == "DNR" ) then
        iLFT = i
        exit
    end if
end do

do i = 2 , mm
    if( basis(i) % residue == "ACP" .AND. basis(i-1) % residue == "RGT" ) then
!    if( basis(i) % residue == "RGT" .AND. basis(i-1) % residue == "DNR" ) then
        iRGT = i - 1
        exit
    end if
end do

z_LFT = ( basis( iLFT - 1 ) % z + basis( iLFT ) % z ) / 2.0d0
z_RGT = ( basis( iRGT + 1 ) % z + basis( iRGT ) % z ) / 2.0d0

print*, z_LFT , z_RGT

CALL Allocate_Brackets( mm , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase )

!   building up the electron and hole wavepackets with expansion coefficients at t = 0  ...
do n = 1 , n_part
    select case( eh_tag(n) )

        case( "el" )

            CALL FMO_analysis ( system , basis , UNI , el_FMO , instance="E" )

            MO_bra( : , n ) = el_FMO%L( orbital(n) , : )
            MO_ket( : , n ) = el_FMO%R( : , orbital(n) )

            Print 591, orbital(n) , el_FMO%erg(orbital(n))

            deallocate( el_FMO%L , el_FMO%R , el_FMO%erg )

        case( "hl" )

            CALL FMO_analysis ( system , basis , UNI , hl_FMO , instance="H" )

            MO_bra( : , n ) = hl_FMO%L( orbital(n) , : )
            MO_ket( : , n ) = hl_FMO%R( : , orbital(n) )

            Print 592, orbital(n) , hl_FMO%erg(orbital(n))
            If( (orbital(n) > hl_FMO%Fermi_State) ) write(*,"(/a)") '>>> warning: hole state above the Fermi level <<<'

            deallocate( hl_FMO%L , hl_FMO%R , hl_FMO%erg )

    end select
end do

! stop here to preview and check input and system info ...
If( preview ) stop

! DUAL representation for efficient calculation of survival probabilities ...

CALL gemm( UNI%L , MO_bra , DUAL_bra , 'T' , 'N' )
CALL gemm( UNI%R , MO_ket , DUAL_ket , 'N' , 'N' )

Sup(1)   = dreal( sum( DUAL_bra(:,1) * DUAL_ket(:,1) * basis(:) % s ) )
Sdown(1) = dreal( sum( DUAL_bra(:,2) * DUAL_ket(:,2) * basis(:) % s ) )

! save populations ...
t  = t_i
it = 1
Pops(it,:,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t_i )
QDyn%dyn(it,:,:,:) = Pops(it,:,:,:)

!If( DensityMatrix ) then
!    If( n_part == 1 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI )
!    If( n_part == 2 ) CALL MO_Occupation( t_i, MO_bra, MO_ket, UNI, UNI )
!End If

! LOCAL representation for film STO production ...
AO_bra = DUAL_bra
AO_ket = dconjg(AO_bra)

CALL probability_flux( AO_bra(:,1) , AO_ket(:,1) , z_LFT , i_up(1,1) , i_down(1,1) , i_SO(1,1) )
CALL probability_flux( AO_bra(:,1) , AO_ket(:,1) , z_RGT , i_up(1,2) , i_down(1,2) , i_SO(1,2) )
CALL probability_flux( AO_bra(:,2) , AO_ket(:,2) , z_LFT , i_up(1,3) , i_down(1,3) , i_SO(1,3) )
CALL probability_flux( AO_bra(:,2) , AO_ket(:,2) , z_RGT , i_up(1,4) , i_down(1,4) , i_SO(1,4) )

i_total(1,1) = i_up(1,1) + i_down(1,1)
i_total(1,2) = i_up(1,2) + i_down(1,2)
i_total(1,3) = i_up(1,3) + i_down(1,3)
i_total(1,4) = i_up(1,4) + i_down(1,4)

!   save the initial GaussianCube file ...
!If( GaussianCube ) then
!
!    do n = 1 , n_part
!    n = 1
!        if( eh_tag(n) == "XX" ) cycle
!        CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it ,t , eh_tag(n) )
!    end do
!
!end If

!-------------------------------------------------------------
!                       Q-DYNAMICS

it_init = it + 1

t_rate = (t_f - t_i) / float(n_t)

DO it = it_init , n_t

    write(33,*) it , n_t

    t = t + t_rate

    phase(:) = cdexp(- zi * UNI%erg(:) * t_rate / h_bar)

    forall( j=1:n_part)
        MO_bra(:,j) = merge( conjg(phase(:)) * MO_bra(:,j) , C_zero , eh_tag(j) /= "XX" )
        MO_ket(:,j) = merge(       phase(:)  * MO_ket(:,j) , C_zero , eh_tag(j) /= "XX" )
    end forall

    ! DUAL representation for efficient calculation of survival probabilities ...
    CALL gemm( UNI%R , MO_ket , DUAL_ket , 'N' , 'N' )
    CALL gemm( UNI%L , MO_bra , DUAL_bra , 'T' , 'N' )

    Sup(it)   = dreal( sum( DUAL_bra(:,1) * DUAL_ket(:,1) * basis(:) % s ) )
    Sdown(it) = dreal( sum( DUAL_bra(:,2) * DUAL_ket(:,2) * basis(:) % s ) )

    ! get populations ...
    Pops(it,:,:,:) = Populations( QDyn%fragments , basis , DUAL_bra , DUAL_ket , t )
    QDyn%dyn(it,:,:,:) = Pops(it,:,:,:)

    ! LOCAL representation for film STO production ...
    AO_bra = DUAL_bra
    AO_ket = dconjg(AO_bra)

    CALL probability_flux( AO_bra(:,1) , AO_ket(:,1) , z_LFT , i_up(it,1) , i_down(it,1) , i_SO(it,1) )
    CALL probability_flux( AO_bra(:,1) , AO_ket(:,1) , z_RGT , i_up(it,2) , i_down(it,2) , i_SO(it,2) )
    CALL probability_flux( AO_bra(:,2) , AO_ket(:,2) , z_LFT , i_up(it,3) , i_down(it,3) , i_SO(it,3) )
    CALL probability_flux( AO_bra(:,2) , AO_ket(:,2) , z_RGT , i_up(it,4) , i_down(it,4) , i_SO(it,4) )

    i_total(it,1) = i_up(it,1) + i_down(it,1)
    i_total(it,2) = i_up(it,2) + i_down(it,2)
    i_total(it,3) = i_up(it,3) + i_down(it,3)
    i_total(it,4) = i_up(it,4) + i_down(it,4)

!    If( GaussianCube .AND. mod(it,GaussianCube_step) == 0 ) then
!
!        do n = 1 , n_part
!            if( eh_tag(n) == "XX" ) cycle
!            CALL Gaussian_Cube_Format( AO_bra(:,n) , AO_ket(:,n) , it , t , eh_tag(n) )
!        end do
!
!    end If

!    if ( DP_Moment ) CALL Dipole_Moment( system , basis , UNI%L , UNI%R , AO_bra , AO_ket , Dual_ket , Total_DP )

END DO

close(35)

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

open( unit=35 , file="Sz.dat" , action="write" , status="unknown" )
write(35,1000) "#                time" , "Sz, Sz(Psi(0))=1" , "Sz, Sz(Psi(0))=-1"
do i = 1 , n_t
    write(35,1001) dfloat(i-1) * t_rate , Sup(i) , Sdown(i)
end do
close(35)

open( unit=40 , file="UpCurrent.dat" , action="write" , status="unknown" )
write(40,297) "# z_LFT = " , z_LFT
write(40,297) "# z_RGT = " , z_RGT
write(40,299) "#             ||          Sz(t=0)=1         ||          Sz(t=0)=-1        |"
write(40,300) "#          time" , "iup_LFT" , "iup_RGT" , "iup_LFT" , "iup_RGT"
do i = 1 , n_t
    write(40,301) dfloat(i-1) * t_rate , i_up(i,1) , i_up(i,2) , i_up(i,3) , i_up(i,4)
end do
close(40)

open( unit=40 , file="DownCurrent.dat" , action="write" , status="unknown" )
write(40,297) "# z_LFT = " , z_LFT
write(40,297) "# z_RGT = " , z_RGT
write(40,299) "#             ||          Sz(t=0)=1         ||          Sz(t=0)=-1        |"
write(40,300) "#          time" , "idown_LFT" , "idown_RGT" , "idown_LFT" , "idown_RGT"
do i = 1 , n_t
    write(40,301) dfloat(i-1) * t_rate , i_down(i,1) , i_down(i,2) , i_down(i,3) , i_down(i,4)
end do
close(40)

open( unit=40 , file="SOCurrent.dat" , action="write" , status="unknown" )
write(40,297) "# z_LFT = " , z_LFT
write(40,297) "# z_RGT = " , z_RGT
write(40,299) "#             ||          Sz(t=0)=1         ||          Sz(t=0)=-1        |"
write(40,300) "#          time" , "idown_LFT" , "idown_RGT" , "idown_LFT" , "idown_RGT"
do i = 1 , n_t
    write(40,301) dfloat(i-1) * t_rate , i_SO(i,1) , i_SO(i,2) , i_SO(i,3) , i_SO(i,4)
end do
close(40)

open( unit=40 , file="TotalCurrent.dat" , action="write" , status="unknown" )
write(40,297) "# z_LFT = " , z_LFT
write(40,297) "# z_RGT = " , z_RGT
write(40,299) "#             ||          Sz(t=0)=1         ||          Sz(t=0)=-1        |"
write(40,300) "#          time" , "i_LFT" , "i_RGT" , "i_LFT" , "i_RGT"
do i = 1 , n_t
    write(40,301) dfloat(i-1) * t_rate , i_total(i,1) , i_total(i,2) , i_total(i,3) , i_total(i,4)
end do
close(40)

open( unit=40 , file="Area.dat" , action="write" , status="unknown" )
write(40,297) "# z_LFT = " , z_LFT
write(40,297) "# z_RGT = " , z_RGT
write(40,299) "#             ||          Sz(t=0)=1         ||          Sz(t=0)=-1        |"
write(40,300) "#          time" , "i_LFT" , "i_RGT" , "i_LFT" , "i_RGT"
a = D_zero
b = D_zero
c = D_zero
d = D_zero
do i = 1 , n_t - 1
    a = a + t_rate * ( i_total(i+1,1) + i_total(i,1) ) / 2.0d0
    b = b + t_rate * ( i_total(i+1,2) + i_total(i,2) ) / 2.0d0
    c = c + t_rate * ( i_total(i+1,3) + i_total(i,3) ) / 2.0d0
    d = d + t_rate * ( i_total(i+1,4) + i_total(i,4) ) / 2.0d0
    write(40,301) dfloat(i-1) * t_rate + t_rate / 2.0d0 , a , b , c , d
end do
close(40)

open( unit=40 , file="SO_Area.dat" , action="write" , status="unknown" )
write(40,297) "# z_LFT = " , z_LFT
write(40,297) "# z_RGT = " , z_RGT
write(40,299) "#             ||          Sz(t=0)=1         ||          Sz(t=0)=-1        |"
write(40,300) "#          time" , "i_LFT" , "i_RGT" , "i_LFT" , "i_RGT"
a = D_zero
b = D_zero
c = D_zero
d = D_zero
do i = 1 , n_t - 1
    a = a + t_rate * ( i_SO(i+1,1) + i_SO(i,1) ) / 2.0d0
    b = b + t_rate * ( i_SO(i+1,2) + i_SO(i,2) ) / 2.0d0
    c = c + t_rate * ( i_SO(i+1,3) + i_SO(i,3) ) / 2.0d0
    d = d + t_rate * ( i_SO(i+1,4) + i_SO(i,4) ) / 2.0d0
    write(40,301) dfloat(i-1) * t_rate + t_rate / 2.0d0 , a , b , c , d
end do
close(40)

deallocate( Pops , MO_bra , MO_ket , AO_bra , AO_ket , DUAL_bra , DUAL_ket , phase , i_up , i_down , i_total , Sup , Sdown )

include 'formats.h'

297 format(a10,f10.3)
298 format(a135)
299 format(a75)
300 format(9a15)
301 format(f15.7,8ES15.5)
1000 format(3a21)
1001 format(3f21.6)

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

! local parameters ...
character(4) :: Sz(2) = ["up","down"]

! local variables ...
logical                       :: erg_done
integer                       :: nf , n , nc , it , spin
complex*16                    :: wp_energy
character(len=:) ,allocatable :: f_name


nc = merge(size(QDyn%fragments)+2,size(QDyn%fragments)+1,SOC)

erg_done = no
do spin = 1 , n_spin
    
    do n = 1 , n_part

        if( eh_tag(n) == "XX" ) cycle

        call FileName( f_name , n , spin , instance="dens" )
        open( unit = 52 , file = f_name , status = "replace" , action = "write" , position = "append" )
        write(52,15) "#" , ( nf+1 , nf=0,nc )  ! <== numbered columns for your eyes only ...

        if( SOC ) then
            write(52,12) "#" , "time" , QDyn%fragments , adjustr(Sz(spin)) , "u+d"
        else
            write(52,12) "#" , "time" , QDyn%fragments , "total"
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
    write(54,12) "#" , "time" , QDyn%fragments , "total"

    do it = 1 , n_t

        write(54,13) QDyn%dyn(it,0,1,1) , ( QDyn%dyn(it,nf,1,1) + QDyn%dyn(it,nf,1,2) , nf=1,nc-1 )

    end do

    close(54)

    if( n_part == 2 ) then

        open( unit = 54 , file = "dyn.trunk/"//"hl"//"_dens.dat" , status = "replace" , action = "write" , position = "append" )

        write(54,15) "#" , ( nf+1 , nf=0,nc-1 )  ! <== numbered columns for your eyes only ...
        write(54,12) "#" , "time" , QDyn%fragments , "total"

        do it = 1 , n_t

            write(54,13) QDyn%dyn(it,0,2,1) , ( QDyn%dyn(it,nf,2,1) + QDyn%dyn(it,nf,2,2) , nf=1,nc-1 )

        end do

        close(54)

    end if
end if

12 FORMAT(/A,15A11)
13 FORMAT(F13.6,14F11.5)
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

 module Data_Output

    use type_m
    use parameters_m        , only  : driver ,      &
                                      n_part ,      &
                                      spectrum ,    &
                                      survival ,    &
                                      NetCharge ,   &
                                      CH_and_DP_step
    use tuning_m            , only  : eh_tag                                      
    use Babel_m             , only  : System_Characteristics    
    use Structure_Builder   , only  : system => Extended_Cell    

    private

    public :: Populations , Dump_stuff , Net_Charge

    interface Populations
        module procedure Populations_vct , Populations_mtx
    end interface Populations

    ! module variables ...
    integer                , save :: counter = 0
    real*8  , allocatable  , save :: Net_Charge(:)

 contains
!
!
!
!==================================================================
 function Populations_vct( QDyn_fragments , basis , bra , ket , t )
!==================================================================
 implicit none
 character(*)    , intent(in)  :: QDyn_fragments(:)
 type(STO_basis) , intent(in)  :: basis(:)
 complex*16      , intent(in)  :: bra(:) , ket(:)
 real*8          , intent(in)  :: t
 real*8                        :: Populations_vct( 0:size(QDyn_fragments)+1 )

! local variables ...
integer             :: nf , N_of_fragments , ati
character(len=1)    :: fragment 

!-------------------------------------------------------------
!              get time-dependent Populations
!
! storage:  Populations_vct( 0:size(QDyn_fragments)+1 )
!
! Populations_vct( [time] + [# of fragments] + [total dens] )
!-------------------------------------------------------------

! time of population ...
Populations_vct(0) = t

! partial populations ...
N_of_fragments = size( QDyn_fragments )

do nf = 1 , N_of_fragments

    fragment = QDyn_fragments (nf)

    Populations_vct(nf) = pop_Slater( basis , bra(:) , ket(:) , fragment )

end do

! total population ...
Populations_vct(N_of_fragments+1) = pop_Slater( basis , bra(:) , ket(:) )

! atomic net-charge ...
do ati = 1 , system%atoms
    Net_Charge(ati) = abs( sum( bra(:)*ket(:) , basis(:)%atom == ati ) )
end do

! dump atomic net-charges for visualization
If ( NetCharge .AND. (mod(counter,CH_and_DP_step)==0) ) CALL dump_NetCharge (t) 

counter = counter + 1

!---------------------------------------------------- 

end function Populations_vct
!
!
!
!==================================================================
 function Populations_mtx( QDyn_fragments , basis , bra , ket , t )
!==================================================================
 implicit none
 character(*)    , intent(in)  :: QDyn_fragments(:)
 type(STO_basis) , intent(in)  :: basis(:)
 complex*16      , intent(in)  :: bra(:,:) , ket(:,:)
 real*8          , intent(in)  :: t
 real*8                        :: Populations_mtx( 0:size(QDyn_fragments)+1 , n_part)

! local variables ...
integer             :: n , nf , N_of_fragments , ati
real*8              :: charge_El , charge_Hl
character(len=1)    :: fragment 

!-----------------------------------------------------------------------
!              get time-dependent Populations
!
! storage:  Populations_mtx( 0:size(QDyn_fragments)+1 , n_part)
!
! Populations_mtx( [time] + [# of fragments] + [total dens] , el_hl )
!-----------------------------------------------------------------------

! time of population ...
Populations_mtx(0,:) = t

! partial populations ...
N_of_fragments = size( QDyn_fragments )

do n = 1 , n_part

    do nf = 1 , N_of_fragments

        fragment = QDyn_fragments (nf)

        Populations_mtx( nf , n ) = pop_Slater( basis , bra(:,n) , ket(:,n) , fragment )

    end do

    ! total populations ...
    Populations_mtx( N_of_fragments+1 , n ) = pop_Slater( basis , bra(:,n) , ket(:,n) )

end do

! atomic net-charge ...
do ati = 1 , system%atoms
    charge_El = abs( sum( bra(:,1)*ket(:,1) , basis(:)%atom == ati ) )
    charge_Hl = abs( sum( bra(:,2)*ket(:,2) , basis(:)%atom == ati ) )

    Net_Charge(ati) = charge_HL - charge_EL
end do

! dump atomic net-charges for visualization
If ( NetCharge .AND. (mod(counter,CH_and_DP_step)==0) ) CALL dump_NetCharge (t) 

counter = counter + 1

!---------------------------------------------------- 

end function Populations_mtx 
!
!
!
!==================================================
 subroutine Dump_stuff( TDOS , PDOS , SPEC , QDyn ) 
!==================================================
implicit none
type(f_grid)  , intent(in)     , optional  :: TDOS
type(f_grid)  , intent(in)     , optional  :: PDOS(:)
type(f_grid)  , intent(in)     , optional  :: SPEC
type(f_time)  , intent(in)     , optional  :: QDyn

! local variables ...
integer         :: i , nr , nf , np , N_of_residues , N_of_fragments
character(22)   :: string

! save TDOS ...
If( present(TDOS) ) then
    OPEN( unit=3 , file='dos.trunk/TDOS.dat' , status='unknown' )
        do i = 1 , size(TDOS%func)
            write(3,10) TDOS%grid(i) , TDOS%average(i) , TDOS%peaks(i) , TDOS%occupation(i)
        end do
    CLOSE(3)
end if

! save PDOS ...
If( present(PDOS) ) then
    N_of_residues = size( PDOS )
    do nr = 1 , N_of_residues
        string = "dos.trunk/PDOS-"//PDOS(nr)%residue//".dat"
        OPEN( unit=3 , file=string , status='unknown' )
            do i = 1 , size(PDOS(nr)%func)
                write(3,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i) , PDOS(nr)%peaks(i) , PDOS(nr)%occupation(i)
            end do
        CLOSE(3)
    end do
end if

! save peak and broadened specs ...
If( present(SPEC) ) then
    OPEN( unit=3 , file='dos.trunk/spectrum.dat' , status='unknown' )
        do i = 1 , size(SPEC%func)
            write(3,11) SPEC%grid(i) , SPEC%average(i) , SPEC%peaks(i)
        end do
    CLOSE(3)
end if

! save time-dependent electron or hole populations ...
If( present(Qdyn) ) then

    N_of_fragments = size( QDyn%fragments )

        do np = 1 , n_part

            if( eh_tag(np) == "XX" ) cycle

            OPEN( unit=3 , file="dyn.trunk/"//eh_tag(np)//"_survival.dat" , status="unknown" )

                write(3,14) "#" ,( nf+1 , nf=0,size(QDyn%fragments)+1 )  ! <== numbered columns for your eyes only ...
                write(3,12) "#" , QDyn%fragments , "total"
                do i = 1 , size( QDyn%dyn(:,1,1) )
                    write(3,13) ( QDyn%dyn(i,nf,np) , nf=0,N_of_fragments+1 )
                end do

            CLOSE(3)

            If( driver == "avrg_confgs" ) then  ! <== also save standard deviations ...

                  OPEN( unit=4 , file="dyn.trunk/"//eh_tag(np)//"_std.dat" , status="unknown" )

                      write(4,14) "#" ,( nf+1 , nf=0,size(QDyn%fragments)+1 )  ! <== numbered columns for your eyes only ...
                      write(4,12) "#" , QDyn%fragments , "total"
                      do i = 1 , size( QDyn%std(:,1,1) )
                          write(4,13) ( QDyn%std(i,nf,np) , nf=0,N_of_fragments+1 )
                      end do

                  CLOSE(4)

            End If

        end do
end if

10   FORMAT(4F12.5)
11   FORMAT(3F13.9)
12   FORMAT(/15A10)
13   FORMAT(F11.6,14F10.5)
14 FORMAT(A,I9,14I10)

end subroutine Dump_stuff
!
!
!
!=================================================
 function pop_Slater( basis , za , zb , fragment )
!=================================================
implicit none
type(STO_basis)              , intent(in) :: basis(:)
complex*16                   , intent(in) :: za(:) , zb(:)
character(*)     , optional  , intent(in) :: fragment

! local variables
real*8       :: pop_Slater
complex*16   :: pop 

pop = C_zero

if( present(fragment) ) then
    pop = sum( za(:) * zb(:) , mask = basis%fragment == fragment )
else
    pop = sum( za(:) * zb(:) )
end if

pop_Slater = real( pop )

end function
!
!
!
!===============================
 subroutine dump_NetCharge ( t )
!===============================
implicit none
real*8  , intent(in) :: t

! local variables ...
integer :: ati , i , j 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!saving net_charge ...
OPEN(unit=21 , file="dyn.trunk/NetCharge.inpt" , status = "unknown", action = "write" , position = "append" )
do ati = 1 , system%atoms
    write(21,'(F9.5)',advance='no') net_charge(ati) 
end do
close(21)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
OPEN(unit=24 , file="dyn.trunk/CH-DP-frames.pdb" , status = "unknown", action = "write" , position = "append" )

If( counter == 0 ) write(24,6) System_Characteristics

write(24,4) 'REMARK' , 'manipulated by charge-transfer'
write(24,5) 'TITLE'  , 'manipulated by charge-transfer     t= ', t
write(24,4) 'REMARK' , 'manipulated by charge-transfer'
write(24,1) 'CRYST1' , system%T_xyz(1) , system%T_xyz(2) , system%T_xyz(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
write(24,3) 'MODEL'  , counter

do i = 1 , system%atoms

             write(24,2)'HETATM'                        ,  &    ! <== non-standard atom
                        i                               ,  &    ! <== global number
                        system%Symbol(i)                ,  &    ! <== atom type
                        system%residue(i)               ,  &    ! <== residue name
                        system%nr(i)                    ,  &    ! <== residue sequence number
                        ( system%coord(i,j) , j=1,3 )   ,  &    ! <== xyz coordinates
                        net_charge(i)                   ,  &    ! <== wavepacket occupancy
                        0.d0                            ,  &    
                        system%Symbol(i)                        ! <== chemical element symbol

end do

write(24,'(a)') 'TER'
write(24,'(a)') 'ENDMDL'

close(24)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3,t56,f8.5,t65,f8.5,t77,a2)
3 FORMAT(a6,i9,11i7)
4 FORMAT(a6,t15,a31)
5 FORMAT(a5,t15,a39,f9.5)
6 FORMAT(a6,a72)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

end subroutine dump_NetCharge
!
!
!
end module Data_Output

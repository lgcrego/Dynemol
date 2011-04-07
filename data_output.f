 module Data_Output

    use type_m
    use parameters_m        , only  : n_part ,      &
                                      spectrum ,    &
                                      survival
    use FMO_m               , only  : eh_tag                                      

    public :: Populations   , Dump_stuff

    private

    interface Populations
        module procedure Populations_vct , Populations_mtx
    end interface Populations

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
integer             :: nf , N_of_fragments
character(len=1)    :: fragment 

!----------------------------------------------------------
!              get time-dependent Populations
!----------------------------------------------------------

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
integer             :: n , nf , N_of_fragments
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

!---------------------------------------------------- 

end function Populations_mtx 
!
!
!
!==========================================================
 subroutine Dump_stuff( TDOS , PDOS , SPEC , QDyn ) 
!==========================================================
implicit none
type(f_grid)  , intent(in)     , optional  :: TDOS
type(f_grid)  , intent(in)     , optional  :: PDOS(:)
type(f_grid)  , intent(in)     , optional  :: SPEC
type(newf_time)  , intent(in)     , optional  :: QDyn

! local variables ...
integer         :: i , nr , nf , np , N_of_residues , N_of_fragments
real*8          :: t 
character(12)   :: string

! save TDOS ...
If( present(TDOS) ) then
    OPEN( unit=3 , file='TDOS.dat' , status='unknown' )
        do i = 1 , size(TDOS%func)
            write(3,10) TDOS%grid(i) , TDOS%average(i) , TDOS%peaks(i) , TDOS%occupation(i)
        end do
    CLOSE(3)
end if

! save PDOS ...
If( present(PDOS) ) then
    N_of_residues = size( PDOS )
    do nr = 1 , N_of_residues
        string = "PDOS-"//PDOS(nr)%residue//".dat" 
        OPEN( unit=3 , file=string , status='unknown' )
            do i = 1 , size(PDOS(nr)%func)
                write(3,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i) , PDOS(nr)%peaks(i) , PDOS(nr)%occupation(i)
            end do
        CLOSE(3)
    end do
end if

! save peak and broadened specs ...
If( spectrum ) then
    OPEN( unit=3 , file='spectrum.dat' , status='unknown' )
        do i = 1 , size(SPEC%func)
            write(3,11) SPEC%grid(i) , SPEC%average(i) , SPEC%peaks(i)
        end do
    CLOSE(3)
end if

! save time-dependent electron or hole populations ...
If( survival ) then

    N_of_fragments = size( QDyn%fragments )

        do np = 1 , n_part

            select case ( eh_tag(np) )
                case ( "el" )
                    OPEN( unit=3 , file="el_survival.dat" , status="unknown" )
                case ( "hl" )
                    OPEN( unit=3 , file="hl_survival.dat" , status="unknown" )
            end select

            write(3,12) "#" , QDyn%fragments , "total"
            do i = 1 , size( QDyn%dyn(:,1,1) )
                write(3,13) ( QDyn%dyn(i,nf,np) , nf=0,N_of_fragments+1 )
            end do

            CLOSE(3)

        end do
end if

10   FORMAT(4F12.5)
11   FORMAT(3F13.9)
12   FORMAT(10A9)
13   FORMAT(10F9.4)

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
end module Data_Output

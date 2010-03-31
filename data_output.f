 module Data_Output

    use type_m
    use projectors

    public :: get_Populations , Dump_stuff

    private

 contains
!
!
!
!==================================================================
 function get_Populations(system,basis,bra,ket) result(Populations)
!==================================================================
 implicit none
 type(structure) , intent(in)       :: system
 type(STO_basis) , intent(in)       :: basis(:)
 complex*16      , intent(in)       :: bra(:,:) , ket(:,:)
 real*8                             :: Populations( size(system%list_of_fragments)+1 )

! local variables ...
integer             :: nf , N_of_fragments
character(len=1)    :: fragment

!----------------------------------------------------------
!              get time-dependent Populations
!----------------------------------------------------------

! partial populations ...

N_of_fragments = size( system % list_of_fragments )

do nf = 1 , N_of_fragments

    fragment = system % list_of_fragments (nf)

    Populations(nf) = pop_Slater(basis,bra,ket,fragment)

end do

! total population ...

 Populations(nf+1) = pop_Slater(basis,bra,ket)

! ---------------------------------------------------- 

 end function get_Populations 
!
!
!
!=========================================================================================
 subroutine Dump_stuff( TDOS , PDOS , SPEC , QDyn , list_of_fragments ) 
!=========================================================================================
implicit none
type(f_grid)  , intent(in)  , optional  :: TDOS
type(f_grid)  , intent(in)  , optional  :: PDOS(:)
type(f_grid)  , intent(in)  , optional  :: SPEC
real*8        , intent(in)  , optional  :: QDyn(:,:)
character(*)  , intent(in)  , optional  ::list_of_fragments(:)

! local variables ...
integer         :: i , nr , nf , N_of_residues , N_of_fragments
real*8          :: t , t_rate
character(12)   :: string

! save TDOS ...
If( present(TDOS) ) then
    OPEN( unit=3 , file='TDOS.dat' , status='unknown' )
        do i = 1 , size(TDOS%func)
            write(3,10) TDOS%grid(i) , TDOS%average(i)
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
                write(3,10) PDOS(nr)%grid(i) , PDOS(nr)%average(i)
            end do
        CLOSE(3)
    end do
end if

! save peak and broadened specs ...
If( spectrum ) then
    OPEN( unit=3 , file='spectrum.dat' , status='unknown' )
        do i = 1 , size(SPEC%func)
            write(3,11) SPEC%grid(i) , SPEC%average(i) 
        end do
    CLOSE(3)
end if

! save time-dependent populations ...
If( survival ) then
    t_rate = (t_f - t_i) / float(n_t)
    N_of_fragments = size( list_of_fragments )
    OPEN( unit=3 , file="survival.dat" , status="unknown" )
    write(3,12) "#" , list_of_fragments , "total"
    do i = 1 , n_t
        t = t_rate * (i-1)
        write(3,13) t , ( QDyn(i,nf) , nf=1,N_of_fragments+1 )
    end do
    CLOSE(3)
end if

10   FORMAT(2F12.5)
11   FORMAT(3F13.9)
12   FORMAT(10A9)
13   FORMAT(10F9.4)

 end subroutine Dump_stuff
!
!
!
end module Data_Output

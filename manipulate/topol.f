module Topology_routines

use types_m
use constants_m

public :: connect , dump_topol

private

contains
!
!
!=========================
 subroutine connect( sys )
!=========================
implicit none
type(universe)                  , intent(inout)   :: sys 

! local varibles ...
character(len=1)                :: yn
character(len=2)                :: S1 , S2
integer                         :: i , j
real*8                          :: cutoff
logical                         :: flag

CALL system( "clear" )

write(*,'(/a)',advance='no') ">>> Connect Specific Bonds ? (y/n) "
read (*,'(a)') yn

if ( yn == "y" ) then

    allocate( sys % topol (sys%N_of_atoms,sys%N_of_atoms) , source = .false. )

    do

        write(*,'(1x,3/a)') "Choose chemical elements whose bonds are to be connected (@ to exit) : "
        read*, S1
        If( S1 == "@" ) exit
        read*, S2

        write(*,'(1x,a)') "cut-off distance for the bond: "
        read*, cutoff

        do j = 1 , sys%N_of_atoms 
            do i = j+1 , sys%N_of_atoms 

                flag = (sys%atom(i)%Symbol==S1 .AND. sys%atom(j)%Symbol==S2)       &
                    .OR.                                                           &  
                       (sys%atom(i)%Symbol==S2 .AND. sys%atom(j)%Symbol==S1) 

                If( flag .AND. sqrt(sum((sys%atom(i)%xyz - sys%atom(j)%xyz)**2)) < cutoff ) then
        
                    sys % topol(i,j) = .true. 
                    sys % topol(j,i) = .true.

                end If

            end do
        end do

    end do

end If

CALL system( "clear" )
    
end subroutine connect

!
!
!========================================
 subroutine dump_topol( sys , file_unit )
!========================================
implicit none
type(universe) , intent(in) :: sys 
integer        , intent(in) :: file_unit

! local varibles ...
integer :: i , j

do j = 1 , sys%N_of_atoms

    if( any(sys%topol(:,j) ) ) then

        write(file_unit,'(/A6,I5)', advance='no') "CONECT" , j 

        do i = 1 , sys % N_of_atoms
            If( sys % topol(i,j) ) write(file_unit,'(I5)', advance='no') i
        end do

    end if

end do
    
write(file_unit,'(/)')
    
end subroutine dump_topol

end module Topology_routines

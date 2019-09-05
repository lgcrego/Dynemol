module Multi_Trajectory_routines

use types_m
use constants_m
use Read_Parms
use RW_routines         , only : view_XYZ , Initialize_System
use Function_routines   , only : ad_hoc_tuning
use diagnosis_m         , only : diagnosis
use EDIT_routines       , only : Eliminate_Fragment 
use Statistics_routines , only : Most_Representative_Configuration

public :: Work_Multiple_Trajectories

private

contains
!
!
!====================================
subroutine Work_Multiple_Trajectories
!====================================
implicit none

! local varibles ...
character(len=1)                :: file_format , YorN , choice , wait
integer                         :: frame , n_frames , i , j
real*8                          :: delta_t
type(universe)  , allocatable   :: trj1(:) , trj2(:)

CALL system( "clear" )

write(*,'(/a)') ' This routine reads the files >frames1< and >frames2< : '
write(*,'(/a)') ' Choose the format : '
write(*,'(/a)') ' (p) = PDB '
write(*,'(/a)') ' (v) = VASP (not implemented)'
write(*,'(/a)') ' (d) = DICE (not implemented)'
write(*,'(/a)') ' (x) = XYZ  (not implemented)'
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') file_format

select case( file_format )

    case( 'p' ) 
        CALL Read_PDB_Trajectories(trj1,"frames1.pdb") 
        CALL Read_PDB_Trajectories(trj2,"frames2.pdb") 

end select

! select frame OR trajectory ...
do
    write(*,'(/a)') ' (1) = edit trajectory with AD-HOC'
    write(*,'(/a)') ' (2) = Nuclear pair-correlations '      
    write(*,'(/a)') ' (3) = DONE                '
    write(*,'(/a)',advance='no') '>>>   '
    read (*,'(a)') choice 

    select case( choice )

        case( '1' ) 

            do i = 1 , size(trj1)

                CALL ad_hoc_tuning( trj1(i) )
                trj1(i)%time = trj1(i)%time * (i-1)

                CALL ad_hoc_tuning( trj2(i) )
                trj2(i)%time = trj2(i)%time * (i-1)

            end do

        case( '2' ) 

            CALL Nuclear_Pair_Correlation( trj1 , trj2 )

        case( '3' ) 

        case default
            exit

    end select

end do    

write(*,'(/a)',advance='no') 'press ENTER '
read (*,'(a)') wait

end subroutine Work_Multiple_Trajectories
!
!
!===================================================
 subroutine Read_PDB_Trajectories( trj , file_name )
!===================================================
type(universe)  , allocatable   , intent(out)   :: trj(:)
character(len=11)               , intent(in)    :: file_name

! local variables ...
integer                         :: openstatus , inputstatus , i , j , k , model , number_of_atoms , n , m , dumb_number
real*8                          :: time_1 , time_2 , delta_t 
character(1)                    :: test
character(4)                    :: keyword
character(5)                    :: MMSymbol_char
character(72)                   :: Surface_Characteristics

open(unit = 31, file = file_name, status = 'old', action = 'read', iostat = openstatus)
if (openstatus > 0) then
    print*, " *** Cannot open the file *** ",file_name
    stop
end if

read( unit = 31 , fmt = 31 ) Surface_Characteristics

! find the number of model frames ...
model=0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( inputstatus /= 0 ) exit
    if ( keyword == 'MODE' ) model = model + 1
end do

! return to the top of the file ...
rewind 31

! read number the atoms and time ...
read(unit = 31, fmt = 35, iostat = inputstatus) keyword
do
    if ( keyword == 'MODE' ) then
        exit
    else
        if ( keyword == 'TITL' ) then
            backspace 31
            do 
                read(unit = 31, fmt = 42,advance='no',iostat = inputstatus) test
                if(test == "=") then
                    read(unit = 31, fmt = 41, iostat = inputstatus) time_1 
                    exit    ! <== time_1 read
                end if
            end do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
        else
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
        end if
    end if
end do
number_of_atoms = 0
do
    read(unit = 31, fmt = 35, iostat = inputstatus) keyword
    if ( keyword == 'ATOM' ) then
        backspace 31
        read(unit = 31, fmt = 32, iostat = inputstatus) dumb_number
        number_of_atoms = number_of_atoms + 1
    else
        if ( keyword == 'TITL' ) then
            backspace 31
            do 
                read(unit = 31, fmt = 42,advance='no',iostat = inputstatus) test
                if(test == "=") then
                    read(unit = 31, fmt = 41, iostat = inputstatus) time_2 
                    exit    ! <== time_2 read
                end if
            end do
            exit    ! <== leave outer do loop
        end if
    end if
end do

delta_t = time_2 - time_1

! return to the top of the file ...
rewind 31

! allocate array of frames ...
allocate( trj(model) )

do j = 1 , model
    if( j == 1 ) then
        trj(j)%time = 0.d0
        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( trj(j)%box(i) , i=1,3 )
            end if
            if( keyword == 'MODE' ) exit
        end do
    
        allocate( trj(j)%atom(number_of_atoms) )
    
        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 33, iostat = inputstatus) MMSymbol_char ,             &
                                                            trj(j)%atom(i)%resid ,      &
                                                            trj(j)%atom(i)%nresid ,     &
                                                            ( trj(j)%atom(i)%xyz(k) , k=1,3 ) 

            trj(j)%atom(i)%MMSymbol = adjustl(MMSymbol_char)

        end do
        CALL MMSymbol_2_Symbol(trj(j)%atom)
        CALL Changing_Fragment(trj(j)%atom%resid,trj(j)%atom)
    else
        trj(j)%time = trj(j-1)%time + delta_t
        do
            read(unit = 31, fmt = 35, iostat = inputstatus) keyword
            if( keyword == 'CRYS' ) then
                backspace 31
                read(unit = 31, fmt = 40, iostat = inputstatus) ( trj(j)%box(i) , i=1,3 )
            end if
            if ( keyword == 'MODE' ) exit
        end do
    
        allocate( trj(j)%atom(number_of_atoms) )
    
        do i = 1 , number_of_atoms
            read(unit = 31, fmt = 37, iostat = inputstatus) ( trj(j)%atom(i)%xyz(k) , k=1,3 )
        end do
    end if
end do

trj%Surface_Characteristics = Surface_Characteristics
trj%N_of_Atoms = number_of_atoms

forall(i = 2:model )
    trj(i)%atom%MMSymbol  = trj(1)%atom%MMSymbol
    trj(i)%atom%resid     = trj(1)%atom%resid
    trj(i)%atom%nresid    = trj(1)%atom%nresid
    trj(i)%atom%Symbol    = trj(1)%atom%Symbol
    trj(i)%atom%fragment  = trj(1)%atom%fragment
end forall

! get list of residues in trj ...
CALL Identify_Residues( trj(1) )

! GROUP residues ...
do i = 1 , size(trj)
    CALL Pack_Residues( trj(i)%atom , trj(1)%list_of_residues ) 
end do

! Defining N_of_atoms and N_of_Solvent_molecules ...
trj%N_of_atoms = number_of_atoms
trj%N_of_Solvent_molecules = maxval( trj(1)%atom%nresid , trj(1)%atom%fragment == 'S' ) - minval( trj(1)%atom%nresid , trj(1)%atom%fragment == 'S' ) + 1

close(31)

! Formats ...
31 format(a72)
32 format(5x, i6)
33  format(t12,a5,t18,a3,t23,i7,t31,f8.3,t39,f8.3,t47,f8.3)                
35 format(a4)
36 format(7x, i7)
37 format(t32, f8.3, t40, f8.3, t48, f8.3)
38 format(10x, a1)
39 format(81x, f7.0)
40 format(6x, 3f9.3)
41 format(f10.5)
42 format(a1)

end subroutine Read_PDB_Trajectories
!
!
!
!==================================
subroutine Changing_Fragment(FGT,a)
!==================================
implicit none
character(3)    , intent(in)  :: FGT(:)
type(atomic)    , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(FGT(i))
        case( 'CCC') 
            a(i)%fragment = 'C' 
        case( 'FFF') 
            a(i)%fragment = 'F' 
        case( 'ACN') 
            a(i)%fragment = 'S' 
        case( 'BLK') 
            a(i)%fragment = 'B' 
        case( 'SRF') 
            a(i)%fragment = 'I' 
    end select

 END DO

end subroutine Changing_Fragment
!
!
!
!
!
!==================================================
 subroutine Nuclear_Pair_Correlation( trj1 , trj2 )
!==================================================
implicit none
type(universe)  , allocatable   :: trj1(:) , trj2(:)

! local varibles ...
integer         :: step 
character(1)    :: choice 
character(2)    :: Symbol
character(3)    :: MMSymbol

real*8  , allocatable :: d1_ij(:) , d2_ij(:) , C_ij(:) , C_ij_mean(:) , C_ij_2mean(:) , Sigma(:) , numerator(:) , denominator(:)
real*8                :: weight, a_mean , b_mean , aa_mean , bb_mean , ab_mean
integer               :: N_of_samples , N_of_atoms , length_trj , i , j, k , n , counter

CALL system( "clear" )

write(*,'(/a)') ' Choose criterium '
write(*,'(/a)') ' (1) = entire structure '
write(*,'(/a)') ' (2) = Chemical element '
write(*,'(/a)') ' (3) = MMSymbol '
write(*,'(/a)') ' (4) = atom index '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') choice

select case( choice )

    case( '1' )
        write(*, '(1x,a)'               ) "Frame step: "
        write(*, '(1x,a)', advance = 'no') 'step = '
        read*, step

    case( '2' )
        write(*, '(1x,a)',advance = 'no') "Chemical Symbol of elements to be analyzed: "
        read*, Symbol
        write(*, '(1x,a)'               ) "Frame step: "
        write(*, '(1x,a)', advance = 'no') 'step = '
        read*, step

    case( '3' )
        write(*, '(1x,a)',advance = 'no') "MMSymbol of elements to be analyzed: "
        read*, MMSymbol
        write(*, '(1x,a)'               ) "Frame step: "
        write(*, '(1x,a)', advance = 'no') 'step = '
        read*, step

end select

If( size(trj1) /= size(trj2) ) stop " >>> trajectories have different lengths <<< "
length_trj = size(trj1)

! determines the number of sample points ...
N_of_samples = sum( [( 1 , j=1,length_trj,step )] )

allocate( d1_ij(N_of_samples) , source = 0.d0 )
allocate( d2_ij(N_of_samples) , source = 0.d0 )

allocate( C_ij       (N_of_samples) , source = 0.d0 )
allocate( C_ij_2mean (N_of_samples) , source = 0.d0 )
allocate( C_ij_mean  (N_of_samples) , source = 0.d0 )

allocate( numerator  (N_of_samples) , source = 0.d0 )
allocate( denominator(N_of_samples) , source = 0.d0 )

allocate( Sigma (N_of_samples) , source = 0.d0 )

! for all atoms in the structure ...
N_of_atoms = trj1(1) % N_of_atoms

counter = N_of_atoms * (N_of_atoms - 1) / 2

do i = 1   , N_of_atoms
do j = i+1 , N_of_atoms

    n = 1
    do k = 1 , length_trj , step

        d1_ij(n) = dsqrt( sum( (trj1(k) % atom(i) % xyz(:) - trj1(k) % atom(j) % xyz(:))**2 ) )
        d2_ij(n) = dsqrt( sum( (trj2(k) % atom(i) % xyz(:) - trj2(k) % atom(j) % xyz(:))**2 ) )
    
        n = n + 1

    end do
  
    do n = 1 , N_of_samples

        a_mean = sum( d1_ij(1:n) ) / float(n)
        b_mean = sum( d2_ij(1:n) ) / float(n)

        aa_mean = sum( d1_ij(1:n) * d1_ij(1:n) ) / float(n)
        bb_mean = sum( d2_ij(1:n) * d2_ij(1:n) ) / float(n)
        ab_mean = sum( d1_ij(1:n) * d2_ij(1:n) ) / float(n)

        denominator(n) = sqrt( (aa_mean - a_mean*a_mean) * (bb_mean - b_mean*b_mean) )

        numerator(n) = (ab_mean - a_mean*b_mean) 
       
    end do

    where( abs(numerator-denominator) < low_prec )
        C_ij = 1.d0
    elsewhere
        C_ij = numerator / denominator
    endwhere


    do n = 1 , N_of_samples

        C_ij_2mean(n) = C_ij_2mean(n) + C_ij(n)*C_ij(n) / counter

        C_ij_mean(n) = C_ij_mean(n) + C_ij(n) / counter

        Sigma(n) = sqrt( C_ij_2mean(n) - C_ij_mean(n)*C_ij_mean(n) ) 

    end do


end do
end do 


do n = 1 , N_of_samples
    write(33,*) n , Sigma(n)
end do
    

deallocate( d1_ij, d2_ij , C_ij , C_ij_mean , C_ij_2mean , Sigma , numerator , denominator )

end subroutine Nuclear_Pair_Correlation
!
!
!
end module Multi_Trajectory_routines

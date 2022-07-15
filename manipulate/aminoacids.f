module Aminoacids

use types_m
use constants_m
use Read_Parms      , only : atomic , atomic_mass 

public :: AminoacidStuff

    ! module variables ...
    type(universe) :: inpt_list

private

contains
!
!
!
!================================
 subroutine AminoacidStuff( sys )
!================================
implicit none
type(universe) , intent(inout)   :: sys 

! local varibles ...
integer :: choice
character(len=1)  :: wait

CALL system('clear')

write(*,'(/a)') ' (1)  = Sort and Trim Aminoacids'
write(*,'(/a)') ' (2)  = Distance between residues '      
write(*,'(/a)') ' (0)  = DONE                '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(I)') choice 

select case( choice )

    case( 1 ) 
        CALL SortAminoacids( sys )

    case( 2 )
        CALL QM_2_MM_distance( sys )

end select

write(*,'(/a)',advance='no') 'press ENTER '
read (*,'(a)') wait

end subroutine AminoacidStuff
!
!
!
!================================
 subroutine SortAminoacids( sys )
!================================
implicit none
type(universe) , intent(inout)   :: sys 

! local varibles ...
integer                        :: i , j , k , L , n , i1 , i2 , max_12 , max_21
integer          , allocatable :: NewOrder(:) , tmp(:)
character(4)     , allocatable :: mismatch(:,:) 
character(7)     , allocatable :: tag(:)
logical                        :: flag 
type(universe)                 :: OPLS
type(universe)                 :: OPLS_list

! Build buffer OPLS ...
CALL Build_OPLS_Aminoacids( OPLS )

! Fix isotope-residues of input.pdb ...
CALL Fix_residues( sys , OPLS )

allocate( mismatch(sys%N_of_atoms,2) , source = 'XXXX'   )
allocate( tag     (sys%N_of_atoms  ) , source = 'XXXXXXX')
allocate( NewOrder(sys%N_of_atoms  ) , source = 0        )

allocate( OPLS_list% atom(sys%N_of_atoms) )

n = 0
i = 1
do while( i <= sys%N_of_atoms )

    ! here, identify the atom's residue as OPLS% amino(k)% residAA ...
    k = sum( [ (j , j=1,OPLS%N_of_aminos) ] , (OPLS% amino(:)% residAA == sys% atom(i)% resid) )
    If( k == 0 ) then  ! <== residue not found ...
         Print*, ">>> residue [",sys% atom(i)% resid,"] not found <<<"
         i = i + 1
         cycle
    end If

    ! here, verify whether atom MMSymbol is found in OPLS% amino(k) ...
    L = sum( [ (j , j=1,size(OPLS%amino(k)%atom)) ] , OPLS% amino(k)% atom(:)% AASymbol == sys% atom(i)% MMSymbol ) 
    If( L == 0 ) then  ! MMSymbol not found in OPLS% amino(k) ...
         write(4,*) k , "[",OPLS%amino(k)%residAA,"]","[",sys%atom(i)%resid,"]", i , " => ","[",sys% atom(i)% MMSymbol,"]" 

         ! if .true. pair is already in mismatch list ...
         flag = any( tag == sys%atom(i)%MMSymbol//sys% atom(i)% resid )
         If( .NOT. flag ) then
              n = n + 1
              mismatch(n,1) = sys% atom(i)% MMSymbol
              mismatch(n,2) = sys% atom(i)% resid

              tag(n) = sys% atom(i)% MMSymbol//sys% atom(i)% resid
         end if
    end If
    NewOrder(i) = L

    i = i + 1
end do

! Build reference inpt ...
i1 = 1
do i = 1 , size(inpt_list%amino)

    ! here, identify the atom's residue as OPLS% amino(k)% residAA ...
    k = sum( [ (j , j=1,OPLS%N_of_aminos) ] , (OPLS% amino(:)% residAA == sys% atom(i1)% resid) )

    i2 = (i1-1) + OPLS% amino(k)% N_of_atoms 

    OPLS_list% atom(i1:i2)% AASymbol = OPLS% amino(k)% atom(:)% AASymbol
    OPLS_list% atom(i1:i2)% indx     = [(j , j=1,OPLS% amino(k)% N_of_atoms)] 

    ! here, identify atoms in OPLS list that do not have a match in input.pdb ...
    forall( j=i1:i2 , .NOT. any( inpt_list% amino(i)% atom(:)% AASymbol == OPLS_list%atom(j)% AASymbol) ) &
    OPLS_list% atom(j)% indx = 0

    do j = i1 , i2

         ! if (input.pbd(j) == OPLS_list(:)) = .true., already matching, no work needed ...
         If(any(OPLS_list% atom(i1:i2)% AASymbol == sys% atom(j)% MMSymbol)) cycle

! for debuginng purposes ...
!        n = 1 
!        do L = i1 , i2
!           print*, n , "[",sys%atom(j)%resid,"]", &
!                   L , "[",OPLS_list%atom(L)%AASymbol,"]", &
!                   j , "[",sys%atom(j)%MMSymbol,"]", &
!                   - (OPLS_list%atom(L)%indx+NewOrder(j)) , &
!                   verify(OPLS_list%atom(L)%AASymbol , sys%atom(j)%MMSymbol), &
!                   verify(sys%atom(j)%MMSymbol , OPLS_list%atom(L)%AASymbol)
!           n= n+1
!        end do

         max_12 = maxval( verify(OPLS_list%atom(i1:i2)%AASymbol , sys%atom(j)%MMSymbol)&
                            - (OPLS_list%atom(i1:i2)%indx+NewOrder(j)) )

         max_21 = maxval( verify(sys%atom(j)%MMSymbol , OPLS_list%atom(i1:i2)%AASymbol)&
                            - (OPLS_list%atom(i1:i2)%indx+NewOrder(j)) )

         If( max_12 > max_21 ) then

              NewOrder(j) = maxloc( verify(OPLS_list%atom(i1:i2)%AASymbol , sys%atom(j)%MMSymbol)&
                                   - (OPLS_list%atom(i1:i2)%indx+NewOrder(j)) , dim=1 )
         else
              NewOrder(j) = maxloc( verify(sys%atom(j)%MMSymbol , OPLS_list%atom(i1:i2)%AASymbol)&
                                   - (OPLS_list%atom(i1:i2)%indx+NewOrder(j)) , dim=1 )
         end If
         ! indx holds the stack position of AASymbol in the corresponding OPLS residue ...
         OPLS_list%atom( (i1-1)+NewOrder(j) )%indx = NewOrder(j)

         ! changes input.pdb MMSymbol by corresponding OPLS AASymbol ...
         sys%atom(j)%MMSymbol = OPLS_list%atom( (i1-1)+NewOrder(j) )%AASymbol 

    end do 

    ! nresid in ascending order, starting from 1, nothing else changed ...
    sys%atom(i1:i2)%nresid = i

    ! turning NewOrder inside-out ...
    allocate( tmp(OPLS% amino(k)% N_of_atoms ), source = NewOrder(i1:i2) )
    forall( j = 1:size(tmp) , tmp(j) /= 0 ) NewOrder( tmp(j)+(i1-1) ) = j
    deallocate( tmp )

    i1 = i2 + 1
end do


i1 = inpt_list% amino(1)% N_of_atoms 
do i = 2 , size(inpt_list%amino)

    i2 = i1 + inpt_list% amino(i)% N_of_atoms 

    where( NewOrder(i1+1:i2) /= 0 ) &
    NewOrder(i1+1:i2) = NewOrder(i1+1:i2) + (inpt_list% amino(i)% indx - 1)

    i1 = i2

end do

write(*,'(/a)') '>>>  Saving seed.pdb  <<<'

CALL Dump_pdb( sys , NewOrder )

stop

 end subroutine SortAminoacids
!
!
!========================================
 subroutine Build_OPLS_Aminoacids( OPLS )
!========================================
implicit none
type(universe) , intent(out) :: OPLS

! local varibles ...
integer                    :: ioerr , inputstatus
integer                    ::  i , nres , indx
character(1)               :: firstChar
character(4) , allocatable :: AAsurrogate(:)
character(4)               :: AASymbol
character(5)               :: keyword1

OPEN(unit=3,file='aminoacids.rtp',form="formatted",access="sequential",status='old',iostat=ioerr,err=10)

allocate( OPLS%amino(100) )

nres = 0
do 
    read( 3 , '(t1,a)' , iostat = inputstatus ) firstChar
    if( inputstatus /=0 ) exit  !! <== end of file ...
    if( firstChar == "[" ) then
         read( 3 , '(t2,a)' , iostat = inputstatus ) firstChar
         if( firstChar == "[" ) then
              backspace 3
              backspace 3
              nres = nres + 1
              read(3,'(t3,a4)') OPLS%amino(nres)%residAA

              read( 3 , '(t4,a5)' , iostat = inputstatus ) keyword1  !! <== [ atoms ] ...
              allocate( AAsurrogate(50) )
              indx = 0
              do
                  read( 3 , * , iostat = inputstatus ) AASymbol
                  if( AASymbol == "[" ) exit
                  indx = indx + 1
                  AAsurrogate(indx) = AASymbol
              end do
              allocate( OPLS%amino(nres)%atom (indx) )
              OPLS% amino(nres)% atom(:)% AASymbol = AAsurrogate(1:indx)  !! <== aminoacid residue is stored ...
              OPLS% amino(nres)% atom(:)% indx     = [( i , i=1,indx)]    !! <== position within residue is stored ...
              OPLS% amino(nres)% N_of_atoms        = indx                 !! <== szie of residue is stored ...
              deallocate( AAsurrogate )
         end if
    end if
end do 
OPLS%N_of_aminos = nres

10 if( ioerr > 0 ) stop 'aminoacid databasis not found; terminating execution'

end subroutine Build_OPLS_Aminoacids
!
!
!
!=====================================
 subroutine Fix_residues( sys , OPLS )
!=====================================
implicit none
type(universe) , intent(inout) :: sys 
type(universe) , intent(in)    :: OPLS

! local varibles ...
integer :: i , j , i1 , indx
logical :: sizes_differ , vis_a_vis

! some preprocessing ...
where( sys% atom% resid == "CYS" ) sys% atom% resid = "CYS2"
where( sys% atom% resid == "HIS" ) sys% atom% resid = "HISH"

sys% atom% indx = [( i , i=1,sys%N_of_atoms)]

!----------------------------------------------------------------------------
! determines the N_of_atoms per residue in input.pdb
allocate( inpt_list% amino( sys%atom(sys%N_of_atoms)% nresid ) )
j  = 0  ; i1 = 1  ; i  = 1
do while( i <= sys%N_of_atoms )
    i = i + 1
    If( (sys%atom(i)%resid /= sys%atom(i1)%resid) .OR. (sys%atom(i)%nresid /= sys%atom(i1)%nresid) ) then
         j = j + 1
         inpt_list% amino(j)% residAA    = sys%atom(i-1)%resid  ! <== name of residue ...
         inpt_list% amino(j)% N_of_atoms = i-i1                 ! <== # of atoms comprising this residue in input.pdb ...
         inpt_list% amino(j)% indx       = i1                   ! <== starting position of this residue in input.pdb ...

         allocate( inpt_list%amino(j)%atom(i-i1) )
         inpt_list% amino(j)% atom(:)% AASymbol = sys%atom(i1:i-1)%MMSymbol   ! <== store corresponding MMSymbols ...

         i1 = i
    end if
end do         
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! look for differences in residue sizes for input.pdb and OPLS buffer ...
! for isotope-residues, correct residue names in input.pdb 
do j = 1 , size(inpt_list%amino)
    do i = 1 , OPLS% N_of_aminos

         ! comparing residues: same name but different sizes ...
         vis_a_vis    = ( OPLS% amino(i)% residAA    == inpt_list% amino(j)% residAA    )
         sizes_differ = ( OPLS% amino(i)% N_of_atoms /= inpt_list% amino(j)% N_of_atoms )

         If( vis_a_vis .AND. sizes_differ ) then

              print*,'[',OPLS%amino(i)%residAA,']','[',inpt_list% amino(j)% residAA,']'  &
                    , OPLS% amino(i)% N_of_atoms , inpt_list% amino(j)% N_of_atoms       &
                    , "input index = ", inpt_list% amino(j)% indx 

              indx = inpt_list% amino(j)% indx 

              select case ( inpt_list% amino(j)% residAA ) 

                   case( "GLU" )
                        select case ( inpt_list% amino(j)% N_of_atoms )
                             case ( 14 )  !<== PGLU
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "PGLU"

                             case ( 16 )  !<== GLUH
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "GLUH"
                        end select

                   case( "LYS" )
                        select case ( inpt_list% amino(j)% N_of_atoms )
                             case ( 22 )  !<== LYSH
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "LYSH"
                        end select

                   case( "ASP" )
                        select case ( inpt_list% amino(j)% N_of_atoms )
                             case ( 13 )  !<== ASPH
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "ASPH"
                        end select

                   case( "ARG" )
                        select case ( inpt_list% amino(j)% N_of_atoms )
                             case ( 23 )  !<== ARGN
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "ARGN"
                        end select

                   case( "CYS2" )
                        select case ( inpt_list% amino(j)% N_of_atoms )
                             case ( 10 )  !<== CYS2
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "CYS2"

                             case ( 11 )  !<== CYSH
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "CYSH"
                        end select

                   case( "HISH" )
                        select case ( inpt_list% amino(j)% N_of_atoms )
                             case ( 17 )  !<== HISD , HIS1 , HISE
                             sys% atom(indx:indx+inpt_list% amino(j)% N_of_atoms-1)% resid = "HISE"
                        end select

              end select
         end if

    end do
end do
!----------------------------------------------------------------------------

end subroutine Fix_residues
!
!
!
!=====================================
 subroutine Dump_pdb( sys , NewOrder )
!=====================================
implicit none 
type(universe) , intent(inout) ::  sys
integer        , intent(in)    :: NewOrder(:)

! local variables ...
integer ::  i , j , k

!----------------------------------------------
!     generate pdb file for GROMACS
!----------------------------------------------

OPEN(unit=4,file='seed.pdb',status='unknown')
write(4,6) sys%Surface_Characteristics

write(4,1) 'CRYST1' , sys%box(1) , sys%box(2) , sys%box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'

do j = 1 , sys%N_of_atoms

    If( NewOrder(j) == 0 ) cycle 

    i = NewOrder(j) 

    write(4,2)  'HETATM'                        ,  &    ! <== non-standard atom
                j                               ,  &    ! <== global number
                sys%atom(i)%MMSymbol            ,  &    ! <== atom type
                ' '                             ,  &    ! <== alternate location indicator
                sys%atom(i)%resid               ,  &    ! <== residue name
                ' '                             ,  &    ! <== chain identifier
                sys%atom(i)%nresid              ,  &    ! <== residue sequence number
                ' '                             ,  &    ! <== code for insertion of residues
                ( sys%atom(i)%xyz(k) , k=1,3 )  ,  &    ! <== xyz coordinates 
                1.00                            ,  &    ! <== occupancy
                0.00                            ,  &    ! <== temperature factor
                ' '                             ,  &    ! <== segment identifier
                ' '                             ,  &    ! <== here only for tabulation purposes
                sys%atom(i)%symbol              ,  &    ! <== chemical element symbol
                sys%atom(i)%charge                      ! <== charge on the atom
end do

! check and print topological connections ...
!If( allocated(sys%topol) ) CALL dump_topol(sys,4)

write(4,3) 'MASTER', 0 , 0 , 0 ,  0 , 0 , 0 , 0 , 0 , count(NewOrder/=0) , 0 , count(NewOrder/=0) , 0
write(4,*) 'END'

close(4)

1 FORMAT(a6,3F9.3,3F7.2,a11,a4)
2 FORMAT(a6,i5,a4,a1,a4,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4)
3 FORMAT(a6,i9,11i5)
6 FORMAT(a6,a72)

end subroutine Dump_pdb
!
!
!
!
!==================================
 subroutine QM_2_MM_distance( sys )
!==================================
implicit none
type(universe) , intent(inout)   :: sys 

!local variables ...
integer                       :: i , j , k , nres_max , ref_res , choice
integer         , allocatable :: aux(:)
real*8                        :: boundary , ref_CofM(3)
real*8          , allocatable :: res_mass(:) , res_CofM_dist(:)
type(R3_vector) , allocatable :: res_CofM(:)
character(4)                  :: residue
character(4)    , allocatable :: tmp(:) , res_name(:)

!-----------------------------------------------------------
! identify the distinct residue numbers ...
allocate( aux(sys%N_of_atoms) , source = 0 )
allocate( tmp(sys%N_of_atoms) )

nres_max =  maxval(sys%atom%nresid)
k = 1
aux(1) = sys% atom(1)% nresid
tmp(1) = sys% atom(1)% resid
do i = 2 , sys% N_of_atoms
    if( .not. any(aux == sys% atom(i)% nresid) ) then
         k = k + 1
         aux(k) = sys% atom(i)% nresid
         tmp(k) = sys% atom(i)% resid
    end if
end do
write(*,*)' Different residue numbers:',k      
write(*,*)' Biggest nresid:', nres_max
!-----------------------------------------------------------

!-----------------------------------------------------------
! find distance between Center of Mass of residues ...

! get the Atomic_Masses ...
sys% atom% mass = Atomic_Mass(sys% atom% AtNo)

allocate( res_CofM(k) , res_CofM_dist(k) , res_mass(k) ) 
allocate( res_name(k) , source = tmp(1:k) )

do i = 1 , k
    do j = 1 , 3
         res_CofM(i)% xyz(j) = sum( sys% atom(:)% xyz(j)*sys% atom(:)% mass , sys% atom% nresid == i )
    end do
    res_mass(i) = sum( sys% atom(:)% mass , sys% atom% nresid == i )
    res_CofM(i)% xyz = res_CofM(i)% xyz / res_mass(i) 
end do

write(*,'(/a)') ' Choose residue NAME (1) or residue NUMBER (2) to fix the origin'
write(*,'(/a)',advance='no') '>>>   '
read (*,'(I)') choice 

select case( choice )

    case( 2 )

         write(*,'(/a)') ' Choose reference residue NUMBER to place the origin of coords'
         write(*,'(/a)',advance='no') '>>>   '
         read (*,'(I)') ref_res 

         ref_CofM = res_CofM(ref_res)% xyz
         do i = 1 , k
             res_CofM(i)% xyz = res_CofM(i)% xyz - ref_CofM
             res_CofM_dist(i) = dsqrt( sum(res_CofM(i)% xyz * res_CofM(i)% xyz) ) 
         end do

    case( 1 )

         write(*,'(/a)') ' Choose reference residue NAME to place the origin of coords'
         write(*,'(/a)',advance='no') '>>>   '
         read (*,*) residue

         forall( j=1:3) ref_CofM(j) = sum( res_CofM(:)% xyz(j) * res_mass(:) , res_name == residue ) / sum( res_mass(:) , res_name == residue )

         do i = 1 , k
             res_CofM(i)% xyz = res_CofM(i)% xyz - ref_CofM
             res_CofM_dist(i) = dsqrt( sum(res_CofM(i)% xyz * res_CofM(i)% xyz) ) 
         end do

end select
!-----------------------------------------------------------
write(*,'(/a)') ' Choose cutoff radius for THERMAL/FROZEN (THR/FRZ) domains'
write(*,'(/a)',advance='no') '>>>   '
read (*,*) boundary

If(.not. any( sys% atom% resid == "AMI" )) stop "THR/FRZ operation acts on AMI residues; no AMI residues found; execution terminated"

where( (res_name == "AMI") .AND. (res_CofM_dist > boundary)  ) res_name = "FRZ"
where( (res_name == "AMI") .AND. (res_CofM_dist <= boundary) ) res_name = "THR"

forall( i=1:sys% N_of_atoms ) sys% atom(i)% resid = res_name( sys% atom(i)% nresid )

end subroutine QM_2_MM_distance
!
!
!
end module Aminoacids

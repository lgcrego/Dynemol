module Read_Parms

use types_m

    type(EHT)                   , public    , protected :: atom(300)
    real*8      , allocatable   , public    , protected :: Atomic_Mass(:)

    public :: Symbol_2_AtNo , AtNo_2_Symbol , MMSymbol_2_Symbol , Pack_Residues , Identify_Residues

!    private

contains
!    
!
!
!==========================
 subroutine read_EHT_params
!==========================
 implicit none

! local variables .... 
 integer          :: ioerr , i , AtNo , Ang , DOS_sum
 character(len=1) :: spdf

 OPEN(unit=3,file='my_eht_parameters.dat',status='old')

 AtNo = 1
 Ang  = 0
 DOS_sum = 0
 do i = 1 , size(atom)

    read(3,*,IOSTAT=ioerr)   &
    atom(AtNo)%symbol      , &
    atom(AtNo)%AtNo        , &
    atom(AtNo)%Nvalen      , &
    atom(AtNo)%Nzeta(Ang)  , &
    atom(AtNo)%Nquant(Ang) , &
    spdf                   , &
    atom(AtNo)%IP(Ang)     , &
    atom(AtNo)%zeta(Ang,1) , &
    atom(AtNo)%zeta(Ang,2) , &
    atom(AtNo)%coef(Ang,1) , &
    atom(AtNo)%coef(Ang,2) , &
    atom(AtNo)%DOS

    if(ioerr < 0) EXIT

    select case (spdf)
        case('s')
            DOS_sum = DOS_sum + 1
        case('p')
            DOS_sum = DOS_sum + 3
        case('d') 
            DOS_sum = DOS_sum + 5
        case('f')
            DOS_sum = DOS_sum + 7
    end select

    if( DOS_sum == atom(AtNo)%DOS) then
        atom(AtNo)%AngMax = Ang
        Ang               = 0
        DOS_sum           = 0
        AtNo              = AtNo + 1
    else
        Ang = Ang + 1
    end if

 end do    

 end subroutine read_EHT_params
! 
!
!==========================
subroutine Read_Atomic_Mass
!==========================
implicit none

! local variables ...
real*8  , allocatable   :: temp(:)
integer                 :: ioerr , i , size_of_array

allocate( temp(300) )

OPEN(unit=3,file='atomic_mass.dat',status='old')
do 
    read(3,*,IOSTAT=ioerr) i , temp(i)  
    if(ioerr < 0) EXIT
end do    
CLOSE(3)

size_of_array = i

allocate( Atomic_Mass(size_of_array) )

Atomic_Mass = temp(1:size_of_array)

deallocate( temp )

end subroutine Read_Atomic_Mass
! 
!
!
!===========================
 subroutine Symbol_2_AtNo(a)
!===========================
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(adjustl(a(i)%symbol))
        case( 'H') 
            a(i)%AtNo = 1 
        case( 'LI','Li') 
            a(i)%AtNo = 3 
        case( 'BE','Be') 
            a(i)%AtNo = 4 
        case( 'B' ) 
            a(i)%AtNo = 5 
        case( 'C') 
            a(i)%AtNo = 6 
        case( 'N') 
            a(i)%AtNo = 7 
        case( 'O') 
            a(i)%AtNo = 8 
        case( 'F') 
            a(i)%AtNo = 9 
        case( 'AL','Al') 
            a(i)%AtNo = 13 
        case( 'SI','Si') 
            a(i)%AtNo = 14 
        case( 'P','p') 
            a(i)%AtNo = 15 
        case( 'S','s') 
            a(i)%AtNo = 16 
        case( 'CL','Cl') 
            a(i)%AtNo = 17 
        case( 'TI','Ti') 
            a(i)%AtNo = 22 
        case( 'CR','Cr') 
            a(i)%AtNo = 24 
        case( 'MN','Mn') 
            a(i)%AtNo = 25 
        case( 'FE','Fe') 
            a(i)%AtNo = 26 
        case( 'ZN','Zn') 
            a(i)%AtNo = 30 
        case( 'RU','Ru') 
            a(i)%AtNo = 44 
        case( 'I') 
            a(i)%AtNo = 53 
        case( 'Cs','CS' ) 
            a(i)%AtNo = 55 
        case( 'Pb','PB') 
            a(i)%AtNo = 82 
        case default
            print*, ' >> Symbol_2_AtNo ; unknown atom found <<' , '[',a(i)%symbol,']' , i
            stop
    end select

 END DO

 end subroutine Symbol_2_AtNo
!
!
!
!===========================
 subroutine AtNo_2_Symbol(a)
!===========================
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer :: i

 DO i = 1 , size(a)

    select case(a(i)%AtNo)
        case( 1) 
            a(i)%Symbol = 'H'
        case( 3) 
            a(i)%Symbol = 'Li'
        case( 4) 
            a(i)%Symbol = 'Be'
        case( 5) 
            a(i)%Symbol = 'B'
        case( 6) 
            a(i)%Symbol = 'C'
        case( 7) 
            a(i)%Symbol = 'N'
        case( 8) 
            a(i)%Symbol = 'O'
        case( 9) 
            a(i)%Symbol = 'F'
        case( 13 ) 
            a(i)%Symbol = 'Al'
        case( 14 ) 
            a(i)%Symbol = 'Si'
        case( 15 ) 
            a(i)%Symbol = 'P' 
        case( 16 ) 
            a(i)%Symbol = 'S '
        case( 17 ) 
            a(i)%Symbol = 'Cl '
        case( 22 ) 
            a(i)%Symbol = 'Ti '
        case( 24 ) 
            a(i)%Symbol = 'Cr '
        case( 25 ) 
            a(i)%Symbol = 'Mn'
        case( 26 ) 
            a(i)%Symbol = 'Fe'
        case( 30 ) 
            a(i)%Symbol = 'Zn'
        case( 44 ) 
            a(i)%Symbol = 'Ru'
        case( 53 ) 
            a(i)%Symbol = 'I'
        case( 55 ) 
            a(i)%Symbol = 'Cs'
        case( 82 ) 
            a(i)%Symbol = 'Pb'
        case( 101:109 )
            a(i)%Symbol = '@@'
        case default
            print*, ' >> AtNo_2_Symbol ; unknown atom found <<' , a(i)%AtNo
            stop
    end select

 END DO

 end subroutine AtNo_2_Symbol
!
!
!
!================================
subroutine Identify_Residues( a )
!================================
implicit none
type(universe)  , intent(inout) :: a

! local variables ...
integer                         :: i , j , counter
character(3)    , allocatable   :: temp(:)
logical                         :: flag

allocate( temp(a%N_of_Atoms) )

temp(1) = a % atom(1) % resid
counter = 1

do i = 1 , a%N_of_Atoms

    flag = .true.
    do j = 1 , counter
        flag = flag .AND. (temp(j) /= a%atom(i)%resid)
    end do

    if( flag ) then
        counter = counter + 1
        temp(counter) = a%atom(i)%resid
    end if

end do

! build list of residues in a ...
allocate( a%list_of_residues(counter) )
a%list_of_residues = temp(1:counter)
deallocate( temp )

end subroutine Identify_Residues
!
!
!
!================================================
 subroutine Pack_Residues( a , list_of_residues )
!================================================
implicit none
type(atomic)  , allocatable  ,  intent(inout) :: a(:)
character(*)                 , intent(in)     :: list_of_residues(:)

! local variables ...
integer                     :: i , N , first , last
type(atomic) , allocatable  :: temp(:) , buffer(:)

allocate( temp(size(a)) )

first = 1
last  = 0

do i = 1 , size(list_of_residues)

    N = count( a%resid == list_of_residues(i) )
    allocate( buffer(N) )

    buffer = pack( a , a%resid == list_of_residues(i) , buffer )

    last = last + N

    temp(first:last) = buffer

    first = last + 1

    deallocate(buffer)

end do

CALL move_alloc( from=temp , to=a )

end subroutine Pack_Residues
!
!
!
!===============================
 subroutine MMSymbol_2_Symbol(a)
!===============================
type(atomic) , intent(inout) :: a(:)

! local variables ...
integer             :: i
character(len=1)    :: element

 DO i = 1 , size(a)

    write( element,'(A1)' ) adjustl( a(i)%MMSymbol )

    select case( element )
        case( 'I' ) 
            a(i)%Symbol = 'I' 
        case( 'S' ) 
            a(i)%Symbol = 'S' 
        case( 'C' ) 
            a(i)%Symbol = 'C' 
        case( 'N' ) 
            a(i)%Symbol = 'N' 
        case( 'O' ) 
            a(i)%Symbol = 'O' 
        case( 'H' ) 
            a(i)%Symbol = 'H' 
        case( 'P' ) 
            a(i)%Symbol = 'P' 
    end select

    select case( adjustl( a(i)%MMSymbol) )
        case( 'Ox' , 'OxH' ) 
            a(i)%Symbol = 'O' 
        case( 'YN' , 'NTr' , 'Nx' ) 
            a(i)%Symbol = 'N' 
        case( 'Al' ) 
            a(i)%Symbol = 'Al' 
        case( 'Ti' , 'TI' ) 
            a(i)%Symbol = 'Ti' 
        case( 'Si' , 'SI' ) 
            a(i)%Symbol = 'Si' 
        case( 'Li' ) 
            a(i)%Symbol = 'Li' 
        case( 'Zn' , 'ZN' )
            a(i)%Symbol = 'Zn' 
        case( 'Ru' ) 
            a(i)%Symbol = 'Ru' 
        case( 'Pb' ) 
            a(i)%Symbol = 'Pb' 
        case( 'HC' ) 
            a(i)%Symbol = 'H'  
        case( 'C=' , 'CTr' , 'CS' , 'CC' , 'CM' , 'YC' ) 
            a(i)%Symbol = 'C'  
        case( 'SS' ) 
            a(i)%Symbol = 'S'  
    end select

 END DO
 
 end subroutine MMSymbol_2_Symbol
!
!
!
end module Read_Parms
    

module Semi_Empirical_Parms

    use type_m
    use parameters_m    , only  : OPT_parms

    type(EHT) , public , protected :: atom(300) 

    public :: read_EHT_parameters , Include_OPT_parameters


    private

    type(EHT)   , allocatable   , save :: EH_atom(:)

    interface Include_OPT_parameters 
        module procedure Basis_OPT_parameters
        module procedure Atomic_OPT_parameters
    end interface

    contains
!
!
! 
!==============================
 subroutine read_EHT_parameters
!==============================
 implicit none

! local variables ... 
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
    atom(AtNo)%zeta(Ang,1) , &      ! <== zetas of my_eht_parameters.dat are given in units of a0^{-1} ...
    atom(AtNo)%zeta(Ang,2) , &      ! <== zetas of my_eht_parameters.dat are given in units of a0^{-1} ...
    atom(AtNo)%coef(Ang,1) , &
    atom(AtNo)%coef(Ang,2) , &
    atom(AtNo)%DOS         , &
    atom(AtNo)%polar                ! <== 10^{-24}*cm^3  
    
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

 close(3)

 ! transform zetas to units of Angs^{-1} ...
 forall( Ang=0:3 , i=1:2 ) atom(:)%zeta(Ang,i) =  atom(:)%zeta(Ang,i) / a_Bohr 

 If( OPT_parms ) CALL read_OPT_parameters

 end subroutine read_EHT_parameters
! 
!
! 
!==============================
 subroutine read_OPT_parameters
!==============================
 implicit none

! local variables ... 
 integer        :: ioerr , nr , i 
 character(4)   :: spdf
 character(2)   :: dumb 
 character(8)   :: Symbol_char
 character(12)  :: EHSymbol_char

OPEN(unit=3,file='OPT_eht_parameters.input.dat',status='old')

! read file heading ...
read(3,*,IOSTAT=ioerr) dumb 

! number of lines of OPT_parameters ...
nr = 0
do 
    read(3,*,IOSTAT=ioerr) dumb 
    if(ioerr < 0) EXIT
    nr = nr + 1
end do    

allocate( EH_atom(nr) )

! rewind and read the OPT-EHT_parameters ...
rewind 3    ;    read(3,*,IOSTAT=ioerr) dumb 

do i = 1 , size(EH_atom)

    read(3,17)  Symbol_char             ,   &
                EHSymbol_char           ,   &
                EH_atom(i)%AtNo         ,   &
                EH_atom(i)%Nvalen       ,   &
                EH_atom(i)%Nzeta(0)     ,   &
                EH_atom(i)%n            ,   &
                spdf                    ,   &
                EH_atom(i)%IP(0)        ,   &
                EH_atom(i)%zeta(0,1)    ,   &      ! <== zetas of OPT_eht_parameters.input.dat are given in units of a0^{-1} ...
                EH_atom(i)%zeta(0,2)    ,   &      ! <== zetas of OPT_eht_parameters.input.dat are given in units of a0^{-1} ...
                EH_atom(i)%coef(0,1)    ,   &
                EH_atom(i)%coef(0,2)    ,   &
                EH_atom(i)%k_WH(0)

    select case ( adjustl(spdf) )
        case('s')
            EH_atom(i)%l = 0
        case('p')
            EH_atom(i)%l = 1
        case('d') 
            EH_atom(i)%l = 2
        case('f')
            EH_atom(i)%l = 3
    end select

    EH_atom(i) % Symbol   =  adjustl( Symbol_char   )
    EH_atom(i) % EHSymbol =  adjustl( EHSymbol_char )
    EH_atom(i) % DOS      =  atom( EH_atom(i)%AtNo ) % DOS
    EH_atom(i) % AngMax   =  atom( EH_atom(i)%AtNo ) % AngMax

end do

close(3)

Print 44
Print 45 , EH_atom%EHSymbol

 ! transform zetas to units of Angs^{-1} ...
 forall( i=1:2 ) EH_atom(:)%zeta(0,i) =  EH_atom(:)%zeta(0,i) / a_Bohr 

include 'formats.h'

17 format(A8,t10,A12,t24,I7,t32,I10,t43,I9,t53,I5,t61,A4,t68,F9.6,t78,F9.6,t88,F9.6,t98,F9.6,t108,F9.6,t118,F9.6)

end subroutine read_OPT_parameters
! 
!
!
!==========================================
 subroutine Basis_OPT_parameters( basis )
!==========================================
implicit none
type(STO_basis) , intent(inout) :: basis(:)

! local variables ...
integer :: i

do i = 1 , size(EH_atom)

    where( (adjustl(basis%EHSymbol) == EH_atom(i)%EHSymbol) .AND. (basis%l == EH_atom(i)%l) ) 
        
        basis%IP        =  EH_atom(i)%IP    (0)
        basis%Nzeta     =  EH_atom(i)%Nzeta (0)
        basis%coef(1)   =  EH_atom(i)%coef  (0,1)
        basis%coef(2)   =  EH_atom(i)%coef  (0,2)
        basis%zeta(1)   =  EH_atom(i)%zeta  (0,1)
        basis%zeta(2)   =  EH_atom(i)%zeta  (0,2)
        basis%k_WH      =  EH_atom(i)%k_WH  (0)

    end where

end do

end subroutine Basis_OPT_parameters
! 
!
!
!==============================================
 subroutine Atomic_OPT_parameters( system )
!==============================================
implicit none
type(structure) , intent(inout) :: system

! local variables ...
integer :: i

do i = 1 , size(EH_atom)

    where( (adjustl(system%MMSymbol) == EH_atom(i)%EHSymbol) ) 
        
        system%Nvalen = EH_atom(i)%Nvalen

    end where

end do

end subroutine Atomic_OPT_parameters
!
!
end module Semi_Empirical_Parms

module Semi_Empirical_Parms

    use type_m

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

    type(EHT)                   , public    , protected :: atom(300)
    real*8      , allocatable   , public    , protected :: Atomic_Mass(:)

    contains
!
!
! 
!=======================================================================
 subroutine Define_EH_Parametrization( Unit_Cell , Characteristics ) 
!=======================================================================
implicit none
type(structure) , intent(inout) :: Unit_Cell
character(*)    , intent(inout) :: Characteristics

! defining the k_WH parameter for the system ...

! TiO2 ...
where( Unit_Cell % residue == "CCC" ) Unit_Cell % k_WH = 1.75d0

! PYR ...
where( Unit_Cell % residue == "PYR" ) unit_cell % k_WH = 1.75d0

! Alq3 ...
!where( Unit_Cell % residue == "ALQ" ) Unit_Cell % k_WH = 2.0d0

! to be compared with structure information ...
Characteristics = "TiO2-pyridine-small"

end subroutine Define_EH_Parametrization 
!    
!
!
!==============================
 subroutine read_EHT_parameters
!==============================

 implicit none

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

 end subroutine read_EHT_parameters
! 
!
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
end module Semi_Empirical_Parms

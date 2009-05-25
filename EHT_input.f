module EHT_parameters

    use type_m

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

    type(EHT) , public :: atom(300)

    contains
!    
!
!
 subroutine read_EHT_parameters

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
 end module EHT_parameters
    

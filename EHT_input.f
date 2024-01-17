module Semi_Empirical_Parms

    use type_m
    use parameters_m    , only  : OPT_parms

    type(EHT) , public , protected :: atom(300) 
    type(EHT) , allocatable , save :: EH_atom(:)

    public :: read_EHT_parameters , Include_OPT_parameters, EH_atom

    private

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

 OPEN(unit=3,file=dynemoldir//'my_eht_parameters.dat',status='old')

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
 forall( Ang=0:3 , i=1:2 ) atom(:)%zeta(Ang,i) = atom(:)%zeta(Ang,i) / a_Bohr 

 ! truncate zeta parameters to 1.d-4 ...
 do concurrent ( Ang=0:3 , i=1:2 )
      atom(:)%zeta(Ang,i) = atom(:)%zeta(Ang,i) * 1.d4
      atom(:)%zeta(Ang,i) = int(atom(:)%zeta(Ang,i))
      atom(:)%zeta(Ang,i) = atom(:)%zeta(Ang,i) * 1.d-4

      atom(:)%coef(Ang,i) = atom(:)%coef(Ang,i) * 1.d4
      atom(:)%coef(Ang,i) = int(atom(:)%coef(Ang,i))
      atom(:)%coef(Ang,i) = atom(:)%coef(Ang,i) * 1.d-4
end do

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
 character(6)   :: spdf
 character(2)   :: dumb 
 character(8)   :: Symbol_char
 character(11)  :: residue_char
 character(12)  :: EHSymbol_char 
 logical        :: flag1 , flag2 , flag3 , flag4

OPEN(unit=3,file=dynemolworkdir//'opt_eht_parms.input',status='old')

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
                residue_char            ,   &
                EH_atom(i)%AtNo         ,   &
                EH_atom(i)%Nvalen       ,   &
                EH_atom(i)%Nzeta(0)     ,   &
                EH_atom(i)%n            ,   &
                spdf                    ,   &
                EH_atom(i)%IP(0)        ,   &
                EH_atom(i)%zeta(0,1)    ,   &      ! <== zetas of opt_eht_parms.input are given in units of Ang^{-1} ...
                EH_atom(i)%zeta(0,2)    ,   &      ! <== zetas of opt_eht_parms.input are given in units of Ang^{-1} ...
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
    EH_atom(i) % residue  =  adjustl( residue_char  )
    EH_atom(i) % DOS      =  atom( EH_atom(i)%AtNo ) % DOS
    EH_atom(i) % AngMax   =  atom( EH_atom(i)%AtNo ) % AngMax

    ! input-error checking ...
    ! "Happy families are all alike; every unhappy family is unhappy in its own way"
    flag1 = ( (EH_atom(i)% zeta(0,2) == 0.0) .AND. (EH_atom(i)% Nzeta(0) == 2) )
    flag2 = ( (EH_atom(i)% zeta(0,2) /= 0.0) .AND. (EH_atom(i)% Nzeta(0) == 1) )

    flag3 = ( (EH_atom(i)% coef(0,2) == 0.0) .AND. (EH_atom(i)% Nzeta(0) == 2) )
    flag4 = ( (EH_atom(i)% coef(0,2) /= 0.0) .AND. (EH_atom(i)% Nzeta(0) == 1) )

    If( flag1 .OR. flag2 .OR. flag3 .OR. flag4 ) then
        CALL warning("error in opt_eht_parms.input ; check Nzeta parameter")
        STOP 
    end If
 
end do

close(3)

Print 44
Print 45 , ( EH_atom(i)% EHSymbol , EH_atom(i)% residue , i = 1,size(EH_atom) )

forall( i=1:2 ) EH_atom(:)%zeta(0,i) =  EH_atom(:)%zeta(0,i) / a_Bohr 

! truncate zeta parameters to 1.d-5 ...
do concurrent ( i=1:2 )
     EH_atom(:)%zeta(0,i) = EH_atom(:)%zeta(0,i) * 1.d5
     EH_atom(:)%zeta(0,i) = anint(EH_atom(:)%zeta(0,i))
     EH_atom(:)%zeta(0,i) = EH_atom(:)%zeta(0,i) * 1.d-5

     EH_atom(:)%coef(0,i) = EH_atom(:)%coef(0,i) * 1.d5
     EH_atom(:)%coef(0,i) = anint(EH_atom(:)%coef(0,i))
     EH_atom(:)%coef(0,i) = EH_atom(:)%coef(0,i) * 1.d-5
end do

include 'formats.h'

17 format(A8,t10,A12,t23,A11,t35,I7,t44,I10,t55,I9,t65,I5,t71,A6,t80,F9.5,t90,F9.6,t100,F9.6,t110,F9.6,t120,F9.6,t130,F9.6)

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

    select case( EH_atom(i)%residue(3:3) )

           case( "*" )
                where( (adjustl(basis%EHSymbol) == EH_atom(i)%EHSymbol) &
                .AND.  (adjustl(basis%residue(1:2) ) == EH_atom(i)%residue(1:2))  &
                .AND. (basis%l == EH_atom(i)%l) ) 
                    
                    basis%IP        =  EH_atom(i)%IP    (0)
                    basis%Nzeta     =  EH_atom(i)%Nzeta (0)
                    basis%coef(1)   =  EH_atom(i)%coef  (0,1)
                    basis%coef(2)   =  EH_atom(i)%coef  (0,2)
                    basis%zeta(1)   =  EH_atom(i)%zeta  (0,1)
                    basis%zeta(2)   =  EH_atom(i)%zeta  (0,2)
                    basis%k_WH      =  EH_atom(i)%k_WH  (0)
                    basis%modified  =  .true.
               
                end where

           case default 
                where( (adjustl(basis%EHSymbol) == EH_atom(i)%EHSymbol) &
                .AND.  (adjustl(basis%residue ) == EH_atom(i)%residue)  &
                .AND. (basis%l == EH_atom(i)%l) ) 
                    
                    basis%IP        =  EH_atom(i)%IP    (0)
                    basis%Nzeta     =  EH_atom(i)%Nzeta (0)
                    basis%coef(1)   =  EH_atom(i)%coef  (0,1)
                    basis%coef(2)   =  EH_atom(i)%coef  (0,2)
                    basis%zeta(1)   =  EH_atom(i)%zeta  (0,1)
                    basis%zeta(2)   =  EH_atom(i)%zeta  (0,2)
                    basis%k_WH      =  EH_atom(i)%k_WH  (0)
                    basis%modified  =  .true.
               
                end where

    end select      

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

    select case( EH_atom(i)%residue(3:3) )

           case( "*" )
                where( (adjustl(system%MMSymbol) == EH_atom(i)%EHSymbol)           &
                .AND.  (adjustl(system%residue(1:2)) == EH_atom(i)%residue(1:2))) 
                    system%Nvalen = EH_atom(i)%Nvalen
                end where

           case default 
                where( (adjustl(system%MMSymbol) == EH_atom(i)%EHSymbol)           &
                .AND.  (adjustl(system%residue) == EH_atom(i)%residue)) 
                    system%Nvalen = EH_atom(i)%Nvalen
                end where

    end select

end do

end subroutine Atomic_OPT_parameters
!
!
end module Semi_Empirical_Parms

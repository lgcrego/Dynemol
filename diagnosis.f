module diagnostic_types_m

    use MM_types

    public :: diagnostic_types

    private

    interface diagnostic_types
        module procedure diagnostic_MM_atomic 
        module procedure diagnostic_MM_molecular
    end interface diagnostic_types

contains
!
!
!
!=======================================
 subroutine diagnostic_MM_molecular( a )
!=======================================
implicit none
type(MM_molecular)  , intent(in)    :: a

! local variables ...
integer :: i , option

do

    write(*,*) ' (0) QUIT           '
    print'("")'      
    write(*,*) ' (1)  atom          '
    write(*,*) ' (2)  bonds         '
    write(*,*) ' (3)  kbond0        '
    write(*,*) ' (4)  angs          '
    write(*,*) ' (5)  kang0         '
    write(*,*) ' (6)  diheds        '
    write(*,*) ' (7)  kdihed0       '
    write(*,*) ' (8)  bonds14       '
    write(*,*) ' (9)  harm          '
    write(*,*) ' (10) fact14        '
    write(*,*) ' (11) funct_bond    '
    write(*,*) ' (12) funct_angle   '
    write(*,*) ' (13) funct_dihed   '
    write(*,*) ' (14) dihedral_type '

    read (*,'(a)') option

    select case( option )

        case(0)
            stop

        case(1)
            CALL diagnostic_MM_atomic( a % atom )

        case(2)
            do i = 1 , size(a%bonds(:,1))
                write(*,*) a % bonds(i,:)
            end do

        case(3)
            do i = 1 , size(a%kbond0(:,1))
                write(*,*) a % kbond0(i,:)
            end do

        case(4)
            do i = 1 , size(a%angs(:,1))
                write(*,*) a % angs(i,:)
            end do

        case(5)
            do i = 1 , size(a%kang0(:,1))
                write(*,*) a % kang0(i,:)
            end do

        case(6)
            do i = 1 , size(a%diheds(:,1))
                write(*,*) a % diheds(i,:)
            end do

        case(7)
            do i = 1 , size(a%kdihed0(:,1))
                write(*,*) a % kdihed0(i,:)
            end do

        case(8)
            do i = 1 , size(a%bonds14(:,1))
                write(*,*) a % bonds14(i,:)
            end do

        case(9)
            write(*,*) a % harm(:)

        case(10)
            write(*,*) a % fact14(:)

        case(11)
            write(*,*) a % funct_bond(:)

        case(12)
            write(*,*) a % funct_angle(:)

        case(13)
            write(*,*) a % funct_dihed(:)

        case(14)
            write(*,*) a % dihedral_type(:)

        case default
            exit

    end select

end do

end subroutine diagnostic_MM_molecular
!
!
!
!====================================
 subroutine diagnostic_MM_atomic( a )
!====================================
implicit none
type(MM_atomic) , intent(in)    :: a(:)

! local variables ...
integer :: i , option

do

    write(*,*) ' (0) QUIT        '
    print'("")'      
    write(*,*) ' (1) AtNo       '
    write(*,*) ' (2) my_id      '
    write(*,*) ' (3) my_species '
    write(*,*) ' (4) nr         '
    write(*,*) ' (5) residue    '
    write(*,*) ' (6) Symbol     '
    write(*,*) ' (7) MMSymbol   '
    write(*,*) ' (8) mass       '
    write(*,*) ' (9) charge     '
    write(*,*) ' (10) MM_charge '
    write(*,*) ' (11) eps       '
    write(*,*) ' (12) sig       '
    write(*,*) ' (13) free      '

    read (*,*) option

    select case( option )

        case(0)
            stop

        case(1)
            write(*,*) a(:) % AtNo

        case(2)
            write(*,*) a(:) % my_id

        case(3)
            write(*,*) a(:) % my_species

        case(4)
            write(*,*) a(:) % nr

        case(5)
            write(*,*) a(:) % residue

        case(6)
            write(*,*) a(:) % Symbol

        case(7)
            write(*,*) a(:) % MMSymbol

        case(8)
            write(*,*) a(:) % mass

        case(9)
            write(*,*) a(:) % charge

        case(10)
            write(*,*) a(:) % MM_charge

        case(11)
            write(*,*) a(:) % eps

        case(12)
            write(*,*) a(:) % sig
            
        case(13)
            write(*,*) a(:) % free

        case default
            exit

    end select

end do

end subroutine diagnostic_MM_atomic

end module diagnostic_types_m

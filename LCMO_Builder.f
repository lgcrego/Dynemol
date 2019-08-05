 module LCMO_m

    use type_m
    use constants_m
    use parameters_m                , only : electron_state , hole_state
                                             
    public :: LCMO_Builder

    private

    ! module variables ...
    integer :: LCMO_nstates
    integer , allocatable   :: LCMO_states_el(:) , LCMO_states_hl(:)
    real*8  , allocatable   :: LCMO_coeff(:) 
    logical                 :: done = .false.
    
 contains


!======================================================
 subroutine LCMO_Builder( wv_FMO , FMO_erg , instance )
!======================================================
 implicit none
 real*8          , allocatable  , intent(inout) :: wv_FMO(:,:) 
 real*8                         , intent(inout) :: FMO_erg(:)
 character(*)    , optional     , intent(in)    :: instance

! local variables ...
 integer                 :: i , rows_wv_FMO , cols_wv_FMO
 real*8                  :: sumsq 
 real*8, allocatable     :: coeff_el(:) , coeff_hl(:) , casida_el(:) , casida_hl(:)

 If ( .NOT. done ) CALL preprocess_Casida

! START Construction of initial el/hl state as LCMO ...

rows_wv_FMO = size( wv_FMO(:,1) )
cols_wv_FMO = size( wv_FMO(1,:) )

! normalize LCMO_el coefficients to 1, for el and hl alike ...
sumsq         = sqrt( sum( sign(D_one,LCMO_coeff) * LCMO_coeff*LCMO_coeff ) )
LCMO_coeff(:) = LCMO_coeff(:) / sumsq

select case (instance)

    case( "E" )

        allocate( coeff_el (rows_wv_FMO) , source = D_zero )

        do i = 1 , rows_wv_FMO
            coeff_el(i) = sqrt( sum( sign(D_one,LCMO_coeff) * LCMO_coeff*LCMO_coeff , LCMO_states_el == i ) )
        end do

        allocate( casida_el (cols_wv_FMO) , source = D_zero)

        do i = 1 , cols_wv_FMO
            casida_el(i) = sum( coeff_el(:) * wv_FMO(:,i) )
        end do

        wv_FMO(electron_state,:) = casida_el(:) 

        FMO_erg(electron_state)  = sum( coeff_el*coeff_el * FMO_erg )

        deallocate( coeff_el , casida_el , LCMO_states_el )

    case( "H" )        

        allocate( coeff_hl (rows_wv_FMO) , source = D_zero )

        do i = 1 , rows_wv_FMO
            coeff_hl(i) = sqrt( sum( sign(D_one,LCMO_coeff) * LCMO_coeff*LCMO_coeff , LCMO_states_hl == i ) )
        end do

        allocate( casida_hl (cols_wv_FMO) , source = D_zero)

        do i = 1 , cols_wv_FMO
            casida_hl(i) = sum( coeff_hl(:) * wv_FMO(:,i) )
        end do

        wv_FMO(hole_state,:) = casida_hl(:)

        FMO_erg(hole_state)  = sum( coeff_hl*coeff_hl * FMO_erg )

        deallocate( coeff_hl , casida_hl , LCMO_states_hl , LCMO_coeff )

end select

end subroutine LCMO_Builder
!
!
!
!=============================
 subroutine preprocess_Casida
!=============================
implicit none

! local variables ...
integer :: j , ioerr , nr 
character(1) :: dumb

OPEN(unit=3,file='casida-input.dat',status='old',iostat=ioerr,err=10)
nr = 0
do 
    read(3,*,IOSTAT=ioerr) dumb
    if(ioerr < 0) EXIT
    nr = nr + 1
end do    

LCMO_nstates = nr - 1

! allocatting module variables ... 
allocate( LCMO_states_el( LCMO_nstates ) )
allocate( LCMO_states_hl( LCMO_nstates ) )
allocate( LCMO_coeff    ( LCMO_nstates ) )

! read the casida-input file ...
rewind 3
read(3,*) dumb

Print 74 
Print 75
do j = 1 , LCMO_nstates

    read(3,*)  LCMO_states_hl(j) ,  LCMO_states_el(j) ,  LCMO_coeff(j)

    write(*,76) LCMO_states_hl(j) ,  LCMO_states_el(j) ,  LCMO_coeff(j) 

end do

CLOSE(3)

Print 43

10 if( ioerr > 0 ) stop "file casida-input.dat not found; terminating execution"

done = .true.

include 'formats.h'

end subroutine preprocess_Casida
!
!
!
end module LCMO_m

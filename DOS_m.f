module DOS_m

    use type_m
    use omp_lib    
    use constants_m
    use parameters_m            , only  : sigma , DOS_range , SOC       
    use Semi_Empirical_Parms    , only  : the_chemical_atom => atom

    public  :: Total_DOS , Partial_DOS

    private

    ! modulo variables ...
    real*8                :: gauss_norm , two_sigma2 , step
    integer , allocatable :: list_of_DOS_states(:)
    integer , allocatable :: atom(:) 
    integer               :: n_of_DOS_states 

 contains
!
!
!
!==========================================================
subroutine  Total_DOS( QM , basis , TDOS , internal_sigma )
!==========================================================
implicit none
type(R_eigen)             , intent(in)    :: QM
type(STO_basis)           , intent(in)    :: basis(:)
type(f_grid)              , intent(inout) :: TDOS
real*8        , OPTIONAL  , intent(in)    :: internal_sigma

! local variables ...
real*8  , allocatable :: erg_MO(:) , peaks(:) , func(:)
real*8                :: sgm , sub_occupation , SpinProject
integer               :: i1 , i2 , npoints, k , j , indx
integer               :: spin , S

if( present(internal_sigma) ) then 
    sgm = internal_sigma
else
    sgm = sigma
end if

npoints = size(TDOS%grid)

gauss_norm = 1.d0  !  1.d0 / (sgm*sqrt2PI)  <== for gauss_norm = 1 the gaussians are not normalized ...
two_sigma2 = 2.d0 * sgm*sgm
step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

! states in the range [DOS_range%inicio,DOS_range%fim]
i1 = maxloc(QM%erg , 1 , QM%erg <  DOS_range%inicio) + 1
i2 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1 + 1

allocate( erg_MO(n_of_DOS_states) )

! find the energies in the range [DOS_range%inicio,DOS_range%fim]
erg_MO = QM%erg( i1 : i2 )

allocate( peaks (npoints) )
allocate( func  (npoints) )

forall(k=1:npoints) TDOS%grid(k) = (k-1)*step + DOS_range%inicio

! the total density of states
TDOS%peaks2 = 0.d0
TDOS%func2  = 0.d0

spin = merge(1,0,SOC)

do S = spin , -spin , -2

    indx = i1-1

    do j = 1 , n_of_DOS_states
    
        indx = indx + 1

        peaks = 0.d0
        where( dabs(TDOS%grid-erg_MO(j)) < (step/two) ) peaks = D_one
    
        func = 0.d0
        where( ((TDOS%grid-erg_MO(j))**2/two_sigma2) < 25.d0 ) func = gauss_norm*exp( -(TDOS%grid-erg_MO(j))**2 / two_sigma2 )

        select case ( S )
           case(0)
               TDOS%peaks2(:,1) = TDOS%peaks2(:,1) + peaks
               TDOS%func2(:,1)  = TDOS%func2 (:,1) + func(:)

           case(1)  
               SpinProject = sum( QM%L(indx,:)*QM%R(:,indx) , basis(:)%s == S )
               TDOS%peaks2(:,1) = TDOS%peaks2(:,1) + SpinProject * peaks
               TDOS%func2(:,1)  = TDOS%func2(:,1)  + SpinProject * func(:)

           case(-1)  
               SpinProject = sum( QM%L(indx,:)*QM%R(:,indx) , basis(:)%s == S )
               TDOS%peaks2(:,2) = TDOS%peaks2(:,2) - SpinProject * peaks
               TDOS%func2(:,2)  = TDOS%func2(:,2)  - SpinProject * func(:)

           end select
    end do

end do

TDOS%average2 = TDOS%average2 + TDOS%func2

! occupation of TDOS ...
select case(spin)

       case(0)
           TDOS%occupation(1) = two * TDOS%peaks2(1,1)
           do k = 2 , npoints 
               TDOS%occupation(k) = TDOS%occupation(k-1) + two*TDOS%peaks2(k,1) 
           end do
           sub_occupation  = two * (i1 - 1) ! <== occupation below DOS window
           
           TDOS%occupation = sub_occupation + TDOS%occupation

       case(1)
           TDOS%occupation(1) = TDOS%peaks2(1,1) + abs(TDOS%peaks2(1,2))
           do k = 2 , npoints 
              TDOS%occupation(k) = TDOS%occupation(k-1) +     TDOS%peaks2(k,1)   ! <== spin up
              TDOS%occupation(k) = TDOS%occupation(k)   + abs(TDOS%peaks2(k,2))  ! <== spin down
           end do
           sub_occupation  = (i1 - 1) ! <== occupation below DOS window
           
           TDOS%occupation = sub_occupation + TDOS%occupation

       end select 

DEALLOCATE( peaks , func )

print*, '>> TDOS done <<'

end subroutine Total_DOS
!
!
!
!==========================================================================
subroutine  Partial_DOS( system , basis , QM , PDOS , nr , internal_sigma )
!==========================================================================
implicit none
type(structure)             , intent(in)    :: system
type(STO_basis)             , intent(in)    :: basis(:)
type(R_eigen)               , intent(in)    :: QM
type(f_grid)  , allocatable , intent(inout) :: PDOS(:)
integer       , OPTIONAL    , intent(in)    :: nr
real*8        , OPTIONAL    , intent(in)    :: internal_sigma

! local variables ...
real*8  , allocatable :: erg_MO(:) , peaks(:) , func(:) 
real*8                :: sgm , sub_occupation , Project
integer               :: i ,  k , i1 , i2 , j , indx , npoints
integer               :: spin , S
logical , allocatable :: mask(:)

if( present(internal_sigma) ) then 
    sgm = internal_sigma
else
    sgm = sigma
end if

npoints = size( PDOS(nr)%grid )

gauss_norm = 1.d0 !  / (sgm*sqrt2PI)    <== for gauss_norm = 1 the gaussians are not normalized ...
two_sigma2 = 2.d0 * sgm*sgm
step = (DOS_range%fim-DOS_range%inicio) / float(npoints-1)

! states in the range [DOS_range%inicio,DOS_range%fim] ...
i1 = maxloc(QM%erg , 1 , QM%erg <  DOS_range%inicio) + 1
i2 = maxloc(QM%erg , 1 , QM%erg <= DOS_range%fim   ) 

n_of_DOS_states = i2 - i1 + 1

allocate( erg_MO(n_of_DOS_states) )

! find the energies in the range [DOS_range%inicio,DOS_range%fim]
erg_MO = QM%erg( i1 : i2 )

allocate( peaks (npoints)     , source = D_zero  )
allocate( func  (npoints)     , source = D_zero  )
allocate( mask  (size(basis)) , source = .false. )

forall(k=1:npoints) PDOS(nr)%grid(k) = (k-1)*step + DOS_range%inicio

! the total density of states
PDOS(nr)%peaks2 = 0.d0
PDOS(nr)%func2  = 0.d0

spin = merge(1,0,SOC)

do S = spin , -spin , -2

    indx = i1-1

    do j = 1 , n_of_DOS_states
    
        indx = indx + 1

        peaks = 0.d0
        where( dabs(PDOS(nr)%grid-erg_MO(j)) < (step/two) ) peaks = D_one
    
        func = 0.d0
        where( ((PDOS(nr)%grid-erg_MO(j))**2/two_sigma2) < 25.d0 ) func = gauss_norm*exp( -(PDOS(nr)%grid-erg_MO(j))**2 / two_sigma2 )

        select case ( S )
           case(0)
               Project = sum( QM%L(indx,:)*QM%R(:,indx) , basis(:)%residue == PDOS(nr)%residue ) 
               PDOS(nr)%peaks2(:,1) = PDOS(nr)%peaks2(:,1) + Project * peaks
               PDOS(nr)%func2(:,1)  = PDOS(nr)%func2 (:,1) + Project * func(:)

           case(1)  
               mask = ( basis(:)%residue == PDOS(nr)%residue ) .AND. ( basis(:)%S == S )
               Project = sum( QM%L(indx,:)*QM%R(:,indx) , mask )
               PDOS(nr)%peaks2(:,1) = PDOS(nr)%peaks2(:,1) + Project * peaks
               PDOS(nr)%func2(:,1)  = PDOS(nr)%func2 (:,1) + Project * func(:)

           case(-1)  
               mask = ( basis(:)%residue == PDOS(nr)%residue ) .AND. ( basis(:)%S == S )
               Project = sum( QM%L(indx,:)*QM%R(:,indx) , mask )
               PDOS(nr)%peaks2(:,2) = PDOS(nr)%peaks2(:,2) - Project * peaks
               PDOS(nr)%func2(:,2)  = PDOS(nr)%func2 (:,2) - Project * func(:)

           end select

    end do

end do

PDOS(nr)%average2 = PDOS(nr)%average2 + PDOS(nr)%func2

deallocate( peaks , func )

! occupation of PDOS(nr) ...
select case(spin)

       case(0)
           PDOS(nr)%occupation(1) = two * PDOS(nr)%peaks2(1,1)
           do k = 2 , npoints 
               PDOS(nr)%occupation(k) = PDOS(nr)%occupation(k-1) + two*PDOS(nr)%peaks2(k,1) 
           end do
           sub_occupation  = two * (i1 - 1) ! <== occupation below DOS window
           
           PDOS(nr)%occupation = sub_occupation + PDOS(nr)%occupation

       case(1)
           PDOS(nr)%occupation(1) = PDOS(nr)%peaks2(1,1) + abs(PDOS(nr)%peaks2(1,2))
           do k = 2 , npoints 
              PDOS(nr)%occupation(k) = PDOS(nr)%occupation(k-1) +     PDOS(nr)%peaks2(k,1)   ! <== spin up
              PDOS(nr)%occupation(k) = PDOS(nr)%occupation(k)   + abs(PDOS(nr)%peaks2(k,2))  ! <== spin down
           end do

           ! occupation below DOS window
           sub_occupation = 0.d0
           do j = 1 , i1-1
              sub_occupation = sub_occupation + sum( QM%L(j,:)*QM%R(:,j) , basis(:)%residue == PDOS(nr)%residue ) 
           end do
           
           PDOS(nr)%occupation = sub_occupation + PDOS(nr)%occupation

       end select 

print*, '>> ',PDOS(nr)%residue,' PDOS done <<'

end subroutine Partial_DOS
!
!
!
end module DOS_m

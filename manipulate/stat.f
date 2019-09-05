module Statistics_routines

use types_m
use constants_m
use Read_parms
use Correlation_m       , only : Correlation_driver

public :: Statistics , Most_Representative_Configuration

private

contains
!
!
!
!===========================
subroutine Statistics( trj )
!===========================
type(universe)  , allocatable   , intent(inout) :: trj(:)

! local varibles ...
integer      :: step
character(1) :: operation
character(2) :: atom
character(3) :: atom_A , atom_B , residue

CALL system( "clear" )

write(*,'(/a)') ' Choose operation analysis '
write(*,'(/a)') ' (l) = Bond Length '
write(*,'(/a)') ' (a) = Bond Angle '
write(*,'(/a)') ' (t) = Bond Torsion '
write(*,'(/a)') ' (r) = Radial Distribution Functions (RDF) '
write(*,'(/a)') ' (z) = Linear Distribution Functions (LDF) '
write(*,'(/a)') ' (c) = Correlation Analysis '
write(*,'(/a)',advance='no') '>>>   '
read (*,'(a)') operation

select case( operation )
    case( 'l' )
        CALL Bond_Length(trj)

    case( 'a' )
        CALL Bond_Angle(trj)

    case( 't' )
        CALL Bond_Torsion(trj)

    case( 'r' )
        write(*, '(1x,a)'               ) "MMSymbol of the first atom (A): "
        write(*, '(a10)', advance = 'no') 'atom A = '
        read*, atom_A
        atom_A = adjustr(atom_A)
        write(*, '(1x,a)'               ) "MMSymbol of the second atom (B): "
        write(*, '(a10)', advance = 'no') 'atom B = '
        read*, atom_B
        atom_B = adjustr(atom_B)
        write(*, '(1x,a)'               ) "Frame step: "
        write(*, '(a10)', advance = 'no') 'step = '
        read*, step

        CALL Radial_Function( trj , atom_A , atom_B , step )

    case( 'z' )
        write(*, '(1x,a  )'                ) "MMSymbol of the atom: "
        write(*, '(1x,a7 )', advance = 'no') 'atom = '
        read*, atom
        write(*, '(1x,a  )'                ) "Residue of the atom: "
        write(*, '(1x,a10)', advance = 'no') 'residue = '
        read*, residue
        write(*, '(1x,a  )'                ) "FRAME step: "
        write(*, '(1x,a7 )', advance = 'no') 'step = '
        read*, step

        CALL Linear_Function( trj , atom , residue , step )

    case( 'c' )
        CALL Correlation_driver

end select
!
!
!
end subroutine Statistics
!
!
!
!============================
subroutine Bond_Length( trj )
!============================
type(universe)  , intent(in) :: trj(:)

! local variables ....
real*8   , allocatable  :: d_AB(:)
real*8                  :: d_AB_avg , d_AB_sigma
integer                 :: i , j , indx1 , indx2

! calcule of the distance between the atoms ...
write (*, '(1x,a)'),  "Enter the index of atoms whose distance (d_AB) is calculated: "
write(*,'(A19)',advance='no') 'index of atom A = '
read (*,'(i3)') indx1
write(*,'(A19)',advance='no') 'index of atom B = '
read (*,'(i3)') indx2

allocate( d_AB(size(trj)) )

forall( j=1:size(trj) ) d_AB(j) = sqrt( sum((trj(j)%atom(indx1)%xyz - trj(j)%atom(indx2)%xyz)**2) )

d_AB_avg   = sum( d_AB ) / size(trj)
d_AB_sigma = sqrt( (sum(d_AB*d_AB - d_AB_avg*d_AB_avg)) / size(trj) )

print*, d_AB_avg , d_AB_sigma

deallocate(d_AB)
!
!
!
end subroutine Bond_Length
!
!
!
!===========================
subroutine Bond_Angle( trj )
!===========================
type(universe)  , intent(in)    :: trj(:)

! local variables ...
real*8  , allocatable   :: angle_ABC(:)
real*8                  :: d2_AB , d2_AC , d2_BC , angle_ABC_sigma , angle_ABC_avg
integer                 :: i , j , indx_A , indx_B , indx_C

! enter data to define the bond angle formed by atoms ABC ...
write (*, '(1x,a)'),  "index of the first atom (A): "
write(*,'(A10)',advance='no') 'atom A = '
read*, indx_A
write (*, '(1x,a)'),  "index of the midle atom (B): "
write(*,'(A10)',advance='no') 'atom B = '
read*, indx_B
write (*, '(1x,a)'),  "index of last atom (C): "
write(*,'(A10)',advance='no') 'atom C = '
read*, indx_C

allocate( angle_ABC(size(trj)) )

! calculates the bond angle ...
! calculates the bond angle ...
do i = 1 , size(trj)
    d2_AB = sum((trj(i)%atom(indx_A)%xyz - trj(i)%atom(indx_B)%xyz)**2)
    d2_AC = sum((trj(i)%atom(indx_A)%xyz - trj(i)%atom(indx_C)%xyz)**2)
    d2_BC = sum((trj(i)%atom(indx_B)%xyz - trj(i)%atom(indx_C)%xyz)**2)

    angle_ABC(i) = dacos(( - d2_AC + d2_AB + d2_BC) / (2 * sqrt(d2_AB) * sqrt(d2_BC) ))
end do

angle_ABC = radian * angle_ABC

! dump dihedro ...
open(unit = 20, file = "angle.dat", status = "unknown", action = "write")
do i = 1 , size(trj)
    write(20,500)  angle_ABC(i)
end do
close(20)

500 format(f17.13)

angle_ABC_avg   = sum( angle_ABC ) / size( trj )
angle_ABC_sigma = sqrt( (sum(angle_ABC*angle_ABC - angle_ABC_avg*angle_ABC_avg)) / size(trj) )

print*, angle_ABC_avg , angle_ABC_sigma

deallocate(angle_ABC)
!
!
!
end subroutine Bond_Angle
!
!
!
!=============================
subroutine Bond_Torsion( trj )
!=============================
type(universe)  , intent(in)    :: trj(:)

! local variables ...
real*8  , allocatable   , dimension(:) :: atom_i , atom_j , atom_k , atom_l , ij , jk , lk , dihedro , m , n
real*8                                 :: angle_dihedro , cosseno , seno  , pico_1 , pico_2 , number_angle_1(50) , number_angle_2(50)
integer                                :: indx_i , indx_j , indx_k , indx_l , frame , t , div
integer , parameter                    :: a = 1

! enter data to define the bond angle formed by atoms ABC ...
write (*, '(1x,a)'),  "index of the atom i: "
write(*,'(A10)',advance='no') 'index i = '
read*, indx_i
write (*, '(1x,a)'),  "index of the atom j: "
write(*,'(A10)',advance='no') 'index j = '
read*, indx_j
write (*, '(1x,a)'),  "index of the atom k: "
write(*,'(A10)',advance='no') 'index k = '
read*, indx_k
write (*, '(1x,a)'),  "index of the atom l "
write(*,'(A10)',advance='no') 'index l = '
read*, indx_l

allocate( atom_i        ( 3 )           )
allocate( atom_j        ( 3 )           )
allocate( atom_k        ( 3 )           )
allocate( atom_l        ( 3 )           )
allocate( ij            ( 3 )           )
allocate( jk            ( 3 )           )
allocate( lk            ( 3 )           )
allocate( dihedro       ( size(trj) )   )
allocate( m             ( 3 )           )
allocate( n             ( 3 )           )

do frame = 1 , size(trj)

    atom_i(:) = trj(frame) % atom(indx_i) % xyz(:)
    atom_j(:) = trj(frame) % atom(indx_j) % xyz(:)
    atom_k(:) = trj(frame) % atom(indx_k) % xyz(:)
    atom_l(:) = trj(frame) % atom(indx_l) % xyz(:)

    forall(t=1:3)
        ij(t) = atom_j(t) - atom_i(t)
        jk(t) = atom_k(t) - atom_j(t)
        lk(t) = atom_k(t) - atom_l(t)
    end forall

!   normal to plane ijk ...
    m(1) = jk(2) * ij(3) - jk(3) * ij(2)
    m(2) = jk(3) * ij(1) - jk(1) * ij(3)
    m(3) = jk(1) * ij(2) - jk(2) * ij(1)

!   normal to plane jkl ...
    n(1) = lk(2) * jk(3) - lk(3) * jk(2)
    n(2) = lk(3) * jk(1) - lk(1) * jk(3)
    n(3) = lk(1) * jk(2) - lk(2) * jk(1)

    cosseno = ( m(1)*n(1) + m(2)*n(2) + m(3)*n(3) ) / &
              ( (sqrt(m(1)**2 + m(2)**2 + m(3)**2)) * (sqrt(n(1)**2 + n(2)**2 + n(3)**2)) )

    seno    = ( (n(1)*ij(1) + n(2)*ij(2) + n(3)*ij(3)) * sqrt(jk(1)**2 + jk(2)**2 + jk(3)**2) ) / &
              ( (sqrt(m(1)**2 + m(2)**2 + m(3)**2)) * (sqrt(n(1)**2 + n(2)**2 + n(3)**2)) )

!   calculate of the angle dihedro ...
    if( seno < 0 ) then
        dihedro(frame) = ( PI - acos(cosseno) ) * radian
    else
        dihedro(frame) = acos(cosseno) * radian
    end if

end do

! compute dihedro ...
angle_dihedro = sum(dihedro(:)) / size(trj)
print'("Angle dihedro =")'
print*, angle_dihedro

! dump dihedro ...
open(unit = 20, file = "dihedro.dat", status = "unknown", action = "write")
do t = 1 , size(trj)
    write(20,500)  dihedro(t)
end do
close(20)

500 format(f18.15)

deallocate( atom_i , atom_j , atom_k , atom_l , ij , jk , lk , dihedro , m , n )

end subroutine Bond_Torsion
!
!
!
!=========================================================
subroutine Radial_Function( trj , atom_A , atom_B , step )
!=========================================================
type(universe)  , intent(in)    :: trj(:)
character(3)    , intent(in)    :: atom_A , atom_B
integer         , intent(in)    :: step

! local variables ...
type(atomic)    , dimension(:)   , allocatable :: trj_PBC
real*8          , dimension(:,:) , allocatable :: vec_A , vec_B , distance_AB , g_AB , g_AB_from_A
real*8          , dimension(:)   , allocatable :: x , g_AB_total 
real*8                                         :: side(3) , radius_max , delta_R , V_local , rho_B_bulk
integer                                        :: l , i , j , n , N_A , N_B , frame , sampling_number , PBC_sys_size
integer         , dimension(:)   , allocatable :: resid_A , resid_B , index_max 
logical         , dimension(:,:) , allocatable :: mask

! local parameters ...
integer                          , parameter   :: N_interval = 500

sampling_number = int( (size(trj) - 1) / step ) + 1

! define the maximum radius
forall( i=1:3 ) side(i) = maxval( trj(1) % atom % xyz(i) ) - minval( trj(1) % atom % xyz(i) )
radius_max  =  minval( side ) / two
delta_R     =  radius_max / N_interval

allocate( x(N_interval-1) )     ;   forall( n=1:N_interval-1 ) x(n) = delta_R * n

! number of atoms types A in the central unit_cell ...
N_A =    count( trj(1) % atom % MMSymbol == adjustl(atom_A) )

! number of atoms types B in central unit_cell + 26 PBC unit_cells ...
N_B = 27*count( trj(1) % atom % MMSymbol == adjustl(atom_B) )

PBC_sys_size = trj(1)%N_of_atoms * 27

allocate( distance_AB       (N_B,N_A)               )
allocate( mask              (N_B,N_A)               )
allocate( vec_A             (N_A,3)                 )
allocate( vec_B             (N_B,3)                 )
allocate( resid_A           (N_A)                   )
allocate( resid_B           (N_B)                   )
allocate( index_max         (N_A)                   )
allocate( g_AB_from_A       (N_interval-1,N_A)      )
allocate( g_AB              (N_interval-1,size(trj)))
allocate( trj_PBC           (PBC_sys_size)          )

do frame = 1 , size(trj) , step

    print*, frame , "/", size(trj)

    forall( i=1:3 ) vec_A(:,i) = pack( trj(frame) % atom(:) % xyz(i) , trj(frame) % atom(:) % MMSymbol == adjustl(atom_A) )

    resid_A(:) = pack( trj(frame) % atom(:) % nresid , trj(frame) % atom(:) % MMSymbol == adjustl(atom_A) )

    CALL Apply_PBC( trj(frame) , trj_PBC )

    forall( i=1:3 ) vec_B(:,i) = pack( trj_PBC(:) % xyz(i) , trj_PBC(:) % MMSymbol == adjustl(atom_B) )

    resid_B(:) = pack( trj_PBC(:) % nresid , trj_PBC(:) % MMSymbol == adjustl(atom_B) )

!   distance of atom_B to the atom_A ...
    forall( i=1:N_A , j=1:N_B )
        distance_AB (j,i)  =  sqrt( sum((vec_B(j,:) - vec_A(i,:))**2) )
        mask        (j,i)  =  ( resid_B(j) == resid_A(i) )
    end forall

!   eliminate statistics between the same residue
    where( mask ) distance_AB = real(ABOVE)

    do i = 1 , N_A
        CALL sort( distance_AB(:,i) )
    end do

    index_max = maxloc( distance_AB , 1 , distance_AB <= radius_max )

!   rho_B_bulk = bulk density of particle type B 
    rho_B_bulk = float( N_B ) / (27*product( trj(frame)%box(:) ))

!   RDF from atom A(i)
    do i = 1 , N_A

        forall( n=1:N_interval-1 ) g_AB_from_A(n,i)     =   &
        count( (delta_R*n < distance_AB(1:index_max(i),i)) .AND. (distance_AB(1:index_max(i),i) <= delta_R*(n+1)) ) / ( 4.0*PI*x(n)**2 * delta_R )

    end do

!   compute final g_AB ...
    forall( n=1:N_interval-1 ) g_AB(n,frame) = (1.d0 / rho_B_bulk) * ( 1.d0 / N_A ) * sum(g_AB_from_A(n,:))

end do

! compute final g_AB ...
allocate( g_AB_total(N_interval-1) , source=0.d0 )

do frame = 1 , size(trj) , step
    forall( n=1:N_interval-1 ) g_AB_total(n) = g_AB_total(n) + g_AB(n,frame) / sampling_number
end do

! dump RDF ...
open(unit = 20, file = "data.dat", status = "unknown", action = "write")
    do n = 1 , N_interval-1
        write(20,500)  x(n) , g_AB_total(n)
    end do
500 format(f10.5, t27, f10.5)
close(20)

deallocate( trj_PBC , vec_A , vec_B , g_AB_from_A , distance_AB , mask , x , g_AB , g_AB_total )

end subroutine Radial_Function
!
!
!
!========================================================
subroutine Linear_Function( trj , atom , residue , step )
!========================================================
type(universe)  , intent(in)    :: trj(:)
character(*)    , intent(in)    :: atom
character(*)    , intent(in)    :: residue
integer         , intent(in)    :: step

! local variables ...
real*8  , dimension(:,:) , allocatable  :: LDF
real*8  , dimension(:)   , allocatable  :: x , location , distance , LDF_final
real*8                                  :: delta_L , Area , surface , norm
integer                                 :: l , i , j , n , N_atom , frame , sampling_number

! local parameters ...
integer , parameter :: N_interval = 500
real*8  , parameter :: L_max = 30.0

! parameters
N_atom          = count( (trj(1)%atom(:)%MMSymbol == atom) .AND. (trj(1)%atom(:)%resid == residue) )
Area            = trj(1)%box(1) * trj(1)%box(2)
delta_L         = L_max / N_interval

allocate( x(N_interval+1) )       ;      forall( n=1:N_interval+1 ) x(n) = delta_L * (n - 1)

allocate( location  ( N_atom               ) )
allocate( distance  ( N_atom               ) )
allocate( LDF       ( size(trj),N_interval ) )
allocate( LDF_final ( N_interval           ) )

sampling_number = 0

do frame = 1 , size(trj) , step

    sampling_number = sampling_number + 1

    location(:) = pack( trj(frame) % atom(:) % xyz(3) , (trj(frame) % atom(:) % MMSymbol == atom) .AND. (trj(frame) % atom(:) % resid == residue) )

    ! normalization factor for producing LDF=1 in the bulk ; surface is the onset of solvent ...
    surface  =  maxval( trj(frame)%atom%xyz(3) , trj(frame)%atom(:)%fragment /= "S" )
    norm     =  float(N_atom) / ((trj(frame)%box(3)-surface) / delta_L)

    ! new surface calculation for measuring the distance from the solid substrate ...
    surface = maxval( trj(frame)%atom%xyz(3) , trj(frame)%atom(:)%resid == "CCC" )

    forall( i=1:N_atom ) distance(i) = location(i) - surface

    forall( n=1:N_interval ) LDF(sampling_number,n) = count( x(n) < distance(:) .and. distance(:) <= x(n+1) )

end do

do n = 1 , N_interval
!    LDF_final(n) = sum( LDF(1:sampling_number,n) ) / sampling_number
    LDF_final(n) = (1.0d0 / norm ) * sum( LDF(1:sampling_number,n) ) / sampling_number
end do

! dump LDF ...
open(unit = 20, file = "data.dat", status = "unknown", action = "write")
do i = 1 , N_interval
    write(20,500)  x(i) , LDF_final(i)
end do
close(20)

deallocate( location , x , LDF , LDF_final , distance )

500 format(f11.6, t27, f11.6)

end subroutine Linear_Function
!
!
!
!====================================
subroutine Apply_PBC( trj , trj_PBC )
!====================================
type(universe)  , intent(in)        :: trj
type(atomic)    , intent(inout)     :: trj_PBC(:)

! local variables ...
integer                 :: i , j , k , n , counter

! local parameters ; number of 3D PBC unit-cells ...
integer ,   parameter   :: Replication_Factor = 27

! replicating the central cell to the surrounding cells ...
forall( i = 1:Replication_Factor ) trj_PBC( trj%N_of_atoms*(i-1)+1 : trj%N_of_atoms*i ) = trj%atom

! defining the coordinates for the surrounding cells ...
counter = 0
do k = -1,+1
do j = -1,+1
do i = -1,+1

    forall( n = 1:trj%N_of_atoms )

        trj_PBC(counter+n) % xyz(1) = trj % atom(n) % xyz(1) + i * trj % box(1)
        trj_PBC(counter+n) % xyz(2) = trj % atom(n) % xyz(2) + j * trj % box(2)
        trj_PBC(counter+n) % xyz(3) = trj % atom(n) % xyz(3) + k * trj % box(3)

    end forall

    counter = counter + trj%N_of_atoms

end do
end do
end do

end subroutine Apply_PBC
!
!
!
!=====================
subroutine  sort( ra )
!=====================
real*8  , intent(inout) :: ra(:)

! local variables
real    :: rra
integer :: l, n, i, j, ir

      n = size( ra(:) )

      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         rra  = ra(l)
      else
         rra  = ra(ir)
         ra(ir)  = ra(1)
         ir = ir - 1
         if(ir .eq. 1) then
             ra(1)  = rra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(ra(j) .lt. ra(j+1)) j = j + 1
        endif
      if(rra .lt. ra(j)) then
        ra(i)  = ra(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      ra(i)  = rra
      goto 10

end subroutine sort
!
!
!
!=========================================================
 subroutine Most_Representative_Configuration( trj , sys )
!=========================================================
type(universe)  , allocatable , intent(inout)  :: trj(:)
type(universe)                , intent(out)    :: sys

! local variables ....
real*8   , allocatable  :: xyz(:,:,:) , cost(:)
real*8                  :: soma
integer                 :: i , j , k , typical
character(1)            :: answer
logical  , allocatable  :: mask(:,:)

! create work matrix to emulate trj%atom%xyz ...
allocate( xyz( size(trj) , trj(1)%N_of_atoms , 3 ) )
do k = 1 , 3
    do j = 1 , trj(1)%N_of_atoms
        do i = 1 , size(trj)
            xyz(i,j,k) = trj(i)%atom(j)%xyz(k)
        end do
    end do
end do

! looking for moving atoms ...
allocate( mask( trj(1)%N_of_atoms , 3 ) , source = .false. )

do j = 1 , trj(1)%N_of_atoms
    do k = 1 , 3
        mask(j,k) = any( xyz(:,j,k) /= trj(1)%atom(j)%xyz(k) ) 
    end do
end do

allocate( cost(size(trj)) , source = 0.d0)
!============================================================================================
! most representative configuration has the lowest cost ... 
!do i1 = 1 , size(trj)
!    do i2 = 1 , size(trj)
!        If( i1 /= i2 ) cost(i1) = cost(i1) + sum( (xyz(i1,:,:)-xyz(i2,:,:)) * (xyz(i1,:,:)-xyz(i2,:,:)) , mask )
!    end do
!end do

!------------------------

! this algorithm does the same thing twice as fast but it is twice as unclear ...
cost = 0.d0
!$omp parallel do schedule(dynamic,300) default(shared) private(i1,i2,j,k,soma)
do i1 = 1 , size(trj)
    do i2 = 1 , size(trj)
    if(i1/=i2) then
        soma = 0.d0
        do j = 1 , trj(1)%N_of_atoms
        do k = 1 , 3

            If( mask(j,k) ) soma = soma + (xyz(i1,j,k) - xyz(i2,j,k))*(xyz(i1,j,k) - xyz(i2,j,k))

        end do
        end do
        cost(i1) = cost(i1) + soma
    end if
    end do
end do    
!$omp end parallel do
!============================================================================================

typical = minloc( cost , dim=1 )

write(*,'(/a,I5)') ' Most Representative Configuration is MODEL = ', typical, ' X frame_step'

write(*,'(/a)') ' >>>  Saving RMSD.dat ' 

OPEN(unit=9,file='RMSD.dat',status='unknown')
do i = 1 , size(trj)
    write(9,*) i , cost(i)
end do

write(*,'(/a)',advance='no') " >>>  Save most typical frame ? (y/n) : "
read (*,'(a)') answer

If( answer == "y" ) then
    allocate( sys%atom(size(trj(1)%atom)) )
    sys = trj(typical)
    sys%Surface_Characteristics = trj(1)%Surface_Characteristics

    deallocate(trj)
end If

deallocate( xyz , mask , cost )

end subroutine Most_Representative_Configuration
!
!
!
end module Statistics_routines

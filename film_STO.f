 module Psi_squared_cube_format

    use type_m
    use constants_m
    use parameters_m            , only : initial_state
    use Babel_m                 , only : System_Characteristics
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : Unit_Cell
    use Slater_Type_Orbitals    , only : s_orb , p_orb , d_x2y2 , d_z2 , d_xyz  

    public :: Gaussian_Cube_Format 

    private

    real*8 , parameter :: fringe = 8.d0

 contains   
!
!
!
!==============================================================
 subroutine  Gaussian_Cube_Format( bra , ket , it , t , el_hl )
!==============================================================
 implicit none
 complex*16                , intent(in) :: bra(:), ket(:)
 real*8                    , intent(in) :: t
 integer                   , intent(in) :: it 
 character(*) , optional   , intent(in) :: el_hl

! local variables ...
 complex*16 , allocatable :: Psi_bra(:) , Psi_ket(:)
 real*8     , allocatable :: xyz(:,:)
 real*8                   :: x , y , z , x0 , y0 , z0 , Psi_2 , dx , dy , dz , a , b , c , r , SlaterOrbital
 integer                  :: AtNo , nx_steps , ny_steps , nz_steps , i , j , ix , iy , iz , k 
 character(len=4)         :: string 
 character(len=23)        :: f_name

 real*8 , parameter :: aB = 0.529d0   ! <== Bohr radius

 allocate(xyz(unit_cell%atoms,3))
 allocate(Psi_bra(size(bra)) , Psi_ket(size(ket)))

 write(string,'(i4.4)') it
 f_name = el_hl//'_dens_shot'//string//'.cube'
 OPEN(unit=4,file=f_name,status='unknown')  

! bounding box for isosurfaces ... 
 CALL BoundingBox( unit_cell ) 

!! ATTENTION !! must have (n_step+1)mod(6) = 0
 nx_steps = 60 
 ny_steps = 71 
 nz_steps = 110 

!  voxel dimensions
 dx = unit_cell%BoxSides(1) / nx_steps 
 dy = unit_cell%BoxSides(2) / ny_steps 
 dz = unit_cell%BoxSides(3) / nz_steps      

!  translation to the center of mass
 forall(i=1:unit_cell%atoms,j=1:3) xyz(i,j) = unit_cell%coord(i,j) - unit_cell%Center_of_Mass(j)

! initial corner of the volume Box 
 a = minval(xyz(:,1)) - fringe / two
 b = minval(xyz(:,2)) - fringe / two
 c = minval(xyz(:,3)) - fringe / two

!  start writing Gaussian cube files 
 write(4,*) System_Characteristics
 write(4,*) 'initial_state = ',initial_state,'  /  time = ', t

 write(4,111) unit_cell%atoms , a/aB , b/aB , c/aB
 write(4,111) nx_steps + 1 , dx/aB   , 0.d0 , 0.d0 
 write(4,111) ny_steps + 1 , 0.d0 , dy/aB   , 0.d0
 write(4,111) nz_steps + 1 , 0.d0 , 0.d0 , dz/aB

 DO i = 1 , unit_cell%atoms
 
    write(4,113) unit_cell%AtNo(i), 0.0 , xyz(i,1)/aB , xyz(i,2)/aB , xyz(i,3)/aB

 END DO

!---------------------------------------------------------    
!  drawing the wavefunction denssity in cube format
!========================================================
!  the order orbitals are stored
! 
!       S      -->  1   --> l = 0  ,  m =  0           
!       Py     -->  2   --> l = 1  ,  m = -1    
!       Pz     -->  3   --> l = 1  ,  m =  0         
!       Px     -->  4   --> l = 1  ,  m = +1
!       Dxy    -->  5   --> l = 2  ,  m = -2      
!       Dyz    -->  6   --> l = 2  ,  m = -1
!       Dz2    -->  7   --> l = 2  ,  m =  0     
!       Dxz    -->  8   --> l = 2  ,  m = +1        
!       Dx2y2  -->  9   --> l = 2  ,  m = +2        
!========================================================
 
 DO ix = 0 , nx_steps
     x = a + ix * dx
 
 DO iy = 0 , ny_steps
     y = b + iy * dy
 
 DO iz = 0 , nz_steps 
     z = c + iz * dz

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

   i = 0 
   DO k = 1 , unit_cell%atoms

     AtNo = unit_cell%AtNo(k)

! distance from the center of the nuclei

     r = dsqrt((x-xyz(k,1))*(x-xyz(k,1)) + (y-xyz(k,2))*(y-xyz(k,2)) + (z-xyz(k,3))*(z-xyz(k,3))) 

! coordinates centered on the nuclei   

    x0 = x - xyz(k,1) 
    y0 = y - xyz(k,2) 
    z0 = z - xyz(k,3) 

     do j = 1 , atom(AtNo)%DOS
        i = i + 1
        select case (j)
            case( 1 ) 
                SlaterOrbital = s_orb(r,AtNo)  
            case( 2 ) 
                SlaterOrbital = p_orb(r,y0,AtNo)
            case( 3 ) 
                SlaterOrbital = p_orb(r,z0,AtNo)
            case( 4 ) 
                SlaterOrbital = p_orb(r,x0,AtNo)
            case( 5 ) 
                SlaterOrbital = d_xyz(r,x0,y0,AtNo)
            case( 6 ) 
                SlaterOrbital = d_xyz(r,y0,z0,AtNo)
            case( 7 ) 
                SlaterOrbital = d_z2(r,x0,y0,z0,AtNo)
            case( 8 ) 
                SlaterOrbital = d_xyz(r,x0,z0,AtNo)
            case( 9 ) 
                SlaterOrbital = d_x2y2(r,x0,y0,AtNo)
        end select
        Psi_bra(i) = bra(i) * SlaterOrbital
        Psi_ket(i) = ket(i) * SlaterOrbital
     end do  ! <== DOS
   end do  ! <== atoms

   Psi_2 = cdabs( cdsqrt( sum(Psi_bra)*sum(Psi_ket) ) )
   
   write(4,112,advance='no') Psi_2

   If( (mod(iz,6) == 5) ) write(4,'(a)',advance='yes') 

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 END DO                                                            ! <==  Z coord
 write(4,'(a)',advance='yes') 
 END DO                                                            ! <==  Y coord
 END DO                                                            ! <==  X coord

 close(4)

 deallocate(xyz, Psi_bra, Psi_ket)

 include 'formats.h'

 end subroutine Gaussian_Cube_Format
!
!
!
!===========================
 subroutine BoundingBox( a )
!===========================
 implicit none
 type(structure) :: a

! local variables ...
 integer :: i

!  size of the box
 forall(i=1:3) &
 a%BoxSides(i) = maxval(a%coord(:,i)) - minval(a%coord(:,i)) + fringe

!  find the center of mass
 forall(i=1:3) &
 a%Center_of_Mass(i) = sum(a%coord(:,i)) / a%atoms

end subroutine BoundingBox
!
!
end module Psi_squared_cube_format

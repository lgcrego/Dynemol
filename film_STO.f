 module Psi_squared_cube_format

    use type_m
    use constants_m
    use parameters_m            , only : electron_state
    use Babel_m                 , only : System_Characteristics
    use Semi_Empirical_Parms    , only : atom
    use Structure_Builder       , only : Extended_Cell
    use Slater_Type_Orbitals    , only : s_orb , p_orb , d_x2y2 , d_z2 , d_xyz  

    public :: Gaussian_Cube_Format 

    private

    real*8 , parameter :: fringe = 8.d0

    interface Gaussian_Cube_Format
        module procedure Gaussian_Cube_Format_Real
        module procedure Gaussian_Cube_Format_Cmplx
    end interface  Gaussian_Cube_Format

 contains   
!
!
!
!===================================================================
 subroutine  gaussian_cube_format_Real( bra , ket , it , t , el_hl )
!===================================================================
 implicit none
 real*8                    , intent(in) :: bra(:), ket(:)
 real*8                    , intent(in) :: t
 integer                   , intent(in) :: it 
 character(*) , optional   , intent(in) :: el_hl

! local variables ...
 integer                  :: many_steps   = 131
 integer                  :: medium_steps = 89
 integer                  :: few_steps    = 59  
 integer                  :: n_xyz_steps(3)
 complex*16               :: TotalPsiKet
 real*8     , allocatable :: xyz(:,:) , Psi(:,:,:)
 real*8                   :: x , y , z , x0 , y0 , z0 , dx , dy , dz , a , b , c , r , SlaterOrbital, dumb
 integer                  :: AtNo , i , j , ix , iy , iz , k , l
 character(len=2)         :: prefix
 character(len=5)         :: string 
 character(len=31)        :: f_name

 ! bra is never used, so to avoid compiler warnings ...
 dumb = bra(1)

 allocate(xyz(extended_cell%atoms,3))

 write(string,'(i5.5)') it
 prefix = merge( "el" , el_hl , .NOT. present(el_hl) )
 f_name = 'MO_trunk/'//prefix//'_MO_shot'//string//'.cube'
 OPEN(unit=4,file=f_name,status='unknown')  

! bounding box for isosurfaces ... 
 CALL BoundingBox( extended_cell ) 

! fix number of steps for each direction according to aspect ratio ... 
 do i = 1 , 3

    ! default case
    n_xyz_steps(i) = medium_steps

    If( extended_cell%BoxSides(i) >= 2.5d1 ) n_xyz_steps(i) = many_steps

    If( extended_cell%BoxSides(i) <= 1.2d1 ) n_xyz_steps(i) = few_steps

 end do

!  voxel dimensions
 dx = extended_cell%BoxSides(1) / n_xyz_steps(1) 
 dy = extended_cell%BoxSides(2) / n_xyz_steps(2) 
 dz = extended_cell%BoxSides(3) / n_xyz_steps(3)      

!  translation to the center of mass
 forall(i=1:extended_cell%atoms,j=1:3) xyz(i,j) = extended_cell%coord(i,j) - extended_cell%Center_of_Mass(j)

!  initial corner of the volume Box 
 a = minval(xyz(:,1)) - fringe / two
 b = minval(xyz(:,2)) - fringe / two
 c = minval(xyz(:,3)) - fringe / two

!  start writing Gaussian cube files 
 write(4,*) System_Characteristics
 write(4,*) 'electron_state = ',electron_state,'  /  time = ', t

 write(4,111) extended_cell%atoms , a/a_Bohr , b/a_Bohr , c/a_Bohr
 write(4,111) n_xyz_steps(1) + 1 , dx/a_Bohr , 0.d0 , 0.d0 
 write(4,111) n_xyz_steps(2) + 1 , 0.d0 , dy/a_Bohr , 0.d0
 write(4,111) n_xyz_steps(3) + 1 , 0.d0 , 0.d0 , dz/a_Bohr

!  coordinates in a.u., because zeta is in units of [a_0^{-1}] ...
 xyz = xyz / a_Bohr
 DO i = 1 , extended_cell%atoms
 
    write(4,113) extended_cell%AtNo(i), 0.0 , xyz(i,1) , xyz(i,2) , xyz(i,3)

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

 allocate( Psi(n_xyz_steps(1)+1,n_xyz_steps(2)+1,n_xyz_steps(3)+1) , source = D_zero ) 

!$OMP parallel private(ix,iy,iz,x,y,z,x0,y0,z0,SlaterOrbital,r,AtNo,i,j,k,l,TotalPsiKet)
!$OMP single
 DO ix = 0 , n_xyz_steps(1)
    x = (a + ix * dx) / a_Bohr
    
    !$OMP task untied
    DO iy = 0 , n_xyz_steps(2)
        y = (b + iy * dy) / a_Bohr
 
    DO iz = 0 , n_xyz_steps(3) 
        z = (c + iz * dz) / a_Bohr

        i = 0 
        TotalPsiKet = C_zero 
        DO k = 1 , extended_cell%atoms

            if( Extended_Cell% QMMM(k) /= "QM" ) cycle

            AtNo = extended_cell%AtNo(k)

            ! distance from the center of the nuclei, in a.u.

            r = dsqrt((x-xyz(k,1))*(x-xyz(k,1)) + (y-xyz(k,2))*(y-xyz(k,2)) + (z-xyz(k,3))*(z-xyz(k,3))) 

            ! coordinates centered on the nuclei, in a.u.   

            x0 = x - xyz(k,1) 
            y0 = y - xyz(k,2) 
            z0 = z - xyz(k,3) 

            do j = 1 , atom(AtNo)%DOS

                i = i + 1
                l = extended_cell%BasisPointer(k) + j

                select case (j)
                    case( 1 ) 
                        SlaterOrbital = s_orb(r,AtNo,l)  
                    case( 2 ) 
                        SlaterOrbital = p_orb(r,y0,AtNo,l)
                    case( 3 ) 
                        SlaterOrbital = p_orb(r,z0,AtNo,l)
                    case( 4 ) 
                        SlaterOrbital = p_orb(r,x0,AtNo,l)
                    case( 5 ) 
                        SlaterOrbital = d_xyz(r,x0,y0,AtNo,l)
                    case( 6 ) 
                        SlaterOrbital = d_xyz(r,y0,z0,AtNo,l)
                    case( 7 ) 
                        SlaterOrbital = d_z2(r,x0,y0,z0,AtNo,l)
                    case( 8 ) 
                        SlaterOrbital = d_xyz(r,x0,z0,AtNo,l)
                    case( 9 ) 
                        SlaterOrbital = d_x2y2(r,x0,y0,AtNo,l)
                end select

                TotalPsiKet = TotalPsiKet + ket(i) * SlaterOrbital

            end do  ! <== DOS
        end do  ! <== atoms

        Psi( ix+1 , iy+1 , iz+1 ) = TotalPsiKet

    END DO          ! <==  Z coord
    END DO          ! <==  Y coord
    !$OMP end task 
 END DO             ! <==  X coord
!$OMP end single 
!$OMP end parallel 

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 DO ix = 0 , n_xyz_steps(1)
 
    DO iy = 0 , n_xyz_steps(2)
 
        DO iz = 0 , n_xyz_steps(3) 

            write(4,112,advance='no') Psi( ix+1 , iy+1 , iz+1 ) 

            If( (mod(iz,6) == 5) ) write(4,'(a)',advance='yes') 

        END DO                                                       
        write(4,'(a)',advance='yes') 
    END DO                                                      
 END DO                                                     

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 close(4)

 deallocate(xyz , Psi)

 include 'formats.h'

 end subroutine Gaussian_Cube_Format_Real
!
!
!
!====================================================================
 subroutine  gaussian_cube_format_Cmplx( bra , ket , it , t , el_hl )
!====================================================================
 implicit none
 complex*16                , intent(in) :: bra(:), ket(:)
 real*8                    , intent(in) :: t
 integer                   , intent(in) :: it 
 character(*) , optional   , intent(in) :: el_hl

! local variables ...
 integer                  :: many_steps   = 131
 integer                  :: medium_steps = 89
 integer                  :: few_steps    = 59  
 integer                  :: n_xyz_steps(3)
 complex*16               :: TotalPsiBra , TotalPsiKet
 real*8     , allocatable :: xyz(:,:) , Psi_2(:,:,:)
 real*8                   :: x , y , z , x0 , y0 , z0 , dx , dy , dz , a , b , c , r , SlaterOrbital
 integer                  :: AtNo , i , j , ix , iy , iz , k , l
 character(len=2)         :: prefix
 character(len=5)         :: string 
 character(len=31)        :: f_name

 allocate(xyz(extended_cell%atoms,3))

 write(string,'(i5.5)') it
 prefix = merge( "el" , el_hl , .NOT. present(el_hl) )
 f_name = 'MO_trunk/'//prefix//'_dens_shot'//string//'.cube'
 OPEN(unit=4,file=f_name,status='unknown')  

! bounding box for isosurfaces ... 
 CALL BoundingBox( extended_cell ) 

! fix number of steps for each direction according to aspect ratio ... 
 do i = 1 , 3

    ! default case
    n_xyz_steps(i) = medium_steps

    If( extended_cell%BoxSides(i) >= 2.5d1 ) n_xyz_steps(i) = many_steps

    If( extended_cell%BoxSides(i) <= 1.2d1 ) n_xyz_steps(i) = few_steps

 end do

!  voxel dimensions
 dx = extended_cell%BoxSides(1) / n_xyz_steps(1) 
 dy = extended_cell%BoxSides(2) / n_xyz_steps(2) 
 dz = extended_cell%BoxSides(3) / n_xyz_steps(3)      

!  translation to the center of mass
 forall(i=1:extended_cell%atoms,j=1:3) xyz(i,j) = extended_cell%coord(i,j) - extended_cell%Center_of_Mass(j)

! initial corner of the volume Box 
 a = minval(xyz(:,1)) - fringe / two
 b = minval(xyz(:,2)) - fringe / two
 c = minval(xyz(:,3)) - fringe / two

!  start writing Gaussian cube files 
 write(4,*) System_Characteristics
 write(4,*) 'electron_state = ',electron_state,'  /  time = ', t

 write(4,111) extended_cell%atoms , a/a_Bohr , b/a_Bohr , c/a_Bohr
 write(4,111) n_xyz_steps(1) + 1  , dx/a_Bohr , 0.d0 , 0.d0 
 write(4,111) n_xyz_steps(2) + 1  , 0.d0 , dy/a_Bohr , 0.d0
 write(4,111) n_xyz_steps(3) + 1  , 0.d0 , 0.d0 , dz/a_Bohr

!  coordinates in a.u., because zeta is in units of [a_0^{-1}] ...
 xyz = xyz / a_Bohr
 DO i = 1 , extended_cell%atoms
 
    write(4,113) extended_cell%AtNo(i), 0.0 , xyz(i,1) , xyz(i,2) , xyz(i,3)

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

 allocate( Psi_2(n_xyz_steps(1)+1,n_xyz_steps(2)+1,n_xyz_steps(3)+1) , source = D_zero ) 

!$OMP parallel private(ix,iy,iz,x,y,z,x0,y0,z0,SlaterOrbital,r,AtNo,i,j,k,l,TotalPsiBra,TotalPsiKet)
!$OMP single
 DO ix = 0 , n_xyz_steps(1)
    x = (a + ix * dx) / a_Bohr
    
    !$OMP task untied
    DO iy = 0 , n_xyz_steps(2)
        y = (b + iy * dy) / a_Bohr
 
    DO iz = 0 , n_xyz_steps(3) 
        z = (c + iz * dz) / a_Bohr

        i = 0 
        TotalPsiBra = C_zero 
        TotalPsiKet = C_zero 
        DO k = 1 , extended_cell%atoms

            if( Extended_Cell% QMMM(k) /= "QM" ) cycle

            AtNo = extended_cell%AtNo(k)

            ! distance from the center of the nuclei

            r = dsqrt((x-xyz(k,1))*(x-xyz(k,1)) + (y-xyz(k,2))*(y-xyz(k,2)) + (z-xyz(k,3))*(z-xyz(k,3))) 

            ! coordinates centered on the nuclei   

            x0 = x - xyz(k,1) 
            y0 = y - xyz(k,2) 
            z0 = z - xyz(k,3) 

            do j = 1 , atom(AtNo)%DOS

                i = i + 1
                l = extended_cell%BasisPointer(k) + j

                select case (j)
                    case( 1 ) 
                        SlaterOrbital = s_orb(r,AtNo,l)  
                    case( 2 ) 
                        SlaterOrbital = p_orb(r,y0,AtNo,l)
                    case( 3 ) 
                        SlaterOrbital = p_orb(r,z0,AtNo,l)
                    case( 4 ) 
                        SlaterOrbital = p_orb(r,x0,AtNo,l)
                    case( 5 ) 
                        SlaterOrbital = d_xyz(r,x0,y0,AtNo,l)
                    case( 6 ) 
                        SlaterOrbital = d_xyz(r,y0,z0,AtNo,l)
                    case( 7 ) 
                        SlaterOrbital = d_z2(r,x0,y0,z0,AtNo,l)
                    case( 8 ) 
                        SlaterOrbital = d_xyz(r,x0,z0,AtNo,l)
                    case( 9 ) 
                        SlaterOrbital = d_x2y2(r,x0,y0,AtNo,l)
                end select

                TotalPsiBra = TotalPsiBra + bra(i) * SlaterOrbital
                TotalPsiKet = TotalPsiKet + ket(i) * SlaterOrbital

            end do  ! <== DOS
        end do  ! <== atoms

        Psi_2( ix+1 , iy+1 , iz+1 ) = cdabs( cdsqrt( TotalPsiBra * TotalPsiKet ) )

    END DO          ! <==  Z coord
    END DO          ! <==  Y coord
    !$OMP end task 
 END DO             ! <==  X coord
!$OMP end single 
!$OMP end parallel 

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 DO ix = 0 , n_xyz_steps(1)
 
    DO iy = 0 , n_xyz_steps(2)
 
        DO iz = 0 , n_xyz_steps(3) 

            write(4,112,advance='no') Psi_2( ix+1 , iy+1 , iz+1 ) 

            If( (mod(iz,6) == 5) ) write(4,'(a)',advance='yes') 

        END DO                                                       
        write(4,'(a)',advance='yes') 
    END DO                                                      
 END DO                                                     

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 close(4)

 deallocate(xyz , Psi_2)

 include 'formats.h'

 end subroutine Gaussian_Cube_Format_Cmplx
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

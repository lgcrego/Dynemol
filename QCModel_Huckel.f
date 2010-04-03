 module QCModel_Huckel

    use type_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Semi_Empirical_Parms
    use Structure_Builder
    use Overlap_Builder
    use dipole_potential_m        , only : DP_phi

    implicit real*8      (a-h,o-y)
    implicit complex*16  (z)

 contains
!
!
!
!==================================================
 subroutine EigenSystem( system, basis, QM , flag )
!==================================================
 type(structure)               , intent(in)    :: system
 type(STO_basis)               , intent(in)    :: basis(:)
 type(C_eigen)                 , intent(out)   :: QM
 integer          , optional   , intent(inout) :: flag          

! . local variables 
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) 
 real*8  , ALLOCATABLE :: h(:,:) , dumb_s(:,:) , S_matrix(:,:)
! . local parameters
 real*8  , parameter   :: one = 1.d0 , zero = 0.d0
 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 CALL Overlap_Matrix(system,basis,S_matrix)

 ALLOCATE(QM%erg(size(basis))) 

 ALLOCATE(h(size(basis),size(basis)),dumb_s(size(basis),size(basis)))

! clone S_matrix because SYGVD will destroy it ... 
 dumb_s = S_matrix
 
 do j = 1 , size(basis)
   do i = 1 , j
     
     h(i,j) = huckel(i,j,S_matrix(i,j),basis)

   end do
 end do  

 CALL SYGVD(h,dumb_s,QM%erg,1,'V','U',info)

 If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '
 If ( present(flag) ) flag = info

 DEALLOCATE(dumb_s)

!     ---------------------------------------------------
!   ROTATES THE HAMILTONIAN:  H --> H*S_inv 
!
!   RIGHT EIGENVECTOR ALSO CHANGE: |C> --> S.|C> 
!
!   Rv = <AO|MO> coefficients
!     ---------------------------------------------------

 ALLOCATE(Lv(size(basis),size(basis)))

    Lv = h

 DEALLOCATE(h)

 ALLOCATE(Rv(size(basis),size(basis)))

    CALL gemm(S_matrix,Lv,Rv,'N','N',one,zero)  

 DEALLOCATE( S_matrix )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

 ALLOCATE(QM%L(size(basis),size(basis))) 
! eigenvectors in the rows of QM%L
    QM%L = transpose(Lv)                 
 DEALLOCATE( Lv )

 ALLOCATE(QM%R(size(basis),size(basis)))
! eigenvectors in the columns of QM%R
    QM%R = Rv             
 DEALLOCATE( Rv )

!  the order of storage is the ascending order of eigenvalues
!----------------------------------------------------------
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

! save energies of the TOTAL system 
 OPEN(unit=9,file='system-ergs.dat',status='unknown')
    do i = 1 , size(basis)
        write(9,*) i , QM%erg(i)
    end do
 CLOSE(9)   

 If( verbose ) Print*, '>> EigenSystem done <<'

 end subroutine EigenSystem
!
!
!
!====================================
 pure function Huckel(i,j,S_ij,basis)
!====================================
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: k_eff , k_WH , Huckel

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / 2.d0

 k_eff = k_WH + c3 + c3 * c3 * (1.d0 - k_WH)

 huckel = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / 2.d0

 IF(i == j) huckel = basis(i)%IP

 huckel = huckel + S_ij*DP_phi(i,j,basis)
   
 end function Huckel
!
!
 end module QCModel_Huckel

 module QCModel_Huckel

    use type_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use Overlap_Builder             , only : Overlap_Matrix
    use dipole_potential_m          , only : DP_phi

 contains
!
!
!
!=============================================================
 subroutine EigenSystem( system , basis , QM , flag1 , flag2 )
!=============================================================
 implicit none
 type(structure)               , intent(in)    :: system
 type(STO_basis)               , intent(in)    :: basis(:)
 type(C_eigen)                 , intent(inout) :: QM
 integer          , optional   , intent(inout) :: flag1          
 integer          , optional   , intent(in)    :: flag2          

! local variables ...
 real*8  , ALLOCATABLE :: Lv(:,:) , Rv(:,:) 
 real*8  , ALLOCATABLE :: h(:,:) , dumb_s(:,:) , S_matrix(:,:)
 integer               :: i , j , info

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 CALL Overlap_Matrix(system,basis,S_matrix)

 If( .NOT. allocated(QM%erg) ) ALLOCATE(QM%erg(size(basis))) 

 ALLOCATE(h(size(basis),size(basis)),dumb_s(size(basis),size(basis)))

! clone S_matrix because SYGVD will destroy it ... 
 dumb_s = S_matrix

 If( DP_field_ ) then
 
    do j = 1 , size(basis)
        do i = 1 , j
     
            h(i,j) = huckel_with_FIELDS(i,j,S_matrix(i,j),basis)

        end do
    end do  

 else

    do j = 1 , size(basis)
        do i = 1 , j
     
            h(i,j) = huckel(i,j,S_matrix(i,j),basis)

        end do
    end do  

 end If

 CALL SYGVD(h,dumb_s,QM%erg,1,'V','U',info)

 If ( info /= 0 ) write(*,*) 'info = ',info,' in SYGVD in EigenSystem '
 If ( present(flag1) ) flag1 = info

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

 ! garantees continuity between basis:  Lv(old)  and  Lv(new) ...
 If( (driver == "eigen_slice") .AND. (flag2 > 1) ) CALL phase_locking( Lv , QM%R )

 ALLOCATE(Rv(size(basis),size(basis)))

 CALL gemm(S_matrix,Lv,Rv,'N','N',D_one,D_zero)  

 DEALLOCATE( S_matrix )

!----------------------------------------------------------
!  normalizes the L&R eigenvectors as < L(i) | R(i) > = 1

 If( .NOT. allocated(QM%L) ) ALLOCATE(QM%L(size(basis),size(basis))) 
! eigenvectors in the rows of QM%L
 QM%L = transpose(Lv)                 
 DEALLOCATE( Lv )

 If( .NOT. ALLOCATED(QM%R) ) ALLOCATE(QM%R(size(basis),size(basis)))
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
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel
 real*8  :: k_eff , k_WH , c1 , c2 , c3

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / 2.d0

 k_eff = k_WH + c3 + c3 * c3 * (1.d0 - k_WH)

 huckel = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / 2.d0

 IF(i == j) huckel = basis(i)%IP

 end function Huckel
!
!
!================================================
 pure function Huckel_with_FIELDS(i,j,S_ij,basis)
!================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: Huckel_with_FIELDS
 real*8  :: k_eff , k_WH , c1 , c2 , c3

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 c1 = basis(i)%IP - basis(j)%IP
 c2 = basis(i)%IP + basis(j)%IP

 c3 = (c1/c2)*(c1/c2)

 k_WH = (basis(i)%k_WH + basis(j)%k_WH) / 2.d0

 k_eff = k_WH + c3 + c3 * c3 * (1.d0 - k_WH)

 huckel_with_FIELDS = k_eff * S_ij * (basis(i)%IP + basis(j)%IP) / 2.d0

 IF(i == j) huckel_with_FIELDS = basis(i)%IP

 huckel_with_FIELDS = huckel_with_FIELDS + S_ij*DP_phi(i,j,basis)
   
 end function Huckel_with_FIELDS
!
!
!
!===================================
 subroutine phase_locking( Lv , CR )
!===================================
implicit none
real*8      , intent(inout) :: Lv(:,:)
complex*16  , intent(in)    :: CR(:,:)

! local parameters ...
integer , parameter    :: n_check = 5  ! <== has to be odd !

! local variables ...
integer                :: J , K , N , indx
real*8                 :: big , signal
real*8                 :: norm(n_check)
real*8  , allocatable  :: temp(:,:) , old_Rv(:,:) , MO_ovlp(:,:)

N = size( CR(:,1) )

allocate( temp    (N,N) )
allocate( old_Rv  (N,N) )
allocate( MO_ovlp (N,N) )

old_Rv = real(CR)
         
CALL gemm( Lv , old_Rv , MO_ovlp , 'T' , 'N' , D_one , D_zero )

DO J = 1 , N
      
    norm = D_zero

    DO K = 1 , n_check
         
        indx = J - (n_check/2) + K - 1
      
        IF( indx <= 0 ) indx = 1
        IF( indx >= N ) indx = N
             
        norm(K) = MO_ovlp(J,indx)
             
    END DO   

    big  = abs(norm(1))  
    indx = J - n_check / 2

    DO K = 2 , n_check
        IF( abs(norm(K)) > big ) then
            big  =  abs(norm(K))
            indx = J - (n_check/2) + K - 1
        END IF   
    END DO   

    IF( indx <= 0 ) indx = 1
    IF( indx >= N ) indx = N

    signal = sign( 1.d0 , MO_ovlp(J,indx) )
      
    temp(:,indx) = signal * Lv(:,indx)

END DO

Lv = temp

deallocate( temp , old_Rv , MO_ovlp )

end subroutine phase_locking
!
!
!
end module QCModel_Huckel

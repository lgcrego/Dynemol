 module Hamiltonians

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use omp_lib
    use constants_m
    use parameters_m          , only : EnvField_ , Induced_ , Environ_type , Environ_step , B_ext
    use Dielectric_Potential  , only : Q_phi
    use DP_potential_m        , only : DP_phi
    use DP_main_m             , only : DP_matrix_AO
    use polarizability_m      , only : Induced_DP_phi
    use Semi_Empirical_Parms  , only : atom

    public :: X_ij , even_more_extended_Huckel , Huckel_with_FIELDS , spin_orbit_h , MF_interaction

    private

    ! module variables ...
    real*8  , allocatable :: DP_4_matrix(:,:,:)
    integer , allocatable :: indx(:)
    logical               :: done = .false.

 contains
!
!
!
!===================================
 pure function X_ij( i , j , basis )
!===================================
 implicit none
 integer         , intent(in) :: i , j
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 real*8  :: X_ij
 real*8  :: k_eff , k_WH , c1 , c2 , c3 , c4

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN
 
 if (i == j) then
    X_ij = basis(i)%IP + basis(i)%V_shift
 else
    c1 = basis(i)%IP - basis(j)%IP
    c2 = basis(i)%IP + basis(j)%IP

    c3 = (c1/c2)*(c1/c2)

    c4 = (basis(i)%V_shift + basis(j)%V_shift) * HALF

    k_WH = (basis(i)%k_WH + basis(j)%k_WH) * HALF

    k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

    X_ij = (k_eff*c2*HALF + c4) 
 endif

 end function X_ij
!
!
!
!
!==============================================================================
 function even_more_extended_huckel( system , basis , S_matrix , it ) result(h)
!==============================================================================
implicit none
type(structure)            , intent(in) :: system
type(STO_basis)            , intent(in) :: basis(:)
real*8                     , intent(in) :: S_matrix(:,:)
integer         , optional , intent(in) :: it

! local variables ...
integer               :: i , j , ia , ib , ja , jb , K , L , N
real*8                :: Rab , DP_4_vector(4)
real*8  , ALLOCATABLE :: h(:,:) 
logical               :: evaluate

! instantiating DP_4_matrix ...
if( .not. done ) CALL allocate_DP_4_matrix

! evaluate or not evaluate DP_phi this time...
If( .not. present(it) ) then
   evaluate = .true.
else If( mod(it-1,Environ_step) == 0 ) then
   evaluate = .true.
else
   evaluate = .false.
end if

! resetting DP_4_matrix before fresh calculation ...
if( evaluate ) DP_4_matrix = D_zero

N = size(basis)
Allocate( h(N,N) , source = D_zero )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!$OMP parallel do &
!$OMP   default(shared) &
!$OMP   schedule(dynamic, 1) &
!$OMP   private(ib, ia, Rab, jb, ja, j, i, K, L, DP_4_vector)
do ib = 1, system%atoms
    do ia = ib+1, system%atoms

        if ((system%QMMM(ib) /= "QM") .OR. (system%QMMM(ia) /= "QM")) then
            cycle
        end if

        Rab = GET_RAB(system%coord(ib,:), system%coord(ia,:))
        if (Rab > cutoff_Angs) then
           cycle
        end if

        K = indx(ia)
        L = indx(ib)

        If( evaluate ) then 
           select case (Environ_Type)
                case('DP_MM','DP_QM')
                    DP_4_vector = DP_phi( system , ia , ib )
                case default
                    DP_4_vector =  Q_phi( system , ia , ib )
           end select
           DP_4_matrix(K,L,:) = DP_4_vector
        else 
           DP_4_vector = DP_4_matrix(K,L,:)
        end if

        do jb = 1, atom(system%AtNo(ib))% DOS
            do ja = 1, atom(system%AtNo(ia))% DOS

                j = system% BasisPointer(ib) + jb
                i = system% BasisPointer(ia) + ja

                h(i,j) = huckel_with_FIELDS(i , j , S_matrix(i,j) , basis , DP_4_vector )

                h(j,i) = h(i,j)

            end do
        end do

    end do
end do  
!$OMP END PARALLEL DO

forall( i=1:N ) h(i,i) = X_ij( i , i , basis ) 

end function even_more_extended_huckel
!
!
!
!=====================================================================
 pure function Huckel_with_FIELDS( i , j , S_ij , basis ,DP_4_vector )
!=====================================================================
 implicit none
 integer         , intent(in) :: i , j
 real*8          , intent(in) :: S_ij
 type(STO_basis) , intent(in) :: basis(:)
 real*8          , intent(in) :: DP_4_vector(:)

! local variables ... 
 real*8  :: DP(4)
 real*8  :: r0(3) , Ri(3) , vector(3)
 real*8  :: Huckel_with_FIELDS
 real*8  :: k_eff , k_WH , c1 , c2 , c3 , c4

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

 DP = D_zero

 IF( i == j ) then

      huckel_with_FIELDS = basis(i)%IP + basis(i)%V_shift

 else

      c1 = basis(i)%IP - basis(j)%IP
      c2 = basis(i)%IP + basis(j)%IP
     
      c3 = (c1/c2)*(c1/c2)
     
      c4 = (basis(i)%V_shift + basis(j)%V_shift)*HALF
     
      k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two
     
      k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)
     
      huckel_with_FIELDS = (k_eff*c2/two + c4)*S_ij
     
      If( abs(S_ij) > low_prec ) then
    
          r0(1) = ( basis(i)%x + basis(j)%x ) / two
          r0(2) = ( basis(i)%y + basis(j)%y ) / two
          r0(3) = ( basis(i)%z + basis(j)%z ) / two
          
          Ri(1) = basis(i)%x
          Ri(2) = basis(i)%y
          Ri(3) = basis(i)%z
          
          vector = DP_matrix_AO(i,j,:) - ( r0 - Ri )*S_ij  ! <== in Angs
          
          If( Induced_  ) DP = Induced_DP_phi( i , j , basis )
          
          If( EnvField_ ) DP = DP + DP_4_vector
          
          huckel_with_FIELDS = huckel_with_FIELDS - ( S_ij*DP(1) + dot_product(vector(1:3),DP(2:4)) )

      end If
      
 end if

end function Huckel_with_FIELDS
!
!
!
!===============================================
 subroutine spin_orbit_h( basis , h , S_matrix )
!===============================================
type(STO_basis)               , intent(in)    :: basis(:)
real*8          , allocatable , intent(inout) :: h(:,:)
real*8                        , intent(in)    :: S_matrix(:,:)

! local variables ...
real*8  :: sum1 , sum2 , eps
integer :: i , j , k , c , l , l1 , l2 , n

!----------------------------------------------------------
!     building the coupling spin-orbit HAMILTONIAN

n = size( basis )

allocate( h ( n , n ) , source = 0.0d0 )

! orbitals i and j in the same sitio ...
do i = 1 , n

    do j = i , n

        if( basis( i ) % atom == basis( j ) % atom .AND. basis( i ) % l == basis( j ) % l ) then

            if( basis( i ) % s == basis( j ) % s .AND. basis( i ) % m == basis( j ) % m ) then
                
                CALL SOC_constant( basis( i ) , eps )
                    
                h( i , i ) = HALF * eps * dfloat( basis( i ) % m ) * dfloat( basis( i ) % s )

            elseif( basis( i ) % s /= basis( j ) % s .AND. basis( j ) % m == basis( i ) % m + 1 ) then
                
                CALL SOC_constant( basis( i ) , eps )

                k = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m + 1 )
                
                h( j , i ) = HALF * eps * dsqrt( dfloat( k ) )
                h( i , j ) = h( j , i ) ! It's unecessary, but is easiest compute the terms with orbitals in different sitios

            end if

        end if

    end do

end do

! orbitals i and j in different sitios ...
do i = 1 , n

    do j = i + 1 , n

        if( basis( j ) % atom /= basis( i ) % atom ) then

            ! first sum ...
            if( basis( j ) % s == 1 ) then

                l1 = 1
                l2 = n / 2

            elseif( basis( j ) % s == - 1 ) then

                l1 = n / 2 + 1
                l2 = n

            end if

            sum1 = 0.0d0
            do k = l1 , l2

                if( basis( k ) % atom == basis( i ) % atom .AND. basis( k ) % l == basis( i ) % l ) then

                    select case( basis( i ) % l )

                        case( 0 )

                            sum1 = S_matrix( k , j ) * h( k , i )

                        case( 1 )

                            sum1 = dot_product( S_matrix( k : k + 2 , j ) , h( k : k + 2 , i ) )

                        case( 2 )

                            sum1 = dot_product( S_matrix( k : k + 4 , j ) , h( k : k + 4 , i ) )

                    end select

                    exit

                end if

            end do

            ! second sum ...
            if( basis( i ) % s == 1 ) then

                l1 = 1
                l2 = n / 2

            elseif( basis( i ) % s == - 1 ) then

                l1 = n / 2 + 1
                l2 = n

            end if

            sum2 = 0.0d0
            do k = l1 , l2

                if( basis( k ) % atom == basis( j ) % atom .AND. basis( k ) % l == basis( j ) % l ) then

                    select case( basis( j ) % l )

                        case( 0 )

                            sum2 = S_matrix( k , i ) * h( k , j )

                        case( 1 )

                            sum2 = dot_product( S_matrix( k : k + 2 , i ) , h( k : k + 2 , j ) )

                        case( 2 )

                            sum2 = dot_product( S_matrix( k : k + 4 , i ) , h( k : k + 4 , j ) )

                    end select

                    exit

                end if

            end do

            h( j , i ) = sum1 + sum2
!            h( j , i ) = 0.4d0 * h( j , i )

        end if

    end do

end do

end subroutine spin_orbit_h
!
!
!
!=================================================
 subroutine MF_interaction( basis , h , S_matrix )
!=================================================
type(STO_basis)               , intent(in)    :: basis(:)
real*8          , allocatable , intent(inout) :: h(:,:)
real*8                        , intent(in)    :: S_matrix(:,:)

! local variables ...
real*8  :: a , tx1 , tx2 , tx , tz
integer :: i , j , k , t , n

if( B_ext( 2 ) /= D_zero ) stop ">> Error: Magnetic field interaction has not been implemented for By/=0 <<"

n = size( basis )
k = n / 2

allocate( h ( n , n ) , source = 0.0d0 )

do i = 1 , n

    do j = i , n

        tx = D_zero
        if( B_ext( 1 ) /= D_zero ) then

            tx1 = D_zero
            tx2 = D_zero

            if( i < n ) then
                t = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m + 1 )
                a = dsqrt( dfloat( t ) )
                tx1 = a * S_matrix( j , i + 1 )
            end if

            if( i > 1 ) then
                t = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m - 1 )
                a = dsqrt( dfloat( t ) )
                tx2 = a * S_matrix( j , i - 1 )
            end if

            t  = basis( i ) % s
            tx = tx1 + tx2 + g_S * S_matrix( j , i + t * k )

        end if

        tz = D_zero
        if( B_ext( 3 ) /= D_zero ) tz = ( 2.0d0 * dfloat( basis( i ) % m ) + g_S * dfloat( basis( i ) % s ) ) * S_matrix( j , i )

        h( j , i ) = HALF * mu_B * ( tx * B_ext( 1 ) + tz * B_ext( 3 ) )

    end do

end do

end subroutine MF_interaction
!
!
!
!======================================
 subroutine SOC_constant( basis , eps )
!======================================
implicit none
type(STO_basis) , intent(in)    :: basis
real*8          , intent(inout) :: eps

! Ref. 1: Koseki, S; et. al; J. Phys. Chem. A, 2019, 123, 2325−2339
!         D'Alessando, D. N; et. al; Inorg. Chem. 2006, 45, 3261−3274
! Ref. 2: Martin, Q. C; JOURNAL OF RESEARCH of toe Notional Burea u of Standards - A. Physics and Chemistry, 1971, 75A, 109-111
!         Equation 1: eps = 0.45 * Z^(2.33) * n_eff ^(-3) , in cm^(-1)
! Ref. 3: Manne, R; et. al; Molecular Physics, 1975, 29, 485−500
! Ref. 4: Wittel, K. e Manne, R; Theoret. Chim. Acta (Berl.), 1974, 33, 347-349

select case( basis % AtNo )

    case( 1 )

        if( basis % l == 0 ) eps = 0.0d0

    case( 6 )

        if( basis % l == 0 ) eps = 0.0d0
        if( basis % l == 1 ) eps = 453.50d-6 ! Ref. 2

    case( 7 )

        if( basis % l == 0 ) eps = 0.0d0
        if( basis % l == 1 ) eps = 649.48d-6 ! Ref. 2

    case( 35 )

!        eps = - 0.35d0
        eps = 0.35d0 ! Ref. 3

    case( 44 )

        if( basis % l == 0 ) eps = 0.0d0
        if( basis % l == 1 ) eps = 13.95d-3  ! Ref. 2 <== Revise it...
        if( basis % l == 2 ) eps = 0.12398d0 ! Ref. 1

    case( 53 )

!        eps = - 0.73d0 ! Ref. 3
        eps = 0.73d0 ! Ref. 3

    case default

        print'("Problem with the spin-orbit coupling constants: Hamiltonians.f --> SCO_constant ")'

end select

end subroutine SOC_constant
!
!
!
!==================================================
pure function GET_RAB(a_coord, b_coord) result(rab)
!==================================================
 implicit none

 ! args
 real*8, intent(in) :: a_coord(:)
 real*8, intent(in) :: b_coord(:)

 ! result
 real*8 :: rab

 rab = SUM((a_coord - b_coord) ** 2)
 rab = SQRT(rab)
end function GET_RAB
!
!
!
!===============================
 subroutine allocate_DP_4_matrix
!===============================
 use Structure_Builder , only : a => Extended_Cell 
 implicit none

! local variables ...
 integer :: N_of_QM , i , j

 N_of_QM = count(a%QMMM == "QM")

 allocate( indx(a% atoms) , source = 0 )

 ! indx(:) is used to access DP_4_matrix data ...
 j = 0
 do i = 1 , a%atoms
    if( a% QMMM(i) == "QM" ) then
       j = j + 1
       indx(i) = j
    end if
 end do

 allocate( DP_4_matrix( N_of_QM , N_of_QM , 4 ) , source = D_zero ) 

 done = .true.

end subroutine allocate_DP_4_matrix
!
!
!
end module Hamiltonians


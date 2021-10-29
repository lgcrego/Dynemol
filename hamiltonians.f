 module Hamiltonians

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use omp_lib
    use constants_m
    use parameters_m          , only : EnvField_ , Induced_ , Environ_type , Environ_step , SOC , B_field
    use Dielectric_Potential  , only : Q_phi
    use DP_potential_m        , only : DP_phi
    use DP_main_m             , only : DP_matrix_AO
    use polarizability_m      , only : Induced_DP_phi
    use Semi_Empirical_Parms  , only : atom

    public :: zX_ij , X_ij , even_more_extended_Huckel , Huckel_with_FIELDS , spin_orbit_h

    private

    ! module variables ...
    real*8  , allocatable :: DP_4_matrix(:,:,:)
    integer , allocatable :: indx(:)
    logical               :: done = .false.

 contains
!
!
!
!========================================
 pure function zX_ij( i , j , basis , N )
!========================================
 implicit none
 integer         , intent(in) :: i , j , N
 type(STO_basis) , intent(in) :: basis(:)

! local variables ... 
 complex*16 :: zX_ij
 complex*16 :: zk_eff , zk_WH , zc1 , zc2 , zc3 , zc4 , IP_i , IP_j
 real*8     :: k_WH , c4

!----------------------------------------------------------
!      building  the  HUCKEL  HAMILTONIAN

if( ( i >= N/2 - 251 .AND. i <= N/2 ) .OR. ( i >= N - 251 .AND. i <= N ) ) then
    IP_i = basis(i) % IP * zi
else
    IP_i = dcmplx(basis(i) % IP , D_zero)
end if
if( ( j >= N/2 - 251 .AND. j <= N/2 ) .OR. ( j >= N - 251 .AND. j <= N ) ) then
    IP_j = basis(j) % IP * zi
else
    IP_j = dcmplx(basis(j) % IP , D_zero)
end if

 if (i == j) then
    zX_ij = IP_i + dcmplx(basis(i)%V_shift,D_zero)
 else
    zc1 = IP_i - IP_j
    zc2 = IP_i + IP_j

    zc3 = (zc1/zc2)*(zc1/zc2)

    c4 = (basis(i)%V_shift + basis(j)%V_shift) * HALF
    zc4 = dcmplx(c4,D_zero)

    k_WH = (basis(i)%k_WH + basis(j)%k_WH) * HALF
    zk_WH = dcmplx(k_WH,D_zero)

    zk_eff = zk_WH + zc3 + zc3 * zc3 * (C_one - zk_WH)

    zX_ij = (zk_eff*zc2*HALF + zc4)
 endif

 end function zX_ij
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
implicit none
type(STO_basis)               , intent(in)    :: basis(:)
complex*16      , allocatable , intent(inout) :: h(:,:)
real*8                        , intent(in)    :: S_matrix(:,:)

! local variables ...
complex*16 , allocatable :: h_sup(:,:)
complex*16 :: sum1 , sum2 , ty1 , ty2 , ty
real*8     :: a , b , eps , tx1 , tx2 , tx , tz , dx , dy , dz , dist_ij , dist_it , dist_jt , dist_maxz , dist_minz 
integer    :: i , j , k , l , l1 , l2 , r , s , c , n , t , atom , atomi , atomj , alpha , beta

!----------------------------------------------------------
!     building the coupling spin-orbit HAMILTONIAN

n = size( basis )

allocate( h_sup ( n , n ) , source = C_zero )

! orbitals i and j in the same sitio ...
atom = 0
j    = n / 2
do i = 1 , j

    if( basis( i ) % l /= 0 .AND. basis( i ) % atom /= atom ) then

        atom = basis( i ) % atom

        CALL SOC_constant( basis( i ) % AtNo , 1 , eps )

        a = HALF * eps

        ! <i,up|SL|j,up> terms, with i,j={p-orbitals}...
        h_sup( i + 2 , i ) = - a * zi
        h_sup( i , i + 2 ) = - h_sup( i + 2 , i )

        ! <i,up|SL|j,down> and <i,down|SL|j,up> terms, with i,j={p-orbitals}...
        h_sup( i + j + 1 , i ) =   a * zi
        h_sup( i , i + j + 1 ) = - h_sup( i + j + 1 , i )
        h_sup( i + j , i + 1 ) = - a * zi
        h_sup( i + 1 , i + j ) = - h_sup( i + j , i + 1 )
        h_sup( i + j + 2 , i + 1 ) = - a
        h_sup( i + 1 , i + j + 2 ) =   h_sup( i + j + 2 , i + 1 )
        h_sup( i + j + 1 , i + 2 ) =   a
        h_sup( i + 2 , i + j + 1 ) =   h_sup( i + j + 1 , i + 2 )

        ! <i,down|SL|j,down> terms, with i,j={p-orbitals}...
        h_sup( i + j + 2 , i + j ) =   a * zi
        h_sup( i + j , i + j + 2 ) = - h_sup( i + j + 2 , i + j )

        if( basis( i + 3 ) % l == 2 ) then

            CALL SOC_constant( basis( i ) % AtNo , 2 , eps )

            a = HALF * eps
            b = dsqrt( THREE )

            ! <i,up|SL|j,up> terms, with i,j={d-orbitals}...
            h_sup( i + 7 , i + 3 ) = - TWO * a * zi
            h_sup( i + 3 , i + 7 ) = - h_sup( i + 7 , i + 3 )
            h_sup( i + 6 , i + 4 ) = - a * zi
            h_sup( i + 4 , i + 6 ) = - h_sup( i + 6 , i + 4 )

            ! <i,up|SL|j,down> and <i,down|SL|j,up> terms, with i,j={d-orbitals}...
            h_sup( i + j + 4 , i + 3 ) =   a
            h_sup( i + 3 , i + j + 4 ) =   h_sup( i + j + 4 , i + 3 )
            h_sup( i + j + 6 , i + 3 ) =   a * zi
            h_sup( i + 3 , i + j + 6 ) = - h_sup( i + j + 6 , i + 3 )
            h_sup( i + j + 3 , i + 4 ) = - a
            h_sup( i + 4 , i + j + 3 ) =   h_sup( i + j + 3 , i + 4 )
            h_sup( i + j + 5 , i + 4 ) =   a * b * zi
            h_sup( i + 4 , i + j + 5 ) = - h_sup( i + j + 5 , i + 4 )
            h_sup( i + j + 7 , i + 4 ) =   a * zi
            h_sup( i + 4 , i + j + 7 ) = - h_sup( i + j + 7 , i + 4 )
            h_sup( i + j + 4 , i + 5 ) = - a * b * zi
            h_sup( i + 5 , i + j + 4 ) = - h_sup( i + j + 4 , i + 5 )
            h_sup( i + j + 6 , i + 5 ) = - a * b
            h_sup( i + 5 , i + j + 6 ) =   h_sup( i + j + 6 , i + 5 )
            h_sup( i + j + 3 , i + 6 ) = - a * zi
            h_sup( i + 6 , i + j + 3 ) = - h_sup( i + j + 3 , i + 6 )
            h_sup( i + j + 5 , i + 6 ) =   a * b
            h_sup( i + 6 , i + j + 5 ) =   h_sup( i + j + 5 , i + 6 )
            h_sup( i + j + 7 , i + 6 ) = - a
            h_sup( i + 6 , i + j + 7 ) =   h_sup( i + j + 7 , i + 6 )
            h_sup( i + j + 4 , i + 7 ) = - a * zi
            h_sup( i + 7 , i + j + 4 ) = - h_sup( i + j + 4 , i + 7 )
            h_sup( i + j + 6 , i + 7 ) =   a
            h_sup( i + 7 , i + j + 6 ) =   h_sup( i + j + 6 , i + 7 )

            ! <i,down|SL|j,down> terms, with i,j={d-orbitals}...
            h_sup( i + j + 7 , i + j + 3 ) =   TWO * a * zi
            h_sup( i + j + 3 , i + j + 7 ) = - h_sup( i + j + 7 , i + j + 3 )
            h_sup( i + j + 6 , i + j + 4 ) =   a * zi
            h_sup( i + j + 4 , i + j + 6 ) = - h_sup( i + j + 6 , i + j + 4 )

        end if
    end if
end do

allocate( h ( n , n ) , source = C_zero )

beta = 0
do i = 1 , N/2

    if( basis( i ) % atom /= beta ) then

        beta = basis( i ) % atom

        alpha = 0
        do j = 1 , N / 2

            if( basis( j ) % atom /= alpha ) then

                alpha = basis( j ) % atom

                do k = j , j + count( basis % atom == alpha ) / 2 - 1
                    do l = j , j + count( basis % atom == alpha ) / 2 - 1

                        do r = 0 , 3
                            do s = r , 3

                                ! <i,up|SL|j,up> terms, with i,j={p-orbitals}...
                                h(i+s,i+r) = h(i+s,i+r) + S_matrix(i+s,k    ) * S_matrix(l    ,i+r) * h_sup(k    ,l    ) + &
                                                          S_matrix(i+s,k+N/2) * S_matrix(l    ,i+r) * h_sup(k+N/2,l    ) + &
                                                          S_matrix(i+s,k    ) * S_matrix(l+N/2,i+r) * h_sup(k    ,l+N/2) + &
                                                          S_matrix(i+s,k+N/2) * S_matrix(l+N/2,i+r) * h_sup(k+N/2,l+N/2)

                                ! <i,down|SL|j,down> terms, with i,j={p-orbitals}...
                                h(i+s+N/2,i+r+N/2) = h(i+s+N/2,i+r+N/2) + &
                                                     S_matrix(i+s+N/2,k    ) * S_matrix(l    ,i+r+N/2) * h_sup(k    ,l    ) + &
                                                     S_matrix(i+s+N/2,k+N/2) * S_matrix(l    ,i+r+N/2) * h_sup(k+N/2,l    ) + &
                                                     S_matrix(i+s+N/2,k    ) * S_matrix(l+N/2,i+r+N/2) * h_sup(k    ,l+N/2) + &
                                                     S_matrix(i+s+N/2,k+N/2) * S_matrix(l+N/2,i+r+N/2) * h_sup(k+N/2,l+N/2)

                            end do
                        end do

                        do r = 0 , 3
                            do s = 0 , 3

                                ! <i,down|SL|j,up> terms, with i,j={p-orbitals}...
                                h(i+s+N/2,i+r) = h(i+s+N/2,i+r) + S_matrix(i+s+N/2,k    ) * S_matrix(l    ,i+r) * h_sup(k    ,l    ) + &
                                                                  S_matrix(i+s+N/2,k+N/2) * S_matrix(l    ,i+r) * h_sup(k+N/2,l    ) + &
                                                                  S_matrix(i+s+N/2,k    ) * S_matrix(l+N/2,i+r) * h_sup(k    ,l+N/2) + &
                                                                  S_matrix(i+s+N/2,k+N/2) * S_matrix(l+N/2,i+r) * h_sup(k+N/2,l+N/2)

                            end do
                        end do
                    end do
                end do
            end if
        end do

        do r = 0 , 2
            do s = r + 1 , 3
                ! <i,up|SL|j,up> terms, with i,j={p-orbitals}...
                h(i+r,i+s) = dconjg(h(i+s,i+r))
                ! <i,down|SL|j,down> terms, with i,j={p-orbitals}...
                h(i+r+N/2,i+s+N/2) = dconjg(h(i+s+N/2,i+r+N/2))
            end do
        end do

        ! <i,down|SL|j,up> terms, with i,j={p-orbitals}...
        do r = 0 , 3
            do s = 0 , 3
                h(i+r,i+s+N/2) = dconjg(h(i+s+N/2,i+r))
            end do
        end do

    end if
end do

! orbitals i and j in different sitios ...
do i = 1 , n
    do j = i + 1 , n

        if( basis(j)% atom /= basis(i)% atom ) then

            ! first sum ...
            if( basis( i ) % s == 1 ) then

                l1 = 1
                l2 = n / 2

            elseif( basis( i ) % s == - 1 ) then

                l1 = n / 2 + 1
                l2 = n

            end if

            sum1 = C_zero
            do k = l1 , l2

                if( basis( j ) % l /= 0 .AND. basis( k ) % atom == basis( j ) % atom .AND. basis( k ) % l == basis( j ) % l ) then

                    select case( basis( j ) % l )

                        case( 1 )

                            sum1 = h_sup( j , k     ) * S_matrix( k     , i ) + h_sup( j , k + 1 ) * S_matrix( k + 1 , i ) + &
                                   h_sup( j , k + 2 ) * S_matrix( k + 2 , i )

                        case( 2 )

                            sum1 = h_sup( j , k     ) * S_matrix( k     , i ) + h_sup( j , k + 1 ) * S_matrix( k + 1 , i ) + &
                                   h_sup( j , k + 2 ) * S_matrix( k + 2 , i ) + h_sup( j , k + 3 ) * S_matrix( k + 3 , i ) + &
                                   h_sup( j , k + 4 ) * S_matrix( k + 4 , i )

                    end select

                    exit

                end if

            end do

            ! second sum ...
            if( basis( j ) % s == 1 ) then

                l1 = 1
                l2 = n / 2

            elseif( basis( j ) % s == - 1 ) then

                l1 = n / 2 + 1
                l2 = n

            end if

            sum2 = C_zero
            do k = l1 , l2

                if( basis( i ) % l /= 0 .AND. basis( k ) % atom == basis( i ) % atom .AND. basis( k ) % l == basis( i ) % l ) then

                    select case( basis( i ) % l )

                        case( 1 )

                            sum2 = S_matrix( j , k     ) * h_sup( k     , i ) + S_matrix( j , k + 1 ) * h_sup( k + 1 , i ) + &
                                   S_matrix( j , k + 2 ) * h_sup( k + 2 , i )

                        case( 2 )

                            sum2 = S_matrix( j , k     ) * h_sup( k     , i ) + S_matrix( j , k + 1 ) * h_sup( k + 1 , i ) + &
                                   S_matrix( j , k + 2 ) * h_sup( k + 2 , i ) + S_matrix( j , k + 3 ) * h_sup( k + 3 , i ) + &
                                   S_matrix( j , k + 4 ) * h_sup( k + 4 , i )

                    end select

                    exit

                end if

            end do

            h( j , i ) = sum1 + sum2

        end if

    end do

end do

deallocate( h_sup )

!if( B_field( 1 ) /= D_zero .OR. B_field( 2 ) /= D_zero .OR. B_field( 3 ) /= D_zero ) then
!
!    k = n / 2
!
!    do i = 1 , n
!
!        do j = i , n
!
!            tx = D_zero
!            if( B_field( 1 ) /= D_zero ) then
!
!                tx1 = D_zero
!                tx2 = D_zero
!
!                if( i < n ) then
!
!                    t   = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m + 1 )
!                    a   = dsqrt( dfloat( t ) )
!                    tx1 = a * S_matrix( j , i + 1 )
!
!                end if
!
!                if( i > 1 ) then
!
!                    t   = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m - 1 )
!                    a   = dsqrt( dfloat( t ) )
!                    tx2 = a * S_matrix( j , i - 1 )
!                    
!                end if
!
!                t  = basis( i ) % s
!                tx = tx1 + tx2 + key * g_S * S_matrix( j , i + t * k )
!
!            end if
!
!            ty = D_zero
!            if( B_field( 2 ) /= D_zero ) then
!
!                ty1 = D_zero
!                ty2 = D_zero
!
!                if( i < n ) then
!
!                    t   = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m + 1 )
!                    a   = dsqrt( dfloat( t ) )
!                    ty1 = - a * S_matrix( j , i + 1 ) * zi
!
!                end if
!
!                if( i > 1 ) then
!
!                    t   = basis( i ) % l * ( basis( i ) % l + 1 ) - basis( i ) % m * ( basis( i ) % m - 1 )
!                    a   = dsqrt( dfloat( t ) )
!                    ty2 = - a * S_matrix( j , i - 1 ) * zi
!                    
!                end if
!
!                t  = basis( i ) % s
!                ty = ty1 - tx2 + key * g_S * S_matrix( j , i + t * k ) * zi
!
!            end if
!
!            tz = D_zero
!            if( B_field( 3 ) /= D_zero ) tz = ( two * dfloat( basis( i ) % m ) + g_S * dfloat( basis( i ) % s ) ) * S_matrix( j , i )
!
!            h( j , i ) = h( j , i ) + HALF * mu_B * ( tx * B_field( 1 ) + ty * B_field( 2 ) + tz * B_field( 3 ) )
!
!        end do
!
!    end do
!
!end if

end subroutine spin_orbit_h
!
!
!
!=========================================
 subroutine SOC_constant( AtNo , l , eps )
!=========================================
implicit none
integer , intent(in)    :: AtNo
integer , intent(in)    :: l
real*8  , intent(inout) :: eps

! Ref. 1: Koseki, S; et. al; J. Phys. Chem. A, 2019, 123, 2325−2339
!         D'Alessando, D. N; et. al; Inorg. Chem. 2006, 45, 3261−3274

! Ref. 2: Martin, Q. C; JOURNAL OF RESEARCH of toe Notional Burea u of
! Standards - A. Physics and Chemistry, 1971, 75A, 109-111 (Table of Spin-Orbit
! Energies for p-Electrons in Neutral Atomic (core) np Configurations)
!         Equation 1: eps = 0.45 * Z^(2.33) * n_eff ^(-3) , in cm^(-1)

! Ref. 3: Manne, R; et. al; Molecular Physics, 1975, 29, 485−500

! Ref. 4: Wittel, K. e Manne, R; Theoret. Chim. Acta (Berl.), 1974, 33, 347-349

! Ref. 5: Geyer, M; et. al; J. Phys. Chem. C, 2019, 123, 27230-27241

! Ref. 6: Cusachs, L. C. and. Aldrich, H. S; Chem.Phys. Let; 1971, 12, 197-198

! Ref. 7: Dulitz, K; et. al; Molecular Physics, 114:19, 2848-2856 (Spin-orbit
! coupling and rovibrational structure in the iododiacetylene radial cation by
! PFI-ZEKE photoelectron spectroscopy)

! Ref. 8: Küfner,S ;et.al; Phys. Rev. B, 2013, 87, 235307 (Structural and
! electronic properties of alpha-tin nanocrystals from first principles)

select case( AtNo )

    case( 4 )

!        if( l == 1 ) eps = 92.988d-6 ! Ref. 2
        if( l == 1 ) eps = 92.988d-3 ! teste

    case( 6 )

!        if( l == 1 ) eps = 6.0d-3    ! Ref. 5
        if( l == 1 ) eps = 6.0d0    ! teste
!        if( l == 1 ) eps = 454.0923d-6 ! Ref. 2
!        if( l == 1 ) eps = 4.541d0 ! teste

    case( 7 )

        if( l == 1 ) eps = 649.3675d-6 ! Ref. 2

    case( 15 )

        if( l == 1 ) eps = - 0.00889d0 ! Ref. 2 <-- Neff = 1.51 (Ref.6)

    case( 35 )

!        if( l == 1 ) eps = - 0.35d0
        if( l == 1 ) eps = 0.35d0 ! Ref. 3

    case( 44 )

        if( l == 1 ) eps = 5.8893d-3 ! Ref. 2
        if( l == 2 ) eps = 0.12398d0 ! Ref. 1

    case( 50 )

        if( l == 1 ) eps = 0.68d0 ! Ref. 8 

    case( 53 )

!        if( l == 1 ) eps = 0.73d0 ! Ref. 3
        if( l == 1 ) eps = - 0.62848d0 ! Ref. 7

    case default

        print'("Problem with the spin-orbit coupling constants: Hamiltonians.f --> SOC_constant ")'

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


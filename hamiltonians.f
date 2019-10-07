 module Hamiltonians

    use f95_precision
    use blas95
    use lapack95
    use type_m
    use omp_lib
    use constants_m
    use parameters_m          , only : EnvField_ , Induced_ , Environ_type , Environ_step
    use Dielectric_Potential  , only : Q_phi
    use DP_potential_m        , only : DP_phi
    use DP_main_m             , only : DP_matrix_AO
    use polarizability_m      , only : Induced_DP_phi
    use Semi_Empirical_Parms  , only : atom

    public :: X_ij , even_more_extended_Huckel , Huckel_with_FIELDS

    private

    ! module variables ...
    real*8 , allocatable :: DP_4_matrix(:,:,:)
    logical              :: done = .false.

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
integer               :: i , j , ia , ib , ja , jb , N
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

! resetting DP_$_matrix before fresh calculation ...
if( evaluate ) DP_4_matrix = D_zero

N = size(basis)
Allocate( h(N,N) , source = D_zero )

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!$OMP parallel do &
!$OMP   default(shared) &
!$OMP   schedule(dynamic, 1) &
!$OMP   private(ib, ia, Rab, jb, ja, j, i, DP_4_vector)
do ib = 1, system%atoms
    do ia = ib+1, system%atoms

        if ((system%QMMM(ib) /= "QM") .OR. (system%QMMM(ia) /= "QM")) then
            cycle
        end if

        Rab = GET_RAB(system%coord(ib,:), system%coord(ia,:))
        if (Rab > cutoff_Angs) then
           cycle
        end if

        If( evaluate ) then 
           select case (Environ_Type)
                case('DP_MM','DP_QM')
                    DP_4_vector = DP_phi( system , ia , ib )
                case default
                    DP_4_vector =  Q_phi( system , ia , ib )
           end select
           DP_4_matrix(ia,ib,:) = DP_4_vector
        else 
           DP_4_vector = DP_4_matrix(ia,ib,:)
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
 integer :: N_of_QM

 N_of_QM = count(a%QMMM == "QM")

 If( minloc( a%QMMM , dim=1 , mask = a%QMMM == "MM" ) < N_of_QM ) then
    stop ">> halting: block of QM atoms must precede solvent atoms if EnvField_ = T_ ; check input data <<"
 end if

 allocate( DP_4_matrix( N_of_QM , N_of_QM , 4 ) , source = D_zero ) 

 done = .true.

end subroutine allocate_DP_4_matrix
!
!
!
end module Hamiltonians

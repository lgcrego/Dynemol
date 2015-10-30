! Program for computing Hellman-Feynman-Pulay forces form Huckel Hamiltonian ...
module HuckelForces_m

    use type_m
    use constants_m
    use parameters_m            , only  : PBC 
    use PBC_m                   , only  : Generate_Periodic_Structure
    use Semi_Empirical_Parms    , only  : ChemAtom => atom
    use Allocation_m            , only  : DeAllocate_Structures    
    use Ehrenfest_Builder       , only  : RotationOverlap, solap, rotar, util_overlap, dlmn

    public :: HuckelForces

    private

    !module variables ...
    real*8  :: delta = 1.d-8

    !module parameters ...
    integer , parameter :: xyz_key(3) = [1,2,3]

contains
!
!
!
!==============================================
 subroutine HuckelForces( system , basis , QM )
!==============================================
 implicit none
 type(structure) , intent(in)  :: system
 type(STO_basis) , intent(in)  :: basis(:)
 type(R_eigen)   , intent(in)  :: QM

! local variables ... 
 integer                       :: i , i1 , i2 , n , n_MO , Fermi_level
 real*8          , allocatable :: bra(:), ket(:), Force(:,:)
 type(structure)               :: pbc_system
 type(STO_basis) , allocatable :: pbc_basis(:)

 CALL util_overlap     

 n_MO = size(QM%erg)
 allocate( bra  ( size(basis)              ) )
 allocate( ket  ( size(basis)              ) )
 allocate( Force( 3*system% atoms , 0:n_MO ) )

 ! if no PBC: pbc_system = system ...
 CALL Generate_Periodic_Structure( system, pbc_system, pbc_basis ) 

 do n = 1 , n_MO
    bra = QM%L(n,:)
    ket = QM%R(:,n)
    do i = 1 , system% atoms

        i1 = (i-1)*3 + 1
        i2 = (i-1)*3 + 3

        Force( i1:i2 ,n ) = Hellman_Feynman_Pulay( system, basis, pbc_system, pbc_basis, bra, ket, QM%erg(n), i )

    end do
    ! center of mass force ...
    Print 200, n , sum( Force(:,n) )
 end do

! Force(:,0) = total force at ground state ...
 Fermi_level = system% N_of_electrons / 2
 forall( i=1:size(Force(:,0)) ) Force(i,0) = sum( Force(i,1:Fermi_level) )

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
OPEN( unit=3 , file='EHM_forces.nmd' , status='unknown' )

write(3,*) "EHM Force Analysis"

write( 3 , '(A6 ,1000A3)'   ) "names "         , (system % Symbol(i)   , i = 1 , system% atoms)
write( 3 , '(A9 ,1000A4)'   ) "resnames "      , (system % residue(i)  , i = 1 , system% atoms)
write( 3 , '(A6 ,1000A2)'   ) "chids "         , [("A"                 , i = 1 , system% atoms)]             
write( 3 , '(A7 ,1000I4)'   ) "resids "        , (system % nr(i)       , i = 1 , system% atoms)
write( 3 , '(A6 ,1000A2)'   ) "betas "         , [("0"                 , i = 1 , system% atoms)]             
write( 3 , '(A12,3000F8.4)' ) "coordinates "   , (system % coord(i,:)  , i = 1 , system% atoms)

do n = 0 , n_MO
    write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , n , Force(:,n) 
end do

close(3)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 CALL Deallocate_Structures(pbc_system)
 deallocate(pbc_basis)

 include 'formats.h'

end subroutine HuckelForces
!
!
!
!========================================================================================================
function Hellman_Feynman_Pulay( system, basis, pbc_system, pbc_basis, bra, ket, erg, site ) result(Force)
!========================================================================================================
use util_m , factorial => fact
implicit none
type(structure)  , intent(in) :: system
type(STO_basis)  , intent(in) :: basis(:)
type(structure)  , intent(in) :: pbc_system
type(STO_basis)  , intent(in) :: pbc_basis(:)
real*8           , intent(in) :: bra(:)
real*8           , intent(in) :: ket(:)
real*8           , intent(in) :: erg
integer          , intent(in) :: site 

! local variables ...
real*8  :: expa , expb , Rab , aux , anor , Overlap_stuff
real*8  :: sux(0:10) , delta_a(3) , delta_b(3) 
integer :: a , b , ia , ib , ja , jb , aa ,xyz
integer :: na , la , ma , nb , lb , mb
integer :: msup , i , j , k , m , n

! local parameters ....
integer , parameter :: mxn = 15 , mxl = 5
real*8  , parameter :: cutoff_Angs = 12.d0

! local arrays ...
real*8  , dimension(0:mxl)                   :: solvec 
real*8  , dimension(0:mxl)                   :: solnorm , sol_partial
real*8  , dimension(-mxl:mxl,-mxl:mxl,0:mxl) :: rl , rl2
real*8                                       :: Force(3)

!force on atom site ...
Force = D_zero
!############################################################################################################

ib = site 
do ia = 1 , pbc_system% atoms  

    ! no self-interaction ...
    If( ia == ib ) cycle

    do xyz = 1 , 3

        do n = 1 , -1 , -2 

        delta_a = D_zero
        delta_b = n * delta * merge(D_one , D_zero , xyz_key == xyz )

        ! calculate rotation matrix for the highest l ...
        call RotationOverlap( system, pbc_system, ia, ib, delta_a, delta_b, Rab, rl, rl2 )

        If(Rab > cutoff_Angs) cycle

        do jb = 1 , ChemAtom(     system%AtNo(ib) )% DOS  ;  b  =      system%BasisPointer (ib) + jb
        do ja = 1 , ChemAtom( pbc_system%AtNo(ia) )% DOS  ;  a  =  pbc_system%BasisPointer (ia) + ja

            nb =     basis(b)% n  ;  lb =     basis(b)% l  ;  mb =     basis(b)% m
            na = pbc_basis(a)% n  ;  la = pbc_basis(a)% l  ;  ma = pbc_basis(a)% m
   
            aux = 1.d0 / ( factorial(na+na) * factorial(nb+nb) )
            msup = min(la,lb)

            solnorm(0:msup) = 0.d0

            do j = 1 ,     basis(b)% Nzeta
            do i = 1 , pbc_basis(a)% Nzeta
  
                expb =     basis(b)% zeta(j)
                expa = pbc_basis(a)% zeta(i)

                ! OVERLAP SUBROUTINE: lined-up frame
                call solap(na, la, expa, nb, lb, expb, Rab, solvec)

                anor = dsqrt((expa+expa)**(na+na+1)*(expb+expb)**(nb+nb+1)*(la+la+1.d0)*(lb+lb+1.d0)*aux) / fourPI

                ! Introduces normalization of the STO in the integrals 
                do m = 0, msup-1
                    sol_partial(m) = solvec(m) * anor
                    anor = anor / dsqrt((la+m+1.d0)*(la-m)*(lb+m+1.d0)*(lb-m))
                    if (m == 0) anor = anor + anor
                enddo
                sol_partial(msup) = solvec(msup) * anor
                forall(k=0:msup) solnorm(k) = solnorm(k) + pbc_basis(a)%coef(i)*basis(b)%coef(j)*sol_partial(k)
  
            end do
            end do

            ! Rotation of overlap integrals to the molecular frame
            sux(0) = solnorm(0) * rl(ma,0,la) * rl(mb,0,lb)
            forall(k=1:msup) sux(k) = solnorm(k) * ( rl(ma,-k,la)*rl(mb,-k,lb) + rl(ma,k,la)*rl(mb,k,lb) )

            Overlap_stuff = sum(sux(0:msup))
           
            ! reduce PBC to real system ...
            aa = ia - (pbc_basis(a)%copy_No) * system%atoms
            a  = a  - (pbc_basis(a)%copy_No) * size(basis)

            ! Force on atom ia due to atom ib ...
            Force(xyz) = Force(xyz) - n * ( Huckel_stuff(a,b,basis) - erg ) * bra(a)*ket(b) * Overlap_stuff 

        enddo 
        enddo 
        enddo 

    end do 

end do

!############################################################################################################

ia = site 
do ib = 1 , system% atoms  

    ! no self_interaction ...
    If( ia == ib ) cycle

    do xyz = 1 , 3

        do n = 1 , -1 , -2 

        delta_a = n * delta * merge(D_one , D_zero , xyz_key == xyz )
        delta_b = D_zero

        ! calculate rotation matrix for the highest l ...
        call RotationOverlap( system, pbc_system, ia, ib, delta_a, delta_b, Rab, rl, rl2 )

        If(Rab > cutoff_Angs) cycle

        do jb = 1 , ChemAtom(     system%AtNo(ib) )% DOS  ;  b  =      system%BasisPointer (ib) + jb
        do ja = 1 , ChemAtom( pbc_system%AtNo(ia) )% DOS  ;  a  =  pbc_system%BasisPointer (ia) + ja

            nb =     basis(b)% n  ;  lb =     basis(b)% l  ;  mb =     basis(b)% m
            na = pbc_basis(a)% n  ;  la = pbc_basis(a)% l  ;  ma = pbc_basis(a)% m
   
            aux = 1.d0 / ( factorial(na+na) * factorial(nb+nb) )
            msup = min(la,lb)

            solnorm(0:msup) = 0.d0

            do j = 1 ,     basis(b)% Nzeta
            do i = 1 , pbc_basis(a)% Nzeta
  
                expb =     basis(b)% zeta(j)
                expa = pbc_basis(a)% zeta(i)

                ! OVERLAP SUBROUTINE: lined-up frame
                call solap(na, la, expa, nb, lb, expb, Rab, solvec)

                anor = dsqrt((expa+expa)**(na+na+1)*(expb+expb)**(nb+nb+1)*(la+la+1.d0)*(lb+lb+1.d0)*aux) / fourPI

                ! Introduces normalization of the STO in the integrals 
                do m = 0, msup-1
                    sol_partial(m) = solvec(m) * anor
                    anor = anor / dsqrt((la+m+1.d0)*(la-m)*(lb+m+1.d0)*(lb-m))
                    if (m == 0) anor = anor + anor
                enddo
                sol_partial(msup) = solvec(msup) * anor
                forall(k=0:msup) solnorm(k) = solnorm(k) + pbc_basis(a)%coef(i)*basis(b)%coef(j)*sol_partial(k)
  
            end do
            end do

            ! Rotation of overlap integrals to the molecular frame
            sux(0) = solnorm(0) * rl(ma,0,la) * rl(mb,0,lb)
            forall(k=1:msup) sux(k) = solnorm(k) * ( rl(ma,-k,la)*rl(mb,-k,lb) + rl(ma,k,la)*rl(mb,k,lb) )

            Overlap_stuff = sum(sux(0:msup))
           
            ! reduce PBC to real system ...
            aa = ia - (pbc_basis(a)%copy_No) * system%atoms
            a  = a  - (pbc_basis(a)%copy_No) * size(basis)

            ! Force on atom ia due to atom ib 
            Force(xyz) = Force(xyz) - n * ( Huckel_stuff(a,b,basis) - erg ) * bra(a)*ket(b) * Overlap_stuff 

        enddo 
        enddo 
        enddo 

    end do 

end do

Force = Force / (TWO*delta)

end function Hellman_Feynman_Pulay
!
!
!
!
!======================================
 function Huckel_stuff( i , j , basis )
!======================================
implicit none
integer         , intent(in) :: i , j
type(STO_basis) , intent(in) :: basis(:)

!local variables ...
real*8 :: Huckel_stuff
real*8 :: k_eff , k_WH , c1 , c2 , c3

!-------------------------------------------------
!    constants for the Huckel Hamiltonian

c1 = basis(i)%IP - basis(j)%IP
c2 = basis(i)%IP + basis(j)%IP

c3 = (c1/c2)*(c1/c2)

k_WH = (basis(i)%k_WH + basis(j)%k_WH) / two

k_eff = k_WH + c3 + c3 * c3 * (D_one - k_WH)

Huckel_stuff = k_eff * c2 * HALF 

end function Huckel_stuff
!
!
!
!
end module HuckelForces_m

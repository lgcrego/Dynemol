module TD_Dipole_m

    use type_m
    use omp_lib
    use constants_m
    use blas95
    use parameters_m                , only : n_part
    use Structure_Builder           , only : Extended_Cell ,                &
                                             ExCell_basis
    use DP_main_m                   , only : DP_Matrix_AO
    use tuning_m                    , only : eh_tag

    public :: wavepacket_DP , DeAllocate_DPs , wavepacket 

    type(dipoles) , save :: wavepacket

    private

contains
!
! this module calculates the contribution to DP vector of solute and solvent molecules(nr) 
! due EXCLUSIVELY to el-hl wavepacket charge ...
! Valence charge contribution will be calculated elsewhere ... 
!
!==================================================================
 subroutine wavepacket_DP( a , basis , AO_bra , AO_ket , Dual_ket )
!==================================================================
implicit none
type(structure) , intent(inout) :: a
type(STO_basis) , intent(in)    :: basis    (:)
complex*16      , intent(in)    :: AO_bra   (:,:)
complex*16      , intent(in)    :: AO_ket   (:,:)
complex*16      , intent(in)    :: Dual_ket (:,:)

! local variables ...
integer       :: nr_max , xyz
real*8        :: xyz_FMO(3) , solvation_shell
type(dipoles) :: solvent , solute 

! center of charge of special DPF ...
forall( xyz=1:3 ) xyz_FMO(xyz) = sum( a%coord(:,xyz) * a%Nvalen(:) , a%DPF ) / sum( a%Nvalen , a%DPF )

! el_hl wavepacket component to DP only for molecules within the solvation_shell centered at xyz_FMO ...
solvation_shell = maxval( a%solvation_hardcore , a%DPF ) + 3.d0

CALL preprocess_wavepacket_DP( a , basis , AO_bra , AO_ket , Dual_ket , solute  , xyz_FMO , solvation_shell , instance = "solute"  )
CALL preprocess_wavepacket_DP( a , basis , AO_bra , AO_ket , Dual_ket , solvent , xyz_FMO , solvation_shell , instance = "solvent" )

nr_max = maxval( a%nr )

If( allocated(wavepacket%el_DP) ) CALL DeAllocate_DPs( wavepacket , flag = "dealloc" )

CALL DeAllocate_DPs( wavepacket , nr_max , flag = "alloc" )

! output of this subroutine: nr , CC , DP of (solute and solvent) molecules 
wavepacket % el_DP = solute % el_DP + solvent % el_DP 
wavepacket % hl_DP = solute % hl_DP + solvent % hl_DP 

CALL DeAllocate_DPs( solvent , flag = "dealloc" )
CALL DeAllocate_DPs( solute  , flag = "dealloc" )

end subroutine wavepacket_DP
!
!
!
!==============================================================================================================================
 subroutine preprocess_wavepacket_DP( a , basis , AO_bra , AO_ket , Dual_ket , molecule , origin , solvation_shell , instance )
!==============================================================================================================================
implicit none
type(structure) , intent(inout) :: a
type(STO_basis) , intent(in)    :: basis   (:)
complex*16      , intent(in)    :: AO_bra  (:,:)
complex*16      , intent(in)    :: AO_ket  (:,:)
complex*16      , intent(in)    :: Dual_ket(:,:)
type(dipoles)   , intent(out)   :: molecule
real*8          , intent(in)    :: origin(3)
real*8          , intent(in)    :: solvation_shell
character(*)    , intent(in)    :: instance

! local variables ...
integer                :: i , j , i1 , i2 , nr , last_nr , first_nr , nr_atoms , nr_max , np
real*8                 :: total_valence 
real*8 , allocatable   :: Qi_Ri(:,:) 

select case( instance )

    case( "solvent" )

        ! find positions of solvent atoms in structure ...
        first_nr = minval( a%nr , a%fragment == "S" )
        last_nr  = maxval( a%nr , a%fragment == "S" )

        a%N_of_Solvent_Molecules = (last_nr - first_nr + 1)

    case( "solute" )

        ! find positions of solute atoms in structure ...
        first_nr = minval( a%nr , a%solute )
        last_nr  = maxval( a%nr , a%solute )

        a%N_of_Solute_Molecules = (last_nr - first_nr + 1)

end select

nr_max = maxval( a%nr )

CALL DeAllocate_DPs( molecule , nr_max , flag = "alloc" )

allocate( Qi_Ri(a%atoms,3) , source = D_zero )     

!$OMP parallel default(shared) private( nr , nr_atoms , a , i1 , i2 , i , j , Qi_Ri , total_valence , np )
!$OMP single
do nr = first_nr , last_nr 

!$OMP task

    ! No of atoms with tag nr ...
    nr_atoms = count( a%nr == nr ) 

    ! position of nr residue in structure ...
    i1 = minloc( a%nr , 1 , a%nr == nr ) 
    i2 = (i1-1) + nr_atoms 

    !-----------------------------------------------------
    ! Center_of_Charge(nr) for molecule(nr) ...

    forall( j=1:3 , i=i1:i2 ) Qi_Ri(i,j) = a%Nvalen(i) * a%coord(i,j)

    total_valence = sum( [ (a%Nvalen(i) , i=i1,i2) ] )

    forall(j=1:3) molecule%CC(nr,j) = sum( Qi_Ri(i1:i2,j) ) / total_valence

    !-----------------------------------------------------

    ! if molecule(nr) is within solvation shell ...
    If( sqrt(dot_product(molecule%CC(nr,:)-origin(:) , molecule%CC(nr,:)-origin(:))) <= solvation_shell ) then

        ! calculate wavepacket component of the dipole vector ...
        do np = 1 , n_part
            If( eh_tag(np) == "el" )    &
            CALL Build_wavepacket_DP( a , basis , AO_bra(:,np) , AO_ket(:,np) , Dual_ket(:,np) , nr , molecule%CC(nr,:) , molecule%el_DP(nr,:) ) 
            
            If( eh_tag(np) == "hl" )    &
            CALL Build_wavepacket_DP( a , basis , AO_bra(:,np) , AO_ket(:,np) , Dual_ket(:,np) , nr , molecule%CC(nr,:) , molecule%hl_DP(nr,:) ) 
        end do

    end If

!$OMP end task    

end do

!$OMP end single
!$OMP end parallel

deallocate( Qi_Ri )

end subroutine preprocess_wavepacket_DP
!
!
!
!
!=================================================================================================
 subroutine Build_wavepacket_DP( system , basis , bra , ket , Dual_ket , nr , Q_center , DP_wave )
!=================================================================================================
implicit none
type(structure) , intent(in)    :: system
type(STO_basis) , intent(in)    :: basis(:)
complex*16      , intent(in)    :: bra(:)
complex*16      , intent(in)    :: ket(:)
complex*16      , intent(in)    :: Dual_ket(:)
integer         , intent(in)    :: nr
real*8          , intent(in)    :: Q_center(3)
real*8          , intent(inout) :: DP_wave(3)

! local variables ...
integer                         :: i , xyz , n_basis 
real*8          , allocatable   :: R_vector(:,:)
complex*16      , allocatable   :: a(:) , b(:,:)
logical         , allocatable   :: nr_mask(:)
type(R3_vector)                 :: origin_Dependent , origin_Independent 

! all atomic positions measured from Center_of_Charge(nr) ...
! if origin = Center_of_Charge ==> Nuclear_DP = (0,0,0) ...
allocate( R_vector(system%atoms,3) )
forall( xyz=1:3 ) R_vector(:,xyz) = system%coord(:,xyz) - Q_center(xyz)

n_basis = size(basis)
 
allocate( a(n_basis        ) )
allocate( b(n_basis,n_basis) )

! mask for molecule(nr) ...
allocate( nr_mask(n_basis) , source = .false. )
where( basis%nr == nr ) nr_mask = .true.

do xyz = 1 , 3

    ! origin dependent DP = sum{bra * vec{R} * S_ij * ket}

    forall( i=1:n_basis ) a(i) = bra(i) * R_vector(basis(i)%atom,xyz)

    origin_Dependent%DP(xyz) = real( sum( a(:)*Dual_ket(:) , nr_mask ) )

    ! origin independent DP = sum{bra * vec{DP_matrix_AO(i,j)} * ket}

    b = DP_matrix_AO(:,:,xyz)

    CALL gemv( b , ket , a , C_one , C_zero )    

    origin_Independent%DP(xyz) = real( sum( bra(:)*a(:) , nr_mask ) )

end do

DP_wave = origin_Dependent%DP + origin_Independent%DP 

deallocate( nr_mask , a , b )

end subroutine Build_wavepacket_DP
!
!
!
!
!===========================================
 subroutine DeAllocate_DPs( DPs , n , flag )
!===========================================
implicit none
type(dipoles)            ,   intent(inout) :: DPs
integer       , optional ,   intent(in)    :: n
character(*)             ,   intent(in)    :: flag

select case( flag )

    case( "alloc" )

        allocate( DPs%CC    ( n , 3 ) , source = D_zero )
        allocate( DPs%el_DP ( n , 3 ) , source = D_zero )
        allocate( DPs%hl_DP ( n , 3 ) , source = D_zero )
        allocate( DPs%nr    ( n )     , source = I_zero )

    case( "dealloc" )

        deallocate( DPs%CC    )
        deallocate( DPs%el_DP )
        deallocate( DPs%hl_DP )
        deallocate( DPs%nr    )

    case( "garbage" )

        ! garbage collection before restart calculations ...
        If( allocated( DPs%CC    ) ) deallocate( DPs%CC    )
        If( allocated( DPs%el_DP ) ) deallocate( DPs%el_DP )
        If( allocated( DPs%hl_DP ) ) deallocate( DPs%hl_DP )
        If( allocated( DPs%nr    ) ) deallocate( DPs%nr    )

end select

end subroutine DeAllocate_DPs
!
!
!
end module TD_Dipole_m

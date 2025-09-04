module F_intra_m

    use type_m   
    use omp_lib
    use constants_m
    use parameters_m      , only: QMMM , driver , n_part
    use Allocation_m      , only: Allocate_Structures
    use for_force         , only: pot_INTER, bdpot, angpot, dihpot, Morspot, LJ_14, LJ_intra, Coul_14, Coul_intra, pot_total, Vself
    use MD_read_m         , only: atom , molecule , MM 
    use Ehrenfest_CSDM    , only: Ehrenfest 
    use Ehrenfest_Builder , only: EhrenfestForce 
    use Surface_Hopping   , only: SH_Force
    use decoherence_m     , only: DecoherenceForce
    use FF_bonds          , only: f_bond
    use FF_angles         , only: f_angle
    use FF_diheds         , only: f_dihed
    use FF_Morse          , only: f_Morse
    use FF_intra_nonbond  , only: f_intra_nonbonding

    use Reax_F_intra      , only: DW_f_intra

    private

    public :: ForceINTRA, pot_INTRA, BcastQMArgs, ForceQMMM

    ! module variables ...
    real*8 :: pot_INTRA

    ! imported module variables ...
    integer         :: PST(2)
    real*8          :: t_rate
    character(2)    :: mode
    type(structure) :: system
    type(R_eigen)   :: QM
    type(STO_basis) , allocatable , dimension(:)   :: basis(:)
    Complex*16      , allocatable , dimension(:,:) :: MO_bra , MO_ket , MO_TDSE_bra , MO_TDSE_ket

    interface BcastQMArgs
        module procedure AllocateQMArgs
        module procedure StoreQMArgs_AO_and_FSSH
        module procedure StoreQMArgs_CSDM
    end interface

contains
!
!
!====================
subroutine ForceINTRA
!====================
implicit none
! local_variables ...
integer :: i

! stretching potential ...
call f_bond()

! Angle potential ...
call f_angle()

! Dihedral Potential Angle ... 
call f_dihed()

! Dissociative Potentials ... 
call f_Morse()

! IntraMolecular Nonbonding Potentials ... 
call f_intra_nonbonding()


!====================================================================
! dissociative forces

if( any(molecule%DWFF) ) then

    call DW_f_intra ()

endif


!
!====================================================================
! factor used to compensate the factor1 and factor2 factors ...
! factor3 = 1.0d-20
pot_INTRA = (bdpot + angpot + dihpot)*factor3 + LJ_14 + LJ_intra + Coul_14 + Coul_intra + Morspot
pot_total = pot_INTER + pot_INTRA - Vself
pot_total = pot_total * (mol*micro/MM % N_of_molecules)

if( QMMM ) then
    select case (driver)

       case( "slice_CSDM" ) 
           CALL Ehrenfest( system, basis, MO_bra, MO_ket, QM )
           CALL DecoherenceForce( system , MO_bra , MO_ket , QM%erg , PST )

       case( "slice_AO")
           CALL EhrenfestForce( system , basis , MO_bra , MO_ket , QM , representation=mode)    

       case("slice_FSSH")
           CALL SH_Force( system , basis , MO_bra , MO_ket , QM , t_rate )
    
    end select
endif

! Get total MM force; force units = J/mts = Newtons ...
do i = 1 , MM % N_of_atoms
    
    atom(i)% f_MM(:) = atom(i)% f_MM(:) + (atom(i) % fbond(:)    +  &
                                           atom(i) % fang(:)     +  &
                                           atom(i) % fdihed(:)   +  &
                                           atom(i) % fnonbd14(:) +  & 
                                           atom(i) % fnonch14(:) +  &
                                           atom(i) % fnonbd(:)   +  & 
                                           atom(i) % fMorse(:)   +  & 
                                           atom(i) % fnonch(:)      &
                                          ) * Angs_2_mts

    atom(i)% ftotal(:) = atom(i)% f_MM(:) 
    enddo

end subroutine ForceINTRA
!
!
!
!=====================
 subroutine ForceQMMM
!=====================
 implicit none

! local variables ...
integer :: i

If( driver == "slice_CSDM" ) then
    do i = 1 , MM % N_of_atoms
         atom(i)% f_QM(:)   = atom(i)% Ehrenfest(:) + atom(i)% f_CSDM(:)
         atom(i)% ftotal(:) = atom(i)% f_MM(:) + atom(i)% f_QM(:)
         end do
else
    do i = 1 , MM % N_of_atoms
         atom(i)% f_QM(:)   = atom(i)% Ehrenfest(:)
         atom(i)% ftotal(:) = atom(i)% f_MM(:) + atom(i)% f_QM(:)
         end do
endif

end subroutine ForceQMMM
!
!
!
!====================================================================
 subroutine StoreQMArgs_CSDM( sys , vec , mtx1 , mtx2 , Eigen , PSE )
!====================================================================
 implicit none
 type(structure), intent(inout):: sys
 type(STO_basis), intent(in)   :: vec(:)
 complex*16     , intent(in)   :: mtx1(:,:)
 complex*16     , intent(in)   :: mtx2(:,:)
 type(R_eigen)  , intent(in)   :: Eigen
 integer        , intent(in)   :: PSE(2)

! local variables ... 

basis  = vec
system = sys
PST    = PSE

MO_bra      = mtx1
MO_ket      = mtx2

QM% erg = Eigen% erg
QM% L   = Eigen% L
QM% R   = Eigen% R
QM% Fermi_state = Eigen% Fermi_state

end subroutine StoreQMArgs_CSDM
!
!
!
!===================================================================================
 subroutine StoreQMArgs_AO_and_FSSH( sys , vec , mtx1 , mtx2 , Eigen , delta , txt )
!===================================================================================
 implicit none
 type(structure)             , intent(inout):: sys
 type(STO_basis)             , intent(in)   :: vec(:)
 complex*16                  , intent(in)   :: mtx1(:,:)
 complex*16                  , intent(in)   :: mtx2(:,:)
 type(R_eigen)               , intent(in)   :: Eigen
 real*8           , optional , intent(in)   :: delta
 character(len=*) , optional , intent(in)   :: txt

! local variables ... 

basis   = vec
system  = sys
MO_bra  = mtx1
MO_ket  = mtx2
if(present(delta)) t_rate  = delta
if(present(txt  )) mode    = txt

QM% erg = Eigen% erg
QM% L   = Eigen% L
QM% R   = Eigen% R
QM% Fermi_state = Eigen% Fermi_state

end subroutine StoreQMArgs_AO_and_FSSH
!
!
!
!===================================================
 subroutine AllocateQMArgs( BasisSize , SystemSize )
!===================================================
implicit none
integer , intent(in) :: BasisSize
integer , intent(in) :: SystemSize

allocate(basis  (BasisSize)        )
allocate(MO_bra (BasisSize,n_part) )
allocate(MO_ket (BasisSize,n_part) )

CALL Allocate_Structures( SystemSize , System )

allocate(QM%erg(BasisSize)          )
allocate(QM%L  (BasisSize,BasisSize))
allocate(QM%R  (BasisSize,BasisSize))

end subroutine AllocateQMArgs
!
!
end module F_intra_m

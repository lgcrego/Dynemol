module good_vibrations_m

    use type_m                  
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use parameters_m            , only : PBC
    use MD_read_m               , only : atom , MM 
    use MM_types                , only : MM_atomic 
    use F_intra_m               , only : ForceIntra
    use MM_ERG_class_m          , only : MM_OPT
    use FF_OPT_class_m          , only : FF_OPT , LogicalKey 
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              

    public :: Optimize_Structure , normal_modes , Optimize_Parameters_Driver

    private 

    ! module variables ...
    type(MM_OPT) :: MM_erg
    type(FF_OPT) :: MM_parms
    type(MM_atomic) , allocatable :: atom0(:)

contains
!
!
!
!======================================
 subroutine Optimize_Parameters_Driver
!======================================
implicit none

! local variables ...
real*8                        :: local_minimum 
type(LogicalKey)              :: key
logical                       :: F_ = .false. , T_ = .true.

! saving reference structure for future use ...
allocate( atom0 , source = atom )

! instantiating MM ...
key%bonds  = [T_,T_,F_]
key%angs   = [T_,T_]
key%diheds = F_

MM_parms = FF_OPT( key , kernel = "energy" )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , local_minimum )

call optimize_structure()
print*, RMSD()

print*, "==> first part of FF_OPT done"

! instantiating MM ...
key%bonds  = [T_,F_,F_]
key%angs   = [T_,F_]
key%diheds = F_

MM_parms = FF_OPT( key , kernel = "NormalModes" )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , local_minimum )

atom = atom0
call optimize_structure()

print*, MM_parms%modes
print*, RMSD()

call normal_modes()
stop











end subroutine Optimize_Parameters_Driver 
!
!
!
!==========================
 subroutine normal_modes( )
!==========================
implicit none

! local variables ...
integer                       :: i , j , k , l , column , size_Hessian , info
type(MM_atomic) , allocatable :: equilibrium(:) , atom_fwd(:) , atom_bwd(:)
real*8          , allocatable :: Hessian(:,:)
type(R_eigen)                 :: Hesse

!local parameters ...
real*8 , parameter :: delta         = 1.d-5             ! <== displacement in Angs.
real*8 , parameter :: eV_2_cm_inv   = 1.d-12*8065.73    ! <== displacement in Angs.

! start the normal mode calculations from an energy minimum ...
CALL Optimize_Structure ( )

allocate( equilibrium ( MM % N_of_atoms ) )

do i = 1 , 3
    equilibrium % xyz(i) = MM_erg % p( (i-1) * MM%N_of_atoms + 1 : i * MM%N_of_atoms ) 
end do

! reset the atomic forces ...
forall(i=1:3) atom % ftotal(i) = D_zero

! start build up of Hessian matrix ...
allocate( atom_fwd (   MM % N_of_atoms ) )
allocate( atom_bwd (   MM % N_of_atoms ) )
allocate( Hessian  ( 3*MM % N_of_atoms , 3*MM % N_of_atoms ) )

column = 1
do k = 1 , MM%N_of_atoms

    do j = 1 , 3

        atom(k) % xyz(j) = equilibrium(k) % xyz(j) + delta
        forall(i=1:3) atom % ftotal(i) = D_zero
        CALL ForceIntra
        forall(i=1:3) atom_fwd % ftotal(i) = atom % ftotal(i)

        atom(k) % xyz(j) = equilibrium(k) % xyz(j) - delta
        forall(i=1:3) atom % ftotal(i) = D_zero
        CALL ForceIntra
        forall(i=1:3) atom_bwd % ftotal(i) = atom % ftotal(i)

        do i = 1 , MM%N_of_atoms
            do l = 1 , 3
                Hessian( (i-1)*3+l , column ) = - (atom_fwd(i) % ftotal(l) - atom_bwd(i) % ftotal(l)) / (two*delta)
            end do
        end do

        column = column + 1

    end do

end do

deallocate( equilibrium , atom_fwd , atom_bwd )

! atom%mass * imol = atomic mass in kg ; this is the unit of mass used in the verlet equations ...
forall( i=1:MM%N_of_atoms , l=1:3 ) Hessian( (i-1)*3 + l , : ) = Hessian( (i-1)*3 + l , : ) / sqrt( atom(i)%mass*imol )
forall( j=1:MM%N_of_atoms , k=1:3 ) Hessian( : , (j-1)*3 + k ) = Hessian( : , (j-1)*3 + k ) / sqrt( atom(j)%mass*imol )

! fixing correct units ...
Hessian = Hessian / Angs_2_mts

size_Hessian = 3 * MM % N_of_atoms

allocate( Hesse%erg(size_Hessian) , source=D_zero )

CALL SYEV( Hessian , Hesse % erg , 'V' , 'U' , info )
If ( info /= 0 ) write(*,*) 'info = ',info,' in SYEV in vibes normal modes '

! transforming back the normal modes: A --> M^{-1/2}*A ...
forall( i=1:MM%N_of_atoms , l=1:3 ) Hessian( (i-1)*3 + l , : ) = Hessian( (i-1)*3 + l , : ) / sqrt(atom(i)%mass)

! convert units to cm^{-1} ...
hesse%erg = sqrt( abs(hesse%erg) ) * h_bar*ev_2_cm_inv

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
OPEN( unit=3 , file='normal_modes.nmd' , status='unknown' )

write(3,*) "Normal Mode Analysis"

write( 3 , '(A6 ,1000A3)'   ) "names "         , atom % Symbol 
write( 3 , '(A9 ,1000A4)'   ) "resnames "      , atom % residue
write( 3 , '(A6 ,1000A2)'   ) "chids "         , [("A" , i=1,MM%N_of_atoms)]             
write( 3 , '(A7 ,1000I4)'   ) "resids "        , atom % nr
write( 3 , '(A6 ,1000A2)'   ) "betas "         , [("0" , i=1,MM%N_of_atoms)]             
write( 3 , '(A12,3000F8.4)' ) "coordinates "   , ( atom(i) % xyz(:) , i = 1 , MM%N_of_atoms )

do i = 7 , size_Hessian
    write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , i ,  Hessian(:,i) 
end do

close(3)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

do i = 1 , size_Hessian
write(30,*) i , hesse%erg(i)
end do

end subroutine normal_modes
!
!
!
!=================================
 subroutine Optimize_Structure( )
!=================================
implicit none

! local variables ...
real*8  :: local_minimum 

! instantiating MM ...
MM_erg = MM_OPT( )

MM_erg % driver = "MM_Optimize"

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_erg , MM_erg%N_of_Freedom , local_minimum )

end subroutine Optimize_Structure
!
!
!=============
 function RMSD 
!=============
implicit none
real*8  :: RMSD

!local variable ...
integer  :: i
real*8   :: rij(3) , distance

distance = D_zero
do i = 1 , size(atom)
    rij(:)   = atom(i) % xyz(:) - atom0(i) % xyz(:)
    rij(:)   = rij(:) - MM % box(:) * DNINT( rij(:) * MM % ibox(:) ) * PBC(:)
    RMSD = RMSD + SQRT( rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3) )
end do

end function RMSD
!
!
!
end module good_vibrations_m

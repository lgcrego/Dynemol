module vibrational_modes_m

    use type_m
    use constants_m
    use mkl95_precision
    use mkl95_blas
    use mkl95_lapack
    use MD_read_m               , only : atom , MM 
    use MM_types                , only : MM_system , MM_atomic , debug_MM    
    use MM_ERG_class_m          , only : CG_OPT
    use F_intra_m               , only : ForceIntra
    use NonlinearMM_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              

    public :: Optimize_Structure , normal_modes

    private 

    ! module variables ...
    type(CG_OPT) :: MM_erg

contains
!
!
!
!=================================
 subroutine normal_modes( )
!=================================
implicit none

! local variables ...
integer                       :: i , j , k , l , column , size_Hessian , info
type(MM_atomic) , allocatable :: equilibrium(:) , atom_fwd(:) , atom_bwd(:)
real*8          , allocatable :: Hessian(:,:)
type(R_eigen)                 :: Hesse

!local parameters ...
real*8 , parameter :: delta = 1.d-5   ! <== displacement in Angs.

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

size_Hessian = 3 * MM % N_of_atoms

allocate( Hesse%erg(size_Hessian) , source=D_zero )

CALL SYEV( Hessian , Hesse % erg , 'V' , 'U' , info )

If ( info /= 0 ) write(*,*) 'info = ',info,' in SYEV in vibes normal modes '

print*, hesse%erg

stop

end subroutine normal_modes
!
!
!
!=================================
 subroutine Optimize_Structure( )
!=================================
implicit none

! local variables ...
integer :: i , GlobalMinimum
real*8  :: local_minimum , InitialERG
real*8  ,  allocatable :: forces(:)

! instantiating MM ...
MM_erg = CG_OPT( )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_erg , MM_erg%N_of_Freedom , local_minimum )

end subroutine Optimize_Structure
!
!
end module vibrational_modes_m

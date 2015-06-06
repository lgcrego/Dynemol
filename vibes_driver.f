module good_vibrations_m

    use type_m                  
    use constants_m
    use f95_precision
    use blas95
    use lapack95
    use parameters_m            , only : PBC
    use MD_read_m               , only : atom , MM , molecule
    use MM_types                , only : MM_atomic 
    use setup_m                 , only : Setup
    use Babel_m                 , only : QMMM_key
    use F_intra_m               , only : ForceIntra
    use cost_MM                 , only : nmd_REF_erg , nmd_NOPT_erg , KeyHolder , LogicalKey 
    use MM_ERG_class_m          , only : MM_OPT
    use FF_OPT_class_m          , only : FF_OPT , atom0
    use NonlinearCG_m           , only : Fletcher_Reeves_Polak_Ribiere_minimization                              
    use GA_m                    , only : Genetic_Algorithm

    public :: Optimize_Structure , Normal_Modes , Optimize_Parameters_Driver

    private 

    ! module variables ...
    logical      :: done = .false.
    type(MM_OPT) :: MM_erg
    type(FF_OPT) :: MM_parms

contains
!
!
!
!======================================
 subroutine Optimize_Parameters_Driver
!======================================
implicit none

! local variables ...
real*8                         :: local_minimum 
real*8           , allocatable :: GA_Selection(:,:)
logical                        :: F_ = .false. , T_ = .true. 
type(LogicalKey)               :: key

print*, "==> kernel = energy"
! instantiating MM ...
key%bonds  = [T_,T_,F_]
key%angs   = [T_,T_]
key%diheds = F_

MM_parms = FF_OPT( key , kernel = "energy" )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , local_minimum )

! preprocess with GA method ...
CALL Genetic_Algorithm( MM_parms , GA_Selection )

CALL CG_driver( GA_Selection ) 

CALL MM_parms % output( 0 ) 

atom = atom0

call normal_modes( )

stop

end subroutine Optimize_Parameters_Driver 
!
!
!
!====================================
 subroutine CG_Driver( GA_Selection )
!====================================
implicit none
real*8           , allocatable , intent(inout) :: GA_Selection(:,:)

! local variables ...
integer                          :: i , k , GlobalMinimum
integer                          :: Top_Selection 
real*8                           :: this_minimum
real*8             , allocatable :: local_minimum(:) , InitialCost(:)
type(LogicalKey)                 :: key 


Top_Selection = size(GA_Selection(1,:)) 

allocate( local_minimum(Top_Selection) , source = real_large )
allocate( InitialCost  (Top_Selection)                       )

do i = 1 , Top_Selection

    atom = atom0

    do k = 1 , size(KeyHolder)

        key = KeyHolder(k)

        MM_parms = FF_OPT( key , kernel = "NormalModes" , weights = "use_no_weights" )

        write(*,190) i , KeyHolder(k) % comment , MM_parms% weights

        If( k == 1 ) then

            MM_parms % p = GA_Selection(:,i)

            InitialCost(i) = MM_parms % cost()

        end If

        CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_parms , MM_parms%N_of_Freedom , this_minimum )

        If( this_minimum == real_large ) goto 21

        ! temporarily stores CG optimized FF parameters here ...
        CALL save_temporary_results( GA_Selection(:,i) )
        local_minimum(i) = this_minimum

    end do

21  print*, "  ==> done"

end do

Print 191, ( InitialCost(i) , local_minimum(i) , i = 1 , Top_Selection )

GlobalMinimum = minloc( local_minimum , dim=1 )
key           = KeyHolder( size(KeyHolder) )
MM_parms      = FF_OPT( key , kernel = "NormalModes" , weights = "use_no_weights" )
MM_parms % p  = GA_Selection(:,GlobalMinimum)

Print 192, GlobalMinimum , MM_parms%cost() , RMSD()
Print 193, MM_parms%nmd_OPT_indx

include 'formats.h'

end subroutine CG_Driver
!
!
!
!==========================
 subroutine normal_modes( )
!==========================
implicit none

! local variables ...
integer                       :: i , j , k , l , column , size_Hessian , info , list_size , N_of_free
type(MM_atomic) , allocatable :: equilibrium(:) , atom_fwd(:) , atom_bwd(:)
real*8          , allocatable :: Hessian(:,:)
real*8                        :: symmetric 
type(R_eigen)                 :: Hesse

!local parameters ...
real*8 , parameter :: delta         = 1.d-8             ! <== displacement in Angs.
real*8 , parameter :: eV_2_cm_inv   = 1.d-12*8065.73    ! <== displacement in Angs.

! start the normal mode calculations from an energy minimum ...
CALL Optimize_Structure ( )

allocate( equilibrium ( MM % N_of_atoms ) )

N_of_free = MM_erg % N_of_freedom / 3
do i = 1 , 3
    equilibrium % xyz(i) = MM_erg % p( (i-1) * N_of_free + 1 : i * N_of_free ) 
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

! just to guarantee that the Hessian is symmetric ...
size_Hessian = 3 * MM % N_of_atoms
do j = 1 , size_Hessian
    do i = j+1 , size_Hessian

         symmetric = ( Hessian(i,j) + Hessian(j,i) ) * HALF
         Hessian(i,j) = symmetric
         Hessian(j,i) = symmetric

    end do
end do

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

write( 3 , '(A6 ,1000A3)'   ) "names "         , (atom(i) % Symbol   , i = 1 , MM%N_of_atoms)
write( 3 , '(A9 ,1000A4)'   ) "resnames "      , (atom(i) % residue  , i = 1 , MM%N_of_atoms)
write( 3 , '(A6 ,1000A2)'   ) "chids "         , [("A"               , i = 1 , MM%N_of_atoms)]             
write( 3 , '(A7 ,1000I4)'   ) "resids "        , (atom(i) % nr       , i = 1 , MM%N_of_atoms)
write( 3 , '(A6 ,1000A2)'   ) "betas "         , [("0"               , i = 1 , MM%N_of_atoms)]             
write( 3 , '(A12,3000F8.4)' ) "coordinates "   , ( atom(i) % xyz(:)  , i = 1 , MM%N_of_atoms )

OPEN( unit=4 , file='nmd_erg.dat' , status='unknown' )

If( MM_parms%driver == "Parametrize" ) then

    write(4,'(a86)',advance="no") "# nmd_indx | OPT nmd ergs | REF nmd ergs | NOPT nmd ergs | OPT %-ERROR  | NOPT %-ERROR "
    write(4,*) " "

    list_size = size(MM_parms%nmd_OPT_indx)

    call sort( MM_parms%nmd_REF_indx )

    do i = 1 , list_size

       write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , MM_parms%nmd_REF_indx(i) , Hessian(:,MM_parms%nmd_OPT_indx(i)) 

       write(4,4)               MM_parms%nmd_REF_indx(i) - 6  , &
                  hesse%erg   ( MM_parms%nmd_OPT_indx(i) )    , &
                  nmd_REF_erg ( i )                           , &
                  nmd_NOPT_erg( i )                           , &
                  abs( hesse%erg(MM_parms%nmd_OPT_indx(i)) - nmd_REF_erg ( i ) ) / nmd_REF_erg(i) * 100.0 , &
                  abs( nmd_NOPT_erg(i)                     - nmd_REF_erg ( i ) ) / nmd_REF_erg(i) * 100.0 
    end do

    CALL system( "mv OPT_nmd_indx.inpt OPT_nmd_indx.old" )

    OPEN( unit=14 , file='OPT_nmd_indx.out' )
        write(14,*) size(MM_parms%nmd_OPT_indx)
        write(14,*) MM_parms % nmd_OPT_indx
    close(14) 

else

    do i = 1 , size_Hessian
        write(4,*) i , hesse%erg(i)
        write( 3 , '(A5 ,I4,3000F8.4)' ) "mode " , i ,  Hessian(:,i) 
    end do

end If

close(3)
close(4)

4 FORMAT(i4,t15,F10.3,t30,F10.3,t45,F10.3,t60,F10.3,t75,F10.3)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

! setting up the MM system ...
If( .not. done ) then

    If( MM%N_of_molecules > I_one ) CALL Setup
    atom( QMMM_key ) % charge = atom( QMMM_key ) % MM_charge
    done = .true. 

end If

! instantiating MM ...
MM_erg = MM_OPT( )

CALL Fletcher_Reeves_Polak_Ribiere_minimization( MM_erg , MM_erg%N_of_Freedom , local_minimum )

end subroutine Optimize_Structure
!
!
!
!=================================
 subroutine Moment_of_Inertia( R )
!=================================
implicit none
type(MM_atomic) , intent(inout) :: R(:) 

! local variables ...
integer              :: i , j , info
real*8               :: MoI(3,3) , eigenvalues(3), delta, R0(3)
real*8 , allocatable :: d2(:) 

R%mass = atom%mass

! translate molecule to CM...
forall(i=1:3) R%xyz(i) =  R%xyz(i) - sum(R%mass*R%xyz(i))/sum(R%mass)

allocate( d2(MM%N_of_atoms) )

forall( i=1:MM%N_of_atoms ) d2(i) = sqrt( sum( R(i)%xyz * R(i)%xyz ) )

do i = 1 , 3
    do j = 1 , 3

         delta = merge( D_one , D_zero , i == j )

         MoI(i,j) = sum( R(:)%mass*( delta*d2(:) - R(:)%xyz(i)*R(:)%xyz(j) ) )

    end do
end do

CALL SYEV( MoI , eigenvalues , 'V' , 'U' , info )
If ( info /= 0 ) write(*,*) 'info = ',info,' in Moment of Inertia in vibes normal modes '

end subroutine Moment_of_Inertia
!
!
!
!===============
 function RMSD() 
!===============
implicit none
real*8  :: RMSD

!local variable ...
integer  :: i
real*8   :: rij(3) , distance

RMSD     = D_zero
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
!======================================
 subroutine save_temporary_results( a )
!======================================
implicit none
real*8 , intent(inout) :: a(:)

! local variables ...
type(LogicalKey) :: key

key = KeyHolder( size(KeyHolder) )

MM_parms = FF_OPT( key , kernel="JustKey" )

a = MM_parms % p

end subroutine save_temporary_results
!
!
!
!===================
 subroutine  sort(a)
!===================
implicit none
integer , intent(inout) :: a(:)

! local variables ...
integer  :: ra, l, n, ir, i, j

!-----------------------------------------------------------
!  SORT IRA(I) , SO THAT THE ELEMENTS IRB(I) FOLLOW TOGETHER
!-----------------------------------------------------------
      n = size(a)
      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         ra  = a(l)
      else
         ra = a(ir)
         a(ir) = a(1)
         ir = ir - 1
         if(ir .eq. 1) then
             a(1) = ra
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(a(j) .lt. a(j+1)) j = j + 1
        endif
      if(ra .lt. a(j)) then
        a(i) = a(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      a(i) = ra
      goto 10

end subroutine sort
!
!
!
!======================
 subroutine  sort2(a,b)
!======================
implicit none
integer , intent(inout) :: a(:)
integer , intent(inout) :: b(:)

! local variables ...
integer  :: ra, rb , l, n, ir, i, j

!-----------------------------------------------------------
!  SORT IRA(I) , SO THAT THE ELEMENTS IRB(I) FOLLOW TOGETHER
!-----------------------------------------------------------
      n = size(a)
      l = n/2+1
      ir = n

10    continue
      if(l .gt. 1) then
         l = l -1
         ra  = a(l)
         rb = b(l)
      else
         ra = a(ir)
         rb = b(ir)
         a(ir) = a(1)
         b(ir) = b(1)
         ir = ir - 1
         if(ir .eq. 1) then
             a(1) = ra
             b(1) = rb
             return
         endif
      endif
      i = l
      j = l + l
20    if(j .le. ir) then
        if(j .lt. ir)then
          if(a(j) .lt. a(j+1)) j = j + 1
        endif
      if(ra .lt. a(j)) then
        a(i) = a(j)
        b(i) = b(j)
        i = j
        j = j + j
      else
      j = ir + 1
      endif
      goto 20
      endif
      a(i) = ra
      b(i) = rb
      goto 10

end subroutine sort2
!
!
!
end module good_vibrations_m

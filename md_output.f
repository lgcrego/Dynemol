module MD_dump_m

    use constants_m
    use MM_types        , only: MM_atomic
    use MM_input        , only: MM_frame_step
    use parameters_m    , only: n_t, restart
    use syst            , only: bath_T, Initial_density
    use MD_read_m       , only: MM , atom , molecule , species

    public :: output , cleanup , saving_MM_frame 

    private 

!   module variables ...
    logical         , save        :: first = .true. , done = .false.  
    type(MM_atomic) , allocatable :: previous(:)

contains    
!
!
!
!========================================
 subroutine OUTPUT( Ttrans , frame , dt ) 
!========================================
use for_force   , only: rcut, pot_INTER, ecoul, eintra, evdw, bdpot, angpot, dihpot,    &
                        LJ_14, LJ_intra, Coul_14, Coul_intra, pot_total, forcefield,    &
                        ryck_dih, proper_dih, harmpot, Morspot
implicit none
real*8  , intent(in)    :: Ttrans
integer , intent(in)    :: frame 
real*8  , intent(in)    :: dt
 
! local variables ... 
 integer :: i

IF( .NOT. done ) then 
    open (10, file='MM_log.out', status='unknown')
    write(10,*)
    write(10,'(''********************************************'')')
    write(10,'(''*                                          *'')')
    write(10,'(''*               MM_Dynamics                *'')')
    write(10,'(''*                                          *'')')
    write(10,'(''********************************************'')')
    write(10,*)
    write(10,*)            
    write(10,*)'Initial Paramenters'
    write(10,*)
    write(10,*)
    write(10,'(''Number of Molecules      :'',I10)') MM % N_of_molecules
    write(10,'(''Initial Bath Temperature :'',F10.2,'' Kelvin'')') bath_T
    write(10,'(''Box dimensions           :'',3F10.2,'' Angstroms'')') MM % box(1:3)
    write(10,'(''Cut off radius           :'',F10.2,'' Angstroms'')') rcut
    write(10,'(''Density                  :'',F10.4,'' g/cm3'')') Initial_density
    write(10,*)
    write(10,'(''System Temperature       :'',F10.2,'' Kelvin'')') Ttrans
    write(10,*)
    write(10,'(''Total Simulation Steps   :'',i10)') n_t
    write(10,'(''Integration time step    :'',F10.3,'' femtosec'')') dt*1.d15
    write(10,'(''frame-output step        :'',i10,'' steps'')') MM_frame_step
    write(10,*)

    select case( forcefield )

        case( 1 )
        write(10,'(''Born Meyer Potential'')')

        case( 2) 
        write(10,'(''Lennard-Jones Potential'')')

    end select

    write(10,*)
    write(10,'(''Box contains'',i3,'' different species'')') MM % N_of_species
    write(10,*)
    do i = 1, MM % N_of_species
        if (species(i) % N_of_atoms > 1) then
            write(10,'(''Species'',i2,'' comprised of '',i4,'' atoms'')') i, species(i) % N_of_atoms
        endif
        write(10,*)
    end do 
    close (10)

    done = .true.

end IF

open (10, file='MM_log.out', status='unknown', access='append')

    write(10,'(''<======  ###############  ==>'')')
    write(10,'(''<====  A V E R A G E S  ====>'')')
    write(10,'(''<==  ###############  ======>'')')
    write(10,*)
    write(10,'(''time :'',F10.4,'' ps'')') frame*dt*1.d12
    write(10,*)'Energies (kJ/mol)'
    write(10,'(''Bond Potential             :'',F12.4)') bdpot     *mol*factor3*1.d-6    
    write(10,'(''Harm Bond Potential        :'',F12.4)') harmpot   *mol*factor3*1.d-6    
    write(10,'(''Morse Bond Potential       :'',F12.4)') Morspot   *mol*factor3*1.d-6    
    write(10,'(''Angle Potential            :'',F12.4)') angpot    *mol*factor3*1.d-6   
    write(10,'(''Dihedral Potential         :'',F12.4)') dihpot    *mol*factor3*1.d-6   
    write(10,'(''Proper Dihedral            :'',F12.4)') proper_dih*mol*factor3*1.d-6   
    write(10,'(''Ryckaert-Bell. Dihedral    :'',F12.4)') ryck_dih  *mol*factor3*1.d-6   
    write(10,'(''Lennard-Jones 1-4          :'',F12.4)') LJ_14     *mol*1.d-6  
    write(10,'(''Lennard-Jones Intra        :'',F12.4)') LJ_Intra  *mol*1.d-6  
    write(10,'(''Coulomb 1-4                :'',F12.4)') Coul_14   *mol*1.d-6  
    write(10,'(''Coulomb Intra              :'',F12.4)') Coul_Intra*mol*1.d-6  
    write(10,'(''Lennard-Jones              :'',F12.4)') evdw      *mol*1.d-6      
    write(10,'(''Coulomb self-interaction   :'',F15.4)') eintra    *mol*1.d-6    
    write(10,'(''Coulomb short-range        :'',F15.4)') ecoul     *mol*1.d-6      
    write(10,'(''Total Coulomb              :'',F15.4)') (-(ecoul + eintra)*mol*1.d-6 ) 
    write(10,'(''Potential (INTER) Energy   :'',F12.4)') pot_INTER*mol*1.d-6 / MM % N_of_molecules
    write(10,'(''Potential Energy(gmx-like) :'',ES16.7E3)') pot_total*mol*1.d-3 / MM % N_of_molecules  
    write(10,*)

close (10)
10 format(A2,3F10.4,2X,F10.4)

end subroutine OUTPUT
!
!
!
!=======================================
subroutine SAVING_MM_frame( frame , dt )
!=======================================
implicit none
integer , intent(in)    :: frame
real*8  , intent(in)    :: dt

! local variables ...
integer :: i , l 

! saving latest atomic velocities for future use ... 
open(11, file='velocity_MM.out', status='unknown')
    do i = 1 , MM % N_of_atoms
        write(11,*) atom(i) % vel(1),  atom(i) % vel(2), atom(i) % vel(3) 
    end do
close(11)
 
! config.pdb ... 

 CALL ReGroupMolecule

 open (14, file='frames-MM.pdb', status='unknown', access='append')
        if( first .AND. (.NOT. restart) ) write(14,*)  "MDFlex , no title"
        write(14,995) 'TITLE'  , 'manipulated by MDFlex     t= ', frame*dt*1.d12
        write(14,991) 'CRYST1' , MM % box(1) , MM % box(2) , MM % box(3) , 90.0 , 90.0 , 90.0 , 'P 1' , '1'
        write(14,993) 'MODEL'  , frame

        do i = 1 , MM % N_of_atoms
             write(14,992)  atom(i) % my_id          ,          &     ! <== global number
                            atom(i) % EHSymbol       ,          &     ! <== atom type of Huckel's 
                            ' '                      ,          &     ! <== alternate location indicator
                            atom(i) % residue        ,          &     ! <== residue name
                            ' '                      ,          &     ! <== chain identifier
                            atom(i) % nr             ,          &     ! <== residue sequence number
                            ' '                      ,          &     ! <== code for insertion of residues
                            ( atom(i) % xyz(l) , l=1,3 )   ,    &     ! <== xyz coordinates
                            1.00                     ,          &     ! <== occupancy
                            0.00                     ,          &     ! <== temperature factor
                            ' '                      ,          &     ! <== segment identifier
                            ' '                      ,          &     ! <== here only for tabulation purposes
                            atom(i) % Symbol         ,          &     ! <== chemical element symbol
                            atom(i) % charge         ,          &     ! <== charge on the atom
                            atom(i) % MMSymbol                        ! <== atom type of Molecular Mechanics
        end do
        write(14,'(''MASTER'')')
        write(14,'(''END'')')
 close(14) 
 first = .false.

991 FORMAT(a6,3F9.3,3F7.2,a11,a4)
992 FORMAT('ATOM  ',i5,a5,a1,a3,a2,i4,a4,3F8.3,2F6.2,a4,a6,a2,F8.4,t90,a3)
993 FORMAT(a5,i8)
994 FORMAT(a6,t15,a21)
995 FORMAT(a5,t15,a29,F12.6)
999 format(' ',1X,I7,2X,4F11.5)
       
end subroutine SAVING_MM_frame
!
!
!
!=========================
subroutine ReGroupMolecule
!=========================
implicit none

! local variables ...
integer :: i , j , j1 , j2 , xyz , nresid 
real*8  :: CartesianDistance
logical :: BigMolecule
 
If( .NOT. allocated(previous) ) allocate( previous(MM%N_of_atoms) , source=atom )

 do i = 1 , MM % N_of_molecules
    nresid = molecule(i) % nr
    j1 = sum(molecule(1:nresid-1) % N_of_atoms) + 1
    j2 = sum(molecule(1:nresid) % N_of_atoms)

    BigMolecule = any( [(maxval(atom(j1:j2)%xyz(j)) - minval(atom(j1:j2)%xyz(j)), j=1,3)] > MM%box(:)*HALF )

    If( BigMolecule ) then

        do xyz = 1 , 3
            CartesianDistance = atom(i) % xyz(xyz) - previous(i) % xyz(xyz)
            If( abs(CartesianDistance) > MM % box(xyz)*HALF ) atom(i) % xyz(xyz) = atom(i) % xyz(xyz) - sign( MM % box(xyz) , CartesianDistance )
        end do

    else

        do j = j1 , j2
            do xyz = 1 , 3
                CartesianDistance = atom(j) % xyz(xyz) - atom(j1) % xyz(xyz)
                If( abs(CartesianDistance) > MM % box(xyz)*HALF ) atom(j) % xyz(xyz) = atom(j) % xyz(xyz) - sign( MM % box(xyz) , CartesianDistance )
            end do
        end do

    end iF

 end do

 previous = atom

end subroutine ReGroupMolecule
!
!
!
!==================
 subroutine CLEANUP 
!==================
 implicit none
 logical :: filefound
 character (len=64) :: cmd, filename

 filename = 'frames-MM.pdb'
 inquire(file=filename, EXIST=filefound)
 if (filefound) then
   write (cmd, '("/bin/rm ", A)' ) TRIM (filename)
   call SYSTEM (cmd)
 endif

end subroutine CLEANUP
!
!
!
end module MD_dump_m

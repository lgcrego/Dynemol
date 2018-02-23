! Identify intramolecular non-bonded atom pairs; general use ...
module NonBondPairs

use constants_m
use MM_types               , only : MM_molecular

private
 
public :: Identify_NonBondPairs


contains
!
!
!
!
!===============================================
 subroutine Identify_NonBondPairs( species , a )
!===============================================
implicit none
type(MM_molecular)  , intent(inout) :: species(:)
integer             , intent(in)    :: a

! local variables ...
integer         , allocatable   :: vector_of_pairs(:,:) , vector_of_pairs14(:,:)
integer                         :: i , j , k , m , n , NN , idx
logical                         :: flagB1, flagB2, flagA1, flagA2, flagD1, flagD2, flagB11, flagB12, flagB21, flagB22
logical                         :: flagB111, flagB112, flagB121, flagB122, flagB211, flagB212, flagB221, flagB222 , flag1, flag2
logical         , allocatable   :: InputRef(:,:) , Input14(:,:)

!==============================================================================================

 If( species(a) % N_of_atoms <= 3 ) return

 NN = species(a)% Nbonds14

 allocate( InputRef      ( NN , NN ) , source = .true. )
 allocate( Input14       ( NN , NN ) , source = .false. )

 ! Intramolecular LJ list generation ... 
 ! Preparing ...
 
 do i = 1 , species(a) % N_of_atoms

     idx = 1
     ! Looking for bonds ...  
     do k = 1 , species(a) % Nbonds
         flagB1 = ( species(a) % atom(i) % my_id == species(a) % bonds(k,1) )
         flagB2 = ( species(a) % atom(i) % my_id == species(a) % bonds(k,2) )

         if ( flagB1 ) then
             InputRef(i,species(a)%bonds(k,2)) = .false.
             idx = idx + 1

             ! looking for angle connections ...  
             do n = 1 , species(a) % Nbonds
                flagB11 = ( species(a) % bonds(k,2) == species(a) % bonds(n,1) )
                flagB12 = ( species(a) % bonds(k,2) == species(a) % bonds(n,2) )

               if ( flagB11 ) then
                   InputRef(i,species(a)%bonds(n,2)) = .false.
                   idx = idx + 1

                   ! looking for dihedral connections ...  
                   do m = 1 , species(a) % Nbonds
                      flagB111 = ( species(a) % bonds(n,2) == species(a) % bonds(m,1) )
                      flagB112 = ( species(a) % bonds(n,2) == species(a) % bonds(m,2) )
                      if ( flagB111 ) then
                         InputRef(i,species(a)%bonds(m,2)) = .false.
                         Input14(i,species(a)%bonds(m,2))  = .true.
                         idx = idx + 1
                      elseif ( flagB112 ) then
                         InputRef(i,species(a)%bonds(m,1)) = .false.
                         Input14(i,species(a)%bonds(m,1))  = .true.
                         idx = idx + 1
                      end if
                   end do

                ! looking for angle connections ...  
                elseif ( flagB12 ) then
                   InputRef(i,species(a)%bonds(n,1)) = .false.
                   idx = idx + 1

                   ! looking for dihedral connections ...  
                   do m = 1 , species(a) % Nbonds
                      flagB121 = ( species(a) % bonds(n,1) == species(a) % bonds(m,1) )
                      flagB122 = ( species(a) % bonds(n,1) == species(a) % bonds(m,2) )
                      if ( flagB121 ) then
                         InputRef(i,species(a)%bonds(m,2)) = .false.
                         Input14(i,species(a)%bonds(m,2))  = .true.
                         idx = idx + 1
                      elseif ( flagB122 ) then
                         InputRef(i,species(a)%bonds(m,1)) = .false.
                         Input14(i,species(a)%bonds(m,1))  = .true.
                         idx = idx + 1

                   end if
                   end do

                end if
             end do

         ! Looking for bonds ...  
         elseif ( flagB2 ) then

             InputRef(i,species(a)%bonds(k,1)) = .false.
             idx = idx + 1

             ! looking for angle connections ... 
             do n = 1 , species(a) % Nbonds
                flagB21 = ( species(a) % bonds(k,1) == species(a) % bonds(n,1) )
                flagB22 = ( species(a) % bonds(k,1) == species(a) % bonds(n,2) )
                if ( flagB21 ) then
                   InputRef(i,species(a)%bonds(n,2)) = .false.
                   idx = idx + 1

                   ! looking for dihedral connections ...  
                   do m = 1 , species(a) % Nbonds
                      flagB211 = ( species(a) % bonds(n,2) == species(a) % bonds(m,1) )
                      flagB212 = ( species(a) % bonds(n,2) == species(a) % bonds(m,2) )
                      if ( flagB211 ) then
                         InputRef(i,species(a)%bonds(m,2)) = .false.
                         Input14(i,species(a)%bonds(m,2))  = .true.
                         idx = idx + 1
                      elseif ( flagB212 ) then
                         InputRef(i,species(a)%bonds(m,1)) = .false.
                         Input14(i,species(a)%bonds(m,1))  = .true.
                         idx = idx + 1
                      end if
                   end do

               ! looking for angle connections ... 
               elseif ( flagB22 ) then
                   InputRef(i,species(a)%bonds(n,1)) = .false.
                   idx = idx + 1

                   ! looking for dihedral connections ...  
                   do m = 1 , species(a) % Nbonds
                      flagB221 = ( species(a) % bonds(n,1) == species(a) % bonds(m,1) )
                      flagB222 = ( species(a) % bonds(n,1) == species(a) % bonds(m,2) )
                      if ( flagB221 ) then
                         InputRef(i,species(a)%bonds(m,2)) = .false.
                         Input14(i,species(a)%bonds(m,2))  = .true.
                         idx = idx + 1
                      elseif ( flagB222 ) then
                         InputRef(i,species(a)%bonds(m,1)) = .false.
                         Input14(i,species(a)%bonds(m,1))  = .true.
                         idx = idx + 1
                      end if
                   end do

                end if
             end do

         end if
     end do

     ! Looking for 'explicit' angles (just in case) ...
     do k = 1 , species(a) % Nangs
         flagA1 = ( species(a) % atom(i) % my_id == species(a) % angs(k,1) )
         flagA2 = ( species(a) % atom(i) % my_id == species(a) % angs(k,3) )

         if ( flagA1 ) then
             InputRef(i,species(a)%angs(k,3)) = .false.
             idx = idx + 1
         elseif ( flagA2 ) then
             InputRef(i,species(a)%angs(k,1)) = .false.
             idx = idx + 1
         end if

     end do

     ! Looking for 'explicit' dihedrals (just in case) ...
     do k = 1 , species(a) % Ndiheds
       if( species(a) %funct_dihed(k) == 3 ) then
         flagD1 = ( species(a) % atom(i) % my_id == species(a) % diheds(k,1) )
         flagD2 = ( species(a) % atom(i) % my_id == species(a) % diheds(k,4) )

         if ( flagD1 ) then
             InputRef(i,species(a)%diheds(k,4)) = .false.
             Input14(i,species(a)%diheds(k,4))  = .true.
             idx = idx + 1
         elseif ( flagD2 ) then
             InputRef(i,species(a)%diheds(k,1)) = .false.
             Input14(i,species(a)%diheds(k,1))  = .true.
             idx = idx + 1
         end if
       end if
     end do
 end do 
 
 ! Intermediate variable ... 
 allocate( vector_of_pairs( species(a) % N_of_Atoms * species(a) % N_of_Atoms , 2 ) , source = I_zero )
 k = 1
 do i = 1 , species(a) % N_of_Atoms - 1

     do j = i , species(a) % N_of_Atoms
         if ( InputRef(i,j) == .true. ) then
             vector_of_pairs(k,1) = i
             vector_of_pairs(k,2) = j
             k = k + 1  
         end if
     end do

 end do


 allocate( species(a) % IntraLJ( size( pack( vector_of_pairs(:,2), vector_of_pairs(:,2) /= 0 ) ), 2 ) )
 species(a)%NintraLJ = size( pack( vector_of_pairs(:,2), vector_of_pairs(:,2) /= 0 ) ) 

 ! Finally associating the nonbonded interactions to species ...
 do i = 1 , size( pack( vector_of_pairs(:,2), vector_of_pairs(:,2) /= 0 ) )
     species(a) % IntraLJ(i,1) = vector_of_pairs(i,1)
     species(a) % IntraLJ(i,2) = vector_of_pairs(i,2)
 end do

! 1--4 Interactions: only if [ pairs ] is not read in the .itp file ...
if ( .NOT. allocated (species(a) % bonds14) ) then 
 ! Eliminating bonding terms in 1--4 interactions ...
 do i = 1 , species(a) % N_of_atoms
   do j = 1 , species(a) % N_of_atoms 
     do k = 1 , species(a) % Nbonds
       flag1 = ( species(a) % bonds(k,1) == i ) .and. ( species(a) % bonds(k,2) == j )
       flag2 = ( species(a) % bonds(k,1) == j ) .and. ( species(a) % bonds(k,2) == i )   
       if( ( flag1 ) .or. ( flag2 ) ) Input14(i,j) = .false.
     end do
   end do
 end do



! Intermediate variable for 1-4... 
allocate( vector_of_pairs14( species(a) % N_of_Atoms * species(a) % N_of_Atoms ,2 ) , source = I_zero )
 k = 1
 do i = 1 , species(a) % N_of_Atoms - 1

     do j = i , species(a) % N_of_Atoms
         if ( Input14(i,j) == .true. ) then
             vector_of_pairs14(k,1) = i
             vector_of_pairs14(k,2) = j
             k = k + 1
         end if
     end do
 end do
species(a) % Nbonds14 = ( size( pack( vector_of_pairs14(:,2), vector_of_pairs14(:,2) /= 0 ) ) ) 
 ! Preparing 1--4 interactions ...

! species(a) % Nbonds14 = ( size( pack( Input14(:,:), Input14(:,:) .eqv. .true. ) ) ) / 2
 
 allocate( species(a) % bonds14(species(a) % Nbonds14,2 ) )
 ! Associating 1--4 interactions to species ...
 k = 1
 do i = 1 , species(a) % N_of_atoms
   do j = i , species(a) % N_of_atoms 
     if( Input14(i,j) .eqv. .true. ) then
        species(a) % bonds14(k,1) = i
        species(a) % bonds14(k,2) = j
        k = k + 1
     end if
   end do
 end do

end if

 deallocate( InputRef , Input14 , vector_of_pairs )

!==============================================================================================
end subroutine Identify_NonBondPairs
!
!
!
end  module NonBondPairs

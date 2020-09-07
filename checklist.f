MODULE setup_checklist

 use type_m
 use parameters_m
 use MM_input        


 public :: checklist , dump_driver_parameters_and_tuning , Checking_Topology

 private

 ! module variables ... 
 logical :: done = .false.

 contains
!
!
!
!====================
subroutine checklist
!====================
implicit none 

! local variables ...

! feel free to add your own dynemol-for-dummies checklist ...
select case( DRIVER )

    case( "avrg_confgs" )

       If( file_type      /= "trajectory"   ) stop ">>> Mind: if DRIVER=^avrg_confgs^: file_type must be ^trajectory^; check parameters.f"
       If( nuclear_matter /= "extended_sys" ) stop ">>> Mind: if DRIVER=^avrg_confgs^: nuclear_matter must be ^extended_sys^; check parameters.f"

    case( "Genetic_Alg" )

       If( nuclear_matter /= "extended_sys" ) stop ">>> Mind: if DRIVER=^Genetic_Alg^: nuclear_matter must be ^extended_sys^; check parameters.f"
       If( EnvField_      == .true.         ) stop ">>> Mind: if DRIVER=^Genetic_Alg^: flag ^EnvField_ = F_^; check parameters.f"

end select


If ( (frame_step /= 1) .AND. (file_type /= "trajectory") ) then
    CALL system("sed '11i >>> halting: frame_step /= 1, only for avrg_confgs or time-slice dynamics <<<' warning.signal |cat")
    stop
End If


end subroutine checklist
!
!
!
!
!============================================
 subroutine dump_driver_parameters_and_tuning
!============================================
implicit none
 
! local variables ... 
 integer :: i
 character(len=3)  :: tag
 character(len=12) :: number_string

open (10, file='log.trunk/driver_parms_and_tuning.log', status='unknown')

    write(10,'(''<======  ###############  ==>'')')
    write(10,'(''<====    PARAMETERS.F   ====>'')')
    write(10,'(''<==  ###############  ======>'')')
    write(10,*)
    write(10,'(" DRIVER          :" , A12  )') DRIVER          

    if( DRIVER /= "MM_Dynamics" ) then

        write(10,'(" QMMM            :" , A10)') merge(".true. <==",".false.   ",QMMM)            
        write(10,'(" OPT_parms       :" , A10)') merge(".true. <==",".false.   ",OPT_parms)       
        write(10,'(" SPECTRUM        :" , A10)') merge(".true. <==",".false.   ",SPECTRUM)        
        write(10,'(" Alpha_Tensor    :" , A10)') merge(".true. <==",".false.   ",Alpha_Tensor)    

        tag = merge("no ","yes",GaussianCube)       
        write(10 , '(" GaussianCube    :" , A10)' , advance=tag) merge(".true. <==",".false.   ",GaussianCube)    
        If( GaussianCube ) &
            write(10,'(" GaussianCube_step = " , I0)') GaussianCube_step      

        tag = merge("no ","yes",NetCharge)       
        write(10 , '(" NetCharge       :" , A10)' , advance=tag) merge(".true. <==",".false.   ",NetCharge)       
        If( NetCharge ) &
            write(10,'(" CH_and_DP_step = " , I0)') CH_and_DP_step      

        write(10,'(" DensityMatrix   :" , A10)') merge(".true. <==",".false.   ",DensityMatrix)   
        write(10,'(" AutoCorrelation :" , A10)') merge(".true. <==",".false.   ",AutoCorrelation) 
        write(10,'(" VDOS_           :" , A10)') merge(".true. <==",".false.   ",VDOS_)           

        tag = merge("no ","yes",EnvField_)       
        write(10 , '(" EnvField_       :" , A10)' , advance=tag) merge(".true. <==",".false.   ",EnvField_)       
        If( EnvField_) then 
            write(10,'(" Environ_Type =" , A6)' , advance=tag) Environ_Type    
            write(10,'("   /   Environ_step = " , I0)') Environ_step    
        end IF
        write(10,'(" Coulomb_        :" , A10)') merge(".true. <==",".false.   ",Coulomb_)        
        write(10,'(" Induced_        :" , A10)') merge(".true. <==",".false.   ",Induced_)        
        write(10,'(" frame_step      : " , I0)') frame_step      

    end if

    write(10,'(" ad_hoc          :" , A10)') merge(".true. <==",".false.   ",ad_hoc)          
    write(10,'(" restart         :" , A10)') merge(".true. <==",".false.   ",restart)         
    write( number_string , '(F12.4)' ) t_f
    write(10,'(" t_f             : " , A12)') adjustl(number_string)
    write(10,'(" n_t             : " ,  I0)') n_t            

    if( (DRIVER(1:5) == "slice") .or.  DRIVER == "q_dynamics" .or. DRIVER == "avrg_confgs") then
        write(10,'(" CT_dump_step    : " , I0)') CT_dump_step   
        write(10,'(" n_part          : " , I0)') n_part         
        write(10,'(" hole_state      : " , I0)') hole_state     
        write(10,'(" electron_state  : " , I0)') electron_state 
    end if

    write(10,'(" nnx , nny       : " , I0,I2)') nnx , nny           
    write(10,'(" PBC             : " , I0,I2,I2)') PBC            
    write(10,*)


 if( nuclear_matter == "MDynamics" ) then

    write(10,'(''<======  ###############  ==>'')')
    write(10,'(''<====    parameters_MM.f   ====>'')')
    write(10,'(''<==  ###############  ======>'')')
    write(10,*)
    tag = merge("no ","yes",driver_MM == "MM_Dynamics")       
    write(10 , '(" driver_MM       : " , A11)' , advance=tag) driver_MM
    If( driver_MM == "MM_Dynamics" ) then 
        write(10,'(" <== thermostat  = "   , A14)' ) thermostat
        If( thermostat /= "Microcanonical" ) then
            write( number_string , '(F6.2)' ) temperature
            write(10,'(t36,"temperature = " , A6)') adjustl(number_string)
            write( number_string , '(F8.5)' ) thermal_relaxation_time
            write(10,'(t36,"relax time  = " , A7)') adjustl(number_string)
        end If
        write(10,'(t36,"read_velocities :" , A10)') merge(".true. <==",".false.   ",read_velocities)           
    end IF
    write( number_string , '(F6.2)' ) cutoff_radius
    write(10,'(" Wolf_cutoff     : " , A6)') adjustl(number_string)
    write( number_string , '(F8.5)' ) damping_Wolf
    write(10,'(" Wolf_damping    : " , A8)') adjustl(number_string)
    write(10,'(" MM_input_format : " , A4)') MM_input_format
    write(10,'(" MM_log_step     : " , I0)') MM_log_step
    write(10,'(" MM_frame_step   : " , I0)') MM_frame_step
    write(10,*)

 end if

close (10)

 if( ad_hoc ) then

    open (10, file='log.trunk/driver_parms_and_tuning.log', status='old', access='append')
    write(10,'(''<======  ###############  ==>'')')
    write(10,'(''<====       TUNING.F       ====>'')')
    write(10,'(''<==  ###############  ======>'')')
    write(10,*)
    close (10)

    if( DRIVER /= "MM_Dynamics" ) then
        open (10, file='log.trunk/driver_parms_and_tuning.log', status='old', access='append')
        write(10,'("==  Electronic Tuning  ==")')
        write(10,*)
        CALL system("awk '/univ/ && NR >= 37 && NR <= 70' < tuning.f | awk '!/ !/' >> log.trunk/driver_parms_and_tuning.log")
        close(10)
    end if


    if( nuclear_matter == "MDynamics" ) then
        open (10, file='log.trunk/driver_parms_and_tuning.log', status='old', access='append')
        write(10,*)
        write(10,'("==  Nuclear Tuning  ==")')
        write(10,*)
        CALL system("awk '/atom/&&/%/&&/=/' < tuning.f | awk '!/univ/' | awk '!/ !/' | awk '{gsub(/^[ \t]+/,x); print}' >> log.trunk/driver_parms_and_tuning.log")
        close(10)
    end if

 end if

 CALL system("echo dyn.trunk/ dos.trunk/ opt.trunk/ | xargs -n 1 cp log.trunk/driver_parms_and_tuning.log ")

end subroutine dump_driver_parameters_and_tuning
!
!
!
!================================================================
 function Checking_Topology( bonds , angs , diheds ) result(TorF)
!================================================================
implicit none
integer , intent(in) :: bonds (:,:)
integer , intent(in) :: angs  (:,:)
integer , intent(in) :: diheds(:,:)
logical              :: TorF
 
! local variables ... 
integer               :: i , j , x , y , z
integer               :: Nbonds , Nangs , Ndiheds , KeyLeft , KeyRight
integer , allocatable :: BondKeys(:) , AngKeys(:)

Nbonds  =  size(bonds (:,1)) 
Nangs   =  size(angs  (:,1))
Ndiheds =  size(diheds(:,1))

! checking bonds topology ...
allocate( BondKeys(Nbonds) )
do i = 1 , Nbonds

     x = bonds(i,1)  ;  y = bonds(i,2) 
     BondKeys(i) = PairingFunction( x , y , verbose = .true. ) 

end do

! checking angs topology ...
do i = 1 , Nangs

     x = angs(i,1)  ;  y = angs(i,2) 
     KeyLeft = PairingFunction(x,y) 
     If( .not. any(KeyLeft == BondKeys) ) call error_message(i,angs,instance="ang")

     x = angs(i,2)  ;  y = angs(i,3) 
     KeyRight = PairingFunction(x,y) 
     If( .not. any(KeyRight == BondKeys) ) call error_message(i,angs,instance="ang")

     If( KeyLeft == KeyRight ) call error_message(i,angs,instance="ang")

end do

! checking diheds topology ...
allocate( AngKeys(Nangs) )
do i = 1 , Nangs

     x = angs(i,1)  ;  y = angs(i,2)   ;  z = angs(i,3) 
     AngKeys(i) = CantorPairing( x , y , z ) 

end do

do i = 1 , Ndiheds

     x = diheds(i,1)  ;  y = diheds(i,2)   ;  z = diheds(i,3) 
     KeyLeft = CantorPairing( x , y , z ) 
     If( .not. any(KeyLeft == AngKeys) ) call error_message(i,diheds,instance="dihed")

     x = diheds(i,2)  ;  y = diheds(i,3)   ;  z = diheds(i,4) 
     KeyRight = CantorPairing( x , y , z ) 
     If( .not. any(KeyRight == AngKeys) ) call error_message(i,diheds,instance="dihed")

end do

! prepare to leave ...
if( done ) then  
    TorF = .true.     ! <==  error detected
    close(10)
else
    TorF = .false.    ! <==  NO error detected
end If

end function Checking_Topology
!
!
!
!
!===============================================
 function CantorPairing(i,j,k) result(R)
! 3-tupling Cantor Function ...
! f(i,j,k) = f(k,j,i)
!===============================================
implicit none
integer            , intent(in) :: i,j,k

! local variables ... 
integer :: R , L , a , b

! Symmetric Pairing for (i,k)-tuple ...
a = max(i,k)  ;  b = min(i,k)
L = a*(a+1)/2 + b 

! Cantor pairing with the center pairing ...
R = (L+j)*(L+j+1)/2 + L 

end function CantorPairing
!
!
!
!===============================================
 function PairingFunction(i,j,verbose) result(k)
!===============================================
implicit none
integer            , intent(in) :: i,j
logical , optional , intent(in) :: verbose

! local variables ... 
integer :: k , a , b

If( (i == j) .and. present(verbose) ) then
    Print 232, i , j
    stop
end If

! the symmetric pairing satisfies f(i,j)=f(j,i) ...

! Symmetric Cantor Pairing ...
!k = (i+j)*(i+j+1)/2 + (i*j) 

! Symmetric Pairing ...
a = max(i,j)  ;  b = min(i,j)
k = a*(a+1)/2 + b 

include 'formats.h'

end function PairingFunction
!
!
!
!===========================================
 subroutine error_message(i , a , instance ) 
!===========================================
implicit none
integer          , intent(in) :: i
integer          , intent(in) :: a(:,:)
character(len=*) , intent(in) :: instance

If( .not. done ) open (10, file='log.trunk/Topology.test.log', status='unknown')

select case (instance)

       case("ang")
       write(10,231) a(i,1) , a(i,2) , a(i,3) 

       case("dihed")
       write(10,233) a(i,1) , a(i,2) , a(i,3)  , a(i,4) 

end select

done = .true.

include 'formats.h'

end subroutine error_message
!
!
!
!
end MODULE setup_checklist

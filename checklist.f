MODULE setup_checklist

 use type_m
 use parameters_m
 use MM_input        


 public :: checklist , dump_driver_parameters_and_tuning

 private

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
        CALL system("awk '/univ/ && NR >= 37 && NR <= 70' < tuning.f | awk '!/\!/' >> log.trunk/driver_parms_and_tuning.log")
        close(10)
    end if


    if( nuclear_matter == "MDynamics" ) then
        open (10, file='log.trunk/driver_parms_and_tuning.log', status='old', access='append')
        write(10,*)
        write(10,'("==  Nuclear Tuning  ==")')
        write(10,*)
        CALL system("awk '/atom/&&/%/&&/=/' < tuning.f | awk '!/univ/' | awk '!/\!/' | awk '{gsub(/^[ \t]+/,x); print}' >> log.trunk/driver_parms_and_tuning.log")
        close(10)
    end if

 end if

 CALL system("echo dyn.trunk/ dos.trunk/ opt.trunk/ | xargs -n 1 cp log.trunk/driver_parms_and_tuning.log ")

end subroutine dump_driver_parameters_and_tuning
!
!
!
end MODULE setup_checklist

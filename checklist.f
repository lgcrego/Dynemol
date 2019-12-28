MODULE setup_checklist

 use type_m
 use parameters_m

 public :: checklist

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
    Print*, " >>> halting: frame_step /= 1, only for avrg_confgs or time-slice dynamics <<<"
    stop
End If


end subroutine checklist

end MODULE setup_checklist

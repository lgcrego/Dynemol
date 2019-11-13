MODULE parameter_checklist

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



end select









end subroutine checklist

end MODULE parameter_checklist

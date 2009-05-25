module interfaces_m

    interface

        subroutine projector(basis,coef_L,coef_R,zR,wv_FMO)
            use type_m
            type(STO_basis)          , intent(in)  :: basis(:)
            complex*16 , ALLOCATABLE , intent(out) :: coef_L(:,:) , coef_R(:,:)
            complex*16               , intent(in)  :: zR(:,:)
            real*8                   , intent(in)  :: wv_FMO(:,:)
        end subroutine projector

        function pop_Slater(za,zb,list_of_atoms,system)
            use type_m
            complex*16      , intent(in) :: za(:,:) , zb(:,:)
            integer         , intent(in) :: list_of_atoms(:)
            type(structure) , intent(in) :: system
        end function pop_Slater

        function Huckel(i,j,basis)
            use type_m
            integer         , intent(in) :: i , j
            type(STO_basis) , intent(in) :: basis(:)
        end function Huckel

        subroutine eigen_FMO(system,wv_FMO)
            use type_m
            type(structure)               , intent(in)  :: system
            real*8          , ALLOCATABLE , intent(out) :: wv_FMO(:,:)
        end subroutine eigen_FMO

        subroutine eigen(zL,zR,erg,basis)
            use type_m
            complex*16 , ALLOCATABLE , intent(out) :: zL(:,:) , zR(:,:)
            real*8     , ALLOCATABLE , intent(out) :: erg(:)
            type(STO_basis)          , intent(in)  :: basis(:)
        end subroutine eigen

        subroutine AllocateStructures
        end subroutine AllocateStructures

    end interface

end module interfaces_m

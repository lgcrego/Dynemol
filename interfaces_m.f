module interfaces_m

    interface Symbol_2_AtNo
        subroutine Sym_2_AtNo_TRJ(a)
            use type_m
            implicit none
            type(atomic) , intent(inout) :: a(:)
        end subroutine Sym_2_AtNo_TRJ

        subroutine Sym_2_AtNo_XYZ(a)
            use type_m
            implicit none
            type(structure) , intent(inout) :: a
        end subroutine Sym_2_AtNo_XYZ
    end interface 

end module interfaces_m

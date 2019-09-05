write(*,'(/,">>> =======  using pre-defined parameters for TiO2  ====== <<<")')

! atoms that belong to the cluster residue ...
where(  system % atom % fragment /= 'S' .AND. &
        system % atom % fragment /= 'F' .AND. &
        system % atom % fragment /= 'D'         ) system%atom%fragment = 'C'

! cluster residue is numbered as 1 (one) ...
where( system % atom % fragment == 'C' .AND. system % atom % nresid == 0 ) system%atom%nresid = 1

! bulk atoms of TiO2 cluster ...
where( system % atom % fragment == 'C' ) fragment = 'B'

! surface atoms of TiO2 cluster ...
where( (system % atom % xyz(3)  > 13.0) .AND. (system % atom % fragment == 'C') ) fragment =  'I'
where( (system % atom % xyz(3)  > 6.0 ) .AND. (system % atom % fragment == 'C') ) fragment =  'I'

! define MM atom types ...

where( (system % atom % symbol == 'Ti') .AND. (fragment == 'I') ) 
    system % atom % MMSymbol = 'Ti'
    system % atom % charge   = +2.196
end where

where( (system % atom % symbol == 'Ti') .AND. (fragment == 'B') ) 
    system % atom % MMSymbol = 'Ti'
    system % atom % charge   = +2.196
end where

where( (system % atom % symbol == 'O' ) .AND. (fragment == 'I') ) 
    system % atom % MMSymbol = 'O' 
    system % atom % charge   = -1.098
end where

where( (system % atom % symbol == 'O' ) .AND. (fragment == 'B') ) 
    system % atom % MMSymbol = 'O' 
    system % atom % charge   = -1.098
end where

where( (system % atom % Symbol == 'Ti') .AND. (system % atom % fragment == 'C') ) system % atom % mass = atomic_mass(22)
where( (system % atom % Symbol == 'O' ) .AND. (system % atom % fragment == 'C') ) system % atom % mass = atomic_mass(8)

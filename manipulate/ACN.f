! J. Comp. Chem. 21 (2000), 901 
! define MM atom types ...
sol_mol % atom(1) % MMSymbol = 'YN'
sol_mol % atom(1) % charge   = -0.490

sol_mol % atom(2) % MMSymbol = 'YC'
sol_mol % atom(2) % charge   = +0.382

sol_mol % atom(3) % MMSymbol = 'CT'
sol_mol % atom(3) % charge   = -0.2376

sol_mol % atom(4:6) % MMSymbol = 'HC'
sol_mol % atom(4:6) % charge   = +0.115 

! define charge group ...
sol_mol % atom(1:6) % nrcg   = 1 

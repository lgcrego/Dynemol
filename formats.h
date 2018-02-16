40    FORMAT(/,"======================   Opitimizing  Parameters  ==========================")

41    FORMAT(/,"MMsymbol    |   s   |   p   |   d   |   IP  |   zeta1   |   zeta2   |   k_WH")

42    FORMAT(A3,t13,I8,t21,I8,t29,I8,t37,I8,t45,I12,t57,I12,t69,I10)

421   FORMAT(A3,t17,I1,t25,I1,t33,I1,t41,I1,t49,I1,t61,I1,t73,I1)

43    FORMAT(/,"============================================================================")

44    FORMAT(/,">>>  Reading Optmized EHT Parameters")

45    FORMAT(A4)

46    FORMAT(/,">>>  Using Ad Hoc tuning ")

47    FORMAT(/,">>>  Saving structure.log ",/)

48    FORMAT("Symbol  |  MMsymbol  |  NoAt  |  Nvalen  |  Nzeta  |  n  |  spdf  |    IP   |  zeta1  |  zeta2  |  coef1  |  coef2  |  k_WH")

50    FORMAT(/,1x,'# of cluster states  = ',I5/,  &
               1x,'# of molecule states = ',I5/,  &
               1x,'# of unit cells = ',I2) 

51    FORMAT(1x,'>>> PBC  in  use  (',3I1,')  <<<')

52    FORMAT(/,1x,'>>> Hamiltonian  ',I4,'    /     ','time =',F8.4)  

53    FORMAT(/,1x,'========= STARTING  S_matrix  =========')

54    FORMAT(1x,' >>> no PBC')

55    FORMAT(1x,'========= FINISHING S_matrix  =========',/) 

56    FORMAT(/,1x,'Initial state of the isolated molecule => ',I3)
 
57    FORMAT(1x,'checking normalization of wave-packet =',F9.6,',',F9.6)

58    FORMAT(1x,'Electrons in the Slater determinant = ',F9.3)

591   FORMAT(/,1x,'Energy of El-packet state(',I3,') = ',F10.5)

592   FORMAT(/,1x,'Energy of Hl-packet state(',I3,') = ',F10.5)

60    FORMAT(1x,'norm of Psi(t) = ',F10.7)     

61    FORMAT(1x,'>> AO_preprocess done <<')

62    FORMAT(/,'>>>  EXCITED State Calculation ')

63    FORMAT(/,'>>>  GROUND State Calculation ')

69    FORMAT(1x,'Sparsity of OVERLAP Mtrx = ',F10.5)

70    FORMAT(/,'>>> System Characteristics: ',A72)

71    FORMAT(1x,' >>> PBC in the X direction')

72    FORMAT(1x,' >>> PBC in the Y direction')

73    FORMAT(1x,'Sparsity of DIPOLE Mtrx =',3F7.4)

74    FORMAT(/,"======================   TDDFT-Casida  Parameters  ==========================")

75    FORMAT(/,"occupied state  |   unoccupied state   |   transition density  ")

76    FORMAT(t6,I4,t26,I4,t49,F10.7)

100   FORMAT(I5,A4,F10.5,F10.5,F10.5)

101   FORMAT(F10.5,F10.5,F10.5,A3,A3,A3)

111   FORMAT(I5,3F12.6)

112   FORMAT(6E13.5)

113   FORMAT(I5,4F12.6)

120   FORMAT(/,1x,'Total number of orbitals  = ',I6)

121   FORMAT(1x,A2,' atoms  = ',I5)

122   FORMAT(/,1x,A3,' residues / atoms  = ',I5,' / ',I5)

123   FORMAT(/,1x,A1,' fragment atoms  = ',I5)

140   FORMAT(/,1x,'Total number of electrons  = ',I6)

141   FORMAT(/,1x,'Total number of atoms  = ',I6,' / # of flex atoms = ',I6)

142   FORMAT(/,1x,'Total number of QM atoms  = ',I6)

143   FORMAT(/,1x,'Total number of MM atoms  = ',I6,/)

152   FORMAT(1x,'REF DP Vector = (',3F8.4,') ==> ',F7.4,' Debye')

153   FORMAT(/,1x,'======== STARTING  DIPOLE MATRIX ANALYSIS =========')

154   FORMAT(1x,'DIPOLE Vector = (',3F8.4,') ==> ',F7.4,' Debye')

155   FORMAT(1x,'========= FINISHING DIPOLE MATRIX ANALYSIS ========',/)

156   FORMAT(1x,' >>> no DIPOLE')

157   FORMAT(/,1x,'N_of_Solvent_Molecules in DP_field = ',I4/)

158   FORMAT(/,1x,'N_of_Solute_Molecules in DP_field = ',I4/)

160   FORMAT(/,1x,'Custo_GA : ',I4,' /',I4)

161   FORMAT(/,1x,'Custo_GA : ',I4,' /',I4 , F15.4 , a25)

165   FORMAT(1x,'========= ',A20,' =========')

166   FORMAT(1x,'truncation error = ',F8.5)

180   FORMAT(/,1x,'step = ',I6,'		time = ',F7.4,a9)

181   FORMAT(1x,'step is = ',a3)

182   FORMAT(1x,'step tried = ',F7.4)

183   FORMAT(1x,'step did = ',F7.4)
 
184   FORMAT(1x,'step to try = ',F7.4)

185   FORMAT(1x,'truncation error = ',F10.7)

186   FORMAT(1x,'time = ',F9.5,' fs')

187   FORMAT(1x,'Sparsity of H_prime Mtrx = ',F10.5)

188   FORMAT(1x,'Polarizability Tensor diagonal elements = (',3F10.3,') Angs^3    ==>  Average = ',F10.3,' Angs^3')

189   FORMAT(1x,'Polarizability Tensor diagonal elements = (',3F10.3,') a_Bohr^3  ==>  Average = ',F10.3,' a_Bohr^3')

190   FORMAT(/,I2,' ', a20 , a30)

191   FORMAT(2F32.6/)

192   FORMAT(1x,'Global Minimun = ',I2,'    /    Cost = ',F12.6,'   /    Geometry RMSD = ',F10.6)  

193   FORMAT(/,1x,"OPT_nmd_indx",/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I,/,6I)  

194   FORMAT(/,1x,"  ==> done with Recursive steps")

200   FORMAT(1x,'Center of Mass Force (MO = ',I3,') = ',F10.7,'   eV/Angs')     


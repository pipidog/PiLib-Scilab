scc.HubU=[1,0;2,0;15,0;16,0]        // U for each state, [state_label, U] or blank
scc.Charge=[1,0,0.5*ones(1,12),0,1] // charge of each state, 1x total state
scc.Mixing=[1]                    	// mixing parameter, 0~1
scc.Iteration=[30]                	// maximal iterations, integer
scc.Converge=[10^-3]              	// convergence criterion, real, at least < 0.1
scc.Mesh=[300,1]                    // k-space mesh for Ef, (1x1,1x2,1x3), large for metal
scc.Temperature=[100]             	// temperature for searching Ef, large for insulator
scc.Memory=['normal']       	  	// 'normal', 'max' or 'HDD' 

============= PiLib Variable =============
scc.E_Fermi, @full, the Fermi level
ORDER=   -7, SIZE=[    1,    1], TYPE=REAL

           1

   -3.141010

============= PiLib Variable =============
scc.E_gap, @full, the band gap
ORDER=   -6, SIZE=[    1,    1], TYPE=REAL

           1

    2.879874

============= PiLib Variable =============
scc.DM_out, @t-sp, the output density matrix
ORDER=    1, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1    0.049998  0.000000
       2       2    0.049998  0.000000
       3       3    0.050001  0.000000
       4       4    0.050001  0.000000
       5       5    0.050001  0.000000
       6       6    0.050001  0.000000
       7       7    0.050001  0.000000
       8       8    0.050001  0.000000
       9       9    0.050001  0.000000
      10      10    0.050001  0.000000
      11      11    0.050001  0.000000
      12      12    0.050001  0.000000
      13      13    0.050001  0.000000
      14      14    0.050001  0.000000
      15      15    0.049998  0.000000
      16      16    0.049998  0.000000

============= PiLib Variable =============
scc.U_mat, @t-sp, the Hubbard potential matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
scc.H_onsite, @t-sp, H_onsite(hop.onsite_E+hop.LS_mat+scc.U_mat)
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

scc.HubU=[]                       // U for each state, [state_label, U] or blank
scc.Charge=[0.5*ones(1:16)]                     // charge of each state, 1x total state
scc.Mixing=[1]                    // mixing parameter, 0~1
scc.Iteration=[30]                // maximal iterations, integer
scc.Converge=[10^-3]              // convergence criterion, real, at least < 0.1
scc.Mesh=[30,30]                       // k-space mesh for Ef, (1x1,1x2,1x3), large for metal
scc.Temperature=[100]             // temperature for searching Ef, large for insulator
scc.Memory=['normal']             // 'normal', 'max' or 'HDD' 

============= PiLib Variable =============
scc.E_Fermi, @full, the Fermi level
ORDER=   -1, SIZE=[    1,    1], TYPE=REAL

           1

   -1.265197

============= PiLib Variable =============
scc.E_gap, @full, the band gap
ORDER=   -1, SIZE=[    1,    1], TYPE=REAL

           1

    3.411688

============= PiLib Variable =============
scc.DM_out, @t-sp, the output density matrix
ORDER=   -1, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1    6.094720  0.000000
       2       2    6.094720  0.000000
       3       3    4.474823  0.000000
       4       4    4.474823  0.000000
       5       5    4.955634  0.000000
       6       6    4.474823  0.000000
       7       7    4.474823  0.000000
       8       8    4.955634  0.000000
       9       9    6.094720  0.000000
      10      10    6.094720  0.000000
      11      11    4.474823  0.000000
      12      12    4.474823  0.000000
      13      13    4.955634  0.000000
      14      14    4.474823  0.000000
      15      15    4.474823  0.000000
      16      16    4.955634  0.000000

============= PiLib Variable =============
scc.U_mat, @t-sp, the Hubbard potential matrix
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
scc.H_onsite, @t-sp, H_onsite(hop.onsite_E+hop.LS_mat+scc.U_mat)
ORDER=    0, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1   -8.868000  0.000000
       2       2   -8.868000  0.000000
       3       4    0.000000  1.500000
       3       8   -1.500000  0.000000
       4       8    0.000000 -1.500000
       5       6    1.500000  0.000000
       5       7    0.000000  1.500000
       6       7    0.000000 -1.500000
       9       9   -8.868000  0.000000
      10      10   -8.868000  0.000000
      11      12    0.000000  1.500000
      11      16   -1.500000  0.000000
      12      16    0.000000 -1.500000
      13      14    1.500000  0.000000
      13      15    0.000000  1.500000
      14      15    0.000000 -1.500000

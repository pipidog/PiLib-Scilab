scc.HubU=...                      // U for each state, [state_label, U] or blank
[1,10;2,10;3,10;4,10;5,10;6,10;7,10;8,10;9,10;10,10]
scc.Charge=...				      // charge of each state, 1x total state
[1,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1]                     
scc.Mixing=[1]                    // mixing parameter, 0~1
scc.Iteration=[30]                // maximal iterations, integer
scc.Converge=[10^-3]              // convergence criterion, real, at least < 0.1
scc.Mesh=[15,15,15]               // k-space mesh for Ef, (1x1,1x2,1x3), large for metal
scc.Temperature=[100]             // temperature for searching Ef, large for insulator
scc.Memory='max'  			      // 'normal', 'max', 'HDD' to store H_k

============= PiLib Variable =============
scc.E_Fermi, @full, the Fermi level
ORDER=    0, SIZE=[    1,    1], TYPE=REAL

           1

    1.107022

============= PiLib Variable =============
scc.E_gap, @full, the band gap
ORDER=    0, SIZE=[    1,    1], TYPE=REAL

           1

    4.678230

============= PiLib Variable =============
scc.DM_out, @t-sp, the output density matrix
ORDER=    1, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1    0.099891  0.000000
       2       2    0.099868  0.000000
       3       3    0.099865  0.000000
       4       4    0.015351  0.000000
       5       5    0.014631  0.000000
       6       6    0.099891  0.000000
       7       7    0.099868  0.000000
       8       8    0.099865  0.000000
       9       9    0.015351  0.000000
      10      10    0.014631  0.000000
      11      11    0.090000  0.000000
      12      12    0.090280  0.000000
      13      13    0.090115  0.000000
      14      14    0.090000  0.000000
      15      15    0.090280  0.000000
      16      16    0.090115  0.000000

============= PiLib Variable =============
scc.U_mat, @t-sp, the Hubbard potential matrix
ORDER=    1, SIZE=[   11,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1   -0.498910  0.000000
       2       2   -0.498680  0.000000
       3       3   -0.498650  0.000000
       4       4    0.346530  0.000000
       5       5    0.353720  0.000000
       6       6   -0.498910  0.000000
       7       7   -0.498680  0.000000
       8       8   -0.498650  0.000000
       9       9    0.346530  0.000000
      10      10    0.353720  0.000000

============= PiLib Variable =============
scc.H_onsite, @t-sp, H_onsite(hop.onsite_E+hop.LS_mat+scc.U_mat)
ORDER=    1, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1   -0.798910  0.000000
       2       2   -0.798680  0.000000
       3       3   -0.798650  0.000000
       4       4    0.346530  0.000000
       5       5    0.353720  0.000000
       6       6   -0.798910  0.000000
       7       7   -0.798680  0.000000
       8       8   -0.798650  0.000000
       9       9    0.346530  0.000000
      10      10    0.353720  0.000000
      11      11   -0.300000  0.000000
      12      12   -0.300000  0.000000
      13      13   -0.300000  0.000000
      14      14   -0.300000  0.000000
      15      15   -0.300000  0.000000
      16      16   -0.300000  0.000000

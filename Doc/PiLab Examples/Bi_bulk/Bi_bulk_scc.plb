scc.HubU=[]                       // U for each state, [state_label, U] or blank
scc.Charge=...                    // charge of each state, 1x total state
[1*ones(1,2),0.5*ones(1,6),1*ones(1,2),0.5*ones(1,6)]   //6s(2)6p(3)
scc.Mixing=[1]                    // mixing parameter, 0~1
scc.Iteration=[30]                // maximal iterations, integer
scc.Converge=[1e-3]              // convergence criterion, real, at least < 0.1
scc.Mesh=[5,5,5]                // k-space mesh for Ef, (1x1,1x2,1x3), large for metal
scc.Temperature=[100]             // temperature for searching Ef, large for insulator
scc.Memory=['max']             // 'normal', 'max' or 'HDD' 

============= PiLib Variable =============
scc.E_Fermi, @full, the Fermi level
ORDER=   -3, SIZE=[    1,    1], TYPE=REAL

           1

   -8.944653

============= PiLib Variable =============
scc.E_gap, @full, the band gap
ORDER=   -1, SIZE=[    1,    1], TYPE=REAL

           1

    2.807355

============= PiLib Variable =============
scc.DM_out, @t-sp, the output density matrix
ORDER=    1, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1    0.098075  0.000000
       2       2    0.098075  0.000000
       3       3    0.050276  0.000000
       4       4    0.050276  0.000000
       5       5    0.051376  0.000000
       6       6    0.050276  0.000000
       7       7    0.050276  0.000000
       8       8    0.051376  0.000000
       9       9    0.098075  0.000000
      10      10    0.098075  0.000000
      11      11    0.050274  0.000000
      12      12    0.050274  0.000000
      13      13    0.051374  0.000000
      14      14    0.050274  0.000000
      15      15    0.050274  0.000000
      16      16    0.051374  0.000000

============= PiLib Variable =============
scc.U_mat, @t-sp, the Hubbard potential matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
scc.H_onsite, @t-sp, H_onsite(hop.onsite_E+hop.LS_mat+scc.U_mat)
ORDER=    1, SIZE=[   29,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1   -1.090610  0.000000
       2       2   -1.090610  0.000000
       3       3   -0.048610  0.000000
       3       4    0.000000  0.075000
       3       8   -0.075000  0.000000
       4       4   -0.048610  0.000000
       4       8    0.000000 -0.075000
       5       5   -0.048610  0.000000
       5       6    0.075000  0.000000
       5       7    0.000000  0.075000
       6       6   -0.048610  0.000000
       6       7    0.000000 -0.075000
       7       7   -0.048610  0.000000
       8       8   -0.048610  0.000000
       9       9   -1.090600  0.000000
      10      10   -1.090600  0.000000
      11      11   -0.048600  0.000000
      11      12    0.000000  0.075000
      11      16   -0.075000  0.000000
      12      12   -0.048600  0.000000
      12      16    0.000000 -0.075000
      13      13   -0.048600  0.000000
      13      14    0.075000  0.000000
      13      15    0.000000  0.075000
      14      14   -0.048600  0.000000
      14      15    0.000000 -0.075000
      15      15   -0.048600  0.000000
      16      16   -0.048600  0.000000

scc.HubU=[]        // U for each state, [state_label, U] or blank
scc.Charge=[0.5*ones(1,64)] // charge of each state, 1x total state
scc.Mixing=[1]                    	// mixing parameter, 0~1
scc.Iteration=[30]                	// maximal iterations, integer
scc.Converge=[10^-3]              	// convergence criterion, real, at least < 0.1
scc.Mesh=[300,1]                    // k-space mesh for Ef, (1x1,1x2,1x3), large for metal
scc.Temperature=[100]             	// temperature for searching Ef, large for insulator
scc.Memory=['max']       	  	// 'normal', 'max' or 'HDD' 

============= PiLib Variable =============
scc.E_Fermi, @full, the Fermi level
ORDER=   -7, SIZE=[    1,    1], TYPE=REAL

           1

    2.564362

============= PiLib Variable =============
scc.E_gap, @full, the band gap
ORDER=   -1, SIZE=[    1,    1], TYPE=REAL

           1

    2.170466

============= PiLib Variable =============
scc.DM_out, @t-sp, the output density matrix
ORDER=    1, SIZE=[   65,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       1       1    0.050000  0.000000
       2       2    0.050000  0.000000
       3       3    0.050000  0.000000
       4       4    0.050000  0.000000
       5       5    0.050000  0.000000
       6       6    0.050000  0.000000
       7       7    0.050000  0.000000
       8       8    0.050000  0.000000
       9       9    0.050000  0.000000
      10      10    0.050000  0.000000
      11      11    0.050000  0.000000
      12      12    0.050000  0.000000
      13      13    0.050000  0.000000
      14      14    0.050000  0.000000
      15      15    0.050000  0.000000
      16      16    0.050000  0.000000
      17      17    0.050000  0.000000
      18      18    0.050000  0.000000
      19      19    0.050000  0.000000
      20      20    0.050000  0.000000
      21      21    0.050000  0.000000
      22      22    0.050000  0.000000
      23      23    0.050000  0.000000
      24      24    0.050000  0.000000
      25      25    0.050000  0.000000
      26      26    0.050000  0.000000
      27      27    0.050000  0.000000
      28      28    0.050000  0.000000
      29      29    0.050000  0.000000
      30      30    0.050000  0.000000
      31      31    0.050000  0.000000
      32      32    0.050000  0.000000
      33      33    0.050000  0.000000
      34      34    0.050000  0.000000
      35      35    0.050000  0.000000
      36      36    0.050000  0.000000
      37      37    0.050000  0.000000
      38      38    0.050000  0.000000
      39      39    0.050000  0.000000
      40      40    0.050000  0.000000
      41      41    0.050000  0.000000
      42      42    0.050000  0.000000
      43      43    0.050000  0.000000
      44      44    0.050000  0.000000
      45      45    0.050000  0.000000
      46      46    0.050000  0.000000
      47      47    0.050000  0.000000
      48      48    0.050000  0.000000
      49      49    0.050000  0.000000
      50      50    0.050000  0.000000
      51      51    0.050000  0.000000
      52      52    0.050000  0.000000
      53      53    0.050000  0.000000
      54      54    0.050000  0.000000
      55      55    0.050000  0.000000
      56      56    0.050000  0.000000
      57      57    0.050000  0.000000
      58      58    0.050000  0.000000
      59      59    0.050000  0.000000
      60      60    0.050000  0.000000
      61      61    0.050000  0.000000
      62      62    0.050000  0.000000
      63      63    0.050000  0.000000
      64      64    0.050000  0.000000

============= PiLib Variable =============
scc.U_mat, @t-sp, the Hubbard potential matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000

============= PiLib Variable =============
scc.H_onsite, @t-sp, H_onsite(hop.onsite_E+hop.LS_mat+scc.U_mat)
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
scc.HubU=[]                       // U for each state, [state_label, U] or blank
scc.Charge=[0.5,0.5]              // charge of each state, 1x total state
scc.Mixing=[1]                    // mixing parameter, 0~1
scc.Iteration=[30]                // maximal iterations, integer
scc.Converge=[10^-3]              // convergence criterion, real, at least < 0.1
scc.Mesh=[15,15]                  // k-space mesh for Ef, (1x1,1x2,1x3), large for metal
scc.Temperature=[100]             // temperature for searching Ef, large for insulator
scc.Memory=['normal']       	  // 'normal', 'max' or 'HDD' 

============= PiLib Variable =============
scc.E_Fermi, @full, the Fermi level
ORDER=   -7, SIZE=[    1,    1], TYPE=REAL

           1

    5.000000

============= PiLib Variable =============
scc.E_gap, @full, the band gap
ORDER=    0, SIZE=[    1,    1], TYPE=REAL

           1

    1.039229

============= PiLib Variable =============
scc.DM_out, @t-sp, the output density matrix
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.500000  0.000000
       2       2    0.500000  0.000000

============= PiLib Variable =============
scc.U_mat, @t-sp, the Hubbard potential matrix
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000

============= PiLib Variable =============
scc.H_onsite, @t-sp, H_onsite(hop.onsite_E+hop.LS_mat+scc.U_mat)
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000

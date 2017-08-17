hop.SiteOrb=[1,0;2,0]                    // specify orbital of each site, nx3, [site, n, l]
hop.Order=[2]                      // order of nearest order coupling, 1x1 integer, must <= lat.Order
hop.SKint=...                      // specify SK parameters, nx7, [Orb1,Orb2,nn_order,ts,tp,td,tf]
[1,2,1,0.3,0,0,0 ; 1,2,2,1,0,0,0]
hop.LS=[0]                        // strength of LS coupling, 1x1, real
hop.Filter=[10^-3]                // fliter of small hopping elements, 1x1, real
hop.Basis=['c']                   // Basis of the hopping matrix, 'c', 's', 'rc', 'rs'
hop.SelState=[]                  // Input state labels to pick states, 1xn, integer
hop.OnsiteE=[]                     // Onsite energy of picked states, 1xn, real

============= PiLib Variable =============
hop.state_info_text, @full, [state_label, site, identifier, l, SubOrb_text]     
ORDER=    0, SIZE=[    4,    5], TYPE=STRING

1 # 1 # 1 # 0 # 1 S  s,d # 
2 # 1 # 1 # 0 # 2 S  s,u # 
3 # 2 # 2 # 0 # 1 S  s,d # 
4 # 2 # 2 # 0 # 2 S  s,u # 

============= PiLib Variable =============
hop.state_info, @full, [state_label, site, identifier, l, SubOrb] 
ORDER=    0, SIZE=[    4,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1       1       0       1
       2       1       1       0       2
       3       2       2       0       1
       4       2       2       0       2

============= PiLib Variable =============
hop.LS_mat, @t-sp, LS coupling matrix
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000

============= PiLib Variable =============
hop.onsite_E, @t-sp, onsite energy matrix
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000

============= PiLib Variable =============
hop.hop_size, @full, size of hop.hop_mat, [sublatt, size(hop.hop_mat(n))]
ORDER=    0, SIZE=[    2,    4], TYPE=INTEGER

       1       2       3       4

       1       4       4       2
       2       4       4       2

============= PiLib Variable =============
hop.hop_mat(1)(:,:,1), @a-sp, hop_mat between site-1 and its 1-th neighbor
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.300000  0.000000
       2       4    0.300000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,2), @a-sp, hop_mat between site-1 and its 2-th neighbor
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    1.000000  0.000000
       2       4    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,1), @a-sp, hop_mat between site-2 and its 1-th neighbor
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       3       1    0.300000  0.000000
       4       2    0.300000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,2), @a-sp, hop_mat between site-2 and its 2-th neighbor
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       3       1    1.000000  0.000000
       4       2    1.000000  0.000000

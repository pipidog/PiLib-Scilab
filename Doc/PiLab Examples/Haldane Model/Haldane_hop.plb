hop.SiteOrb=[1,1;2,1]           // specify orbital of each site, nx3, [site, n, l]
hop.Order=[2]                       // order of nearest order coupling, 1x1 integer, must <= lat.Order
hop.SKint=...                       // specify SK parameters, nx7, [Orb1,Orb2,nn_order,ts,tp,td,tf]
[1,2,1,0,1,0,0;...
1,1,2,0,0.1*exp(1*%i*%pi/2),0,0;...
2,2,2,0,0.1*exp(1*%i*%pi/2),0,0]
hop.LS=[0]                          // strength of LS coupling, 1x1, real
hop.Filter=[10^-3]                  // fliter of small hopping elements, 1x1, real
hop.Basis=['c']                     // Basis of the hopping matrix, 'c', 's', 'rc', 'rs'
hop.SelState=[3,9]            // Input state labels to pick states, 1xn, integer
hop.OnsiteE=[0,0]               // Onsite energy of picked states, 1xn, real
// Comment:
// To simulate Haldane, manually complex conjugate the value of hop.hop_mat
// of neighbor 5,6,9 for sublattice 1 and 4,7,8 for sublattice 2 
// To break parity symmetry, manually put hop.onsite_E with +E and -E
// on the two matrix elements.

============= PiLib Variable =============
hop.state_info_text, @full, [state_label, site, identifier, l, SubOrb_text] 
ORDER=    0, SIZE=[    2,    5], TYPE=STRING

1 # 1 # 1 # 1 # 5 P  z,d # 
2 # 2 # 2 # 1 # 5 P  z,d # 

============= PiLib Variable =============
hop.state_info, @full, [state_label, site, identifier, l, SubOrb] 
ORDER=    0, SIZE=[    2,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1       1       1       5
       2       2       2       1       5

============= PiLib Variable =============
hop.LS_mat, @t-sp, LS coupling matrix
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

	   2       2     0.000001  0.000000

============= PiLib Variable =============
hop.onsite_E, @t-sp, onsite energy matrix
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       1       1    0.000000  0.000000
	   2       2    0.000000  0.000000

============= PiLib Variable =============
hop.hop_size, @full, size of hop.hop_mat, [sublatt, size(hop.hop_mat(n))]
ORDER=    0, SIZE=[    2,    4], TYPE=INTEGER

       1       2       3       4

       1       2       2       9
       2       2       2       9

============= PiLib Variable =============
hop.hop_mat(1)(:,:,1), @a-sp, hop_mat between site-1 and its 1-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,2), @a-sp, hop_mat between site-1 and its 2-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,3), @a-sp, hop_mat between site-1 and its 3-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,4), @a-sp, hop_mat between site-1 and its 4-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  0.100000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,5), @a-sp, hop_mat between site-1 and its 5-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  -0.100000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,6), @a-sp, hop_mat between site-1 and its 6-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  -0.100000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,7), @a-sp, hop_mat between site-1 and its 7-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  0.100000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,8), @a-sp, hop_mat between site-1 and its 8-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  0.100000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,9), @a-sp, hop_mat between site-1 and its 9-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  -0.100000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,1), @a-sp, hop_mat between site-2 and its 1-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,2), @a-sp, hop_mat between site-2 and its 2-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,3), @a-sp, hop_mat between site-2 and its 3-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    1.000000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,4), @a-sp, hop_mat between site-2 and its 4-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  -0.100000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,5), @a-sp, hop_mat between site-2 and its 5-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  0.100000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,6), @a-sp, hop_mat between site-2 and its 6-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  0.100000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,7), @a-sp, hop_mat between site-2 and its 7-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  -0.100000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,8), @a-sp, hop_mat between site-2 and its 8-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  -0.100000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,9), @a-sp, hop_mat between site-2 and its 9-th neighbor
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  0.100000

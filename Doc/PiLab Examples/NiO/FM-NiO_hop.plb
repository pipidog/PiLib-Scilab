hop.SiteOrb=[1,2;2,1]         // specify orbital of each site, nx3, [site, l]
hop.Order=[1]                     // order of nearest order coupling, 1x1 integer, must <= lat.Order
hop.SKint=[1,2,1,2,1,0,0]         // specify SK parameters, nx7, [Orb1,Orb2,nn_order,ts,tp,td,tf]
hop.LS=[0]                        // strength of LS coupling, 1x1, real
hop.Filter=[10^-3]                // fliter of small hopping elements, 1x1, real
hop.Basis=['c']                   // Basis of the hopping matrix, 'c', 's', 'rc', 'rs'
hop.SelState=[1:16]              // Input state labels to pick states, 1xn, integer
hop.OnsiteE=...                   // Onsite energy of picked states, 1xn, real
[-3,-3,-3,0,0,-3,-3,-3,0,0,-3,-3,-3,-3,-3,-3]

============= PiLib Variable =============
hop.state_info_text, @full, [state_label, site, identifier, l, SubOrb_text]     
ORDER=    0, SIZE=[   16,    5], TYPE=STRING

1 # 1 # 1 # 2 # 9 D  xy,d # 
2 # 1 # 1 # 2 # 10 D  yz,d # 
3 # 1 # 1 # 2 # 11 D  zx,d # 
4 # 1 # 1 # 2 # 12 D  x2-y2,d # 
5 # 1 # 1 # 2 # 13 D  3z2-r2,d # 
6 # 1 # 1 # 2 # 14 D  xy,u # 
7 # 1 # 1 # 2 # 15 D  yz,u # 
8 # 1 # 1 # 2 # 16 D  zx,u # 
9 # 1 # 1 # 2 # 17 D  x2-y2,u # 
10 # 1 # 1 # 2 # 18 D  3z2-r2,u # 
11 # 2 # 2 # 1 # 3 P  x,d # 
12 # 2 # 2 # 1 # 4 P  y,d # 
13 # 2 # 2 # 1 # 5 P  z,d # 
14 # 2 # 2 # 1 # 6 P  x,u # 
15 # 2 # 2 # 1 # 7 P  y,u # 
16 # 2 # 2 # 1 # 8 P  z,u # 

============= PiLib Variable =============
hop.state_info, @full, [state_label, site, identifier, l, SubOrb] 
ORDER=    0, SIZE=[   16,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1       1       2       9
       2       1       1       2      10
       3       1       1       2      11
       4       1       1       2      12
       5       1       1       2      13
       6       1       1       2      14
       7       1       1       2      15
       8       1       1       2      16
       9       1       1       2      17
      10       1       1       2      18
      11       2       2       1       3
      12       2       2       1       4
      13       2       2       1       5
      14       2       2       1       6
      15       2       2       1       7
      16       2       2       1       8

============= PiLib Variable =============
hop.LS_mat, @t-sp, LS coupling matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
hop.onsite_E, @t-sp, onsite energy matrix
ORDER=    1, SIZE=[   13,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1   -0.300000  0.000000
       2       2   -0.300000  0.000000
       3       3   -0.300000  0.000000
       6       6   -0.300000  0.000000
       7       7   -0.300000  0.000000
       8       8   -0.300000  0.000000
      11      11   -0.300000  0.000000
      12      12   -0.300000  0.000000
      13      13   -0.300000  0.000000
      14      14   -0.300000  0.000000
      15      15   -0.300000  0.000000
      16      16   -0.300000  0.000000

============= PiLib Variable =============
hop.hop_size, @full, size of hop.hop_mat, [sublatt, size(hop.hop_mat(n))]
ORDER=    0, SIZE=[    2,    4], TYPE=INTEGER

       1       2       3       4

       1      16      16       6
       2      16      16       6

============= PiLib Variable =============
hop.hop_mat(1)(:,:,1), @a-sp, hop_mat between site-1 and its 1-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3      13    0.100000  0.000000
       4      11    0.173205  0.000000
       5      11   -0.100000  0.000000
       8      16    0.100000  0.000000
       9      14    0.173205  0.000000
      10      14   -0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,2), @a-sp, hop_mat between site-1 and its 2-th neighbor
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1      11    0.100000  0.000000
       2      13    0.100000  0.000000
       4      12   -0.173205  0.000000
       5      12   -0.100000  0.000000
       6      14    0.100000  0.000000
       7      16    0.100000  0.000000
       9      15   -0.173205  0.000000
      10      15   -0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,3), @a-sp, hop_mat between site-1 and its 3-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       2      12   -0.100000  0.000000
       3      11   -0.100000  0.000000
       5      13   -0.200000  0.000000
       7      15   -0.100000  0.000000
       8      14   -0.100000  0.000000
      10      16   -0.200000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,4), @a-sp, hop_mat between site-1 and its 4-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       2      12    0.100000  0.000000
       3      11    0.100000  0.000000
       5      13    0.200000  0.000000
       7      15    0.100000  0.000000
       8      14    0.100000  0.000000
      10      16    0.200000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,5), @a-sp, hop_mat between site-1 and its 5-th neighbor
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1      11   -0.100000  0.000000
       2      13   -0.100000  0.000000
       4      12    0.173205  0.000000
       5      12    0.100000  0.000000
       6      14   -0.100000  0.000000
       7      16   -0.100000  0.000000
       9      15    0.173205  0.000000
      10      15    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,6), @a-sp, hop_mat between site-1 and its 6-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3      13   -0.100000  0.000000
       4      11   -0.173205  0.000000
       5      11    0.100000  0.000000
       8      16   -0.100000  0.000000
       9      14   -0.173205  0.000000
      10      14    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,1), @a-sp, hop_mat between site-2 and its 1-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       4   -0.173205  0.000000
      11       5    0.100000  0.000000
      13       3   -0.100000  0.000000
      14       9   -0.173205  0.000000
      14      10    0.100000  0.000000
      16       8   -0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,2), @a-sp, hop_mat between site-2 and its 2-th neighbor
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       1   -0.100000  0.000000
      12       4    0.173205  0.000000
      12       5    0.100000  0.000000
      13       2   -0.100000  0.000000
      14       6   -0.100000  0.000000
      15       9    0.173205  0.000000
      15      10    0.100000  0.000000
      16       7   -0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,3), @a-sp, hop_mat between site-2 and its 3-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       3    0.100000  0.000000
      12       2    0.100000  0.000000
      13       5    0.200000  0.000000
      14       8    0.100000  0.000000
      15       7    0.100000  0.000000
      16      10    0.200000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,4), @a-sp, hop_mat between site-2 and its 4-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       3   -0.100000  0.000000
      12       2   -0.100000  0.000000
      13       5   -0.200000  0.000000
      14       8   -0.100000  0.000000
      15       7   -0.100000  0.000000
      16      10   -0.200000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,5), @a-sp, hop_mat between site-2 and its 5-th neighbor
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       1    0.100000  0.000000
      12       4   -0.173205  0.000000
      12       5   -0.100000  0.000000
      13       2    0.100000  0.000000
      14       6    0.100000  0.000000
      15       9   -0.173205  0.000000
      15      10   -0.100000  0.000000
      16       7    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,6), @a-sp, hop_mat between site-2 and its 6-th neighbor
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       4    0.173205  0.000000
      11       5   -0.100000  0.000000
      13       3    0.100000  0.000000
      14       9    0.173205  0.000000
      14      10   -0.100000  0.000000
      16       8    0.100000  0.000000

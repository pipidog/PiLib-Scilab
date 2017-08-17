hop.SiteOrb=...                    // specify orbital of each site, nx3, [site, n, l]
[1,0;2,0;3,0;4,0;...
5,0;6,0;7,0;8,0]
hop.Order=[1]                     // order of nearest order coupling, 1x1 integer, must <= lat.Order
hop.SKint=...                      // specify SK parameters, nx7, [Orb1,Orb2,nn_order,ts,tp,td,tf]
[ 1,2,1,1,0,0,0;...       		
  2,3,1,1,0,0,0;...
  3,4,1,1,0,0,0;...
  4,5,1,1,0,0,0;...
  5,6,1,1,0,0,0;...
  6,7,1,1,0,0,0;...
  7,8,1,1,0,0,0]
hop.LS=[0]                        // strength of LS coupling, 1x1, real
hop.Filiter=[10^-3]                // fliter of small hopping elements, 1x1, real
hop.Basis=['c']                   // Basis of the hopping matrix, 'c', 's', 'rc', 'rs'
hop.SelState=[]                  // Input state labels to pick states, 1xn, integer
hop.OnsiteE=[]                     // Onsite energy of picked states, 1xn, real

============= PiLib Variable =============
hop.state_info_text, @full, [state_label, site, identifier, l, SubOrb_text]     
ORDER=    0, SIZE=[   16,    5], TYPE=STRING

1 # 1 # 1 # 0 # 1 S  s,d # 
2 # 1 # 1 # 0 # 2 S  s,u # 
3 # 2 # 2 # 0 # 1 S  s,d # 
4 # 2 # 2 # 0 # 2 S  s,u # 
5 # 3 # 3 # 0 # 1 S  s,d # 
6 # 3 # 3 # 0 # 2 S  s,u # 
7 # 4 # 4 # 0 # 1 S  s,d # 
8 # 4 # 4 # 0 # 2 S  s,u # 
9 # 5 # 5 # 0 # 1 S  s,d # 
10 # 5 # 5 # 0 # 2 S  s,u # 
11 # 6 # 6 # 0 # 1 S  s,d # 
12 # 6 # 6 # 0 # 2 S  s,u # 
13 # 7 # 7 # 0 # 1 S  s,d # 
14 # 7 # 7 # 0 # 2 S  s,u # 
15 # 8 # 8 # 0 # 1 S  s,d # 
16 # 8 # 8 # 0 # 2 S  s,u # 

============= PiLib Variable =============
hop.state_info, @full, [state_label, site, identifier, l, SubOrb] 
ORDER=    0, SIZE=[   16,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1       1       0       1
       2       1       1       0       2
       3       2       2       0       1
       4       2       2       0       2
       5       3       3       0       1
       6       3       3       0       2
       7       4       4       0       1
       8       4       4       0       2
       9       5       5       0       1
      10       5       5       0       2
      11       6       6       0       1
      12       6       6       0       2
      13       7       7       0       1
      14       7       7       0       2
      15       8       8       0       1
      16       8       8       0       2

============= PiLib Variable =============
hop.LS_mat, @t-sp, LS coupling matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
hop.onsite_E, @t-sp, onsite energy matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
hop.hop_size, @full, size of hop.hop_mat, [sublatt, size(hop.hop_mat(n))]
ORDER=    0, SIZE=[    8,    4], TYPE=INTEGER

       1       2       3       4

       1      16      16       2
       2      16      16       3
       3      16      16       3
       4      16      16       3
       5      16      16       3
       6      16      16       3
       7      16      16       3
       8      16      16       2

============= PiLib Variable =============
hop.hop_mat(1)(:,:,1), @a-sp, hop_mat between site-1 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       3    0.100000  0.000000
       2       4    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,2), @a-sp, hop_mat between site-1 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       3    0.100000  0.000000
       2       4    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,1), @a-sp, hop_mat between site-2 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       1    0.100000  0.000000
       4       2    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,2), @a-sp, hop_mat between site-2 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       1    0.100000  0.000000
       4       2    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,3), @a-sp, hop_mat between site-2 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       5    0.100000  0.000000
       4       6    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(3)(:,:,1), @a-sp, hop_mat between site-3 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       3    0.100000  0.000000
       6       4    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(3)(:,:,2), @a-sp, hop_mat between site-3 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       7    0.100000  0.000000
       6       8    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(3)(:,:,3), @a-sp, hop_mat between site-3 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       7    0.100000  0.000000
       6       8    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(4)(:,:,1), @a-sp, hop_mat between site-4 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       5    0.100000  0.000000
       8       6    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(4)(:,:,2), @a-sp, hop_mat between site-4 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       5    0.100000  0.000000
       8       6    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(4)(:,:,3), @a-sp, hop_mat between site-4 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       9    0.100000  0.000000
       8      10    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(5)(:,:,1), @a-sp, hop_mat between site-5 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9       7    0.100000  0.000000
      10       8    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(5)(:,:,2), @a-sp, hop_mat between site-5 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9      11    0.100000  0.000000
      10      12    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(5)(:,:,3), @a-sp, hop_mat between site-5 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9      11    0.100000  0.000000
      10      12    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(6)(:,:,1), @a-sp, hop_mat between site-6 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       9    0.100000  0.000000
      12      10    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(6)(:,:,2), @a-sp, hop_mat between site-6 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       9    0.100000  0.000000
      12      10    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(6)(:,:,3), @a-sp, hop_mat between site-6 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11      13    0.100000  0.000000
      12      14    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(7)(:,:,1), @a-sp, hop_mat between site-7 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      11    0.100000  0.000000
      14      12    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(7)(:,:,2), @a-sp, hop_mat between site-7 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      15    0.100000  0.000000
      14      16    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(7)(:,:,3), @a-sp, hop_mat between site-7 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      15    0.100000  0.000000
      14      16    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(8)(:,:,1), @a-sp, hop_mat between site-8 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      15      13    0.100000  0.000000
      16      14    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(8)(:,:,2), @a-sp, hop_mat between site-8 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      15      13    0.100000  0.000000
      16      14    0.100000  0.000000

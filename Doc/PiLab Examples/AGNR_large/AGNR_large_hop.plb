hop.SiteOrb=...                    // specify orbital of each site, nx3, [site, n, l]
[1,0;2,0;3,0;4,0;...
5,0;6,0;7,0;8,0;...
9,0;10,0;11,0;12,0;...
13,0;14,0;15,0;16,0;...
17,0;18,0;19,0;20,0;...
21,0;22,0;23,0;24,0;...
25,0;26,0;27,0;28,0;...
29,0;30,0;31,0;32,0]
hop.Order=[1]                     // order of nearest order coupling, 1x1 integer, must <= lat.Order
hop.SKint=...                      // specify SK parameters, nx7, [Orb1,Orb2,nn_order,ts,tp,td,tf]
[ 1,2,1,1,0,0,0;...
  1,16,1,1,0,0,0;... // y PBC
  1,17,1,1,0,0,0;...
  2,3,1,1,0,0,0;...
  2,18,1,1,0,0,0;...
  3,4,1,1,0,0,0;...
  3,19,1,1,0,0,0;...
  4,5,1,1,0,0,0;...
  4,20,1,1,0,0,0;...
  5,6,1,1,0,0,0;...  
  5,21,1,1,0,0,0;...
  6,7,1,1,0,0,0;...
  6,22,1,1,0,0,0;...
  7,8,1,1,0,0,0;...
  7,23,1,1,0,0,0;...
  8,9,1,1,0,0,0;...
  8,24,1,1,0,0,0;...
  9,10,1,1,0,0,0;...
  9,25,1,1,0,0,0;...
  10,11,1,1,0,0,0;...  
  10,26,1,1,0,0,0;...
  11,12,1,1,0,0,0;...
  11,27,1,1,0,0,0;...
  12,13,1,1,0,0,0;...
  12,28,1,1,0,0,0;...
  13,14,1,1,0,0,0;...
  13,29,1,1,0,0,0;...
  14,15,1,1,0,0,0;...
  14,30,1,1,0,0,0;...
  15,16,1,1,0,0,0;...
  15,31,1,1,0,0,0;...
  16,32,1,1,0,0,0;...
  17,18,1,1,0,0,0;...
  17,32,1,1,0,0,0;... // y PBC
  18,19,1,1,0,0,0;...
  19,20,1,1,0,0,0;...
  20,21,1,1,0,0,0;...
  21,22,1,1,0,0,0;...
  22,23,1,1,0,0,0;...
  23,24,1,1,0,0,0;...  
  24,25,1,1,0,0,0;...
  25,26,1,1,0,0,0;...
  26,27,1,1,0,0,0;...
  27,28,1,1,0,0,0;...
  28,29,1,1,0,0,0;...
  29,30,1,1,0,0,0;...
  30,31,1,1,0,0,0;...
  31,32,1,1,0,0,0...
]
hop.LS=[0]                        // strength of LS coupling, 1x1, real
hop.Filiter=[10^-3]                // fliter of small hopping elements, 1x1, real
hop.Basis=['c']                   // Basis of the hopping matrix, 'c', 's', 'rc', 'rs'
hop.SelState=[]                  // Input state labels to pick states, 1xn, integer
hop.OnsiteE=[]                     // Onsite energy of picked states, 1xn, real

============= PiLib Variable =============
hop.state_info_text, @full, [state_label, site, identifier, l, SubOrb_text]     
ORDER=    0, SIZE=[   64,    5], TYPE=STRING

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
17 # 9 # 9 # 0 # 1 S  s,d # 
18 # 9 # 9 # 0 # 2 S  s,u # 
19 # 10 # 10 # 0 # 1 S  s,d # 
20 # 10 # 10 # 0 # 2 S  s,u # 
21 # 11 # 11 # 0 # 1 S  s,d # 
22 # 11 # 11 # 0 # 2 S  s,u # 
23 # 12 # 12 # 0 # 1 S  s,d # 
24 # 12 # 12 # 0 # 2 S  s,u # 
25 # 13 # 13 # 0 # 1 S  s,d # 
26 # 13 # 13 # 0 # 2 S  s,u # 
27 # 14 # 14 # 0 # 1 S  s,d # 
28 # 14 # 14 # 0 # 2 S  s,u # 
29 # 15 # 15 # 0 # 1 S  s,d # 
30 # 15 # 15 # 0 # 2 S  s,u # 
31 # 16 # 16 # 0 # 1 S  s,d # 
32 # 16 # 16 # 0 # 2 S  s,u # 
33 # 17 # 17 # 0 # 1 S  s,d # 
34 # 17 # 17 # 0 # 2 S  s,u # 
35 # 18 # 18 # 0 # 1 S  s,d # 
36 # 18 # 18 # 0 # 2 S  s,u # 
37 # 19 # 19 # 0 # 1 S  s,d # 
38 # 19 # 19 # 0 # 2 S  s,u # 
39 # 20 # 20 # 0 # 1 S  s,d # 
40 # 20 # 20 # 0 # 2 S  s,u # 
41 # 21 # 21 # 0 # 1 S  s,d # 
42 # 21 # 21 # 0 # 2 S  s,u # 
43 # 22 # 22 # 0 # 1 S  s,d # 
44 # 22 # 22 # 0 # 2 S  s,u # 
45 # 23 # 23 # 0 # 1 S  s,d # 
46 # 23 # 23 # 0 # 2 S  s,u # 
47 # 24 # 24 # 0 # 1 S  s,d # 
48 # 24 # 24 # 0 # 2 S  s,u # 
49 # 25 # 25 # 0 # 1 S  s,d # 
50 # 25 # 25 # 0 # 2 S  s,u # 
51 # 26 # 26 # 0 # 1 S  s,d # 
52 # 26 # 26 # 0 # 2 S  s,u # 
53 # 27 # 27 # 0 # 1 S  s,d # 
54 # 27 # 27 # 0 # 2 S  s,u # 
55 # 28 # 28 # 0 # 1 S  s,d # 
56 # 28 # 28 # 0 # 2 S  s,u # 
57 # 29 # 29 # 0 # 1 S  s,d # 
58 # 29 # 29 # 0 # 2 S  s,u # 
59 # 30 # 30 # 0 # 1 S  s,d # 
60 # 30 # 30 # 0 # 2 S  s,u # 
61 # 31 # 31 # 0 # 1 S  s,d # 
62 # 31 # 31 # 0 # 2 S  s,u # 
63 # 32 # 32 # 0 # 1 S  s,d # 
64 # 32 # 32 # 0 # 2 S  s,u # 

============= PiLib Variable =============
hop.state_info, @full, [state_label, site, identifier, l, SubOrb] 
ORDER=    0, SIZE=[   64,    5], TYPE=INTEGER

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
      17       9       9       0       1
      18       9       9       0       2
      19      10      10       0       1
      20      10      10       0       2
      21      11      11       0       1
      22      11      11       0       2
      23      12      12       0       1
      24      12      12       0       2
      25      13      13       0       1
      26      13      13       0       2
      27      14      14       0       1
      28      14      14       0       2
      29      15      15       0       1
      30      15      15       0       2
      31      16      16       0       1
      32      16      16       0       2
      33      17      17       0       1
      34      17      17       0       2
      35      18      18       0       1
      36      18      18       0       2
      37      19      19       0       1
      38      19      19       0       2
      39      20      20       0       1
      40      20      20       0       2
      41      21      21       0       1
      42      21      21       0       2
      43      22      22       0       1
      44      22      22       0       2
      45      23      23       0       1
      46      23      23       0       2
      47      24      24       0       1
      48      24      24       0       2
      49      25      25       0       1
      50      25      25       0       2
      51      26      26       0       1
      52      26      26       0       2
      53      27      27       0       1
      54      27      27       0       2
      55      28      28       0       1
      56      28      28       0       2
      57      29      29       0       1
      58      29      29       0       2
      59      30      30       0       1
      60      30      30       0       2
      61      31      31       0       1
      62      31      31       0       2
      63      32      32       0       1
      64      32      32       0       2

============= PiLib Variable =============
hop.LS_mat, @t-sp, LS coupling matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000

============= PiLib Variable =============
hop.onsite_E, @t-sp, onsite energy matrix
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000

============= PiLib Variable =============
hop.hop_size, @full, size of hop.hop_mat, [sublatt, size(hop.hop_mat(n))]
ORDER=    0, SIZE=[   32,    4], TYPE=INTEGER

       1       2       3       4

       1      64      64       2
       2      64      64       3
       3      64      64       3
       4      64      64       3
       5      64      64       3
       6      64      64       3
       7      64      64       3
       8      64      64       3
       9      64      64       3
      10      64      64       3
      11      64      64       3
      12      64      64       3
      13      64      64       3
      14      64      64       3
      15      64      64       3
      16      64      64       2
      17      64      64       2
      18      64      64       3
      19      64      64       3
      20      64      64       3
      21      64      64       3
      22      64      64       3
      23      64      64       3
      24      64      64       3
      25      64      64       3
      26      64      64       3
      27      64      64       3
      28      64      64       3
      29      64      64       3
      30      64      64       3
      31      64      64       3
      32      64      64       2

============= PiLib Variable =============
hop.hop_mat(1)(:,:,1), @a-sp, hop_mat between site-1 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       1       3    0.100000  0.000000
       2       4    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(1)(:,:,2), @a-sp, hop_mat between site-1 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       1      33    0.100000  0.000000
       2      34    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,1), @a-sp, hop_mat between site-2 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       3       1    0.100000  0.000000
       4       2    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,2), @a-sp, hop_mat between site-2 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       3       5    0.100000  0.000000
       4       6    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(2)(:,:,3), @a-sp, hop_mat between site-2 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       3      35    0.100000  0.000000
       4      36    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(3)(:,:,1), @a-sp, hop_mat between site-3 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       5       3    0.100000  0.000000
       6       4    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(3)(:,:,2), @a-sp, hop_mat between site-3 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       5      37    0.100000  0.000000
       6      38    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(3)(:,:,3), @a-sp, hop_mat between site-3 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       5       7    0.100000  0.000000
       6       8    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(4)(:,:,1), @a-sp, hop_mat between site-4 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       7       9    0.100000  0.000000
       8      10    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(4)(:,:,2), @a-sp, hop_mat between site-4 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       7      39    0.100000  0.000000
       8      40    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(4)(:,:,3), @a-sp, hop_mat between site-4 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       7       5    0.100000  0.000000
       8       6    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(5)(:,:,1), @a-sp, hop_mat between site-5 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       9       7    0.100000  0.000000
      10       8    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(5)(:,:,2), @a-sp, hop_mat between site-5 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       9      11    0.100000  0.000000
      10      12    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(5)(:,:,3), @a-sp, hop_mat between site-5 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
       9      41    0.100000  0.000000
      10      42    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(6)(:,:,1), @a-sp, hop_mat between site-6 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      11       9    0.100000  0.000000
      12      10    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(6)(:,:,2), @a-sp, hop_mat between site-6 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      11      43    0.100000  0.000000
      12      44    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(6)(:,:,3), @a-sp, hop_mat between site-6 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      11      13    0.100000  0.000000
      12      14    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(7)(:,:,1), @a-sp, hop_mat between site-7 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      13      15    0.100000  0.000000
      14      16    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(7)(:,:,2), @a-sp, hop_mat between site-7 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      13      45    0.100000  0.000000
      14      46    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(7)(:,:,3), @a-sp, hop_mat between site-7 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      13      11    0.100000  0.000000
      14      12    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(8)(:,:,1), @a-sp, hop_mat between site-8 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      15      13    0.100000  0.000000
      16      14    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(8)(:,:,2), @a-sp, hop_mat between site-8 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      15      17    0.100000  0.000000
      16      18    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(8)(:,:,3), @a-sp, hop_mat between site-8 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      15      47    0.100000  0.000000
      16      48    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(9)(:,:,1), @a-sp, hop_mat between site-9 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      17      15    0.100000  0.000000
      18      16    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(9)(:,:,2), @a-sp, hop_mat between site-9 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      17      19    0.100000  0.000000
      18      20    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(9)(:,:,3), @a-sp, hop_mat between site-9 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      17      49    0.100000  0.000000
      18      50    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(10)(:,:,1), @a-sp, hop_mat between site-10 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      19      17    0.100000  0.000000
      20      18    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(10)(:,:,2), @a-sp, hop_mat between site-10 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      19      21    0.100000  0.000000
      20      22    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(10)(:,:,3), @a-sp, hop_mat between site-10 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      19      51    0.100000  0.000000
      20      52    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(11)(:,:,1), @a-sp, hop_mat between site-11 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      21      19    0.100000  0.000000
      22      20    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(11)(:,:,2), @a-sp, hop_mat between site-11 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      21      53    0.100000  0.000000
      22      54    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(11)(:,:,3), @a-sp, hop_mat between site-11 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      21      23    0.100000  0.000000
      22      24    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(12)(:,:,1), @a-sp, hop_mat between site-12 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      23      55    0.100000  0.000000
      24      56    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(12)(:,:,2), @a-sp, hop_mat between site-12 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      23      21    0.100000  0.000000
      24      22    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(12)(:,:,3), @a-sp, hop_mat between site-12 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      23      25    0.100000  0.000000
      24      26    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(13)(:,:,1), @a-sp, hop_mat between site-13 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      25      27    0.100000  0.000000
      26      28    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(13)(:,:,2), @a-sp, hop_mat between site-13 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      25      57    0.100000  0.000000
      26      58    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(13)(:,:,3), @a-sp, hop_mat between site-13 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      25      23    0.100000  0.000000
      26      24    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(14)(:,:,1), @a-sp, hop_mat between site-14 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      27      25    0.100000  0.000000
      28      26    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(14)(:,:,2), @a-sp, hop_mat between site-14 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      27      59    0.100000  0.000000
      28      60    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(14)(:,:,3), @a-sp, hop_mat between site-14 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      27      29    0.100000  0.000000
      28      30    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(15)(:,:,1), @a-sp, hop_mat between site-15 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      29      31    0.100000  0.000000
      30      32    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(15)(:,:,2), @a-sp, hop_mat between site-15 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      29      61    0.100000  0.000000
      30      62    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(15)(:,:,3), @a-sp, hop_mat between site-15 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      29      27    0.100000  0.000000
      30      28    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(16)(:,:,1), @a-sp, hop_mat between site-16 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      31      29    0.100000  0.000000
      32      30    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(16)(:,:,2), @a-sp, hop_mat between site-16 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      31      63    0.100000  0.000000
      32      64    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(17)(:,:,1), @a-sp, hop_mat between site-17 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      33       1    0.100000  0.000000
      34       2    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(17)(:,:,2), @a-sp, hop_mat between site-17 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      33      35    0.100000  0.000000
      34      36    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(18)(:,:,1), @a-sp, hop_mat between site-18 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      35       3    0.100000  0.000000
      36       4    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(18)(:,:,2), @a-sp, hop_mat between site-18 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      35      33    0.100000  0.000000
      36      34    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(18)(:,:,3), @a-sp, hop_mat between site-18 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      35      37    0.100000  0.000000
      36      38    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(19)(:,:,1), @a-sp, hop_mat between site-19 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      37       5    0.100000  0.000000
      38       6    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(19)(:,:,2), @a-sp, hop_mat between site-19 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      37      35    0.100000  0.000000
      38      36    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(19)(:,:,3), @a-sp, hop_mat between site-19 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      37      39    0.100000  0.000000
      38      40    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(20)(:,:,1), @a-sp, hop_mat between site-20 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      39      41    0.100000  0.000000
      40      42    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(20)(:,:,2), @a-sp, hop_mat between site-20 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      39       7    0.100000  0.000000
      40       8    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(20)(:,:,3), @a-sp, hop_mat between site-20 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      39      37    0.100000  0.000000
      40      38    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(21)(:,:,1), @a-sp, hop_mat between site-21 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      41      39    0.100000  0.000000
      42      40    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(21)(:,:,2), @a-sp, hop_mat between site-21 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      41      43    0.100000  0.000000
      42      44    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(21)(:,:,3), @a-sp, hop_mat between site-21 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      41       9    0.100000  0.000000
      42      10    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(22)(:,:,1), @a-sp, hop_mat between site-22 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      43      41    0.100000  0.000000
      44      42    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(22)(:,:,2), @a-sp, hop_mat between site-22 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      43      11    0.100000  0.000000
      44      12    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(22)(:,:,3), @a-sp, hop_mat between site-22 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      43      45    0.100000  0.000000
      44      46    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(23)(:,:,1), @a-sp, hop_mat between site-23 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      45      47    0.100000  0.000000
      46      48    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(23)(:,:,2), @a-sp, hop_mat between site-23 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      45      13    0.100000  0.000000
      46      14    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(23)(:,:,3), @a-sp, hop_mat between site-23 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      45      43    0.100000  0.000000
      46      44    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(24)(:,:,1), @a-sp, hop_mat between site-24 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      47      45    0.100000  0.000000
      48      46    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(24)(:,:,2), @a-sp, hop_mat between site-24 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      47      49    0.100000  0.000000
      48      50    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(24)(:,:,3), @a-sp, hop_mat between site-24 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      47      15    0.100000  0.000000
      48      16    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(25)(:,:,1), @a-sp, hop_mat between site-25 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      49      47    0.100000  0.000000
      50      48    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(25)(:,:,2), @a-sp, hop_mat between site-25 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      49      51    0.100000  0.000000
      50      52    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(25)(:,:,3), @a-sp, hop_mat between site-25 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      49      17    0.100000  0.000000
      50      18    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(26)(:,:,1), @a-sp, hop_mat between site-26 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      51      49    0.100000  0.000000
      52      50    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(26)(:,:,2), @a-sp, hop_mat between site-26 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      51      53    0.100000  0.000000
      52      54    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(26)(:,:,3), @a-sp, hop_mat between site-26 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      51      19    0.100000  0.000000
      52      20    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(27)(:,:,1), @a-sp, hop_mat between site-27 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      53      51    0.100000  0.000000
      54      52    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(27)(:,:,2), @a-sp, hop_mat between site-27 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      53      21    0.100000  0.000000
      54      22    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(27)(:,:,3), @a-sp, hop_mat between site-27 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      53      55    0.100000  0.000000
      54      56    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(28)(:,:,1), @a-sp, hop_mat between site-28 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      55      23    0.100000  0.000000
      56      24    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(28)(:,:,2), @a-sp, hop_mat between site-28 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      55      53    0.100000  0.000000
      56      54    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(28)(:,:,3), @a-sp, hop_mat between site-28 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      55      57    0.100000  0.000000
      56      58    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(29)(:,:,1), @a-sp, hop_mat between site-29 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      57      59    0.100000  0.000000
      58      60    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(29)(:,:,2), @a-sp, hop_mat between site-29 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      57      25    0.100000  0.000000
      58      26    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(29)(:,:,3), @a-sp, hop_mat between site-29 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      57      55    0.100000  0.000000
      58      56    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(30)(:,:,1), @a-sp, hop_mat between site-30 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      59      57    0.100000  0.000000
      60      58    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(30)(:,:,2), @a-sp, hop_mat between site-30 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      59      27    0.100000  0.000000
      60      28    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(30)(:,:,3), @a-sp, hop_mat between site-30 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      59      61    0.100000  0.000000
      60      62    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(31)(:,:,1), @a-sp, hop_mat between site-31 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      61      63    0.100000  0.000000
      62      64    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(31)(:,:,2), @a-sp, hop_mat between site-31 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      61      29    0.100000  0.000000
      62      30    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(31)(:,:,3), @a-sp, hop_mat between site-31 and its 3-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      61      59    0.100000  0.000000
      62      60    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(32)(:,:,1), @a-sp, hop_mat between site-32 and its 1-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      63      61    0.100000  0.000000
      64      62    0.100000  0.000000

============= PiLib Variable =============
hop.hop_mat(32)(:,:,2), @a-sp, hop_mat between site-32 and its 2-th neighbor
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      64      64    0.000000  0.000000
      63      31    0.100000  0.000000
      64      32    0.100000  0.000000
